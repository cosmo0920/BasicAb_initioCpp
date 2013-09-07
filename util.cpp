#include "param.h"
#include "util.h"
#include <cmath>
#if defined (__linux__)
#  include <f2c.h>
#  include "include/clapack.h"
#  undef abs
#endif
#include <cstdio>
#include <array>

namespace Quantum 
{
  double f0(double arg) 
  {    
    // calculates the f-functions
    using namespace Quantum::Const;
    double f0;
    if (arg<1.0e-6) {
      f0 = 1.0-arg/3.0;
#ifdef DEBUG
      std::printf("[%s] *arg=0*\n", __FILE__);
      std::printf("[%s] arg:%f f0:%f\n",__FILE__, arg, f0);
#endif      
      return f0;
    }
    f0 = sqrt(pi/arg)*derf(sqrt(arg))/2.0;
#ifdef DEBUG
    std::printf("[%s] arg:%f f0:%f\n",__FILE__, arg, f0);
#endif
    return f0;
  }
  
#ifdef USE_DERF
  double derf(double arg) 
  {
    //calculate error function
    using namespace Quantum::Const;
    const double p = 0.3275911;
    std::array<double, 5> a =
        {{ 0.254829592,-0.284496736,1.421413741,-1.453152027,1.061405429 }};
    double t, tn, poly, derf;
    t = 1.0/(1.0+p*arg);
    tn = t;
    poly = a.at(0)*tn;
    for (int i = 1, size = a.size() ; i < size; ++i) {
      tn = tn*t;
      poly = poly + a.at(i)*tn;
    } 
    derf = 1.0-poly*exp(-arg*arg);
    return derf;
  }
#else //normally use here
  double derf(double arg) 
  {   
    return std::erf(arg);
  }
#endif
  //template matout function
  template <typename T>
  void matout(T a, const char *label)
  {
    int size = a.dimension();
    printf("\n  The %.5s matrix\n", label);
    for(int i = 0; i < size; ++i) {
      printf("\t\t\t\%d", i+1);
    }
    printf("\n");
    for(int i = 0; i < size; i++) {
      printf("  %5.5d\t\t ", i+1);
      for(int j = 0; j < size; j++) {
        printf("%12.10E\t\t ",a.at(i,j));
      }
      printf("\n");
    }
  }
  //T = FixMatrix2d version matout
  void matout(FixMatrix2d a, const char *label) 
  {
    matout<FixMatrix2d>(a,label);
  }
  
  void formg()
  {
    using namespace Quantum::Variable;
    using namespace Quantum::Array;
    int size = g.dimension();
    for(int i = 0; i < size; ++i) {
      for(int j = 0; j < size; ++j) {
        g.at(i,j) = 0.0;
        for(int k = 0; k < size; ++k) {
          for(int l = 0; l < size; ++l) {
            g.at(i,j) = g.at(i,j)+p.at(k,l)*(tt.at(i,j,k,l)-0.5*tt.at(i,l,k,j));
          }
        }
      }
    }
  }
  //template mult function
  template <typename T>
  void mult(T a, T b, T &c)
  {
    int size = a.dimension();
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
        c.at(i,j) = 0.0;
        for(int k = 0; k < size; ++k) {
          c.at(i,j) = c.at(i,j) + a.at(i,k) * b.at(k,j);
        }
      }
    }
  }
  //T = FixMatrix2d version mult
  void mult(FixMatrix2d a, FixMatrix2d b, FixMatrix2d &c)
  {
    mult<FixMatrix2d>(a,b,c);
  }
  
  /*
   * diagolize Matrix F to give eigenvectors in Matrix C
   * and eigenvalues in Matrix E
   */
  void diag(FixMatrix2d f, FixMatrix2d &c, FixMatrix2d &e)
  {
    using namespace Quantum::Variable;
    using namespace Quantum::Array;
    using namespace Quantum::Const;
    double theta = 0.0, temp = 0.0;
    if(std::abs(f.at(0,0) - f.at(1,1)) > 1.0e-20) {
      //solution for heteronuclear diatomic
      theta = 0.50*atan(2.0*f.at(0,1)/(f.at(0,0)-f.at(1,1)));
    } else {
      //symmetry determined solution
      theta = pi/4.0;
    }
    c.at(0,0) =  cos(theta);
    c.at(1,0) =  sin(theta);
    c.at(0,1) = -sin(theta);
    c.at(1,1) =  cos(theta);
    e.at(0,0) = f.at(0,0)*pow(cos(theta),2)+f.at(1,1)*pow(sin(theta),2)
        + f.at(0,1)*sin(2.0*theta);
    e.at(1,1) = f.at(1,1)*pow(cos(theta),2)+f.at(0,0)*pow(sin(theta),2)
        - f.at(0,1)*sin(2.0*theta);
    e.at(1,0) = e.at(0,1) = 0.0;
    if(e.at(1,1) > e.at(0,0)) {
      return;
    }
        
    //order eigenvalues and eigenvectors
    temp = e.at(1,1);
    e.at(1,1) = e.at(0,0);
    e.at(0,0) = temp;
    temp = c.at(0,1);
    c.at(0,1) = c.at(0,0);
    c.at(0,0) = temp;
    temp = c.at(1,1);
    c.at(1,1) = c.at(1,0);
    c.at(1,0) = temp;
  }
}
