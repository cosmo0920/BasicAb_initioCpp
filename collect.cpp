#include <cmath>
#include "collect.h"
#include "param.h"
#include "util.h"
#include <cstdio>
namespace Quantum
{
  void colect(int iop,int n, double r, double& zeta1, double& zeta2,
              double za, double zb)
  {
    using namespace Quantum::Variable;
    using namespace Quantum::Array;
    int MatDim = 0;
    //assemble core hamiltonian
    h.at(0,0) = t11+v11a+v11b;
    h.at(0,1) = t12+v12a+v12b;
    h.at(1,0) = h.at(0,1);
    h.at(1,1) = t22+v22a+v22b;
    //assemble overlap matrix
    s.at(0,0) = 1.0;
    s.at(0,1) = s12;
    s.at(1,0) = s.at(0,1);
    s.at(1,1) = 1.0;
    //canonical orthogonalization
    x.at(0,0) = 1.0/sqrt(2.0*(1.0+s12));
    x.at(1,0) = x.at(0,0);
    x.at(0,1) = 1.0/sqrt(2.0*(1.0-s12));
    x.at(1,1) = -x.at(0,1);
    //transpose of transformation matrix
    xt.at(0,0) = x.at(0,0);
    xt.at(0,1) = x.at(1,0);
    xt.at(1,0) = x.at(0,1);
    xt.at(1,1) = x.at(1,1);
    //matrix of two-electron integrals
    tt.at(0,0,0,0) = v1111;
    tt.at(1,0,0,0) = v2111;
    tt.at(0,1,0,0) = v2111;
    tt.at(0,0,1,0) = v2111;
    tt.at(0,0,0,1) = v2111;
    tt.at(1,0,1,0) = v2121;
    tt.at(0,1,1,0) = v2121;
    tt.at(1,0,0,1) = v2121;
    tt.at(0,1,0,1) = v2121;
    tt.at(1,1,0,0) = v2211;
    tt.at(0,0,1,1) = v2211;
    tt.at(1,1,1,0) = v2221;
    tt.at(1,1,0,1) = v2221;
    tt.at(1,0,1,1) = v2221;
    tt.at(0,1,1,1) = v2221;
    tt.at(1,1,1,1) = v2222;
    if(iop == 0) {
      return;
    }
    matout(s,"S");
    matout(x,"X");
    matout(h,"H");
    MatDim = h.dimension();
    for(int i = 0; i < MatDim; ++i) {
      for(int j = 0; j < MatDim; ++j) {
        for(int k = 0; k < MatDim; ++k) {
          for(int l = 0; l < MatDim; ++l) {
            printf("  ( %d %d %d %d),  %f\n", i+1, j+1, k+1, l+1, tt.at(i,j,k,l));
          }
        }
      }
    }
  }
}
