#include "param.h"
#include "util.h"
#include "one_integral.h"
#include <cmath>
#ifdef DEBUG
#  include <cstdio>
#endif
namespace Quantum 
{
  /*
   *  T = (ab/a+b)(3 - 2ab/(a+b) * (Ra-Rb)^2)
   *      *  (pi/a+b)^(3/2)exp(-ab/(a+b) * (Ra-Rb)^2)
   */
  double t(double a,double b,double rab2) 
  {
    using namespace Quantum::Const;
    double t;
    t = a*b/(a+b)*(3.0-2.0*a*b*rab2/(a+b))*(pow((pi/(a+b)),(3.0/2.0)))
                                            * exp(-a*b*rab2/(a+b));
#ifdef DEBUG
    std::printf("[%s]:t a:%f, b:%f, rab2:%f ,t:%f\n",__FILE__, a, b, rab2, t);
#endif
    return t;
  }
  
  /*
   * V = (-2pi*Zc/a+b)exp(-ab/(a+b) * (Ra-Rb)^2)*F0((a+b)(Rp-Rc)^2)
   * Rp = (a*Ra+b*Rb)/(a+b)
   */
  double v(double a, double b, double rab2, double rcp2, double zc)
  {
    using namespace Quantum::Const;
    double v;
    v = (-2.0*pi*zc)/(a+b)*f0((a+b)*rcp2)*exp(-a*b*rab2/(a+b));
 #ifdef DEBUG
    std::printf("[%s]:v a:%f, b:%f v:%f rab2:%f rcp2:%f\n",
                __FILE__,a, b, v, rab2, rcp2);
#endif
    return v;
  }

  /*
   * S = (pi/a+b)^(3/2)exp(-ab/(a+b) * (Ra-Rb)^2)
   */
  double s(double a, double b, double rab2)
  {
    //calculate overlaps for un-normalized primitives
    using namespace Quantum::Const;
    double s;
    s = (pow((pi/(a+b)),(3.0/2.0)))*exp(-a*b*rab2/(a+b));
#ifdef DEBUG
    std::printf("[%s]:s a:%f b:%f rab2:%f s:%f\n",__FILE__, a, b, rab2, s);
#endif
    return s;
  }  
}
