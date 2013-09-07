#include "util.h"
#include "param.h"
#include <cmath>
#ifdef DEBUG
#  include <cstdio>
#endif
namespace Quantum 
{
  double twoe(double a, double b, double c, double d,
              double rab2, double rcd2, double rpq2)
  {
    //calculate two-electron integrals     
    using namespace Quantum::Const;
    double twoe;
    twoe = 2.0*pow(pi,(5.0/2.0))/((a+b)*(c+d)*sqrt(a+b+c+d))
        * f0((a+b)*(c+d)*rpq2/(a+b+c+d))
        * exp(-a*b*rab2/(a+b)-c*d*rcd2/(c+d));
#ifdef DEBUG
    std::printf("[%s] a:%f b:%f c:%f d:%f rab2:%f rcd2:%f rpq2:%f\n",
                __FILE__, a, b, c, d, rab2, rcd2, rpq2);
#endif
    return twoe;
  }
}
