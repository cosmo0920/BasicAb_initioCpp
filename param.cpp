#include <cmath>
#include "param.h"

namespace Quantum 
{
  namespace SI
  {
    namespace prefix
    {
      const double micro = 1.0e-6, nano = 1.0e-9, pico = 1.0e-12;
    }
  }
  namespace Const 
  {
    const double pi        = 4.0 * atan(1.0);
    const double conver    = SI::prefix::micro;
    const int    maxcycle  = 100;
  }
  namespace Variable 
  {
    double r2, s12, t11, t12, t22, v11a, v12a, v22a, v11b, v12b, v22b,
        v1111, v2111, v2121, v2211, v2221, v2222;
  }
  namespace Array
  {
    FixMatrix2d s, x, xt, h, f, g, c,
        fprime ,cprime ,p ,oldp,e;
    FixMatrix2x4d tt;
  }
}
