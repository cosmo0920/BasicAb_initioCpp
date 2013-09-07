#ifndef __PARAM_H_
#define __PARAM_H_
#include <cmath>
#include "include/FixMatrix.hpp"
typedef FixMatrix<double, 2, 2>       FixMatrix2d;
typedef FixMatrix<double, 3, 3>       FixMatrix3d;
typedef FixMatrix<double, 2, 2, 2, 2> FixMatrix2x4d;

namespace Quantum 
{
  namespace Const 
  {
    extern const double pi;
    extern const double conver;
    extern const double micro, nano, pico;
    extern const int    maxcycle;
  }
  namespace Variable 
  {
    extern double r2, s12, t11, t12, t22, v11a, v12a, v22a, v11b, v12b, v22b,
        v1111, v2111, v2121, v2211, v2221, v2222;
  }
  namespace Array
  {
    extern FixMatrix2d s, x, xt, h, f, g, c,
        fprime ,cprime ,p ,oldp,e;
    extern FixMatrix2x4d tt;
  }
}
#endif //__PARAM_H_
