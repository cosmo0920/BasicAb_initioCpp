#ifndef __UTIL_H_
#define __UTIL_H_
#include "param.h"
namespace Quantum 
{
  const int N = 2;
  double f0(double arg);
  double derf(double arg);
  void matout(FixMatrix2d a, const char *label);
  //void matout(FixMatrix2d a,int im,int in, const char *label);
  void formg();
  void mult(FixMatrix2d a, FixMatrix2d b, FixMatrix2d &c);
  void diag(FixMatrix2d f, FixMatrix2d &c, FixMatrix2d &e);
}
#endif //__UTIL_H_
