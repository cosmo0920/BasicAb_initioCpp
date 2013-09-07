#include <iostream>
#include "hfcalc.h" 
#if defined (__linux__)
#  include <fpu_control.h>
#endif
using namespace Quantum;
int main(void) 
{
#if defined (__linux__)
  // この2行を追加
  fpu_control_t control_word = _FPU_RC_NEAREST | _FPU_IEEE;
  _FPU_SETCW(control_word);
#endif
  const int iop = 2, numG = 3;
  double r = 1.4632, zeta1 = 2.0925, zeta2 = 1.24, za = 2.0, zb = 1.0;

  hfcalc(iop, numG, r, zeta1, zeta2, za ,zb);
  return 0;
}
