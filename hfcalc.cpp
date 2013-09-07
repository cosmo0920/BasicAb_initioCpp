#include <cstdio>
#include "hfcalc.h"
#include "integral.h"
#include "collect.h"
#include "scf.h"
namespace Quantum 
{
  void hfcalc(int iop, int n, double r,
              double& zeta1, double& zeta2, double za, double zb)
  {
    if(iop != 0) {
      printf("1  STO-%dG For Atomic numbers %5.2f and %5.2f\n", n, za, zb);
    }
    intgrl(iop,n,r,zeta1,zeta2,za,zb);
    colect(iop,n,r,zeta1,zeta2,za,zb);
    // do scf iterations
    scf(iop,n,r,zeta1,zeta2,za,zb);
  }
}

