#include "util.h"
#include <cstdio>
#include "param.h"
#include <utility> 
namespace Quantum 
{
  void scf(int iop, int n, double r, double zeta1, double zeta2,
           double za, double zb)
  {
    using namespace Quantum::Const;
    using namespace Quantum::Array;
    int cycle = 0;
    double q[2] = {0.0, 0.0};
    double en = 0.0, delta, ent = 0.0, N = 0.0;
    for (int i = 0, size = p.dimension(); i < size; ++i) {
      for(int j = 0, size = p.dimension(); j < size; ++j) {
        p.at(i,j) = 0.0;
      }
    }
    if(iop < 2) {
      return;//?
    }
    matout(p,"P");
    while(cycle <= maxcycle) {
      
      //start iteration loop
      cycle = cycle + 1;
      printf("\nCycle%3d:\n",cycle);
      //form two-electron part of fock matrix from P
      formg();
      matout(g,"G");
      //add core hamiltonian to get Fock matrix
      for(int i = 0, size = f.dimension(); i < size; ++i) {
        for(int j = 0, size = f.dimension(); j < size; ++j) {
          f.at(i,j) = h.at(i,j) + g.at(i,j);
        }
      }
      //calculate electronic energy
      en = 0.0;
      for(int i = 0, size = h.dimension(); i < size; ++i) {
        for(int j = 0, size = h.dimension(); j < size; ++j) {
          en = en + 0.5*p.at(i,j) * (h.at(i,j) + f.at(i,j));
        }
      }
      matout(f,"F");
      printf("\n  electronic energy = %20.12E\n", en);
      //F * x = G 
      mult(f,x,g);
      //xt * G = F'
      mult(xt,g,fprime);
      //diagonalize transformed Fock matrix
      // F' = λC' 
      diag(fprime,cprime,e);
      //transform eigenvectors to get matrix C
      // x * C' = C
      mult(x,cprime,c);
      //form new density matrix
      for(int i = 0, size = p.dimension(); i < size; ++i) {
        for(int j = 0, size = p.dimension(); j < size; ++j) {
          //save present density matrix before creating new one
          oldp.at(i,j) = p.at(i,j);
          
          p.at(i,j) = 0.0;
          for(int k = 0, size = c.dimension(); k < size/2; ++k) {
            p.at(i,j) = p.at(i,j) + 2.0*c.at(i,k)*c.at(j,k);
#ifdef DEBUG
            std::printf("[%s] oldp:%f p:%f\n", __FILE__, oldp.at(i,j), p.at(i,j));
#endif
          }
        }
      }
      matout(fprime,"F'");
      matout(cprime,"C'");
      matout(e,"E");
      matout(c,"C");
      matout(p,"P");
      //initialize every cycle: delta=.0f
      delta=.0f;
      for(int i = 0, size = p.dimension(); i < size; ++i) {
        for(int j = 0, size = p.dimension(); j < size; ++j) {
          delta = delta + pow(p.at(i,j)-oldp.at(i,j),2);
        }
      }
      delta = sqrt(delta/4.0);
      printf("delta(convergence of denstiy matrix) = %12.6E\n", delta);
    
      if(delta < conver) {
        ent = en + za*zb/r;
        printf("\nscf calculation is converged\n");
        printf("iteration number  = %4d\n", cycle);
        printf("electronic energy = %20.12E\n", en);
        printf("total energy      = %20.12E\n", ent);
        if(iop == 1) {
          matout(g,"G");
          matout(f,"F");
          matout(e,"E");
          matout(c,"C");
          matout(p,"P");
        }
        mult(p,s,oldp);
        if(iop == 0) {
          return;
        }
        matout(oldp,"PS");
        for(int i = 0, size = oldp.dimension() ; i < size; ++i) {
          N += oldp.at(i,i);
        }
        q[0] = za - oldp.at(0,0);
        q[1] = zb - oldp.at(1,1);
        printf("\nelectron     = %20.15f\n", N);
        printf("charge[%2.1f]  = %20.15f\n", za, q[0]);
        printf("charge[%2.1f]  = %20.15f\n", zb, q[1]);
        return;
      }
      if(cycle > maxcycle) {
        printf("scf claculation is not converged...\n");
        ent = en + za*zb/r;
        printf("iteration number  = %4d\n", cycle);
        printf("electronic energy = %20.12E\n", en);
        printf("total energy      = %20.12E\n", ent);
        mult(p,s,oldp);
        matout(oldp,"PS");
        break;
      }
    }
  }
}
