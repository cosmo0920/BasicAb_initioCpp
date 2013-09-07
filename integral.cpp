#include "param.h"
#include "one_integral.h"
#include "two_integral.h"
#include "integral.h"
#include <cstdio>

namespace Quantum 
{
  void intgrl(int iop, int n, double r,
              double& zeta1, double& zeta2, double za, double zb)
  {
    using namespace Quantum::Const;
    using namespace Quantum::Variable;
    FixMatrix3d coef, expon;
    double d1[3], a1[3], d2[3], a2[3];
    double rap, rap2, raq, raq2, rbp, rbp2, rbq, rbq2, rpq, rpq2;
    //initialize coef and expon
    coef.at(0,0) = 1.0;      coef.at(1,0) = coef.at(2,0) = 0.0;
    coef.at(0,1) = 0.678914; coef.at(1,1) = 0.430120; coef.at(2,1) = 0.0;
    coef.at(0,2) = 0.444635; coef.at(1,2) = 0.535328; coef.at(2,2) = 0.154329;
    
    expon.at(0,0) = 0.270950; expon.at(1,0) = expon.at(2,0) = 0.0;
    expon.at(0,1) = 0.151623; expon.at(1,1) = 0.851819; expon.at(2,1) = 0.0;
    expon.at(0,2) = 0.109818; expon.at(1,2) = 0.405771; expon.at(2,2) = 2.22766;

    r2 = r * r;
    for(int i = 0; i < n; ++i) {
      a1[i] = expon.at(i,n-1)*(zeta1*zeta1);
      d1[i] =  coef.at(i,n-1)*(pow((2.0*a1[i]/pi),0.75));
      a2[i] = expon.at(i,n-1)*(zeta2*zeta2);
      d2[i] =  coef.at(i,n-1)*(pow((2.0*a2[i]/pi),0.75));
#ifdef DEBUG
      std::printf("[%s] a1:%f d1:%f a2:%f d2:%f\n",
                  __FILE__, a1[i],d1[i], a2[i], d2[i]);
#endif
    }
    s12   = 0.0;
    t11   = 0.0; t12   = 0.0; t22   = 0.0;
    v11a  = 0.0; v12a  = 0.0; v22a  = 0.0;
    v11b  = 0.0; v12b  = 0.0; v22b  = 0.0;
    v1111 = 0.0; v2111 = 0.0; v2121 = 0.0;
    v2211 = 0.0; v2221 = 0.0; v2222 = 0.0;
    /* calculate 1-electron integrals
     * center a is first atom, center b is second atom
     * origin on center a
     * v12a = off-diagonal nuclear attraction to center a */
    for (int i = 0; i < n; ++i) {
      for(int j = 0; j < n; ++j) {
        // rap2 is squared distance between center a and center p
        rap  = a2[j]*r/(a1[i]+a2[j]);
        rap2 = rap*rap;
        rbp2 = (r-rap)*(r-rap);
        s12  = s12 + s(a1[i],a2[j],r2)*d1[i]*d2[j];
        t11  = t11+t(a1[i],a1[j],0.0)*d1[i]*d1[j];
        t12  = t12+t(a1[i],a2[j], r2)*d1[i]*d2[j];
        t22  = t22+t(a2[i],a2[j],0.0)*d2[i]*d2[j];
        v11a = v11a+v(a1[i],a1[j],0.0,0.0,za)*d1[i]*d1[j];
        v12a = v12a+v(a1[i],a2[j],r2,rap2,za)*d1[i]*d2[j];
        //v12a = -1.102913;
        v22a = v22a+v(a2[i],a2[j],0.0, r2,za)*d2[i]*d2[j];
        v11b = v11b+v(a1[i],a1[j],0.0, r2,zb)*d1[i]*d1[j];
        v12b = v12b+v(a1[i],a2[j],r2,rbp2,zb)*d1[i]*d2[j];
        v22b = v22b+v(a2[i],a2[j],0.0,0.0,zb)*d2[i]*d2[j];
#ifdef DEBUG
        std::printf("[%s] a1[%d]:%f, a2[%d]:%f\n",__FILE__, i, a1[i], j, a2[j]);
        std::printf("[%s] rap2:%f\n",__FILE__, rap2);
        std::printf("[%s] rbp2:%f\n",__FILE__, rbp2);
        std::printf("[%s] r2:%f za:%f d1[%d]:%f d2[%d]:%f v12a=%f\n",
                    __FILE__, r2, za, i, d1[i], j, d2[j], v12a);
#endif
      }
    }
    //calculate 2-electron integrals
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
          for(int l = 0; l < n; ++l) {
            rap = a2[i]*r/(a2[i]+a1[j]);
            rbp = r-rap;
            raq = a2[k]*r/(a2[k]+a1[l]);
            rbq = r-raq;
            rpq = rap-raq;
            rap2 = rap*rap;
            rbp2 = rbp*rbp;
            raq2 = raq*raq;
            rbq2 = rbq*rbq;
            rpq2 = rpq*rpq;
            v1111 = v1111 + twoe(a1[i],a1[j],a1[k],a1[l],0.0,0.0,0.0)*
                d1[i]*d1[j]*d1[k]*d1[l];
            v2111 = v2111 + twoe(a2[i],a1[j],a1[k],a1[l],r2,0.0,rap2)*
                d2[i]*d1[j]*d1[k]*d1[l];
            v2121 = v2121 + twoe(a2[i],a1[j],a2[k],a1[l],r2,r2,rpq2)*
                d2[i]*d1[j]*d2[k]*d1[l];
            v2211 = v2211 + twoe(a2[i],a2[j],a1[k],a1[l],0.0,0.0,r2)*
                d2[i]*d2[j]*d1[k]*d1[l];
            v2221 = v2221 + twoe(a2[i],a2[j],a2[k],a1[l],0.0,r2,rbq2)*
                d2[i]*d2[j]*d2[k]*d1[l];
            v2222 = v2222 + twoe(a2[i],a2[j],a2[k],a2[l],0.0,0.0,0.0)*
                d2[i]*d2[j]*d2[k]*d2[l];
#ifdef DEBUG
            std::printf("[%s] v2121:%f r2:%f twoe:%f a2[%d]:%f a1[%d]:%f d2[%d]:%f d1[%d]:%f d2[%d]:%f d1[%d]%f\n",
                        __FILE__, v2121, r2, twoe(a2[i],a1[j],a2[k],a1[l],r2,r2,rpq2),
                   i, a2[i], j, a1[j], i, d2[i], j, d1[j], k, d2[k],
                   l, d1[l]);
            std::printf("[%s] v2221:%f d2[%d]:%f d2[%d]:%f d2[%d]:%f d1[%d]%f\n",
                        __FILE__, v2221, i, d2[i], j, d2[j], k, d2[k], l, d1[l]);
            std::printf("[%s] rap:%.16E raq:%.16E rpq:%.16E\n",
                        __FILE__, rap, raq, rpq);
            std::printf("[%s] r2:%.16E rpq2:%20.12E\n",__FILE__, r2, rpq2);
            std::printf("[%s] [twoe]:%.16E\n",
                        __FILE__, twoe(a2[i],a1[j],a2[k],a1[l],r2,r2,rpq2));
            std::printf("[%s] d2[%d]:%.16E d1[%d]:%.16E\n",
                        __FILE__, i, d2[i], j, d1[j]);
            std::printf("[%s] d2[%d]:%.16E d1[%d]:%.16E\n",
                        __FILE__, k, d2[k], l, d1[l]);
            std::printf("[%s] a2[%d]:%.16E a1[%d]:%.16E\n",
                        __FILE__, i, a2[i], j, a1[j]);
            std::printf("[%s] a2[%d]:%.16E a1[%d]:%.16E\n",
                        __FILE__, k, a2[k], l, a1[l]);
            std::printf("[%s] v2121:%.16E\n",__FILE__, v2121);
#endif
          }
        }
      }
    }
    if( iop == 0 ) return;
    printf("   r          zeta1      zeta2      s12        t11\n");
    printf("   %f   %f   %f   %f   %f\n",
           r, zeta1, zeta2, s12, t11);
    printf("   t12        t22        v11a       v12a       v22a\n");
    printf("   %f   %f   %f   %f   %f\n",
           t12, t22, v11a, v12a, v22a);
    printf("   v11b       v12b        v22b      v1111      v2111\n");
    printf("   %f   %f   %f   %f   %f\n",
           v11b, v12b, v22b, v1111, v2111);
    printf("   v2121      v2211       v2221     v2222\n");
    printf("   %f   %f   %f   %f\n",
           v2121, v2211, v2221, v2222);
  }
}
