/* kinetic energy, pressure tensor, and final a[0], a[1] update */
/* VERLET = 3 (LEAPFROG): 2v(t)^2 = ([r(t+h)-r(t)]/h)^2 + ([r(t)-r(t-h)]/h)^2 */
    VVO(v_mol1,=v_mol2,=0)
    m_mol=0;

    loop (i,0,ns) {
      vector V2; /* [r(t+h)-r(t)]/h */
      VVV(V2,=p[i],-r[i])

      /* DEVELOPMENT NOTE:
         VV(V2,/=scale) and VV(V2,*=scale)
         NPT: fluctuations in total energy increased, drift not affected
         tau.rho<0: wrong! */

      mi=si[i].mass;
      m_mol+=mi;

      VV(v_mol1,+=mi*V1)
      VV(v_mol2,+=mi*V2)
      En.kin+=(SQR(V1)+SQR(V2))*mi;
#if PRESSURETENSOR&PT_KIN
      En.Pkin[0]+=(V1[0]*V1[0]+V2[0]*V2[0])*mi;
      En.Pkin[1]+=(V1[1]*V1[1]+V2[1]*V2[1])*mi;
      En.Pkin[2]+=(V1[2]*V1[2]+V2[2]*V2[2])*mi;
#  if PRESSURETENSOR&PT_OFF
      En.Pkin[3]+=(V1[1]*V1[2]+V2[1]*V2[2])*mi;
      En.Pkin[4]+=(V1[2]*V1[0]+V2[2]*V2[0])*mi;
      En.Pkin[5]+=(V1[0]*V1[1]+V2[0]*V2[1])*mi;
#  endif /*# PRESSURETENSOR&PT_OFF */
#endif /*# PRESSURETENSOR&PT_KIN */

      VV(r1[i],=V2)
      VV(r[i],=p[i]) }

    En.kin_tr+=(SQR(v_mol1)+SQR(v_mol2))/m_mol;

#if PRESSURETENSOR&PT_KIN && PRESSURETENSOR&PT_MOL

    En.PKin[0]+=(v_mol1[0]*v_mol1[0]+v_mol2[0]*v_mol2[0])/m_mol;
    En.PKin[1]+=(v_mol1[1]*v_mol1[1]+v_mol2[1]*v_mol2[1])/m_mol;
    En.PKin[2]+=(v_mol1[2]*v_mol1[2]+v_mol2[2]*v_mol2[2])/m_mol;

#if PRESSURETENSOR&PT_MOM
    En.PKin2[0]+=Sqr((v_mol1[0]*v_mol1[0]+v_mol2[0]*v_mol2[0])/m_mol);
    En.PKin2[1]+=Sqr((v_mol1[1]*v_mol1[1]+v_mol2[1]*v_mol2[1])/m_mol);
    En.PKin2[2]+=Sqr((v_mol1[2]*v_mol1[2]+v_mol2[2]*v_mol2[2])/m_mol);
    En.PKin3[0]+=Cub((v_mol1[0]*v_mol1[0]+v_mol2[0]*v_mol2[0])/m_mol);
    En.PKin3[1]+=Cub((v_mol1[1]*v_mol1[1]+v_mol2[1]*v_mol2[1])/m_mol);
    En.PKin3[2]+=Cub((v_mol1[2]*v_mol1[2]+v_mol2[2]*v_mol2[2])/m_mol);
#endif

#  if PRESSURETENSOR&PT_OFF
    En.PKin[3]+=(v_mol1[1]*v_mol1[2]+v_mol2[1]*v_mol2[2])/m_mol;
    En.PKin[4]+=(v_mol1[2]*v_mol1[0]+v_mol2[2]*v_mol2[0])/m_mol;
    En.PKin[5]+=(v_mol1[0]*v_mol1[1]+v_mol2[0]*v_mol2[1])/m_mol;
#  endif /*# PRESSURETENSOR&PT_OFF */

#endif /*# PRESSURETENSOR&PT_KIN && PRESSURETENSOR&PT_MOL */
