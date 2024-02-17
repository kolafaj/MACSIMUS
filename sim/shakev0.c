/* kinetic energy, pressure tensor, and final a[0], a[1] update */
/* VERLET = 0 (SHIFTED): v(t) = [r(t+h)-r(t)]/h */
/* development version - do not use for real simulations */
    VO(v_mol1,=0)
    m_mol=0;

    loop (i,0,ns) {
      mi=si[i].mass;
      m_mol+=mi;

      VV(v_mol1,+=mi*V1)
      En.kin+=SQR(V1)*mi; /* 2*kin.E*/
#if PRESSURETENSOR&PT_KIN
      En.Pkin[0]+=V1[0]*V1[0]*mi;
      En.Pkin[1]+=V1[1]*V1[1]*mi;
      En.Pkin[2]+=V1[2]*V1[2]*mi;
#  if PRESSURETENSOR&PT_OFF
      En.Pkin[3]+=V1[1]*V1[2]*mi;
      En.Pkin[4]+=V1[2]*V1[0]*mi;
      En.Pkin[5]+=V1[0]*V1[1]*mi;
#  endif /*# PRESSURETENSOR&PT_OFF */
#endif /*# PRESSURETENSOR&PT_KIN */

      VVV(r1[i],=p[i],-r[i])
      VV(r[i],=p[i]) }

    En.kin_tr+=SQR(v_mol1)/m_mol;

#if PRESSURETENSOR&PT_KIN && PRESSURETENSOR&PT_MOL
    En.PKin[0]+=v_mol1[0]*v_mol1[0]/m_mol;
    En.PKin[1]+=v_mol1[1]*v_mol1[1]/m_mol;
    En.PKin[2]+=v_mol1[2]*v_mol1[2]/m_mol;

#  if PRESSURETENSOR&PT_MOM
    En.PKin2[0]+=Sqr(v_mol1[0]*v_mol1[0]/m_mol);
    En.PKin2[1]+=Sqr(v_mol1[1]*v_mol1[1]/m_mol);
    En.PKin2[2]+=Sqr(v_mol1[2]*v_mol1[2]/m_mol);
    En.PKin3[0]+=Cub(v_mol1[0]*v_mol1[0]/m_mol);
    En.PKin3[1]+=Cub(v_mol1[1]*v_mol1[1]/m_mol);
    En.PKin3[2]+=Cub(v_mol1[2]*v_mol1[2]/m_mol);
#  endif /*# PRESSURETENSOR&PT_MOM */

#  if PRESSURETENSOR&PT_OFF
    En.PKin[3]+=v_mol1[1]*v_mol1[2]/m_mol;
    En.PKin[4]+=v_mol1[2]*v_mol1[0]/m_mol;
    En.PKin[5]+=v_mol1[0]*v_mol1[1]/m_mol;
#  endif /*# PRESSURETENSOR&PT_OFF */

#endif /*# PRESSURETENSOR&PT_KIN && PRESSURETENSOR&PT_MOL */
