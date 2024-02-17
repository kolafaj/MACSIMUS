/* kinetic energy, pressure tensor, and final a[0], a[1] update */
/* VERLET=1 (standard Verlet): v(t) = [r(t+h)-r(t-h)]/2h */
/* equivalent to velocity Verlet */
    VO(v_mol1,=0)
    m_mol=0;

    loop (i,0,ns) {
      vector V2; /* [r(t+h)-r(t)]/h */
      vector V;  /* [r(t+h)-r(t-h)]/h (doubled velocity) */
      VVV(V2,=p[i],-r[i])

      mi=si[i].mass;
      m_mol+=mi;
      
      VVV(V,=V1,+V2)
      VV(v_mol1,+=mi*V)
      En.kin+=SQR(V)*mi; /* 8*kin.E: resolved later */
#if PRESSURETENSOR&PT_KIN
      En.Pkin[0]+=V[0]*V[0]*mi;
      En.Pkin[1]+=V[1]*V[1]*mi;
      En.Pkin[2]+=V[2]*V[2]*mi;
#  if PRESSURETENSOR&PT_OFF
      En.Pkin[3]+=V[1]*V[2]*mi;
      En.Pkin[4]+=V[2]*V[0]*mi;
      En.Pkin[5]+=V[0]*V[1]*mi;
#  endif /*# PRESSURETENSOR&PT_OFF */
#endif /*# PRESSURETENSOR&PT_KIN */

      VV(r1[i],=V2)
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
