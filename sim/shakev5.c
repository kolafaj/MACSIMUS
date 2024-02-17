/* kinetic energy, pressure tensor, and final a[0], a[1] update */
/* VERLET = 5 (EQUIPARTITION): 12v(t)^2 = 5V1^2 + 5V2^2 + 2V1V2
   V2 = [r(t+h)-r(t)]/h, V1 = [r(t)-r(t-h)]/h
   best equipartition for a vibrationg diatomic */

    VVO(v_mol1,=v_mol2,=0)
    m_mol=0;

    loop (i,0,ns) {
      vector V2; /* [r(t+h)-r(t)]/h */
      VVV(V2,=p[i],-r[i])

      mi=si[i].mass;
      m_mol+=mi;

      VV(v_mol1,+=mi*V1)
      VV(v_mol2,+=mi*V2)
      En.kin+=(5./12*(SQR(V1)+SQR(V2))+1./6*SCAL(V1,V2))*mi;
#if PRESSURETENSOR&PT_KIN
      En.Pkin[0]+=(5./12*(V1[0]*V1[0]+V2[0]*V2[0])+1./6*V1[0]*V2[0])*mi;
      En.Pkin[1]+=(5./12*(V1[1]*V1[1]+V2[1]*V2[1])+1./6*V1[1]*V2[1])*mi;
      En.Pkin[2]+=(5./12*(V1[2]*V1[2]+V2[2]*V2[2])+1./6*V1[2]*V2[2])*mi;
#  if PRESSURETENSOR&PT_OFF
      En.Pkin[3]+=(5./12*(V1[1]*V1[2]+V2[1]*V2[2])+1./6*V1[1]*V2[2])*mi;
      En.Pkin[4]+=(5./12*(V1[2]*V1[0]+V2[2]*V2[0])+1./6*V1[2]*V2[0])*mi;
      En.Pkin[5]+=(5./12*(V1[0]*V1[1]+V2[0]*V2[1])+1./6*V1[0]*V2[1])*mi;
#  endif /*# PRESSURETENSOR&PT_OFF */
#endif /*# PRESSURETENSOR&PT_KIN */

      VV(r1[i],=V2)
      VV(r[i],=p[i]) }

    En.kin_tr+=(5./12*(SQR(v_mol1)+SQR(v_mol2))+1./6*SCAL(v_mol1,v_mol2))/m_mol;

#if PRESSURETENSOR&PT_KIN && PRESSURETENSOR&PT_MOL

    En.PKin[0]+=(5./12*(v_mol1[0]*v_mol1[0]+v_mol2[0]*v_mol2[0])+1./6*v_mol1[0]*v_mol2[0])/m_mol;
    En.PKin[1]+=(5./12*(v_mol1[1]*v_mol1[1]+v_mol2[1]*v_mol2[1])+1./6*v_mol1[1]*v_mol2[1])/m_mol;
    En.PKin[2]+=(5./12*(v_mol1[2]*v_mol1[2]+v_mol2[2]*v_mol2[2])+1./6*v_mol1[2]*v_mol2[2])/m_mol;

#  if PRESSURETENSOR&PT_MOM
#    error not implemented
#  endif /*# PRESSURETENSOR&PT_MOM */

#  if PRESSURETENSOR&PT_OFF
    En.PKin[3]+=(5./12*(v_mol1[1]*v_mol1[2]+v_mol2[1]*v_mol2[2])+1./6*v_mol1[1]*v_mol2[2])*mi;
    En.PKin[4]+=(5./12*(v_mol1[2]*v_mol1[0]+v_mol2[2]*v_mol2[0])+1./6*v_mol1[2]*v_mol2[0])*mi;
    En.PKin[5]+=(5./12*(v_mol1[0]*v_mol1[1]+v_mol2[0]*v_mol2[1])+1./6*v_mol1[0]*v_mol2[1])*mi;
#  endif /*# PRESSURETENSOR&PT_OFF */

#endif /*# PRESSURETENSOR&PT_KIN && PRESSURETENSOR&PT_MOL */
