/* 
   Summing up forces, their virial, the electrostatic potential including the
   elst field due to the dipole moment of the cell, Q->M, and the Yeh-Berkowitz
   slab corrections.
   #included twice by ewald.c with different macro ADDFORCES:
   1) Standard evaluation for neutral systems.
   2) With background version (el.bg) of the Yeh-Berkowitz correction: 
      doi:10.1063/1.3216473, eq. (30) 
      [V. Ballenegger, A. Arnold and J. J. Cerda, JCP 131, 094107 (2009)]
*/

      loop (n,0,No.N) {
        mn=molec+n;
        si=spec[mn->sp]->si;
        ns=mn->ns;
#if PRESSURETENSOR&PT_VIR
        r=rof(mn,rp);
#endif /*# PRESSURETENSOR&PT_VIR */
        fr=rof(mn,frp);
        loop (i,0,ns) {
#if defined(POLAR) && POLAR&32
          if (si[i].qtype&FQ) {
            aux=polarrof(mn,rp)[i][0]; /* fluctuating charge */
            ptr=fr[i];
            ADDFORCES(r[i])
            polptr=(real*)((char*)fr[i]+polar_off);
            /* M-dependent term missing (epsinf=infty only)? - to re-consider! */
            polptr[0]+=Phi[iq]/(2*PI)-2*alsPI*aux; /* elst. potential */
            iq++;
          } else
#endif /*# defined(POLAR) && POLAR&32 */
          if ((aux=si[i].charge)) {
            ptr=fr[i];
            ADDFORCES(r[i])
            iq++; }
#ifdef POLAR
          // should not test if si[i].qtype&FQ...
          if ((aux=si[i].chargepol)) {
            vector rboth;
            real *dr=polarrof(mn,rp)[i];

            VVV(rboth,=r[i],+dr)
            ptr=(real*)((char*)fr[i]+polar_off);
            ADDFORCES(rboth)
            iq++; }
#endif /*# POLAR */
        } }
      if (iq!=Nq) ERROR(("Ewald internal iq=%d Nq=%d",iq,Nq))
