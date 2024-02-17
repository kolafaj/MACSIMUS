/*
   Ewald pass 1 (calculate Fourier transform of charges)
   #included from ewald.c
*/
  loopto (kx,0,kxmax) {
    if (kx) iqeloop(ex) CxE(qtab[iqt].x)
    iqloop qtab[iqt].xy = qtab[iqt].xY = qtab[iqt].x;

    /* kymax=(int)(L[1]*sqrt(Kk-Sqr(kx/L[0]))); */
    kymax=ktab[kx].kymax;

    loopto (ky,0,kymax) {
      if (ky) iqeloop(ey) { CxE(qtab[iqt].xy) CxiE(qtab[iqt].xY) }

      iqloop {
        qtab[iqt].xyz = qtab[iqt].xyZ = qtab[iqt].xy;
        qtab[iqt].xYz = qtab[iqt].xYZ = qtab[iqt].xY; }

      /* kzmax=(int)(L[2]*sqrt(Kk-Sqr(ky/L[1])-Sqr(kx/L[0]))); */
      kzmax=ktab[kx].kzmax[ky];
      loopto (kz,!(kx+ky),kzmax) {
        if (kz)	iqeloop(ez) {
  	  CxE(qtab[iqt].xyz) CxiE(qtab[iqt].xyZ)
          CxE(qtab[iqt].xYz) CxiE(qtab[iqt].xYZ) }

#ifdef GAUSSIANCHARGES
        A.re=B.re=C.re=D.re=A.im=B.im=C.im=D.im=0;
        AT.re=BT.re=CT.re=DT.re=AT.im=BT.im=CT.im=DT.im=0;
        AS.re=BS.re=CS.re=DS.re=AS.im=BS.im=CS.im=DS.im=0;

        switch ((kx>0) + (ky>0) + (kz>0)) {

/* +++ */ case 3:
            iqloop {
              Cadd(A,qtab[iqt].xyz) /* +++ */
              Cadd(AT,sekk[ikk][qtab[iqt].site]*qtab[iqt].xyz)
              Cadd(AS,ssekk[ikk][qtab[iqt].site]*qtab[iqt].xyz)
              Cadd(B,qtab[iqt].xyZ) /* ++- */
              Cadd(BT,sekk[ikk][qtab[iqt].site]*qtab[iqt].xyZ)
              Cadd(BS,ssekk[ikk][qtab[iqt].site]*qtab[iqt].xyZ)
              Cadd(C,qtab[iqt].xYz) /* +-+ */
              Cadd(CT,sekk[ikk][qtab[iqt].site]*qtab[iqt].xYz)
              Cadd(CS,ssekk[ikk][qtab[iqt].site]*qtab[iqt].xYz)
              Cadd(D,qtab[iqt].xYZ) /* +-- */
              Cadd(DT,sekk[ikk][qtab[iqt].site]*qtab[iqt].xYZ)
              Cadd(DS,ssekk[ikk][qtab[iqt].site]*qtab[iqt].xYZ)
            }
            Qk[iQ].Q=A; Qk[iQ].T=AT; Qk[iQ].S=AS; iQ++;
            Qk[iQ].Q=B; Qk[iQ].T=BT; Qk[iQ].S=BS; iQ++;
            Qk[iQ].Q=C; Qk[iQ].T=CT; Qk[iQ].S=CS; iQ++;
            Qk[iQ].Q=D; Qk[iQ].T=DT; Qk[iQ].S=DS; iQ++;
          break;

/* 0++ */
/* +0+ */
/* ++0 */ case 2:
            if (kz==0) iqloop {
              Cadd(A,qtab[iqt].xyz) /* ++0 */
              Cadd(AT,sekk[ikk][qtab[iqt].site]*qtab[iqt].xyz)
              Cadd(AS,ssekk[ikk][qtab[iqt].site]*qtab[iqt].xyz)
              Cadd(B,qtab[iqt].xYz) /* +-0 */
              Cadd(BT,sekk[ikk][qtab[iqt].site]*qtab[iqt].xYz)
              Cadd(BS,ssekk[ikk][qtab[iqt].site]*qtab[iqt].xYz)
            }
            else iqloop {
              Cadd(A,qtab[iqt].xyz) /* ..+ */
              Cadd(AT,sekk[ikk][qtab[iqt].site]*qtab[iqt].xyz)
              Cadd(AS,ssekk[ikk][qtab[iqt].site]*qtab[iqt].xyz)
              Cadd(B,qtab[iqt].xyZ) /* ..- */
              Cadd(BT,sekk[ikk][qtab[iqt].site]*qtab[iqt].xyZ)
              Cadd(BS,ssekk[ikk][qtab[iqt].site]*qtab[iqt].xyZ)
            }
            Qk[iQ].Q=A; Qk[iQ].T=AT; Qk[iQ].S=AS; iQ++;
            Qk[iQ].Q=B; Qk[iQ].T=BT; Qk[iQ].S=BS; iQ++;
          break;

/* +00 */
/* 0+0 */
/* 00+ */ case 1:
            iqloop {
              Cadd(A,qtab[iqt].xyz)
              Cadd(AT,sekk[ikk][qtab[iqt].site]*qtab[iqt].xyz)
              Cadd(AS,ssekk[ikk][qtab[iqt].site]*qtab[iqt].xyz)
            }
          Qk[iQ].Q=A; Qk[iQ].T=AT; Qk[iQ].S=AS; iQ++;
          break;

/* ??? */ default: ERROR(("Ewald switch (k=%d %d %d)",kx,ky,kz)) }

#else

        A.re=B.re=C.re=D.re=A.im=B.im=C.im=D.im=0;

        switch ((kx>0) + (ky>0) + (kz>0)) {

/* +++ */ case 3:
            iqloop {
              Cadd(A,qtab[iqt].xyz) /* +++ */
              Cadd(B,qtab[iqt].xyZ) /* ++- */
              Cadd(C,qtab[iqt].xYz) /* +-+ */
              Cadd(D,qtab[iqt].xYZ) /* +-- */ }
            Qk[iQ++]=A; Qk[iQ++]=B; Qk[iQ++]=C; Qk[iQ++]=D;
          break;

/* 0++ */
/* +0+ */
/* ++0 */ case 2:
            if (kz==0) iqloop {
              Cadd(A,qtab[iqt].xyz) /* ++0 */
              Cadd(B,qtab[iqt].xYz) /* +-0 */ }
            else iqloop {
              Cadd(A,qtab[iqt].xyz) /* ..+ */
              Cadd(B,qtab[iqt].xyZ) /* ..- */ }
            Qk[iQ++]=A; Qk[iQ++]=B;
          break;

/* +00 */
/* 0+0 */
/* 00+ */ case 1:
            iqloop {
              Cadd(A,qtab[iqt].xyz) }
            Qk[iQ++]=A;
          break;

/* ??? */ default: ERROR(("Ewald switch (k=%d %d %d)",kx,ky,kz)) }
#endif
      ikk++;
      } } } /* kz,ky,kx */
