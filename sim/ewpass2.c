/* 
   Ewald pass 2 (calculate forces - no measurement version)
   #included from ewald.c 
   see ewpass2m.c for more info
*/
  loopto (kx,0,kxmax) {
    /* NOTE: 
       calculations of qtab are unnecessarily repeated in the 2nd pass,
       however, storing all the factors would require a lot of memory,
       namely [ikk] indexed array of qtab's */
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
        if (kz) iqeloop(ez) {
          CxE(qtab[iqt].xyz) CxiE(qtab[iqt].xyZ)
          CxE(qtab[iqt].xYz) CxiE(qtab[iqt].xYZ) }

        e=ekk[ikk]; 

        ek[0]=e*kx/box.L[0]; ek[1]=e*ky/box.L[1]; ek[2]=e*kz/box.L[2];

        switch ((kx>0) + (ky>0) + (kz>0)) {

#ifdef GAUSSIANCHARGES
        /********************************************/
        /**  pass 2 - Ewald with Gaussian charges  **/
        /********************************************/
/* +++ */ case 3:
            A=Qk[iQ].Q; AT=Qk[iQ].T; AS=Qk[iQ].S; iQ++;
            B=Qk[iQ].Q; BT=Qk[iQ].T; BS=Qk[iQ].S; iQ++;
            C=Qk[iQ].Q; CT=Qk[iQ].T; CS=Qk[iQ].S; iQ++;
            D=Qk[iQ].Q; DT=Qk[iQ].T; DS=Qk[iQ].S; iQ++;
            iqfloop {
              FA = AT.re*qtab[iqt].xyz.im - AT.im*qtab[iqt].xyz.re /* +++ */
                + (A.re*qtab[iqt].xyz.im -  A.im*qtab[iqt].xyz.re)*sekk[ikk][qtab[iqt].site];
              FB = BT.re*qtab[iqt].xyZ.im - BT.im*qtab[iqt].xyZ.re /* ++- */
                +  (B.re*qtab[iqt].xyZ.im -  B.im*qtab[iqt].xyZ.re)*sekk[ikk][qtab[iqt].site];
              FC = CT.re*qtab[iqt].xYz.im - CT.im*qtab[iqt].xYz.re /* +-+ */
                +  (C.re*qtab[iqt].xYz.im -  C.im*qtab[iqt].xYz.re)*sekk[ikk][qtab[iqt].site];
              FD = DT.re*qtab[iqt].xYZ.im - DT.im*qtab[iqt].xYZ.re /* +-- */
                +  (D.re*qtab[iqt].xYZ.im -  D.im*qtab[iqt].xYZ.re)*sekk[ikk][qtab[iqt].site];
#  if defined(POLAR) && POLAR&32
#    error "POLAR&32 + GAUSSIANCHARGES not implemented"
              Phi[iqt] += e*(A.re*qtab[iqt].xyz.re + A.im*qtab[iqt].xyz.im
                             + B.re*qtab[iqt].xyZ.re + B.im*qtab[iqt].xyZ.im
                             + C.re*qtab[iqt].xYz.re + C.im*qtab[iqt].xYz.im
                             + D.re*qtab[iqt].xYZ.re + D.im*qtab[iqt].xYZ.im);
#  endif /*# defined(POLAR) && POLAR&32 */
              F[iqt][0] += QTIMES ek[0]*((FA+FB)+(FC+FD));
              F[iqt][1] += QTIMES ek[1]*((FA+FB)-(FC+FD));
              F[iqt][2] += QTIMES ek[2]*(FA-FB+FC-FD); }
            break;
/* 0++ */
/* +0+ */
/* ++0 */ case 2:
            A=Qk[iQ].Q; AT=Qk[iQ].T; AS=Qk[iQ].S; iQ++;
            B=Qk[iQ].Q; BT=Qk[iQ].T; BS=Qk[iQ].S; iQ++;
/* ..+ */   if (kz) {
              iqfloop {
                FA = AT.re*qtab[iqt].xyz.im - AT.im*qtab[iqt].xyz.re /* ..+ */
                  + (A.re*qtab[iqt].xyz.im  - A.im*qtab[iqt].xyz.re)*sekk[ikk][qtab[iqt].site];
                FB = BT.re*qtab[iqt].xyZ.im - BT.im*qtab[iqt].xyZ.re /* ..- */
                  + (B.re*qtab[iqt].xyZ.im  - B.im*qtab[iqt].xyZ.re)*sekk[ikk][qtab[iqt].site];
#  if defined(POLAR) && POLAR&32
                Phi[iqt] += e*(A.re*qtab[iqt].xyz.re + A.im*qtab[iqt].xyz.im
                               + B.re*qtab[iqt].xyZ.re + B.im*qtab[iqt].xyZ.im);
#  endif /*# defined(POLAR) && POLAR&32 */
                if (kx) F[iqt][0] += QTIMES ek[0]*(FA+FB);
                else    F[iqt][1] += QTIMES ek[1]*(FA+FB);
                F[iqt][2] += QTIMES ek[2]*(FA-FB); } }
/* ++0 */   else {
              iqfloop {
                FA = AT.re*qtab[iqt].xyz.im - AT.im*qtab[iqt].xyz.re /* .+0 */
                  + (A.re*qtab[iqt].xyz.im -  A.im*qtab[iqt].xyz.re)*sekk[ikk][qtab[iqt].site];
                FB = BT.re*qtab[iqt].xYz.im - BT.im*qtab[iqt].xYz.re /* .-0 */
                  + ( B.re*qtab[iqt].xYz.im -  B.im*qtab[iqt].xYz.re)*sekk[ikk][qtab[iqt].site];
#  if defined(POLAR) && POLAR&32
                Phi[iqt] += e*(A.re*qtab[iqt].xyz.re + A.im*qtab[iqt].xyz.im
                               + B.re*qtab[iqt].xYz.re + B.im*qtab[iqt].xYz.im);
#  endif /*# defined(POLAR) && POLAR&32 */
                F[iqt][0] += QTIMES ek[0]*(FA+FB);
                F[iqt][1] += QTIMES ek[1]*(FA-FB); } }
            break;

/* +00 */
/* 0+0 */
/* 00+ */ case 1:
            A=Qk[iQ].Q; AT=Qk[iQ].T; AS=Qk[iQ].S; iQ++;
            iqfloop {
              FA = AT.re*qtab[iqt].xyz.im - AT.im*qtab[iqt].xyz.re
                + ( A.re*qtab[iqt].xyz.im - A.im*qtab[iqt].xyz.re)*sekk[ikk][qtab[iqt].site];
#  if defined(POLAR) && POLAR&32
              Phi[iqt] += e*(A.re*qtab[iqt].xyz.re + A.im*qtab[iqt].xyz.im);
#  endif /*# defined(POLAR) && POLAR&32 */
              if (kx)      F[iqt][0] += QTIMES ek[0]*FA;
              else if (ky) F[iqt][1] += QTIMES ek[1]*FA;
              else         F[iqt][2] += QTIMES ek[2]*FA; }
            break;

#else /*# GAUSSIANCHARGES */
          /********************************/
          /**    pass 2 - normal Ewald   **/
          /********************************/

/* +++ */ case 3:
            A=Qk[iQ++]; B=Qk[iQ++]; C=Qk[iQ++]; D=Qk[iQ++];
            iqfloop {
              FA = A.re*qtab[iqt].xyz.im - A.im*qtab[iqt].xyz.re; /* +++ */
              FB = B.re*qtab[iqt].xyZ.im - B.im*qtab[iqt].xyZ.re; /* ++- */
              FC = C.re*qtab[iqt].xYz.im - C.im*qtab[iqt].xYz.re; /* +-+ */
              FD = D.re*qtab[iqt].xYZ.im - D.im*qtab[iqt].xYZ.re; /* +-- */
#  if defined(POLAR) && POLAR&32
              Phi[iqt] += e*(A.re*qtab[iqt].xyz.re + A.im*qtab[iqt].xyz.im
                           + B.re*qtab[iqt].xyZ.re + B.im*qtab[iqt].xyZ.im
                           + C.re*qtab[iqt].xYz.re + C.im*qtab[iqt].xYz.im
                           + D.re*qtab[iqt].xYZ.re + D.im*qtab[iqt].xYZ.im);
#  endif /*# defined(POLAR) && POLAR&32 */
              F[iqt][0] += QTIMES ek[0]*((FA+FB)+(FC+FD));
              F[iqt][1] += QTIMES ek[1]*((FA+FB)-(FC+FD));
              F[iqt][2] += QTIMES ek[2]*(FA-FB+FC-FD); }
              break;

/* 0++ */
/* +0+ */
/* ++0 */ case 2:
            A=Qk[iQ++]; B=Qk[iQ++];
/* ..+ */   if (kz) {
              iqfloop {  
                FA = A.re*qtab[iqt].xyz.im - A.im*qtab[iqt].xyz.re; /* ..+ */
                FB = B.re*qtab[iqt].xyZ.im - B.im*qtab[iqt].xyZ.re; /* ..- */
#  if defined(POLAR) && POLAR&32
                Phi[iqt] += e*(A.re*qtab[iqt].xyz.re + A.im*qtab[iqt].xyz.im
                             + B.re*qtab[iqt].xyZ.re + B.im*qtab[iqt].xyZ.im);
#  endif /*# defined(POLAR) && POLAR&32 */
                if (kx) F[iqt][0] += QTIMES ek[0]*(FA+FB);
                else    F[iqt][1] += QTIMES ek[1]*(FA+FB);
                        F[iqt][2] += QTIMES ek[2]*(FA-FB); } }
/* ++0 */   else {
              iqfloop {
                FA = A.re*qtab[iqt].xyz.im - A.im*qtab[iqt].xyz.re; /* .+0 */
                FB = B.re*qtab[iqt].xYz.im - B.im*qtab[iqt].xYz.re; /* .-0 */
#  if defined(POLAR) && POLAR&32
                Phi[iqt] += e*(A.re*qtab[iqt].xyz.re + A.im*qtab[iqt].xyz.im
                             + B.re*qtab[iqt].xYz.re + B.im*qtab[iqt].xYz.im);
#  endif /*# defined(POLAR) && POLAR&32 */
                F[iqt][0] += QTIMES ek[0]*(FA+FB);
                F[iqt][1] += QTIMES ek[1]*(FA-FB); } }
            break;

/* +00 */
/* 0+0 */
/* 00+ */ case 1:
            A=Qk[iQ++];
            iqfloop {
              FA = A.re*qtab[iqt].xyz.im - A.im*qtab[iqt].xyz.re;
#  if defined(POLAR) && POLAR&32
              Phi[iqt] += e*(A.re*qtab[iqt].xyz.re + A.im*qtab[iqt].xyz.im);
#  endif /*# defined(POLAR) && POLAR&32 */
              if (kx)      F[iqt][0] += QTIMES ek[0]*FA;
              else if (ky) F[iqt][1] += QTIMES ek[1]*FA;
              else         F[iqt][2] += QTIMES ek[2]*FA; }
            break;

#endif /*#!GAUSSIANCHARGES */
/* ??? */   default: ERROR(("Ewald switch (k=%d %d %d)",kx,ky,kz)) }
        ikk++; } } } /* kz,ky,kx */

  if (iQ!=Q->N) ERROR(("Ewald: iQ=%d != Q->N=%d; L=[%g %g %g]",
                               iQ,      Q->N,    VARG(box.L)))
