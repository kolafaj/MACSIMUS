/*
  include file for simpot.c - intermolecular interactions
  atom-atom version
  cookmeas.c or cooknmea.c must be #included first
*/

#ifndef POLAR
ss=sstab[si1[0].st];
if ( (q=si1[0].charge) ) {
  if ( (qq=q*si2[0].charge) ) LJQQX(&ss[si2[0].st],
                                    r1[0],r2[0],f1[0],f2[0],qq GLOB);
  else LJX(&ss[si2[0].st], r1[0],r2[0],f1[0],f2[0] GLOB); }
else
  LJX(&ss[si2[0].st], r1[0],r2[0],f1[0],f2[0] GLOB);
#else /* POLAR */
ss=sstab[si1[0].st];
polarLJQQX(&ss[si2[0].st],
           r1[0],r2[0],f1[0],f2[0],
           si1,si2 GLOB);
#endif
