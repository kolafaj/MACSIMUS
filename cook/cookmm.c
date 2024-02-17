/*
  include file for cookpot.c - intermolecular interactions
  cookmeas.c or cooknmea.c must be #included first
*/

#ifndef POLAR

loop (i,0,ns1) {
  ss=sstab[si1[i].st];
  r1i=r1[i]; f1i=f1[i];
  if ( (q=si1[i].charge) )
    loop (j,0,ns2)
      if ( (qq=q*si2[j].charge) ) LJQQX(&ss[si2[j].st],
					r1i,r2[j],f1i,f2[j],qq GLOB);
      else LJX(&ss[si2[j].st], r1i,r2[j],f1i,f2[j] GLOB);
  else
    loop (j,0,ns2) LJX(&ss[si2[j].st], r1i,r2[j],f1i,f2[j] GLOB); }

#else /* POLAR */

loop (i,0,ns1) {
  ss=sstab[si1[i].st];
  r1i=r1[i]; f1i=f1[i];
  loop (j,0,ns2) polarLJQQX(&ss[si2[j].st],
			    r1i,r2[j],f1i,f2[j],
			    si1+i,si2+j GLOB); }

#endif
