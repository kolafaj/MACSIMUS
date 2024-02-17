/*
  include file for cookpot.c - intramolecular interactions
  cookmeas.c or cooknmea.c must be #included first
  standard version: 
    I0=0, I1=ns1, GLOB= (empty)
  PARALLEL==3 version: 
    I0=i0, I1=i1, GLOB=,glob (w. measure) or empty (not meas)
*/

#ifndef POLAR

loop (i,I0,I1) {
  ss=sstab[si1[i].st];
  r1i=r1[i]; f1i=f1[i];
  if ( (q=si1[i].charge) ) {

    j0=0; exc=si1[i].exc;
    do {
      j1=exc->indx;
      if ( (qq=q*si2[j1].charge) ) {
        if (exc->type==ONEFOUR) LJQQ14X(&sstab14[si1[i].st][si2[j1].st],
					r1i,r2[j1],f1i,f2[j1],qq GLOB);
        else if (exc->type==EXCL) XQQX(r1i,r2[j1],f1i,f2[j1],qq GLOB); }
      else
        if (exc->type==ONEFOUR) LJX(&sstab14[si1[i].st][si2[j1].st],
				    r1i,r2[j1],f1i,f2[j1] GLOB);

      loop (j,j0,j1) {
        if ( (qq=q*si2[j].charge) ) LJQQX(&ss[si2[j].st],
					  r1i,r2[j],f1i,f2[j],qq GLOB);
        else LJX(&ss[si2[j].st],
		 r1i,r2[j],f1i,f2[j] GLOB); }

      exc++;
      } while ( (j0=j1+1)<i );
    }

  else { /* zero charge */

    j0=0; exc=si1[i].exc;
    do {
      j1=exc->indx;
      if (exc->type==ONEFOUR) LJX(&sstab14[si1[i].st][si2[j1].st],
				  r1i,r2[j1],f1i,f2[j1] GLOB);

      loop (j,j0,j1) LJX(&ss[si2[j].st],
			 r1i,r2[j],f1i,f2[j] GLOB);

      exc++;
      } while ( (j0=j1+1)<i );
    }
  }

#else /* POLAR */

loop (i,0,ns1) {
  ss=sstab[si1[i].st];
  r1i=r1[i]; f1i=f1[i];

  j0=0; exc=si1[i].exc;
  do {
    j1=exc->indx;
    if (exc->type==ONEFOUR) polarLJQQ14X(&sstab14[si1[i].st][si2[j1].st],
					 r1i,r2[j1],f1i,f2[j1],
					 si1+i,si2+j1,1,factor14_1 GLOB);
#if 1
    else if (exc->type==EXCL) polarXQQX(&ss[si2[j1].st],
					   r1i,r2[j1],f1i,f2[j1],
					   si1+i,si2+j1 GLOB);
#else
    else if (exc->type==EXCL) polarLJQQ14X(&ss[si2[j1].st],
					   r1i,r2[j1],f1i,f2[j1],
					   si1+i,si2+j1,0,-1 GLOB);
#endif

    loop (j,j0,j1) polarLJQQX(&ss[si2[j].st],
			      r1i,r2[j],f1i,f2[j],
			      si1+i,si2+j GLOB);

    exc++;
    } while ( (j0=j1+1)<i );
  }

#endif /* POLAR */
