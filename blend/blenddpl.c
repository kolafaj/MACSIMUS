/* debug/control/visualization dump of convergence */
    { 
      static int pass;
      static FILE *f=NULL;
      float h[3],d[3],z[3],rr;

      if (!f) {
	int i;

	z[0]=z[1]=0; z[2]=1;

	f=fopen("poldump.mol","wt");
	fprintf(f,"poldump\n\
\n\
number_of_atoms = %d\n\
\n\
atoms\n",4*ns);
	loop (i,0,ns) fprintf(f,"\
%d FROM C 0 0 1 %d\n\
%d TO   C 0 0 3 %d %d %d\n\
%d A1   C 0 0 1 %d\n\
%d A2   C 0 0 1 %d\n\
",4*i,4*i+1,
4*i+1,4*i,4*i+2,4*i+3,
4*i+2,4*i+1,
4*i+3,4*i+1);
	fclose(f);

	f=fopen("poldump.plb","wb");
	h[0]=ns*4;
	h[1]=0;
	fwrite(h,4,2,f); }

#define SCALE 8

      VV(h,=site[i].r)
      fwrite(h,4,3,f);
      VVV(h,=site[i].r,+SCALE*mu[i])
      fwrite(h,4,3,f);

      VECT(d,mu[i],z)
      rr=SQR(d);
      if (rr) rr=.05/sqrt(rr);
      VO(d,*=rr)

      rr=SQR(mu[i]);
      if (rr) rr=SCALE-.14/sqrt(rr);
      VVVV(h,=site[i].r,+rr*mu[i],+d)
      fwrite(h,4,3,f);
      VVVV(h,=site[i].r,+rr*mu[i],-d)
      fwrite(h,4,3,f);
	
      if (++pass==ns*10) {
	fclose(f); exit(0); } }
