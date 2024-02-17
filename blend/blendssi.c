/* all site-site interactions */
/* the list of exceptions must be ordered! */
loop (i,0,ns) {

#if 0 /* debug print */
  if (site[i].exc)
    prt("r[%d]=%6.3f %6.3f %6.3f %d/%d",
	i,r[i][0],r[i][1],r[i][2],site[i].exc->indx,site[i].exc->type);
  else
    prt("r[%d]=%6.3f %6.3f %6.3f", i,r[i][0],r[i][1],r[i][2]);
#endif

  PREPARE
  for (j0=0, exc=site[i].exc; ; exc=exc->next,j0=j1+1) {
    if (exc) {
      j1=exc->indx;
      if (exc->type==ONEFOUR) {
        IF_IN_CUTOFF_CUBE(j1)
	  SSPOT14; }
#ifdef POLAR
      else {
	if (polar&4) 
	  IF_IN_CUTOFF_CUBE(j1) /* test is void (turn on in future?) */
	    SSPOT12; }  
#endif
      }
    else
      j1=i;

    loop (j,j0,j1)
      IF_IN_CUTOFF_CUBE(j) SSPOT;

    if (!exc) break; }
  }
