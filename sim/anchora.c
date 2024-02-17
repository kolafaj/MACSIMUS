/* #include code for Shake
   (called after SHAKE iterations have been finished)

1/ constrain center-of-mass of a molecule
   Note that SHAKE preserves center-of-mass so that this code need not 
   be inside the SHAKE iteration loop.

2/ printing site or center-of-mass constraint force except for a group

3/ constrain selected axis
*/

#if SHAKE!=1
#  error SHAKE=1 required for (this version of) ANCHOR
#endif /*# SHAKE!=1 */

/* anchor center-of-mass of whole molecule */
if (mn->anchor & ANCHOR_c) {
  vector cm;
  double m=spec[mn->sp]->mass,rr;

  VO(cm,=0)
  loop (i,0,ns) VV(cm,+=si[i].mass*p[i])
  VO(cm,/=m)
  VV(cm,-=mn->r0)
  normalizer(cm);
  /* NOTES: called only once (=> no +=)
     -m here because we want to have a constraint force */
  if (mn->xyz==7) {
    VV(cf,=(-m/(h*h))*cm)
    rr=SQR(cm); }
  else {
    /* for selected coordinates only 
       NB: other parts of the code is inefficiently calculated for whole vectors) */
    int k,b=1;
    
    rr=0;
    loop (k,0,3) {
      if (mn->xyz&b) {
        rr+=Sqr(cm[k]);
        cf[k]=(-m/(h*h))*cm[k]; }
      else
        cf[k]=cm[k]=0;
      b<<=1; } }

  if (rr>anchor.rr) {
    notinplace++;
    rr=sqrt(anchor.rr/rr);
    VO(cm,*=rr) }
  loop (i,0,ns) VV(p[i],-=cm)
}


/* print any constraint force (but groups) */
if (mn->anchor & ANCHOR_pos)
  printanchor(FORCEUNIT,cf,mn->xyz);


/* axes */
if (mn->anchor & ANCHOR_axis) {
  vector from,aux1,aux2,axis0,axis1,axis2,rot;
  double ca,sa;
  double m;

  /* determining the axis of rotation */
  switch (mn->anchor & ANCHOR_axis) {
    case ANCHOR_i: {
      static double **IM,**R;
      int j,k,l;

      if (!IM) {
        allocarray(IM,3);
        allocarray(R,3);
        loop (j,0,3) {
          allocarray(IM[j],3);
          allocarray(R[j],3); } }

      /* principal axis of the inertia tensor */
      VO(from,=0)
      loop (i,0,ns) VV(from,+=si[i].mass*p[i])
      m=spec[mn->sp]->mass;
      VO(from,/=m)

      loop (j,0,3) loop (k,0,3) IM[j][k]=0;
      loop (i,0,ns) {
        vector dr;
        double m=si[i].mass;

        VVV(dr,=p[i],-from)
        loop (j,0,3) loop (k,0,3) IM[j][k]+=m*dr[j]*dr[k]; }

      m=Jacobi(-3,IM,R,1e-13); /* error reached returned */
      if (m>1e-9) ERROR(("inertia tensor: diagonalization failed, err=%g",m))

      /* ordering */
      i=0,j=1,k=2;
    again:
      if (IM[i][i]>IM[j][j]) l=i,i=j,j=l;
      if (IM[j][j]>IM[k][k]) { l=k,k=j,j=l; goto again; }

      if (IM[i][i]/IM[j][j]>0.99)
        ERROR(("inertia tensor: ratio of 2 largest axes closer than 0.99"))

      VV(axis0,=R[i]) }

      break;

    case ANCHOR_a:
      /* (CM,site) */
      VO(from,=0)
      loop (i,0,ns) VV(from,+=si[i].mass*p[i])
      m=spec[mn->sp]->mass;
      VO(from,/=m)
      VVV(axis0,=p[mn->iaxis[0]],-from)
      break;

    case ANCHOR_p:
      /* (site,site) */
      VV(from,=p[mn->iaxis[0]])
      VVV(axis0,=p[mn->iaxis[1]],-from)
      break;

    case ANCHOR_t: {
      /* perpendicular to plane (site,site,site) */
      VV(from,=p[mn->iaxis[0]])
      VVV(aux1,=p[mn->iaxis[1]],-from)
      VVV(aux2,=p[mn->iaxis[2]],-from)
      VECT(axis0,aux1,aux2) }
      break;
    default: ERROR(("bad anchor")); }

  ca=sqrt(SQR(axis0)); VO(axis0,/=ca);
  ca=SCAL(axis0,mn->axis);

  VECT(rot,axis0,mn->axis)
  sa=sqrt(SQR(rot)); VVO(axis1,=rot,/sa);
  sa=SCAL(rot,axis1);

/*.....  prt_("%g %g  ",ca,sa);*/
  if (ca<anchor.cos) {
    notinplace++;
    sa=anchor.sin;
    ca=anchor.cos; }
/*.....  prt("  %g %g",ca,sa);*/

  VECT(axis2,axis0,axis1) /* axis0,axis1,axis2 = basis */

  /* rotate around axis1 so that axis0 becomes mn->axis
     and calculate the moment of inertia (around axis1) */
  m=0;
  loop (i,0,ns) {
    vector dr,b,bb;

    /* centering */
    VVV(dr,=p[i],-from)

    /* coordinates in the basis */
    b[0]=SCAL(dr,axis0);
    b[1]=SCAL(dr,axis1); /* parallel to rot */
    b[2]=SCAL(dr,axis2);

    m+=si[i].mass*(Sqr(b[0])+Sqr(b[2]));

    /* rotation in the basis */
    bb[0]= b[0]*ca+b[2]*sa;
    bb[1]= b[1];
    bb[2]=-b[0]*sa+b[2]*ca;

    /* back */
    dr[0]=axis0[0]*bb[0]+axis1[0]*bb[1]+axis2[0]*bb[2];
    dr[1]=axis0[1]*bb[0]+axis1[1]*bb[1]+axis2[1]*bb[2];
    dr[2]=axis0[2]*bb[0]+axis1[2]*bb[1]+axis2[2]*bb[2];

    VVV(p[i],=from,+dr) }

  /* calculate the moment of constraint forces
     (here, m=moment of inertia with respoct to axis1) */
  VV(aux1,=(m/(h*h))*rot)

    printanchor(TORQUEUNIT,aux1,mn->xyz); }

/* groups: print constraint forces (collected in anchorg.c) */
if (mn->anchor & ANCHOR_g) {
  struct anchorgroup_s *group;

  looplist (group,mn->group)
    printanchor(FORCEUNIT,group->cf,mn->xyz); }
