/*
  returns the dihedral angle (including the correct sign)
  for molecule      r1--r2
  (or part of       /    \
  molecule)       r0      r3

  #included by:
    blend/blendmin.c
    blend/ramachan.c
    show/show.c 

   see also dihedralpot() and improperpot() (sign!) in sim/intrapot.c
*/

double dihedral(fvector r0,fvector r1,fvector r2,fvector r3)
{
  vector v01,v12,v13,v23, a,b,ap,cp;
  double s,q,sinphi,aa,bb,norm,ab,df,abbb,cf,bcbb,apq,cpq,bc;

  VVV(v01,=r1,-r0)
  VVV(v12,=r2,-r1)
  VVV(v13,=r3,-r1)

  VECT(a,v01,v12) /* vector a is perpendicular to 012 plane */

  q=SQR(v12);
  s=SCAL(v12,v13)/q;
  VVV(b,=v13,-s*v12) /* vector b lies in 123 plane an is perpendicular to 12 */

  norm=sqrt( (aa=SQR(a))*(bb=SQR(b)) );
  ab=SCAL(a,b);
  sinphi=ab/norm;

  /* to fix improbable numerical things like |sinphi|=1+1e-16 */
  if (sinphi<-1) sinphi=-1;
  if (sinphi>1) sinphi=1;
  df=asin(sinphi);

  /* ? */
  VECT(ap,v12,v13)
  if (SCAL(ap,a)<0) {
    /* Error("|improper torsion angle|>PI/2"); */
    /* fix for |phi|>PI/2 (should not normally occur!) */
    if (df>0) df=PI-df;
    else df=-PI-df;
    norm=-norm; }

  /* fool-proof check from dihedralpot() */
  VVV(v01,=r1,-r0) VVV(v12,=r2,-r1) VVV(v23,=r3,-r2)
  ab=SCAL(v01,v12); bc=SCAL(v12,v23); bb=SQR(v12);
  abbb=ab/bb; VVV(ap,= -v01,+abbb*v12)
  bcbb=bc/bb; VVV(cp,=  v23,-bcbb*v12)
  /* vectors ap and cp are perpendicular to v12 */

  norm=1/sqrt( (apq=SQR(ap)) * (cpq=SQR(cp)) );

  /* cos phi */
  cf=SCAL(ap,cp)*norm;

  if (fabs(fabs(df)-acos(cf))>1e-6)
    fprintf(stderr,"imprecise dihedral (almost singular?): %8.4f %8.4f deg (dif=%.2g rad)\n",
            180/PI*fabs(df),180/PI*acos(cf),fabs(df)-acos(cf));

  df*=180/PI;
  if (fullangle) df=fmod(df+360,360);

  return df;
}

