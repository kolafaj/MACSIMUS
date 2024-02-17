/* 
  #included by ground, or can be #included directly
  see also fortran.h 
*/

double fstrtod(const char *nptr, char **endptr) /******************* fstrtod */
/* the same as strtod but accepting also FORTRAN numbers as 1.2D+3 */
{
  const char *b=nptr;
  char *end=*endptr;
  static int pass;
  long int fexp;
  double x=strtod(b,&end);

  *endptr=(char *)nptr;
  if (b==end) return 0;
  
  if (*end && strchr("dD",*end)!=NULL && end[1] && strchr("+-0123456789",end[1])) {
    b=end+1;
    fexp=strtol(b,&end,10);
    if (b!=end) {
      if (!pass) {
        pass=1;
        fprintf(stderr,"WARNING: FORTRAN number with D-exponent parsed (more warnings suppressed)\n"); }

      *endptr=end; 
        
      return x*powi(10,fexp); } }

  *endptr=end;
  return x;
}

double fatof(const char *nptr) /************************************** fatof */
/* the same as atof but accepting also FORTRAN numbers as 1.2D+3 */
{
  char *end;

  return fstrtod(nptr,&end);
}
