/* the scaling potential is:

      E              E
  ----------  +  ----------     ( bounding in <xi0,xi1> )
  (xi-xi0)^n     (xi1-xi)^n

  +  A xi  +  B xi^2            ( shape )

    xi (1-xi) K0 K1 (r1-r0)^2
  + -------------------------   ( fixup for changing bond )
        xi K0 + (1-xi) K1
*/

extern struct scale_s {
/* user data */
  real mass;
  real E,A,B;
  real xi0,xi1;
  real K0,K1,r0,r1;
  real xi; /* if set then a[0]->xi set (for minim.) */
  int n;
/* program data */
  vector *f1,*f2,*v,*r;
  } scale;
  
void scaleinit(void);
void scalebegin(ToIntPtr A,ToIntPtr V);
double scaleforce(ToIntPtr B,ToIntPtr A);

pot2_t scale2;
void printscaledpot(void);

