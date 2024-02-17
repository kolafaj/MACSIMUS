#define VELOCITY 
/* 
selects the rhs with velocities and may emit optimized code
required by MACSIMUS
cf. Lagrange and Hamilton in old MACSIMUS (or peo*.c project)
! not checked without VELOCITY recently !
*/

/* MAXGEARORDER is #defined in simglob.h */

extern struct gear_s {
  double C[MAXGEARORDER]; /* corrector coefficients to be changed for init>0 
                             -9e9 denotes coefficients to remain the same */
  int init;               /* version: 1=Gear no velocity, 2=velocity, 3=Janek */
  int order;              /* option('m') copied; former GearOrder */
  int lastorder;          /* to check (order change not allowed) */
} gear;

#ifndef POLAR
void Gear(int neq,ToIntPtr *a);
#else
void Gear2pol(int neq,int nscf,int M,int K,ToIntPtr *a);
#endif
