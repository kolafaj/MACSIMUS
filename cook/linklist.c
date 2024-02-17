/* C2c linklist.C

  !!! DO NOT EDIT linklist.c !!!
      * edit linklist.C and run C2c linklist.C
      * optionally, run cppnest linklist.c
*/

#include "ground.h"
#include "sds.h"
#include "simglob.h"
#include "norm.h"
#include "linklist.h"

/* returns double number by 1 lsb smaller that given positive double number */
double minuseps(double x) /**************************************** minuseps */
{
  double y,z;
  int i;

  if (x<=0) ERROR(("not positive"))
  y=x*0.999999999999999667; /* -3 of the lsb */

  loop (i,0,5) {
    z=casttodouble((x+y)/2);
    if (z>=x) return y;
    y=z; }

  return y;
}

/* ultimate error limit for box.r0 */
#define R0LIMIT 3
#define OUTOFBOX 1e-5

linklist_t *site_storage=NULL;

#if PARALLEL==1
	linklist_t ****linklist(vector *f1p,vector *f2p,vector *rp)
#else
	linklist_t ****linklist(vector *fp,vector *rp)
#endif
/*
 * A link cell list {No.cell[0] * No.cell[1] * No.cell[2]} of the whole
   configuration is built according to the x,y,z coordinates of sites.
 * Memory needed is allocated in the 1st pass and re-used later.
 * x,y,z are normalized to [0,L),
 * type (linklist_t****) is interpreted as
   array[No.cell[0]][No.cell[1]][No.cell[2]] of type (linklist_t*)
 */
{
  int i,sp,ix,iy,iz,n,ns;
  static int oldncell[3],oldns=0;
  static linklist_t ****list=NULL;
  linklist_t *l;
  vector rL;
#if PARALLEL==1
	vector *r,*f1,*f2;
#else
	vector *r,*f;
#endif
#ifdef POLAR
  vector *rpol;
#  if PARALLEL==1
	vector *f1pol,*f2pol;
#  else
	vector *fpol;
#  endif
#endif /*# POLAR */
  molecule_t *mn;

  ix=0;
  loop (i,0,DIM) {
    if (No.cell[i]<2) ERROR(("invalid No.cell[%d]=%d",i,No.cell[i]))
    if (No.cell[i]!=oldncell[i]) ix++; }

  /*
     Trick so that int(L[k]*rL) never gives No.cell[k]; one test is thus
     spared.  Because of possible optimization in the numerical unit, rL is
     cast to double to avoid keeping in 10-byte registers.
  */
  loop (i,0,3) {
    double x;

    rL[i]=casttodouble(No.cell[i]/box.L[i]);
    x=minuseps(box.L[i]);
    while ((int)(x*rL[i])>=No.cell[i]) rL[i]=minuseps(rL[i]); }

  if (!site_storage) {
    int s=sizeof(linklist_t);

    prt("sizeof(struct linklist_s)=%d, CACHELINE=%d (#defined in gen/ground.h)",s,CACHELINE);
    if (s%CACHELINE)
      prt("\
NOTE: The linked list entries do not fit entire cache lines.\n\
Depending on the architecture and size, speed may be gained by padding.\n\
Check CACHELINE and/or add the following line to simopt.h:\n\
#define LINKLIST_PADDING %d\n\
see also cook/linklist.h (cook/linklist.H)\n",(CACHELINE-s%CACHELINE)/sizeof(double));

    /* all memory for sites allocated in advance */
    allocarray(site_storage,No.s);
    oldns=No.s; }
  else
    if (oldns!=No.s)
      ERROR(("linklist: No.s=%d has changed (old=%d)",No.s,oldns))

  if (ix && list) {
    prt("number of cells changed => link-cell list reallocated");
    loop (ix,0,oldncell[0]) {
      loop (iy,0,oldncell[1]) free(list[ix][iy]);
      free(list[ix]); }
    free(list); /* note that my free clears the pointer to NULL ! */ }

  if (!list) {
    /* 3D array of cells allocated */
    allocarray(list,No.cell[0]);
    loop (ix,0,No.cell[0]) {
      allocarray(list[ix],No.cell[1]);
      loop (iy,0,No.cell[1])
	allocarray(list[ix][iy],No.cell[2]); }
    VV(oldncell,=No.cell) }

  /* clear pre-allocated 3D array of cells */
  loop (ix,0,No.cell[0])
    loop (iy,0,No.cell[1])
      loop (iz,0,No.cell[2]) list[ix][iy][iz]=NULL;

  l=site_storage;

  loop (n,0,No.N) if ((sp=(mn=molec+n)->sp)>=0) {
    ns=mn->ns;
    r=rof(mn,rp);
#if PARALLEL==1
    f1=rof(mn,f1p);
    f2=rof(mn,f2p);
#else /*# PARALLEL==1 */
    f=rof(mn,fp);
#endif /*#!PARALLEL==1 */
#ifdef POLAR
    rpol=(vector*)((char*)r+polar_off);
#  if PARALLEL==1
    f1pol=(vector*)((char*)f1+polar_off);
    f2pol=(vector*)((char*)f2+polar_off);
#  else /*# PARALLEL==1 */
    fpol=(vector*)((char*)f+polar_off);
#  endif /*#!PARALLEL==1 */
#endif /*# POLAR */
    loop (i,0,ns) {
      /* store r in the list item, normalize to [0,L) */
#ifdef WORM
#  define IF while
#else /*# WORM */
#  define IF if
#endif /*#!WORM */
      l->r[0]=r[i][0]; IF (l->r[0]>=box.L[0]) l->r[0]-=box.L[0];
      l->r[1]=r[i][1]; IF (l->r[1]>=box.L[1]) l->r[1]-=box.L[1];
      l->r[2]=r[i][2];
#ifndef SLIT
                       IF (l->r[2]>=box.L[2]) l->r[2]-=box.L[2];
#endif /*# SLIT */

#define CELL(IX,I) \
      if (box.r0>0) \
        if (l->r[I]<-box.r0) { \
          WARNING(("linkcell: loss of precision in this MD step because site out of box\n\
*** mol.site=%d.%d %c=%g < -box.r0\n*** box.r0 reset: old=%g, new=%g\n\
*** (box.r0 is a growable copy of user-given box.rmin)\n\
*** to avoid precision loss, specify box.rmin in the input data",\
                   n,i,'x'+I,l->r[I],box.r0,-l->r[I]*1.1)) \
          box.r0=-l->r[I]*1.1; } \
        IX=(int)(l->r[I]*rL[I]); \
      if (IX>=No.cell[I] || IX<0) ERROR(("mol.site=%d.%d %c=%g out of box",n,i,'x'+I,l->r[I]))

      CELL(ix,0)
      CELL(iy,1)
      CELL(iz,2)

#ifdef POLAR
      VV(l->rpol,=rpol[i])
#endif /*# POLAR */

      /* if LINKCELL&1,2 then pointers to indirected forces assigned,
         else physical forces are cleared */
#if PARALLEL==1

#  if LINKCELL&1
	l->f1=f1[i];
#  else
	VO(l->f1,=0)
#  endif
#  if LINKCELL&2
	l->f2=f2[i];
#  else
	VO(l->f2,=0)
#  endif
#  if LINKCELL==0
	l->fptr=f1[i];
#  endif
#  ifdef POLAR
#    if LINKCELL&1
	l->f1pol=f1pol[i];
#    else
	VO(l->f1pol,=0)
#    endif
#    if LINKCELL&2
	l->f2pol=f2pol[i];
#    else
	VO(l->f2pol,=0)
#    endif
#  endif /*# POLAR */

#else /*# PARALLEL==1  */

#  if LINKCELL&1
	l->f=f[i];
#  else
	VO(l->f,=0) l->fptr=f[i];
#  endif
#  ifdef POLAR
#    if LINKCELL&1
	l->fpol=fpol[i];
#    else
	VO(l->fpol,=0)
#    endif
#  endif /*# POLAR */

#endif /*#!PARALLEL==1  */

      l->sp=sp;
      l->n=n;
      l->si=&(spec[sp]->si[i]);
      l->next=list[ix][iy][iz];
      list[ix][iy][iz]=l;
      l++; } }

  if (box.r0>box.r0limit) {
    if (box.r0>R0LIMIT) ERROR(("\
limit box.r0=%g for sites in the basic cell is too large\n\
*** check box.rmin, tau.rho, tau.P, norm, cutoff vs. box, integrator",box.r0))
    else WARNING(("\
limit box.r0=%g for sites in the basic cell > %g\n\
*** check box.rmin, tau.rho, tau.P, norm, cutoff vs. box, integrator\n\
*** (warning limit increased)",box.r0,box.r0limit))
           box.r0limit*=2;
    Min(box.r0limit,R0LIMIT) }

  if (No.occup) {
    static int *hist,nhist;
    int n;

    if (!nhist) {
      nhist=No.occup;
      if (nhist<2) nhist=No.s*3/PROD(No.cell)+3;
      allocarray(hist,nhist+1); }

    arrayzero(hist,nhist+1);

    loop (ix,0,No.cell[0])
      loop (iy,0,No.cell[1])
        loop (iz,0,No.cell[2]) {
          n=0;
	  looplist (l,list[ix][iy][iz]) n++;
	  Min(n,nhist)
	  hist[n]++; }
    prt("cell occupancy:");
    loopto (n,0,nhist) if (hist[n]) prt("%3d  %d",n,hist[n]); }

  return list;
} /* linklist */


#if 0
int testlist(char *msg,linklist_t ***list) /*********************** testlist */
/* # of sites in the list is returned */
{
  int i,j,s=0;
  linklist_t *l;

  loop (i,0,No.cell[0])
    loop (j,0,No.cell[1])
      looplist (l,list[i][j]) {
	s++;
	if (l->m->st<0 || l->m->st>NSITES) {
	  debugf("error %s %i %i %i %6.3f %6.3f %6.3f",
		 msg,l->m->st,l->m->indx,l->m->desttr,
		 l->m->r[0],l->m->r[1],l->m->r[2]);
	  Error(""); }
        }

  return s;
} /* testlist */
#endif /*# 0 */

static int logint(double x)
{
  int i=x+0.9999;

  if (i>1)
    if (x/(i-1)<i/x) i--;

  return i;
}

void lcsetup(vector L) /******************************************** lcsetup */
/* called from main once at start */
{
  int i,prdc=1,nundef=0,xyz=0;
  double prdL=1;
  vector cell;
  double x,m;

  underline("linked-cell list setup");
  /* used to be: 2.5+log(log(No.s)) */
  /* approximately optimum number of sites / simulation cell */
  if (No.percell<=0) No.percell=1+log(No.s);

  loop (i,0,DIM)
    if (No.cell[i]==0) {
      xyz |= 1<<i; nundef++;
      prdL*=L[i]; }
    else
      prdc*=No.cell[i];

  if (nundef) {
    loop (i,0,DIM) if (No.cell[i]==0) {
      x=pow(No.s/(No.percell*prdc*prdL),1./nundef);
#if PARALLEL==1
      if (!No.th) ERROR(("unknown number of threads - why?"))
      if (i==0)
        No.cell[i]=No.th*logint(x*L[i]/No.th);
      else
#endif
        No.cell[i]=logint(x*L[i]);
      Max(No.cell[i],2)
      prdL/=L[i];
      prdc*=No.cell[i];
      nundef--; } }

  prt("No.percell=%g -> %g, used to set No.cell[] in coordinates: %s",
      No.percell,(double)No.s/prdc,xyz?prtxyz(xyz):"none");
  prt("No.cell=(%d,%d,%d), (%.3g,%.3g,%.3g) in cutoff, (%.3g,%.3g,%.3g) final",
      No.cell[0],No.cell[1],No.cell[2],
      box.cutoff*No.cell[0]/box.L[0],box.cutoff*No.cell[1]/box.L[1],box.cutoff*No.cell[2]/box.L[2],
      box.cutoff*No.cell[0]/    L[0],box.cutoff*No.cell[1]/    L[1],box.cutoff*No.cell[2]/    L[2]);

  m=0;
  loop (i,0,DIM) cell[i]=L[i]/No.cell[i];
  x=cell[0]/cell[1]; if (x<1) x=1/x;
  m=x;
  x=cell[1]/cell[2]; if (x<1) x=1/x;
  Max(m,x)
  x=cell[0]/cell[2]; if (x<1) x=1/x;
  Max(m,x)
  if (m>2.5) WARNING(("Cell too elongated (max.aspect=%g)\n\
*** Check the linked-cell list setup!",m))
  else if (m>1.6) prt("WARNING: cell probably too elongated (max.aspect=%g), check the linked-cell list setup",m);
  else prt("maximum cell aspect=%g (should not differ too much from 1)",m);

  prt("parameter to predict 1-4 interaction radius box.over14=%g",box.over14);

  if (box.r0limit==0) box.r0limit=0.5;

  if (box.rmin<0) {
    box.r0=0;
    prt("box.rmin<0 specified => x,y,z>=0 for all atoms assumed by the linked-cell list");
    if (tau.P!=0 || tau.rho!=0 || option('m')>2 || No.depend[DEP_R] ||  No.depend[DEP_L])
      WARNING(("You specified box.rmin<0 with variable box, Gear, or out-of-plane dependants,\n\
*** which may cause negative coordinates. Some non-bonded interactions\n\
*** for larger distances my be omitted.")) }
  else if (box.rmin>0) {
    box.r0=box.rmin;
    prt("box.rmin=%g specified => x,y,z>=0 for all atoms assumed by the linked-cell list",box.rmin); }
  else {
    if (tau.P==0 && tau.rho==0 && option('m')==2) box.r0=0;
    if (tau.P!=0 || tau.rho!=0) box.r0=0.05;
    if (option('m') || No.depend[DEP_R] ||  No.depend[DEP_L]) box.r0=0.15;
    prt("NOTE: box.rmin has not been specified, I guess box.rmin=%g",box.r0); }

  prt("Explanation: the linked-cell list method requires that all coordinates are\n\
positive, but some methods (variable box, Gear prediction) may slightly violate\n\
this condition. box.rmin (internally box.r0) is the allowed margin.\n");
}
