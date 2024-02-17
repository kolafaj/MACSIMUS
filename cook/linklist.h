/* C2c linklist.H

  !!! DO NOT EDIT linklist.h !!!
      * edit linklist.H and run C2c linklist.H
      * optionally, run cppnest linklist.h
*/

#if PARALLEL==3
#  error PARALLEL==3 in linklist.h
#endif /*# PARALLEL==3 */

#ifndef LINKCELL
#  error "LINKCELL not #defined"
#endif /*# LINKCELL */

#if LINKCELL!=0 && LINKCELL!=1 && LINKCELL!=3
#  error wrong LINKCELL, should be one of {0,1,3}
#endif /*# LINKCELL!=0 && LINKCELL!=1 && LINKCELL!=3 */

typedef struct linklist_s {
  vector r;      /* site position physically here (periodic b.c. optimizing) */
#  ifdef POLAR
	vector rpol;   /* also Drude physically here */
#  endif

#if PARALLEL==1 
  /*****
   * PARALLEL==1 requires two copies of forces: f1 and f2 should 
   * point to separate arrays when LJqqm() are called in parallel 
   * flags LINKCELL&1,2 control whether these forces are indirected (flag set)
   * or copied to the structure (flag cleared)
   *****/
#  if LINKCELL&1
	real *f1;           /* 1=1st of pair indirected */
#  else
	vector f1;          /* 0=copied here */
#  endif
#  if LINKCELL&2
	real *f2;           /* 2nd of pair indirected */
#  else
	vector f2;          /* 0=copied here */
#  endif

#  ifdef POLAR
  /* once more for forces to Drude particles */
#    if LINKCELL&1
	real *f1pol;        /* 1st of pair indirected */
#    else
	vector f1pol;       /* 0=copied here */
#    endif
#    if LINKCELL&2
	real *f2pol;        /* 2nd of pair indirected */
#    else
	vector f2pol;       /* 0=copied here */
#    endif
#  endif /*# POLAR */

#else /*# PARALLEL==1  */
  /***** serial version
   * flag LINKCELL&1 controls whether the forces are indirected (flag set)
   * or copied to the structure (flag cleared)
   * flag LINKCELL&2 is ignored 
   *****/
#  if LINKCELL&1
	real *f;            /* forces indirected */
#  else
	vector f;           /* forces copied here */
#  endif
#  ifdef POLAR
#    if LINKCELL&1
	real *fpol;         /* forces indirected */
#    else
	vector fpol;        /* forces copied here */
#    endif
#  endif /*# POLAR */
#endif /*#!PARALLEL==1  */

  /* if all forces are copied, we must know a pointer where they originally reside */
#  if LINKCELL==0
	real *fptr;
#  endif

  /* possible cache line padding: */
  //  double dummy[XXX];
#ifdef LINKLIST_PADDING 
  double dummy[LINKLIST_PADDING];
#endif

  siteinfo_t *si;/* = spec[sp]->si[atom_#] */
  struct linklist_s *next;
  int sp;        /* species */
  int n;         /* number of the molecule containing the site */
} linklist_t;

extern linklist_t *site_storage;

#if PARALLEL==1
linklist_t ****linklist(vector *f1p,vector *f2p,vector *rp);
#else /*# PARALLEL==1 */
linklist_t ****linklist(vector *fp,vector *rp);
#endif /*#!PARALLEL==1 */

int testlist(char *msg,linklist_t ***list);

void lcsetup(vector L);
