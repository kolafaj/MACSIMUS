/* for site-site radial distribution function for one pair of sites */
/* see rdfconv for the old version (not good for > 92681 atoms) */
/* #included by simglob.h and prtrdf.h */

#include "int4.h"

typedef struct /* SDS */ {
  int4 size;        /* sds (compatible w. 16 bit DOS because big endian) */
  char name[2][8];  /* site names */
  int4 indx[2];     /* original site number */
  int4 ns[2];       /* # of sites */
  double npair;     /* # of pairs of given type in the system CHANGED */
  double grid;      /* 1/grid is step in the adopted units */
  double V;         /* sum of volumes (to work also for NPT) */
  int4 nmeas;       /* # of measurements on the system */
  int4 nhist;       /* length of array hist */
  unsigned4 hist[2];/* [nhist] ([2] is here because of purifying alignments) */
} rdf_t;

