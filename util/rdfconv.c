/* make rdfconv
*/

#include "ground.h"
#include "sds.h"
#include "varfile.h"

#include "rdf.h"

typedef struct /* SDS */ {
  int4 size;        /* sds (compatible w. 16 bit DOS because big endian) */
  char name[2][8];  /* site names */
  int4 indx[2];     /* original site number */
  int4 ns[2];       /* # of sites */
  unsigned4 npair;  /* # of pairs of given type in the system CHANGED */
  double grid;      /* 1/grid is step in the adopted units */
  double V;         /* sum of volumes (to work also for NPT) */
  int4 nmeas;       /* # of measurements on the system */
  int4 nhist;       /* length of array hist */
  unsigned4 hist[2];/* [nhist] ([2] is here because of purifying alignments) */
} rdfold_t;

struct list_s {
  struct list_s *next;
  rdf_t *rdf; 
} *head,*list;

int main(int narg, char **arg)
{
  rdfold_t *rdfold;
  rdf_t *rdf;
  int4 i4,nrdf=0;
  int off=(char*)(&rdf->grid)-(char*)rdf;
  int offold=(char*)(&rdfold->grid)-(char*)rdfold;
  double d8;

  initscroll(0);

  if (narg<2) {
    fprintf(stderr,"\
Convert old .rdf (produced by cook prior V2.7f) to the new version. Call by:\n\
  %s SIMNAME.rdf\n",arg[0]);
  exit(1); }

  put2(offold,off)

  VarOpen(arg[1],"r");
  while (VarFile.size) {
    int4 vs=VarFile.size;

    if (vs==4) {
      VarRead(&i4,4);
      fprintf(stderr,"nsites=%d\n",i4); }
    if (vs==8) {
      VarRead(&d8,8);
      fprintf(stderr,"grid=%g\n",d8); }

    if (vs>8) {
      nrdf++;
      rdfold=VarGetSds();
      if (rdfold->size<sizeof(rdfold_t)-2*sizeof(unsigned4) || rdfold->size>8000000)
	ERROR(("rdf->size=%d is suspicious - endian?",rdfold->size))

      if (rdfold->nhist*sizeof(rdfold->hist[0])+((char*)rdfold->hist-(char*)rdfold) != rdfold->size)
        ERROR(("This does not look like an RDF file of cook prior V2.7f"))

      fprintf(stderr,"%s-%s\n",rdfold->name[0],rdfold->name[1]);

      alloc(rdf,rdfold->size+(off-offold));
      memcpy(rdf,rdfold,rdfold->size);
      rdf->size+=(off-offold);
      memcpy(&rdf->grid,&rdfold->grid,rdfold->size-offold);
      rdf->npair=rdfold->npair; 

      if (head) {
        alloconezero(list->next);
        list=list->next; }
      else {
        alloconezero(list);
        head=list; }
      list->rdf=rdf;

      free(rdfold); } }

  VarClose();

  if (!nrdf) ERROR(("no RDF data found"))

  system(string("mv \"%s\" \"%s~\"",arg[1],arg[1]));

  VarOpen(arg[1],"w");
  VarPut(&i4,4);
  VarPut(&d8,8);

  looplist (list,head) VarPutSds(list->rdf);

  VarClose();
  prt("%s converted",arg[1]);

  return 0;
}
