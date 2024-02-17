#include "cpmark.h"

typedef struct rec_s {
  struct rec_s *next;
  float r[1]; /* [NCP] */
} rec_t;

float *getcpheader(FILE *f,int *NCP,int verbose);

int4 appendpakcp(
  FILE *f,    /* file must be opened for writing in binary mode */
  int NCP,    /* # of CP items recorded */
  int nbit,   /* resolution is 2^nbit in the min-max range */
  void *data, /* if norec=0 then the head of rec_t* list 
                 if norec>0 then pointer to r[norec][NCP] (one contig.array) */
  int4 norec);

int getpakcp(FILE *f); 
int nextcprec(float *r);
void endgetpakcp(FILE *f);

void endian(char *a);
