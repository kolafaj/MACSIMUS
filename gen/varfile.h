#include "int4.h"

/*
  NOTE: quite obsolete piece of software kept here for compatibility reasons
  04/2016: #ifdef sdsalloc used instead of #ifdef SDS (sdsalloc must be macro)
           Reverse removed (see old+misc/varfile.h)
  02/2009: extern variables merged to one structure
  04/2004: VarRead macro: added (void*) to suppress overcorrect warning of gcc
*/

extern struct varfile_s {
/* global: */
  FILE *file;
  int4 size;   /* size of last written/read variable */
  int eof;     /* EOF flag */
} VarFile;  

void VarOpen(char *fn, char *mode);
void VarPut(void *v, int4 size);
void VarReadN(char *name,void *v, int4 size);
#define VarRead(N,S) VarReadN(#N,(void*)(N),S)
void VarClose(void);
int4 KeyClose(int4 key);

#ifdef sdsalloc
void VarPutSds(void *s);
void *VarGetSds(void); /* allocates and reads sds */
void VarReadSds(void *s); /* reads sds that has been allocated - size must match */
#endif
