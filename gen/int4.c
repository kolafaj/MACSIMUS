/* cc -O0 -o int4 int4.c
 */
#include "int4.h"
#include <stdlib.h>
#include <stdio.h>

int main(int narg,char **arg)
{
  int Int=1,i;
  char *end=(char*)&Int;

  printf("***** LEGACY CHECK OF INTEGER SIZES *****\n\
%s endian  short=%d int=%d long=%d pointer=%d  int4=%d unsigned4=%d\n",
        *end==1?"little":*end==0?"BIG":"?",
         (int)sizeof(short),(int)sizeof(int),(int)sizeof(long),(int)sizeof(void*),
         (int)sizeof(int4),(int)sizeof(unsigned4));

  if (sizeof(int4)!=4 || sizeof(unsigned4)!=4) {
    printf("\
***** WARNING *****\n\
The sizes of int4 and unsigned4 should be 4. Check file gen/int4.h!\n");
    return 1; }
  else
    printf("OK\n");                                                                   

  return 0;
}
