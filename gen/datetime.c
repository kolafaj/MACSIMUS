#include <time.h>
/*
  gets time and returns time string as "Wed Jun 30 21:49:08 1993"
  (without that stupid \n at the end)
  from,to = character range (0,0=whole string)
  to=25 leaves \n
Mon Sep 10 18:05:36 2001
          11111111112222
012345678901234567890123
*/

char *datetime(int from,int to)
{
 static char lt[26],*c;
 time_t t;

 if (to<=0 || to>25) to=24;
 if (from<0 || from>25) from=0;

 time(&t);
 c=ctime(&t);
 strcpy(lt,c);
 lt[to]=0; 
 return lt+from;
}
