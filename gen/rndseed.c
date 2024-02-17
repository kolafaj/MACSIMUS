#ifndef RNDSEED_C
#  define RNDSEED_C
#  include <time.h>

#  ifdef CHEAPTIME
/* standard method with resolution of 1 s */
unsigned seed_from_time(void)
{
  time_t t;
  
  time(&t);
  
  return t;
}  
#  else /*# CHEAPTIME */
/* 
  better initial rnd seed using microseconds 
  it is assumed that initializations are less frequent than 0.0001 s
  the same seed may occur in 5 days again
*/
#    include <sys/time.h>
unsigned seed_from_time(void)
{
  struct timeval tv;
  struct timezone tz;
  
  gettimeofday(&tv,&tz);
  // fprintf(stderr,"%d\n",(int)tv.tv_sec);
  return (unsigned)tv.tv_usec/100+10000U*(unsigned)tv.tv_sec;
}  
#  endif /*#!CHEAPTIME */
#endif /*# RNDSEED_C */
