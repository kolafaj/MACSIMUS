#include "ground.h"
#include "sds.h"
#include "simopt.h"
#include "options.h"

#include <time.h>
#include "cputime.h"

#ifdef CHEAPTIME
/* real time in resolution by 1 s */
double mytime(void) /************************************************ mytime */
{
  time_t t0;
  time(&t0);
  return (double)t0;
}
#else
#include <sys/time.h>
/* real time in better resolution (10^-6 s if provided by the system) */
double mytime(void) /************************************************ mytime */
{    
  struct timeval tv;
  struct timezone tz;

  gettimeofday(&tv,&tz);

  return (unsigned)tv.tv_usec/1e6+tv.tv_sec;
}
#endif

char *myctime(double t) /******************************************* myctime */
{
  time_t it=t;
  return ctime(&it);
}

#define CHARPERCOL 15

typedef struct {
  double t;
  char n[CHARPERCOL+1]; 
} titem_t;

typedef struct {
  int size;
  titem_t item[1/*ntime*/]; 
} timetab_t;

int ntimetabs, ntime, timetablen;
timetab_t *timetab=NULL;
double time0;

double inittime(int maxitems) /************************************ inittime */
{
  if ( (ntime=maxitems) ) {
    sdsalloc(timetab,sizeof(timetab_t)+ntime*sizeof(titem_t));
    sdszero(timetab); }

  time0=mytime();
  return (double)time0;
} /* inittime */


void CPUtime(char *n) /********************************************* CPUtime */
/*
  time since previous call (with any parameter) is recorded
*/
{
  int i;
  double time1=mytime();

  if (timetab) {
    loop (i,0,ntime) {
      if (!timetab->item[i].n[0]) strcpy(timetab->item[i].n,n);
      if (!strcmp(timetab->item[i].n,n)) {
        timetab->item[i].t += time1-time0; goto ex; } }
    Error("time table overflow");
  ex:; }
  time0=time1;
} /* CPUtime */

void printtime(void) /******************************************* printtime */
{
  titem_t *ti;
#ifdef CHEAPTIME
#define FMT "%15s %10.0f %6.2f"
  int sum=0;
#else
#define FMT "%15s %10.3f %6.2f"
  double sum=0;
#endif

  header("       chunk        time/s    %  ");
  loop (ti,timetab->item,timetab->item+ntime) sum+=ti->t;
  if (sum) loop (ti,timetab->item,timetab->item+ntime)
  if (ti->n[0]) prt(FMT,ti->n,ti->t,100.0*ti->t/sum);
  header("");
}

#if PARALLEL
struct partimes_s partimes;
/* note: for PARALLEL==1 partimes[0]=k-space, partimes[1]=r-space */

void printpartimes(void) /************************************ printpartimes */
{
  if (option('t')) {
    double sum,rsum,ksum=0;
    int i;

    if (!partimes.rspace || !partimes.kspace)
      ERROR(("printpartimes: option -t but control structures not initialized"))

    underline("sampled execution times [s]");
    header("thread    r-space    k-space    sum");
    loop (i,0,max(partimes.nrspace,partimes.nkspace)) {
      prt_(" %3d ",i);
      sum=0;
      if (i<partimes.nrspace) {
        rsum+=partimes.rspace[i];
        prt_(" %9.2f",partimes.rspace[i]);
        sum+=partimes.rspace[i]; }
      if (i<partimes.nkspace) {
        ksum+=partimes.kspace[i];
        prt_(" %9.2f",partimes.kspace[i]);
        sum+=partimes.kspace[i]; }
      prt(" %9.2f",sum); }
    prt("-----------------------------------");
    prt(" sum  %9.2f %9.2f %9.2f",rsum,ksum,rsum+ksum);
    header(""); }

  free(partimes.rspace); partimes.nrspace=0;
  free(partimes.kspace); partimes.nkspace=0;
}
#endif
