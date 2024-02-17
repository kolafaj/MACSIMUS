/* interrupt handler */

#include <signal.h>
#include "ground.h"
#include "simglob.h"

volatile int sig=0;

static char *action(int sig) /*************************************** action */
{
  return 
    sig==1 ? "next data" :
    sig==2 ? "finish cycle" :
    sig==5 ? "break and save" : "?";
}

void asksig(int signr) /********************************************* asksig */
{
#if PARALLEL==0
  static volatile int countdown=3;

  countdown--;
  if (countdown<=0) {
    fprintf(stderr,"Too many interrupts. Stop.\n");
    exit(-1); }
#endif /*# PARALLEL==0 */

  signal(signr,asksig);

 badsig:

  if (option('s') && option('b')>2) fputc('\a',stderr);
  fprintf(stderr,"\nINTERRUPTED(SIG=%d) - ",signr);

  if (option('i')) {
    sig=option('i');
    fprintf(stderr,"action = %d = %s\n",sig,action(sig));
    prt("\nINTERRUPTED(%d) - action = %d = %s\n",signr,sig,action(sig));
    return; }

  /* BUG: this does not work with PARALLEL (never use -i0) */

  fprintf(stderr,"select number or character\n\
 0 = (r)esume = (c)ontinue calculations\n\
 1 = (i)nterrupt: finish cycle then save all and read next data\n\
 2 = (.)stop: finish cycle then save all and stop\n\
 5 = (j)ump to stop: finish step then save all and stop (incomplete cycle!)\n\
-1 = (e)(x)it immediately = (q)uit (nothing saved!)\n");
#ifdef SCR
  if (option('s')) fprintf(stderr," 9 = (s)croll (then type ? for help)\n");
#endif /*# SCR */

  {
    char s[16];
    
    do {
      if (!fgets(s,16,stdin)) { sig=2; return; }
    } while (s[0]=='\n');
    sig=atoi(s);
    switch (tolower(s[0])) {
      case 'c': case 'r': sig=0; break;
      case 'i': sig=1; break;
      case '.': sig=2; break;
      case 'j': sig=5; break;
      case 'x': case 'q': case 'e': sig=-1; break;
      case 's': sig=9; break; }
  }

  switch (sig) {
#ifdef SCR
    case 9: scroll(); sig=0;
#endif /*# SCR */
    case 0: case 1: case 2: case 5: break;
    case -1: exit(-1);
    default: goto badsig; }

  if (option('s') && option('b')>2) fputc('\a',stderr);
  fprintf(stderr,"%i selected\n",sig);

#if PARALLEL==0
  countdown=3; /* refresh ^C countdown */
#endif
}
