/* make cppak

  see macsimus/util/cpztest

  This is a version using linked-cell list - to be changed (cppaknolist-unfinished.c)

  02/2024: bug list fixed
  12/2019: reverse endian removed, update to new time stamps
  6/2017: -i added/fixed, default nbit changed
  8/2010: options improved, multiple files supported

  Compress/decompress convergence profile files (cppak0.c, 1995)
*/

#include "ground.h"
#include "bitfile.h"
#include "pakcp.h"

rec_t *rechead;

int main(int narg,char **arg) /*************************************** main */
{
  char nm[256];
  FILE *cp,*pakcp;
  char *c,*dot;
  int NCP,i,iarg;
  int nbit=15;
  int verbose=1,erase=0,overwrite=0,ignore=0;
  rec_t *rec,*last=NULL;
  float *r;
  enum mode_e { AUTO,COMPRESS,DECOMPRESS } mode=AUTO,dir=AUTO;
  char line[128];

  //  AllocTrace++; AllocRange=1024;

  if (narg<2) {
    fprintf(stderr,"\
Lossy compression of MACSIMUS convergence profile.  Call by:\n\
  cppak [OPTIONS] FILE.EXT [[OPTIONS] FILE.EXT ..]\n\
FILES:\n\
  FILE.cp  = MACSIMUS convergence profile file\n\
  FILE.cpz = packed convergence profile file\n\
OPTIONS (apply to the rest of arguments)\n\
  -nNBITS  number of bits in [min,max] (default=%d, max=24)\n\
  -a       compress FILE.cp, decompress FILE.cpz (default mode)\n\
  -d       force decompress\n\
  -c       force compress\n\
  -i       ignore error \"unexpected EOF\", not recommended with -e\n\
  -o       overwrite the target without warning (default = ask)\n\
  -e       erase the original after (de)compression\n\
  -q       quiet\n\
  -y       shorthand for -o -e -q -i\n\
See also:\n\
  cook* with CPnbit=NBITS\n\
  showcp cp2cp plbpak\n",nbit);
  exit(0); }

  initscroll(0);

  loop (iarg,1,narg)
    if (arg[iarg][0]=='-')
      switch (arg[iarg][1]) {
        case 'n':
          nbit=atoi(arg[iarg]+2);
          if (nbit>24 || nbit<6) {
            ERROR(("number of bits (-n) must be in range 6..24\n*** (15 will be used)"))
            nbit=15; }
          break;
        case 'a': mode=AUTO; break;
        case 'c': mode=COMPRESS; break;
        case 'd': mode=DECOMPRESS; break;
        case 'i': ignore=1; break;
        case 'q': verbose=0; break;
        case 'e': erase=1; break;
        case 'y': erase=1; verbose=0; ignore=1;
        case 'o': overwrite=1; break;
        default: ERROR(("%s: unnown option",arg[iarg])) }
    else {
      strcpy(nm,arg[iarg]);
      /* determine extension */
      dot=NULL;
      for (c=nm; *c; c++) if (*c=='.') dot=c;
      dir=mode;
      if (mode==AUTO) {
        if (!dot) ERROR(("%s: no extension and none of -c -d specified",arg[iarg]))
        if (!strcmp(dot,".cp")) dir=COMPRESS;
        if (!strcmp(dot,".cpz")) dir=DECOMPRESS;
        if (dir==AUTO) ERROR(("%s: unknown extension and none of -c -d specified",arg[iarg]))  }

      if (dir==COMPRESS) {

        /*** compression ***/

        prt("compressing %s...",nm);
        cp=fopen(nm,"rb");
        if (!cp) Error("no such file");
        if (dot) strcpy(dot,".cpz"); else strcat(nm,".plz");
        if (!overwrite) {
          pakcp=fopen(nm,"rb");
          if (pakcp) {
            prt("%s exists. Overwrite (y/n)? ",nm);
            gets(line);
            if (toupper(line[0])!='Y') continue;
            fclose(pakcp); } }
        pakcp=fopen(nm,"wb");

        r=getcpheader(cp,&NCP,verbose);
        if (fwrite(r,sizeof(float),NCP,pakcp)!=NCP)
          ERROR(("%s: cannot write",arg[iarg]))

        while (!feof(cp)) {
          rechead=last=NULL;

          for (;;) {
            i=fread(r,sizeof(float),NCP,cp);
            if (i!=NCP) {
              if (i) {
                if (ignore) prt("%s: unexpected EOF ignored (output may be damaged)",arg[iarg]);
                else ERROR(("%s: unexpected EOF",arg[iarg])) }
              if (rechead) break;
              goto done; }

            if (r[0]<=CPmark) break;
            alloc(rec,sizeof(rec_t)-sizeof(float)+NCP*sizeof(float));
            copy(rec->r,r,NCP*sizeof(float));
            if (rechead) last=last->next=rec;
            else rechead=last=rec;
            last->next=NULL; }

          appendpakcp(pakcp,NCP,nbit,rechead,0);

          if (r[0]<=CPmark)
            if (fwrite(r,sizeof(float),NCP,pakcp)!=NCP) Error("cannot write");

          /* free all */
          while (rechead) {
            rec=rechead->next;
            free(rechead);
            rechead=rec; } } }

      else {

        /*** decompression ***/

        prt("decompressing %s...",nm);
        pakcp=fopen(nm,"rb");
        if (!pakcp) Error("no such file");
        if (dot) strcpy(dot,".cp"); else strcat(nm,".cp");
        if (!overwrite) {
          cp=fopen(nm,"rb");
          if (cp) {
            prt("%s exists. Overwrite (y/n)? ",nm);
            gets(line);
            if (toupper(line[0])!='Y') continue;
            fclose(cp); } }
        cp=fopen(nm,"wb");

        r=getcpheader(pakcp,&NCP,verbose);
        if (fwrite(r,sizeof(float),NCP,cp)!=NCP) Error("cannot write");

        while (getpakcp(pakcp)) {
          while (nextcprec(r))
            if (fwrite(r,sizeof(float),NCP,cp)!=NCP) Error("cannot write");
          endgetpakcp(pakcp); } }

    done:
      fclose(pakcp);
      fclose(cp);
      if (erase) {
        fflush(stdout);
        remove(arg[iarg]);
        if (verbose) prt("%s erased.\n",arg[iarg]); } }

  return 0;
}
