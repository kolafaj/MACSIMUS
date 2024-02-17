/* cc -O2 -o lc lc.c
  rename file(s) to lowercase, (c) J.Kolafa 1996-2008
  2022: -m,-t,-s added
  2014  -z added
  2008: -8 changed to -7, -8 added (utf-8)
  2008: -f, -x, -L added
  2005: rename via other file (vfat)
  2001: also Capitalize and UPPERCASE
  1999: check of file overwrite
*/
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "../gen/mygets.c"
#include <time.h>
#include <sys/types.h>
#include <utime.h>
int utimes(const char *filename, const struct timeval times[2]); // missing header

enum {NONE,LOWER,CAPITALIZE,UPPER,EXT,UPPER2LOWER,PHOTO,REPRINT} style = NONE;
enum {DIRECT,BIT7,CP1250,ISO88592,CP852,ZIP,UTF8,HASHU} code = DIRECT;

#define N 1024

void iso2asc(char *c) /********************************************* iso2asc */
{
  char asc[128]="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~A~L$LS~~SSTZ~ZZ~a~l~ls~~sstz~zzRAAAALCCCEEEEIIDDNNOOOOxRUUUUYTSraaaalccceeeeiiddnnoooo-ruuuuyt";

  for (; *c; c++) *c = *c&128 ? asc[*c&127] : *c;
}

void win2asc(char *c) /********************************************* win2asc */
{
  char asc[128]="~~~~~~~~~~S~STZZ~~~~~~--~~s~stzz~~~L~A~S~~S~~~~Z~~~l~~~~~as~L~l~RAAAALCCCEEEEIIDDNNOOOO~RUUUUYTSraaaalccceeeeiiddnnoooo~ruuuuyt";

  for (; *c; c++) *c = *c&128 ? asc[*c&127] : *c;
}

void dos2asc(char *c) /********************************************* dos2asc */
{
  char asc[128]="CueaauccleOoiZACELlooLlSsOUTtL caiouAaZzEe zCs~~~~~~~AAES~~~~Zz~~~~~~~Aa~~~~~~~ dDDEdNIIe~~~~TU~O~ONnnSsRUrUyYt      S~~~~~uRr~~";

  for (; *c; c++) *c = *c&128 ? asc[*c&127] : *c;
}

void zip2asc(char *c) /********************************************* zip2asc */
{
  /*             20      21      22      23      24      25      26      27      30      31      32      33      34      35      36      37
                 01234567012345670123456701234567012345670123456701234567012345670123456701234567012345670123456701234567012345670123456701234567*/
  char asc[128]="~~~c~~~~~~S~STZZ~~~~~~--~~s~stzz~~~L~A~S~~S~~~~Z~~rl~~~~~az~L~l~RAAAALCCCEEEEIIeDNNOOOO~RUUUUYTSuaaaalccceeeeiiddnnoooo~ruuuuys ";

  for (; *c; c++) *c = *c&128 ? asc[*c&127] : *c;
}

int tenpow[10]= {1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000};

#include "4lc.c"

char *mediaprefix;
int removedigits;
struct timeval mediatimeval[2];
struct tm mediatm;
int touch=-9999;
int summer=0;

void media(char *s,char *end) /*************************************** media */
{
  char fn[N*2],*f,*f0,*c;;
  int ndig=0,suffix=0,maxndig;

  for (c=s; *c; c++ ) {
    if (*c=='.') break;
    if (isdigit(*c)) ndig++; }

  if (ndig<12) {
    /* not renamed */
    strcpy(fn,s);
    return; }

  memset(fn,0,N*2);
  memset(&mediatm,0,sizeof(mediatm));
  strcpy(fn,mediaprefix);
  f0=fn+strlen(mediaprefix);

  if ( ndig==14 || ndig<17 && strchr(s+19,'-') )
    /* 2022-07-19 18.16.58-2.jpg
       01234567890123456789      */
    maxndig=14;
  else
    /* PXL_20220723_083731281.jpg */
    maxndig=ndig;

  ndig=0;
  for (f=f0,c=s; *c; c++ ) {
    if (ndig<maxndig && isdigit(*c)) {
      if (ndig++>=removedigits) *f++=*c;
      if (ndig<=4)       mediatm.tm_year=10*mediatm.tm_year+(*c-'0');
      else if (ndig<=6)  mediatm.tm_mon =10*mediatm.tm_mon +(*c-'0');
      else if (ndig<=8)  mediatm.tm_mday=10*mediatm.tm_mday+(*c-'0');
      else if (ndig<=10) mediatm.tm_hour=10*mediatm.tm_hour+(*c-'0');
      else if (ndig<=12) mediatm.tm_min =10*mediatm.tm_min +(*c-'0');
      else if (ndig<=14) mediatm.tm_sec =10*mediatm.tm_sec +(*c-'0');
      else               mediatm.tm_yday=10*mediatm.tm_yday+(*c-'0');
      if (ndig==8) *f++='_'; /* to separate DATE_TIME */ }
    if (ndig==maxndig) break;
  }

  /* all maxndig digits have been copied */
  c++;

  if (maxndig==14 && *c=='-') {
    /* dropbox-style with suffix -#*/
    suffix=atoi(++c);
    if (suffix<1 || suffix>26) {
      fprintf(stderr,"lc: Dropbox-style name with suffix %d out of range\n",suffix);
      exit(4); }
    c+=1+(suffix>9);
    *f++='`'+suffix; }

  for (; *c; c++ ) *f++=*c;

  if (s+strlen(fn)>=end) {
    fprintf(stderr,"lc: re-coded name too long");
    exit(4); }

  if (ndig>10 && (mediatm.tm_year<2000 || mediatm.tm_year>2060))  {
    fprintf(stderr,"lc: year %d is out of expected range\n",mediatm.tm_year);
    exit(2); }

  if (touch>-9999) {
    mediatm.tm_year-=1900;
    mediatm.tm_mon--;
    mediatm.tm_isdst=summer; /* no summer time, will be solved below */

    if (ndig>14 && ndig<=20) mediatimeval[0].tv_usec=mediatm.tm_yday*tenpow[20-ndig];
    mediatm.tm_yday=0;
    mediatimeval[0].tv_sec=mktime(&mediatm)+3600*touch;

    if (mediatimeval[0].tv_sec<0)   {
      fprintf(stderr,"lc: bad time (cannot convert timeval structure)\n");
      exit(2); }

    mediatimeval[1]=mediatimeval[0]; }

  strcpy(s,fn);
}

void utf82asc(char *s,char *end) /********************************* utf82asc */
{
  char fn[N*2],*f=fn,*c;
  struct utf8_s *u;

  memset(fn,0,N*2);

  for (c=s; *c; ) {
    for (u=utftab; u->utf; u++)
      if (u->eq[0] && !memcmp(c,u->utf,strlen(u->utf))) {
        strcat(fn,u->eq);
        f+=strlen(u->eq);
        c+=strlen(u->utf);
        goto cont; }
    *f++=*c++;
    cont:; }

  if (s+strlen(fn)>=end) {
    fprintf(stderr,"recoded name too long");
    exit(4); }

  strcpy(s,fn);
}

int ncmemcmp(char const *a, char const *b,int len) /*************** ncmemcmp */
{
  int i;

  for (i=0; i<len; i++) {
    int d = tolower((unsigned)a[i]) - tolower((unsigned)b[i]);
    if (d) return d; }

  return 0;
}

void hashu2asc(char *s,char *end) /******************************* hashu2asc */
{
  char fn[N*2],*f=fn,*c;
  struct utf8_s *u;

  memset(fn,0,N*2);

  for (c=s; *c; ) {
    for (u=utftab; u->utf; u++)
      if (!ncmemcmp(c,u->U,6)) {
        strcat(fn,u->utf);
        f+=strlen(u->utf);
        c+=6;
        goto cont; }
    *f++=*c++;
    cont:; }

        fprintf(stderr,"===%s===\n",fn);
  if (s+strlen(fn)>=end) {
    fprintf(stderr,"recoded name too long");
    exit(4); }

  strcpy(s,fn);
}

int main(int narg,char **arg) /**************************************** main */
{
  int i,dry=0,vfat=0,haslc,repstat;
  char *reprint="--";
  char fn0[N],fn[N],*a,*f,*ar,xfrom=0,xto='-',yes=0;
  FILE *list,*over;
  char space=' ',*dot,*from;

  if (narg<2) {
    fprintf(stderr,"\
Rename files to lowercase, re-code names, etc.  Call by:\n\
  lc [OPTIONS] [@][PATH/]FILE [OPTIONS] [[@][PATH/FILE...]\n\
Re-coding (performed first):\n\
  -7 = remove bit 8 (-> ~ if not possible)\n\
  -8 = remove accents or try most similar from UTF-8 coding\n\
  -i = remove accents from ISO-8859-2 characters\n\
  -p = remove accents from CP852 (PC Latin II) characters\n\
  -U = change html-like codes (as #U011B,#u011b) to accented charactes\n\
  -w = remove accents from CP1250 (windows) characters\n\
  -z = remove accents from some unzipped files (undocumented)\n\
Character change:\n\
  -xAB = change char A into B (only one -x option)\n\
Case change (-nX etc.: change space into X, performed after -xAB)\n\
  -c = Capitalize\n\
  -e = lowercase extension (Aa.DOC->Aa.doc)\n\
  -f = lowercase first 6 chars and extension (IMG123-Az.JPG->img123-Az.jpg)\n\
  -l = lowercase\n\
  -L = Lowercase only names with all letters uppercase (A.SH->a.sh, Aa.c->Aa.c)\n\
  -n = no case change (default)\n\
  -rXAB = uppercase between chars AB (default AB=--), otherwise lowercase\n\
  -u = UPPERCASE\n\
Media rename+touch:\n\
  -mX = rename names containing date with separators, add/replace prefix X\n\
  -m-X = as above + remove year (4 digits)\n\
  -t#  = touch (-t = -t0 = UTC) or with time shift # h (1 for CET, 2 for CEST)\n\
         good for Pixel 6 with UTC time\n\
  -s   = set summer time for -t (-t2 -s = -t1)\n\
  NB: both access and modification times are changed BEFORE rename\n\
More options:\n\
  -d = dry run (show what would be renamed but do not actually rename)\n\
  -v = rename via ./vfat.aux: for case-insensitive filesystem (vfat)\n\
  -y = replace target without asking (default=asks for confirmation)\n\
File args:\n\
  PATH/ = path to file (not renamed)\n\
  FILE = file to rename\n\
  @FILE = list of [PATH/]FILEs to rename\n\
Examples:\n\
Lowercase all files beginning with uppercase on a Windows-compatible disk:\n\
  LC_ALL=C.UTF-8 lc -v -l_ /media/MyDisk/[A-Z]*\n\
Remove accents from UTF8 file names, change spaces to _, do not change case:\n\
  lc -n_ -8 *\n\
See also:\n\
  mvs ren kostnice kostnice2 utfpine renwinzip.sh\n\
");
    exit(0); }

  for (i=1; i<narg; i++) {

    if (arg[i][0]=='-') {
      switch (arg[i][1]) {
        case '7': code=BIT7; break;
        case '8': code=UTF8; break;
        case 'c': style=CAPITALIZE; break;
        case 'd': dry++; break;
        case 'e': style=EXT; break;
        case 'f': style=PHOTO; break;
        case 'i': code=ISO88592; break;
        case 'l': style=LOWER; break;
        case 'L': style=UPPER2LOWER; break;
        case 'm':
          if (arg[i][2]=='-') {
            mediaprefix=arg[i]+3;
            removedigits=4; }
          else
            mediaprefix=arg[i]+2;
          break;
        case 'n': style=NONE; break;
        case 'p': code=CP852; break;
        case 'r': style=REPRINT;
          if (arg[i][3]) reprint=arg[i]+3;
          break;
        case 't': touch=atoi(arg[i]+2); break;
        case 's': summer=1; break;
        case 'u': style=UPPER; break;
        case 'U': code=HASHU; break;
        case 'v': vfat++; break;
        case 'w': code=CP1250; break;
        case 'x': xfrom=arg[i][2],xto=arg[i][3]; break;
        case 'y': yes++; break;
        case 'z': code=ZIP; break;
        default: fprintf(stderr,"%s is bad option\n",arg[i]); exit(0); }
      if (strchr("nclLuefr",arg[i][1]) && arg[i][2]) space=arg[i][2];
      continue; }

    else if (arg[i][0]=='@') {
      list=fopen(arg[i]+1,"rt");
      if (!list) {
        fprintf(stderr,"%s not found\n",arg[i]+1);
        exit(-1); } }
    else
      list=NULL;

    if (xfrom && !xto) {
      fprintf(stderr,"bad -x\n");
      exit(-1); }

    if (touch>-9999) {
      /* means do not rename but get time from fn */
      if (!mediaprefix) mediaprefix="-~"; }

    if (mediaprefix && strlen(mediaprefix)>=N) {
      fprintf(stderr,"too long prefix in -m\n");
      exit(-1); }

    ar=arg[i];

  again:
    if (list) {
      if (fscanf(list,"%s",fn0)!=1) goto done;
      ar=fn0; }

    repstat=0;

    if (strlen(ar)>=N) {
      fprintf(stderr,"too long arg\n");
      exit(0); }
    strcpy(fn,ar);

    dot=NULL;
    from=fn;
    haslc=0;
    for (a=fn; *a; a++) {
      if (islower(*a)) haslc++;
      if (*a=='.') dot=a;
      if (*a=='/') from=a+1,dot=NULL,haslc=0; }

    /* coding change */
    switch (code) {
      case DIRECT: break;
      case BIT7:
        for (f=from; *f; f++) if (*(unsigned char*)f&128) {
          *(unsigned char*)f&=127;
          if (!isalnum(*f)) *f='~'; }
        break;
      case ISO88592: iso2asc(from); break;
      case CP1250:   win2asc(from); break;
      case ZIP:      zip2asc(from); break;
      case CP852:    dos2asc(from); break;
      case UTF8:     utf82asc(from,fn+N); break;
      case HASHU:    hashu2asc(from,fn+N); break;
    }

    /* lc, uc, capitalize */
    for (f=from; *f; f++) {
      switch (style) {
        case NONE: break;
        case UPPER: *f=toupper(*f); break;
        case LOWER: *f=tolower(*f); break;
        case REPRINT:
          if (repstat==1) *f=toupper(*f);
          else *f=tolower(*f);
          if (*f==reprint[0] && !repstat) repstat=1;
          else if (*f==reprint[1] && repstat) repstat++;
          break;
        case PHOTO:
          if (f<from+6 || f>dot) *f=tolower(*f);
          break;
        case EXT:
          if (f>dot) *f=tolower(*f);
          break;
        case UPPER2LOWER:
          if (!haslc) *f=tolower(*f);
          break;
        case CAPITALIZE:
          if (f==from) *f=toupper(*f); }
      if (xfrom) if (*f==xfrom) *f=xto;
      if (*f==' ') *f=space; }

    /* normalize media files */
    if (mediaprefix) media(from,fn+N);

    if (touch>-9999) {
      if (dry)
        printf("%s would be touched by %ld.%06ld s (since 1970)\n",
               ar,mediatimeval[0].tv_sec,mediatimeval[0].tv_usec);
      else
        utimes(ar, mediatimeval); }

    if (strcmp(ar,fn) && (!mediaprefix || memcmp("-~",mediaprefix,2))) {
      /* rename */
      printf("%s -> %s",ar,fn);
      if (!yes) {
        over=fopen(fn,"r");
        if (over) {
          char line[16];

          fclose(over);
          printf(" : target exists. Overwrite (y/n)? ");
          getsbufsize=16;
          gets(line);
          if (!strchr("yY",line[0])) goto cont; } }

      if (dry)
        printf(" would be renamed\n");
      else {
        if (vfat) {
          if (rename(ar,"vfat.aux")) {
            printf(" : rename failed (no such file, access denied...)\n");
            exit(0); }
          if (rename("vfat.aux",fn)) {
            printf(" : rename failed (no such file, access denied...)\n");
            exit(0); } }
        else if (rename(ar,fn))
          printf(" : rename failed (no such file, access denied...)\n");
        else printf(" renamed\n"); }
    cont:; }

    if (list) goto again;

  done:; } /* i */

  return 0;
}
