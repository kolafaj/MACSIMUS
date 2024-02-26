/* cc -O2 -o ren ren.c
 */
#include "../gen/include.h"

int main(int narg,char **arg) /**************************************** main */
{
  char *opt="",*op="mv";
  int from=3,i,ok=0,ko=0;
  int Ask=0,ask=1,dry=0,verb=1,all=0,cp=0,never=0;
  enum mode_e { ANYWHERE=0,BEGIN,END,PREFIX } mode = ANYWHERE;

  if (narg<4) {
    fprintf(stderr,"\
Rename files, replacement for the old good \"rename\". Call by:\n\
  ren OLD NEW [-MERGED_OPTION_STRING] FILE [FILE ...]\n\
Substring OLD (no regexp) is replaced by NEW (but option -n).\n\
OPTION_STRING (merged as e.g. -fq):\n\
  --    empty option (good if FILE starts with '-')\n\
  -a    replace all substring occurrences (default=1st only)\n\
  -A    always ask (even if the target does not exist)\n\
  -b    substring OLD can appear at the beginning only\n\
  -c    copy (cp) instead of rename (mv)\n\
  -d    dry run (do not use -q -Q)\n\
  -e    substring OLD can appear at the end only\n\
  -f -y force rename even if the target exists (default=ask)\n\
  -n    never rename if the target exists (ask ignored)\n\
  -p    flags-preserving copy (cp -p) instead of rename (mv)\n\
  -q -Q quiet (only final statistics printed), very quiet\n\
  -x    replace prefix (before OLD) to NEW (useful in scripts)\n\
Examples:\n\
  ren -2- -3- abc-2-x.zip       # the same as \"mv abc-2-x.zip abc-3-x.zip\"\n\
  ren .bak \'\' -e *              # strip off all suffixes .bak\n\
  ren \'\' X- -yb sys-*.txt       # prepend \"X-\", overwrite X-sys-?.txt if any\n\
  ren sys- X-sys- -yb *.txt     # the same as above\n\
  ren . goodname -n badname.jpg # will rename badname.jpg -> goodname.jpg\n\
See also:\n\
  lc mvs mv. difdir mvmess\n\
");
    exit(0); }

  if (arg[3][0]=='-') {
    opt=arg[3],from=4;
    if (strchr(opt,'a')) all=1;
    if (strchr(opt,'A')) Ask=1;
    if (strchr(opt,'c')) cp=1,op="cp";
    if (strchr(opt,'p')) cp=2,op="cp -p";
    if (strchr(opt,'q')) verb=0;
    if (strchr(opt,'Q')) verb=-1;
    if (strchr(opt,'b')) {
      if (mode==END || mode==PREFIX)
        Error("options -ebx are mutually exclusive");
      mode=BEGIN; }
    if (strchr(opt,'d')) dry=1;
    if (strchr(opt,'e')) {
      if (mode==BEGIN || mode==PREFIX)
        Error("options -ebx are mutually exclusive");
      mode=END; }
    if (strchr(opt,'n')) never=1;
    if (strchr(opt,'x')) {
      if (mode==BEGIN || mode==END)
        Error("options -ebx are mutually exclusive");
      mode=PREFIX; }
    if (strpbrk(opt,"fy")) ask=0; }

  if (never) ask=1;

  if (mode && all) Error("-a with -ebx: only one replacement allowed");

  loop (i,from,narg) if (strlen(arg[1])<=strlen(arg[i])) {
    char *match=NULL;
    FILE *f;
    char newname[1024];
    char line[1024];

    switch (mode) {
      case PREFIX:
        match=strstr(arg[i],arg[1]);
        break;
      case ANYWHERE:
        match=strstr(arg[i],arg[1]);
        break;
      case BEGIN:
        if (!memcmp(arg[i],arg[1],strlen(arg[1]))) match=arg[i];
        break;
      case END:
        match=arg[i]+strlen(arg[i])-strlen(arg[1]);
        if (memcmp(match,arg[1],strlen(arg[1]))) match=NULL;
        break;
      default:
        Error("internal"); }

    if (strlen(arg[i])+strlen(arg[2])>512) Error("too long name");

    if (match) {
      int nrepl=0;

      if (mode==PREFIX) {
        strcpy(newname,arg[2]);
        strcat(newname,match); }
      else {
        strcpy(newname,arg[i]);
        strcpy(newname+(match-arg[i]),arg[2]);
        strcat(newname,match+strlen(arg[1])); }

      if (all)
        while ( (match=strstr(newname,arg[1])) ) {
          memmove(match+strlen(arg[2]),match+strlen(arg[1]),strlen(match+strlen(arg[1]))+1);
          memcpy(match,arg[2],strlen(arg[2]));
          if (strlen(newname)>256) Error("too long name while replacement");
          if (nrepl++>256) Error("too many replacements (recursion?)"); }

      if (!strcmp(arg[i],newname))
        Error("the new name is identical to the old one");
      if (strlen(newname)>255)
        Error("the final name is longer than 255 B (ext2 limit)");

      if (Ask) {
        printf("\n%s \"%s\" \"%s\"\n",op,arg[i],newname);
        printf("Should I run the above command (y/N/a)? ");
        gets(line); line[0]=tolower(line[0]);
        if (line[0]=='a') Ask=0,line[0]='y';
        if (line[0]!='y') goto ex; }

      if (ask) {
        f=fopen(newname,"r");
        if (f) {
          fclose(f);
          if (never) {
            if (verb>0)
              printf("\"%s\" not renamed to existing \"%s\"\n",arg[i],newname);
            goto ex; }
          printf("%s exists. Overwrite (y/N/a)? ",newname);
          gets(line); line[0]=tolower(line[0]);
          if (line[0]=='a') ask=0,line[0]='y';
          if (line[0]!='y') goto ex; } }

      ko++;
      if (dry) {
        if (verb>0) printf("%s \"%s\" \"%s\"\n",op,arg[i],newname); }
      else if (cp) {
        static char sys[2048];

        sprintf(sys,"%s \"%s\" \"%s\"",op,arg[i],newname);
        if (system(sys))
          fprintf(stderr,"%s %s %s failed\n",op,arg[i],newname);
        else {
          ok++;
          if (verb>0)
            printf("\"%s\" copied to \"%s\"\n",arg[i],newname); } }
      else {
        if (rename(arg[i],newname))
          fprintf(stderr,"rename %s %s failed\n",arg[i],newname);
        else {
          ok++;
          if (verb>0)
            printf("\"%s\" renamed to \"%s\"\n",arg[i],newname); } }

      ex:; } }

  if (verb<0) return 0;

  ko-=ok;
  if (dry)
    printf("DRY RUN: %d file%c to be renamed or copied\n",ko,"s"[ko==1]);
  else if (ko)
    printf("%d file%c %s and %d failed\n",ok,"s"[ok==1],cp?"copied":"renamed",ko);
  else
    printf("%d file%c %s\n",ok,"s"[ok==1],cp?"copied":"renamed");

  return 0;
}
