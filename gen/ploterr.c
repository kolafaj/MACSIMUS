/*
  2017: slowdraw+key added, also non-key or comment interrupts line
  2009: x error bars added (PLOTERRBAR)
  2006: bug fixed (asymmetric error bars)
  This module draws a file, using formula, 2004.
  Error bars supported.
  Used by plot.
*/
#include "ground.h" /* #define CALC 3 required */
#include "ploterr.h"
#include "xdraw.h"
#include "mydraw.h"

#define LINELEN 16384

/* separator characters; must not include: ! # (and any number components) */
#define SEP " \t\n\r,;"

double *count; /* needed in tabinc.c */
struct plotopt_s plotopt; /* passed from plot */
struct _Idlist_s *id;

int drawerrbar(char *fn,
               char *colx,char *coly,
               char *coldy1,char *coldy2,
               int maxcol,int color,int style,int errstyle)
/*
   fn:
     file name
   colx,coly:
     expressions for x,y
   coldy1,coldy2:
     expressions for error bars (see errstyle)
   maxcol:
     maximum column needed in the expressions = max ID found in :EXPR
   color:
     color&255 = one of 14 standard colors
     color&256 = denotes mode new line in data changes color
   style:
      ...
     -2 = medium point
     -1 = small point
      1 = dotted line
      0 = line, uses previous line style which should be thick=1
      2 = line, thick=3
      3 = line, thick=5
   errstyle:
     0=no error bars
     1=YY
     2=XX
     3=XY
*/
{
  FILE *f;
  char line[LINELEN];
  int i,OK,ncolormode=color&256;
  double xy[4];

  if (!fn || fn[0]==0) f=stdin;
  else f=fopen(fn,"rt");

  if (!f) {
    fprintf(stderr,"ploterr: `%s\': No such file\n",fn);
    return -1; }

  color &= 255;
  mysetcolor(color);

  pointsize=0.798; /* ? */

  if (errstyle) {
    mysetlinestyle(0,0,1);
    if (style<0) pointsize=(-1-style)*0.798; }
  else {
    /* style==0: "setlinestyle(0,0,1)" should have been set */
    if (style==0) mysetlinestyle(2,0x0301,style);
    if (style>1) mysetlinestyle(0,0,2*style-1);
    else if (style<0) my_style.pointsize=pointsize=(-1-style)*0.798; }
  my_style.pointsize=pointsize; /* exported to psstring?  */
  myup();
  *count=0;

  while (fgets(line,LINELEN,f)) {
    char *e,*tok;
    int n; /* number of columns found in the line +1 */

    //    fprintf(stderr,"DEBUG line=%s",line);

    /* not KEY -> blank line (even if comment) */
    if (plotopt.key && !strstr(line,plotopt.key)) tok=NULL;
    else tok=strtok(line,SEP);

    if (!tok || strchr("!#",tok[0])) {
      /* comment or blank -> pen up */
      myup();
      if (!tok) {
        /* blank line resets counter n */
        *count=0;
        if (ncolormode) {
          /* if option -n, blank line also changes color */
          static int plotcolor[14]={
            WHITE,YELLOW,LIGHTCYAN,LIGHTMAGENTA,LIGHTGREEN,LIGHTRED,LIGHTBLUE,
            BROWN,CYAN,MAGENTA,LIGHTGRAY,GREEN,RED,BLUE};
          int c;

          loop (c,0,14) if (color==plotcolor[c]) break;
          c++;
          color=plotcolor[c%14];
          mysetcolor(color); } } }

    else {
      /* _Id.head->val = count, _Id.head->next->val = 1st number, etc. */
      n=1;
      id=_Id.head->next;
      // fprintf(stderr,"  %g %d %s",_Id.head->val,_Id.head->used,_Id.head->id);
      do {
        id->val=strtod(tok,&e);
        if (e==tok) id->val=FP_NAN;
        //        fprintf(stderr,"  %g %d %s",id->val,id->used,id->id);
        id=id->next;
        n++;
        tok=strtok(NULL,SEP);
      } while (tok && id);

      //      fprintf(stderr," more: %g %d %s",id->val,id->used,id->id);

      OK=1;
      loop (i,0,errstyle?4:2) {
        double z;
        char *expr =
          i==0 ? colx :
          i==1 ? coly :
          i==2 ? coldy1 : coldy2;
        e=expr; /* err */
        if (maxcol<n) z=Calc(expr,&e); else OK=0;
        
        if (e==expr) OK=0;
        //        fprintf(stderr,"\n  expr=\"%s\"=%g maxcol=%d n=%d OK=%d",expr,z,maxcol,n,OK);
        if (OK) xy[i]=z; }
      
      //      fprintf(stderr,"\n");

      *count+=1.0;

      if (isfinite(xy[0]) && isfinite(xy[1])) {
        if (errstyle)
          switch (errstyle) {
            case 1: myerrybar(xy[0],xy[1]-xy[2],xy[1]+xy[3]); break;
            case 2: myerrxbar(xy[0]-xy[2],xy[0]+xy[3],xy[1]); break;
            case 3: myerrxbar(xy[0]-xy[2],xy[0]+xy[2],xy[1]);
                    myerrybar(xy[0],xy[1]-xy[3],xy[1]+xy[3]); }
        else
          if (OK) {
            if (style>=0) mydraw(xy[0],xy[1]);
            else mysmartdot(xy[0],xy[1]); }
          else
            myup(); } }
    if (plotopt.slowdraw>0) {
      usleep(plotopt.slowdraw);
      XFlush(display); }
  }

  fclose(f);

  return 0;
}

int drawfile(char *fn,char *colx,char *coly,int maxcol,int color,int style)
{
  return drawerrbar(fn,colx,coly,NULL,NULL,maxcol,color,style,0);
}

int minmaxfile(char *fn,char *colx,char *coly,int maxcol,
               double *X0,double *X1,double *Y0,double *Y1)
{
  FILE *f;
  char line[LINELEN];
  int i,OK;
  double xy[2];

  *count=0;

  if (fn==NULL || (f=fopen(fn,"rt"))==NULL) {
    fprintf(stderr,"ploterr: `%s\': No such file\n",fn);
    return -1; }

  while (fgets(line,LINELEN,f)) {
    char *e,*tok;
    int n;

    /* not KEY -> blank line (even if comment) */
    if (plotopt.key && !strstr(line,plotopt.key)) tok=NULL;
    else tok=strtok(line,SEP);

    if (!tok || strchr("!#",tok[0])) {
      /* comment or blank -> pen up */
      myup();
      /* blank line resets counter n */
      if (!tok) *count=0; }

    else {
      /* _Id.head->val = count, _Id.head->next->val = 1st number, etc. */
      n=1;
      id=_Id.head->next;

      do {
        id->val=strtod(tok,&e);
        if (e==tok) id->val=FP_NAN;
        id=id->next;
        n++;
        tok=strtok(NULL,SEP);
      } while (tok && id);


      OK=1;
      loop (i,0,2) {
        double z;
        char *expr=i?coly:colx;

        e=expr; /* err */
        if (maxcol<n) z=Calc(expr,&e); else OK=0;
        if (e==expr) OK=0;
        if (OK) xy[i]=z; }

      *count+=1.0;

      if (OK) {
        if (isfinite(xy[0])) { Min(*X0,xy[0]) Max(*X1,xy[0]) }
        if (isfinite(xy[1])) { Min(*Y0,xy[1]) Max(*Y1,xy[1]) } } } }

  fclose(f);

  return 0;
}
