/* \make psstring 
 */
#include "ground.h"
#include "draw.h"
/* EXTSTRING required in mydraw.c */
#include "mydraw.h"

int utf8(char *u,char *c,char *x,char *a,char *d)
{
  if (u[0]==c[0] && u[1]==c[1]) {
    *a=x[0],*d=x[1];
    return 1; }
  else
    return 0;
} 

int main(int narg,char **arg)
{
  int i, align=-1, isutf8=0;
  double x=0,y=0;
  char *c,*asc,*dia,*a,*d;

  if (narg<2) {
    fprintf(stderr,"Call by:\n\
  psstring [OPTIONS] STRING [[OPTIONS] STRING ...]\n\
OPTIONS\n\
  position: -xX, -yY\n\
  alignment: -l[eft] (default), -r[ight], -c[enter]\n\
  font size: -sSIZE (default=14)\n\
  font: -fFONT (default=Helvetica)\n\
STRING\n\
  Symbol (Greek etc.): \\a\n\
  Index/exponent: A_i, x^2\n\
Czech accented characters (in UTF-8):\n\
  see tex/makra.eps; cannot precede sub/super-script, bad iIdt, no umlaut...\n\
See also:\n\
  plot textppm tex2pbm tex2ppm tex/makra.eps\n");
   exit(0); }

  my_style.ps=stdout;
  my_style.font="Helvetica";

  loop (i,1,narg) 
    if (arg[i][0]=='-') switch (arg[i][1]) {
      case 'c': align=0; break;
      case 'r': align=1; break;
      case 'l': align=-1; break;
      case 's': my_style.fontsize=atof(arg[i]+2); break;
      case 'f': my_style.font=arg[i]+2; break;
      case 'x': x=atof(arg[i]+2); break;
      case 'y': y=atof(arg[i]+2); break;
      default: fprintf(stderr,"%s is bad option\n",arg[i]); }
    else {
      isutf8=0;
      for (c=arg[i]; *c; c++) if (*c&128) isutf8=1;
      if (isutf8) {
        asc=strdup(arg[i]);
        dia=strdup(arg[i]);
        for (c=arg[i],a=asc,d=dia; *c; c++,a++,d++) {
          *a=*c,*d=' ';
          if (utf8("ě",c,"ev",a,d)) c++; 
          if (utf8("š",c,"sv",a,d)) c++; 
          if (utf8("č",c,"cv",a,d)) c++; 
          if (utf8("ř",c,"rv",a,d)) c++; 
          if (utf8("ž",c,"zv",a,d)) c++; 
          if (utf8("ý",c,"y,",a,d)) c++; 
          if (utf8("á",c,"a,",a,d)) c++; 
          if (utf8("í",c,"ii",a,d)) c++; 
          if (utf8("ů",c,"uo",a,d)) c++; 
          if (utf8("ú",c,"u,",a,d)) c++; 
          if (utf8("ó",c,"o,",a,d)) c++; 
          if (utf8("ň",c,"nv",a,d)) c++; 
          if (utf8("ď",c,"dV",a,d)) c++; 
          if (utf8("ť",c,"tV",a,d)) c++; 
          
          if (utf8("Ě",c,"EV",a,d)) c++; 
          if (utf8("Š",c,"SV",a,d)) c++; 
          if (utf8("Č",c,"CV",a,d)) c++; 
          if (utf8("Ř",c,"RV",a,d)) c++; 
          if (utf8("Ž",c,"ZV",a,d)) c++; 
          if (utf8("Ý",c,"Y,",a,d)) c++; 
          if (utf8("Á",c,"A/",a,d)) c++; 
          if (utf8("Í",c,"I/",a,d)) c++; 
          if (utf8("Ů",c,"UO",a,d)) c++; 
          if (utf8("Ú",c,"U/",a,d)) c++; 
          if (utf8("Ó",c,"O/",a,d)) c++; 
          if (utf8("Ň",c,"NV",a,d)) c++; 
          if (utf8("Ď",c,"DV",a,d)) c++; 
          if (utf8("Ť",c,"TV",a,d)) c++; 
        }
        *a=*d=0; }
      else 
        asc=arg[i];
      
      printf("%% %s\n",arg[i]);
      printf("/%s findfont %.2f scalefont setfont\n",
             my_style.font,my_style.fontsize);
      printf("%.2f %.2f moveto\n",x,y);
      psstring(asc,align);
      putchar('\n'); 
      if (isutf8) {
        for (a=asc,d=dia; *a; a++,d++)  {
          char aa=*a;
          *a=0;
          if (*d!=' ') printf("%.2f (%s) %.2f %.2f %s\n",
            my_style.fontsize,asc,x,y,
            *d==','?"carka":
            *d=='i'?"icarka":
            *d=='/'?"CARKA":
            *d=='v'?"hacek":
            *d=='V'?"HACEK":
            *d=='o'?"krouzek":
            *d=='O'?"KROUZEK":"UNKNOWN" );
          *a=aa; } }

      y-=my_style.fontsize*1.5; }

  return 0;
}
