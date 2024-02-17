#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "draw.h"
#include "plotfile.h"
#include "mydraw.h"

#define Max(A,B) { if (B>A) A=B; }
#define Min(A,B) { if (B<A) A=B; }

static double myatof(char *c)
{
if (strchr("0123456789+-.",c[0])) return atof(c);
else return -4.1e33;
}

int drawfile(char *fn,int colx,int coly,int color,int style)
{
FILE *f;
char line[1024];
int l=0;

if (!fn || fn[0]==0) f=stdin;
else f=fopen(fn,"rt");

if (!f) {
  fprintf(stderr,"\"%s\" does not exist\n",fn);
  return -1; }

mysetcolor(color);
if (style>1) mysetlinestyle(0,0,2*style+1);
else if (style>0) mysetlinestyle(style,0,style);
else pointsize=(-1-style)*0.798;

myup();

while (fgets(line,1024,f)) {
  char *c;

  for (c=line; *c!=0 && *c<=' '; c++);
  if (c[0]=='#' || c[0]=='!' || c[0]==0) {
    myup();
    if (c[0]==0) l=0; }
  else {
    char *t=strtok(line," ,\t\n");
    int col=1;
    double x=-3.1e33,y=-3.1e33;
  
    if (colx==0) x=l;
    if (coly==0) y=l;
    l++;

    while (t) {
      if (col==colx) x=myatof(t);
      if (col==coly) y=myatof(t);
      t=strtok(NULL," ,\t\n");
      col++; }
  
    if (x<=-4e33 || y<=-4e33);
    else if (x<=-3e33 || y<=-3e33)
      myup();
    else
      if (style>=0) mydraw(x,y);
      else mysmartdot(x,y); } }

fclose(f);
return 0;
}

int minmaxfile(char *fn,int cont,int colx,int coly,
	       double *X0,double *X1,double *Y0,double *Y1)
{
FILE *f;
char line[1024];
int l=0;

if (cont==0) {
  *X0=*Y0=3e33;
  *X1=*Y1=-3e33; }
  
if (fn==NULL || (f=fopen(fn,"rt"))==NULL) {
  fprintf(stderr,"\"%s\" does not exist\n",fn);
  return -1; }

while (fgets(line,1024,f))
  if (!strchr("!#",line[0])) {
    char *t=strtok(line," ,\t\n");
    int col=1;
    double x=-3.1e33,y=-3.1e33;

    if (colx==0) x=l;
    if (coly==0) y=l;
    l++;

    while (t) {
      if (col==colx) x=myatof(t);
      if (col==coly) y=myatof(t);
      t=strtok(NULL," ,\t\n");
      col++; }
  
    if (x<=-4e33 || y<=-4e33);
    else if (x>-3e33 && y>-3e33) {
      Min(*X0,x) Max(*X1,x)
      Min(*Y0,y) Max(*Y1,y) }
 }

fclose(f);
return 0;
}
