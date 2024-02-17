/* permutations

  1. generate all n! permutations of integers ren[0]..ren[n-1]
     int n;
     int ren[n];
     int stat[n+1];
     initpermute(ren,stat,n);
     do {
       ...
       } while (permute(ren,stat,n));


  2. generate sets of permutations, see main() below

  speed is about 4e6 permutations / s on 333MHz Celeron

*/

#ifndef PERMUTETYPE /* test, otherwise to be #included */

#include "include.h"

#define PERMUTETYPE int 
/* PERMUTETYPE=void* is possible but some warnings... */

int *initpermutes(int nspec,int *stoich);
int permutes(PERMUTETYPE *ren);


main()
{
 int *ren;
 int stoich[3]={4,3,2};
 int i,n=0;
 
 ren=initpermutes(3,stoich);
 /* see also initpermutes0 without allocating ren */
 do {
   n++;
   loop (i,0,2+3+4) printf("%d",ren[i]); _n
   } while (permutes(ren));
 put(n)
}

#endif /* test */

void initpermute(PERMUTETYPE *ren,int *stat,int n)
{
 int i;
 loop (i,0,n) ren[i]=i,stat[i+1]=0;
}

int permute(PERMUTETYPE *ren,int *stat,int n)
{
 int i,j;
 PERMUTETYPE aux;

 loopto (i,2,n) {
   aux=ren[0];
   loop (j,1,i) ren[j-1]=ren[j];
   ren[i-1]=aux;

   if (++stat[i]==i) stat[i]=0;
   else return 1; }

 return 0;
}

static int *ren_stat,*ren_stat0,ren_nm,*ren_stoich,ren_nspec;

int *initpermutes(int nspec,int *stoich)
{
 int i;
 PERMUTETYPE *ren;

 ren_nm=0;
 if (ren_stat0) free(ren_stat0);
 ren_stoich=stoich;
 ren_nspec=nspec;

 loop (i,0,nspec) ren_nm+=stoich[i];
 allocarrayzero(ren_stat0,ren_nm);
 ren_stat=ren_stat0-1;
 alloc(ren,ren_nm*sizeof(ren[0]));

 loop (i,0,ren_nm) ren[i]=(PERMUTETYPE)i;
 return ren;
}

void initpermutes0(int nm,PERMUTETYPE *ren,int nspec,int *stoich)
/* as above, ren not allocated */
{
 int i,check=0;

 ren_nm=nm;
 ren_stoich=stoich;
 ren_nspec=nspec;

 loop (i,0,nspec) check+=stoich[i];
 if (check!=nm) Error("initpermutes0: check!=nm");
 if (ren_stat0) free(ren_stat0);
 allocarrayzero(ren_stat0,ren_nm);
 ren_stat=ren_stat0-1;

 loop (i,0,ren_nm) ren[i]=(PERMUTETYPE)i;
}

int permutes(PERMUTETYPE *ren)
{
 int off=0;
 int i;

 loop (i,0,ren_nspec) {
   if (permute(ren+off,ren_stat+off,ren_stoich[i])) return 1;
   off+=ren_stoich[i]; }

 free(ren_stat0);

 return 0;
}
