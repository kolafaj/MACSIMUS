/***********************************************************************
 * $Author: markv $
 * $Revision: 1.3 $
 * $Date: 88/10/04 14:33:37 $
 * $Log:        screen.c,v $
 * Revision 1.3  88/10/04  14:33:37  markv
 * Revamped most of this code.  Added code for no antialiasing, jittering, 
 * filtered, and statistically optimized antialiasing.   Statistically
 * optimized anti aliasing still doesn't work properly.
 * 
 * Revision 1.2  88/09/12  10:00:54  markv
 * Fixed lingering (possible) bug with curbuf memory allocation.  
 * Size returned by malloc is now correct in both spots.
 * Thanks tim@ben.
 * 
 * Revision 1.1  88/09/11  11:00:43  markv
 * Initial revision
 * 
 ***********************************************************************/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "defs.h"
#include "pic.h"
#include "extern.h"

Flt
rndcos ()
{
/* JK/random number in (-1,1) */
static double seed;
long x;

if (seed==0) seed=time(&x);
seed=fmod(seed*78125e0,2147483647e0);
return (Flt)(seed/1073741823.5-1.0);
}

int
Screen(view, picfile, xres, yres)
 Viewpoint *view ;
 char *picfile ;
 int xres, yres ;
/* JK WARNING:
 * aspect and scale are globals (should be parameters),
 * defined in data.c, declared in extern.h
*/
{
/*      Pixel   *buf, *oldbuf, *curbuf, *tmp ;
        int red, green, blue ;
        Point   dir ;
        Ray     ray ;
        Color   color;
        int     x, y ; */

        Flt     frustrumwidth ;

        Pic     *pic, *PicOpen();
        Point   viewvec, leftvec, upvec ;

        pic = PicOpen(picfile, xres, yres) ;

        /*
         * calculate the "up" vector...
         */

        VecCopy(view -> view_up, upvec) ;
        VecNormalize(upvec) ;

        /*
         * and the "view" vector....
         */

        VecSub(view -> view_at, view-> view_from, viewvec);
        VecNormalize(viewvec);

        /*
         * And the "left" vector...
         */

        VecCross(view -> view_up, viewvec, leftvec);
        VecNormalize(leftvec);

        /*
         * Calculate the width of the view frustrum in world coordinates.
         * and then scale the left and up vectors appropriately.
         */

        frustrumwidth = (view -> view_dist) * ((Flt) tan(view -> view_angle)) ;
        VecScale(-frustrumwidth, upvec) ;
        VecScale(-frustrumwidth, leftvec) ;

        if (filtflag) {
                FilterScan( pic,
                        view -> view_from,
                        viewvec,
                        upvec,
                        leftvec, 
                        xres, yres) ;
        } else if (jitterflag) {
                JitterScan(   pic,
                        view -> view_from,
                        viewvec,
                        upvec,
                        leftvec, 
                        xres, yres) ;
        } else {
                Scan(   pic,
                        view -> view_from,
                        viewvec,
                        upvec,
                        leftvec, 
                        xres, yres) ;
        }

        PicClose(pic) ;
return 0;
}


int
Scan(pic, eye, viewvec, upvec, leftvec, xres, yres)
 Pic * pic ;
 Vec eye ;
 Vec viewvec ;
 Vec upvec ;
 Vec leftvec ;
 int xres, yres ;
{
        Ray ray ;
        int x, y ;
        Flt xlen, ylen ;
        Color color ;

        /*
         * Generate the image...
         */

        VecCopy(eye, ray.P) ;
        for (y = 0 ; y < yres ; y++) {
                for (x = 0 ; x < xres ; x++) {
                        xlen = ( ((Flt) (2 * x) / (Flt) xres) - 1.0 ) / scale ;
                        ylen = ( ((Flt) (2 * y) / (Flt) yres) - 1.0 )
                               * (aspect/scale) ;
                        VecComb(xlen, leftvec, ylen, upvec, ray.D) ;
                        VecAdd(ray.D, viewvec, ray.D) ;
                        VecNormalize(ray.D);
                        Trace(0, 1.0, &ray, color);
                        color[0] = min(1.0, color[0]) ;
                        color[1] = min(1.0, color[1]) ;
                        color[2] = min(1.0, color[2]) ;
                        PicWritePixel(pic, color) ;
                }
                if (tickflag) fprintf(stderr, "%4.1f%%\r",100.0*y/yres) ;
        }
        if (tickflag) fprintf(stderr, "\a\r") ;

return 0;
}


int
FilterScan(pic, eye, viewvec, upvec, leftvec, xres, yres)
 Pic * pic ;
 Vec eye ;
 Vec viewvec ;
 Vec upvec ;
 Vec leftvec ;
 int xres, yres ;
{
        Ray ray ;
        int x, y, i ;
        Flt xlen, ylen ;
        Color color ;
        Color * nbuf, *obuf, *tmp  ;    /* new and old buffers, resp.   */
        Color avg ;

        /*
         * allocate enough memory for the filter buffer
         */

        nbuf = (Color *) calloc ((xres + 1), sizeof(Color)) ;
        obuf = NULL ;

        VecCopy(eye, ray.P) ;

        for (y = 0 ; y <= yres ; y++) {
                for (x = 0 ; x <= xres ; x++) {
                        xlen = ( ((Flt) (2 * x) / (Flt) xres) - 1.0 ) / scale ;
                        ylen = ( ((Flt) (2 * y) / (Flt) yres) - 1.0 )
                               * (aspect/scale) ;
                        VecComb(xlen, leftvec, ylen, upvec, ray.D) ;
                        VecAdd(ray.D, viewvec, ray.D) ;
                        VecNormalize(ray.D);
                        Trace(0, 1.0, &ray, color);
                        color[0] = min(1.0, color[0]) ;
                        color[1] = min(1.0, color[1]) ;
                        color[2] = min(1.0, color[2]) ;
                        VecCopy(color, nbuf[x]) ;
                }
                if (obuf) {
                        for (i = 0 ; i < xres ; i++) {
                                VecAdd(nbuf[i], nbuf[i+1], avg) ;
                                VecAdd(obuf[i], avg, avg) ;
                                VecAdd(obuf[i+1], avg, avg) ;
                                VecScale(0.25, avg) ;
                                PicWritePixel(pic, avg) ;
                        }
                        tmp = obuf ;
                        obuf = nbuf ;
                        nbuf = tmp ;
                } else {
                        /* 
                         * first scan line, set it up wierdly...
                         */
                        obuf = nbuf ;
                        nbuf = (Color *) calloc ((xres + 1), sizeof(Color)) ;
                }
                if (tickflag) fprintf(stderr, "%4.1f%%\r",100.0*y/yres) ;
        }
        if (tickflag) fprintf(stderr, "\a\r") ;
return 0;
}

int alreadyodd(int x,int y)
{
return (x*y)&1;
}

int alreadyeven(int x,int y)
{
return (x%3==1) && (y%3==1);
}

int
JitterScan(pic, eye, viewvec, upvec, leftvec, xres, yres)
 Pic * pic ;
 Vec eye ; 
 Vec viewvec ; 
 Vec upvec ; 
 Vec leftvec ;
 int xres, yres ;
{
        Ray ray ;
        int x, y, i ;
        Flt xlen, ylen ;
        Flt xwidth, ywidth ;
        Color color, avg ;

        /* JK changes : */
        double smartcoeff;
        int nnest=0;
        int ms=abs(maxsamples);
        double thresholdms=threshold*ms;
        int square=(int)sqrt((double)ms);
        int smartsquare=0;
        int smartantialiasing=0;
        int (*already)(int,int);

        if (maxsamples==-1) {
          fprintf(stderr,"-j-1 is nonsense\n");
          exit(-1); }
        else if (square*square != ms) {
          if (maxsamples<0) {
            fprintf(stderr,"-j-# supported for #=square only\n");
            exit(-1); }
          square=0;
          printf("random jitter pattern %d points\n",ms); }
        else if (maxsamples<0) {
          if (diffuselight!=0) {
            if (square&1) {
              smartsquare=square*2-1;
              already=alreadyodd; }
            else {
              already=alreadyeven;
              smartsquare=square*3; }

            smartcoeff=(double)square/smartsquare;
            smartcoeff*=smartcoeff; 
            printf("square %dx%d antialiasing pattern, extendable to %dx%d with error threshold %g\n",
                   square,square,smartsquare,smartsquare,threshold); }
          else {
            smartantialiasing++;
            printf("square %dx%d antialiasing pattern for contrast>%g\n",
                   square,square,threshold); } }
        else
          if (square==1)
	    printf("no supersampling/antialiasing\n");
	  else
	    printf("square %dx%d supersampling/antialiasing pattern\n",square,square);

        xwidth = 1.0 / (Flt) xres / scale;
        ywidth = aspect / (Flt) yres / scale; /* JK */

        VecCopy(eye, ray.P) ;
if (smartantialiasing) 
#include "antialia.c"
else {
        for (y = 0 ; y < yres ; y++) {
                ylen = ( ((Flt) (2 * y) / (Flt) yres) - 1.0 )
                       * (aspect/scale) ;
                for (x = 0 ; x < xres ; x++) {
                        xlen = ( ((Flt) (2 * x) / (Flt) xres) - 1.0 ) / scale ;

                        avg[0] = avg[1] = avg[2] = 0.0 ;

/* JK */
if (square) 
  if (smartsquare) {
    /* square jitter pattern */
    int x,y;
    Vec fluct;
    double fl;
    static double w[3]={0.3,0.59,0.11};
    
    fluct[0] = fluct[1] = fluct[2] = 0.0 ;

    fl=0;
    for (x=0; x<square; x++)
      for (y=0; y<square; y++) {
        VecComb(xlen + (Flt)(2*x+1-square)/(Flt)(square) * xwidth, leftvec,
                ylen + (Flt)(2*y+1-square)/(Flt)(square) * ywidth, upvec, ray.D) ;
        VecAdd(ray.D, viewvec, ray.D) ;
        VecNormalize(ray.D);
        Trace(0, 1.0, &ray, color);
        fluct[0]+=color[0]*color[0];
        fluct[1]+=color[1]*color[1];
        fluct[2]+=color[2]*color[2];
        VecAdd(color, avg, avg) ; } 
    
    for (x=0; x<3; x++) fl+=sqrt(fluct[x]-avg[x]*avg[x]/ms)*w[x];
    if (fl>thresholdms) {
      nnest++;
      for (x=0; x<smartsquare; x++)
        for (y=0; y<smartsquare; y++) if (!already(x,y)) {
          VecComb(xlen + (Flt)(2*x+1-smartsquare)/(Flt)(smartsquare) * xwidth, leftvec,
                  ylen + (Flt)(2*y+1-smartsquare)/(Flt)(smartsquare) * ywidth, upvec, ray.D) ;
          VecAdd(ray.D, viewvec, ray.D) ;
          VecNormalize(ray.D);
          Trace(0, 1.0, &ray, color);
          VecAdd(color, avg, avg) ; }

      avg[0]*=smartcoeff;
      avg[1]*=smartcoeff;
      avg[2]*=smartcoeff;
    }
  }

  else {
    /* square jitter pattern */
    int x,y;
    for (x=0; x<square; x++)
      for (y=0; y<square; y++) {
        VecComb(xlen + (Flt)(2*x+1-square)/(Flt)(square) * xwidth, leftvec,
                ylen + (Flt)(2*y+1-square)/(Flt)(square) * ywidth, upvec, ray.D) ;
        VecAdd(ray.D, viewvec, ray.D) ;
        VecNormalize(ray.D);
        Trace(0, 1.0, &ray, color);
        VecAdd(color, avg, avg) ; } }

                        /* original: random jitter pattern */
else                    for (i = 0 ; i < ms ; i++) {
                                VecComb(xlen + rndcos() * xwidth, leftvec,
                                        ylen + rndcos() * ywidth, upvec, ray.D) ;
                                VecAdd(ray.D, viewvec, ray.D) ;
                                VecNormalize(ray.D);
                                Trace(0, 1.0, &ray, color);
                                VecAdd(color, avg, avg) ;
                        }

                        VecScale(1.0 / (Flt) ms, avg) ;
                        avg[0] = min(1.0, avg[0]) ;
                        avg[1] = min(1.0, avg[1]) ;
                        avg[2] = min(1.0, avg[2]) ;
                        PicWritePixel(pic, avg) ;
                }
                if (tickflag) fprintf(stderr, "%4.1f%%\r",100.0*y/yres) ;
        }
        if (tickflag) fprintf(stderr, "\a\r") ;

if (smartsquare) printf("%.2f%% smart squares\n",100.0*nnest/(xres*yres));
}

return 0;
}
