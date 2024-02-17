/***********************************************************************

 % changes by J.Kolafa 9/1997:
 %
 % added flags:
 %   -x x-resolution (this corresponds to `angle' parameter in the nff file)
 %   -y y-resolution
 %   -a pixel aspect ratio
 %   -s scale (effectively changes angle); scale>1 enlarges the central
 %            part of the picture and cuts borders
 %   example: to make picture for VGA screen so that `angle' is the
 %            angle of the y-size of the picture, use the following flags:
 %            -x 320 -y 200 -a 1.2 -s 0.75
 % -t prints progress in %, not dots
 % -b FILE will use FILE (in ppm format) as background
 %    BUG: works only for from = (0,0,value), as generated e.g. by show
 % -B BGMODE: 0=use background color outside picture (=when reflected)
 %            1,2=tile the background image
 %            -1,-2=mirror-tile the background image (mirror to be continuous)
 %            1,-1=for front, use the background color
 %            2,-2=for front, calc. the background color as average over
 %                 the picture
 % -X -Y = scaling for the background picture
 % -X0 -Y0 = background pixel --> pixel
 % -R shift of the BG image right, in its x-size
 % -U shift of the BG image up, in its y-size
 %    (e.g., 0.3--0.6 for the sky in the up half of the picture)
 % rnd() changed into rndcos() and made portable
 % jittering/antialiasing improved: if -j# is a square, a square pattern 
 %   is used, otherwise random
 % cone.c: bug fixed
 % small bugs and non-portabilities fixed
 % names longer than DOS standard shortened
 %
 % 1999 more improvements:
 %   -d diffuse Gaussian light, arg is Gaussian sigma of random light shift
 %      in real length units: should use jittering, e.g. -j225
 %   -l light multiplication factor (default=1, useful for more lights)
 %   -A -N ambient light added
 %   -v added (launches your favorite viewer)
 %   jittering/antialiasing improved
 % 2003 fog added, option -f/-F; old -f renamed to -u
 %   -f fog from (z coordinate)
 %   -F fog to (z coordinate)
 %      fog attenuates exponentially to the background color from -f
 %      there is 100% visibility in front of -f and 1/e at -F
 %   -u 2x2 blur (old option -f)

 * $Author: markv $
 * $Revision: 1.6 $
 * $Date: 88/10/31 14:33:44 $
 * $Log:        main.c,v $
 * Revision 1.6  88/10/31  14:33:44  markv
 * Removed non-functioning antialiasing...
 * 
 * Revision 1.5  88/10/04  14:29:16  markv
 * Added flags having to do with antialiasing, both statistically
 * optimized and jittered.
 * 
 * Revision 1.4  88/09/13  16:08:21  markv
 * Moved InitSlabs, must occur before ReadSceneFile
 * 
 * Revision 1.3  88/09/13  16:05:44  markv
 * Usage message now prints -t flag as well.
 * 
 * Revision 1.2  88/09/12  12:53:11  markv
 * Added printout of new statistic on shadow cache hits.
 * 
 * Revision 1.1  88/09/11  11:02:23  markv
 * Initial revision
 * 
 ***********************************************************************/

#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include "defs.h"
#include "extern.h"
#include "pic.h"

extern double atof(const char *);
char **mainargv;
int mainargc;

extern FILE *yyin, *yyout ;

int
main(argc, argv)
 int argc ;
 char * argv[] ;
{
        extern int opt_ind ;
        extern char *opt_arg ;
        int c ;
        char * infilename = NULL ;
        char * outfilename = "out.ppm" ;
        char * bgfilename = NULL;
        struct tms pbuf, tbuf ;
        Flt pixelaspect=1,bgxscale=1,bgyscale=1,bgxshift=0,bgyshift=0;
        char *viewer=NULL;

        yyin=stdin, yyout=stdout ;

        mainargc=argc;
        mainargv=argv;
        if ((Progname = rindex(argv[0], '/')) == NULL) 
                Progname = argv[0] ;
        else 
                Progname ++ ;

        /* was: "r:a:o:i:tfj:" */
        while ((c = get_opt(argc, argv,
          "r:a:o:i:b:htuj:x:y:s:S:f:F:B:X:Y:n:U:R:l:N:I:d:v:c:G:")) != EOF) {
                switch (c) {
                /* JK: */
                case 'f':
                        fog=atof(opt_arg);
			break;
                case 'F':
                        fogthick=atof(opt_arg);
			break;
                case 'n':
                        {
                        int l=strlen(opt_arg);
                        infilename=(char*)malloc(l+5);
                        outfilename=(char*)malloc(l+5);
                        sprintf(infilename,"%s.nff",opt_arg);
                        sprintf(outfilename,"%s.ppm",opt_arg);
                        }
                        break;
                case 'd': 
                        diffuselight=atof(opt_arg);                        
                        break;
                case 'l':
                        lightfactor=atof(opt_arg);
                        break;
                case 'c':
                        threshold=atof(opt_arg);
                        break;
                case 'I':
                        isotropiclight=atof(opt_arg);
                        break;
                case 'N':
                        normaldirlight=atof(opt_arg);
                        break;
                case 'a':
                        aspectflag=atof(opt_arg);
                        break;
                case 'b':
                        bgfilename = opt_arg;
                        break;
                case 'B':
                        bgmode=atoi(opt_arg);
                        break;
                case 'X':
                        bgxscale=atof(opt_arg);
                        break;
                case 'Y':
                        bgyscale=atof(opt_arg);
                        break;
                case 'U':
                        bgyshift=atof(opt_arg);
                        break;
                case 'R':
                        bgxshift=atof(opt_arg);
                        break;
                case 'x':
                        xresolutionflag = atoi(opt_arg) ;
                        break ;
                case 'y':
                        yresolutionflag = atoi(opt_arg) ;
                        break ;
                case 's':
                        scale =  atof(opt_arg);
                        break;

                case 'r':
                        xresolutionflag = yresolutionflag = atoi(opt_arg) ;
                        break ;
                case 'S':
                        relresolutionflag = atof(opt_arg) ;
                        break ;
                case 'u': /* old -f */
                        filtflag = 1 ;
                        break ;
                case 'j':
                        jitterflag = 1 ;
                        maxsamples = atoi(opt_arg) ;
#if 0
                        if (maxsamples <= 0) {
                                fprintf(stderr, "%s: samples must be > 0\n",
                                        Progname) ;
                                exit(1) ;
                        }
#endif
                        break ;
                case 'o':
                        outfilename = opt_arg ;
                        break ;
                case 'i':
                        infilename = opt_arg ;
                        break ;
                case 't':
                        tickflag = 1 ;
                        break ;
                case 'v':
                        viewer=strdup(opt_arg);
                        break;
		case 'G':
		        Gamma=atof(opt_arg);
		        break;
                case '?':
                case 'h':
                        printf("Reasonably Intelligent Raytracer by Mark VandeWettering & Jiri Kolafa\n\
Usage: %s { -Option Argument | -OptionArgument } ...\n\
-h      : this help\n\
-i FILE : input scene file (Neutral File Format)\n\
-o FILE : output file (P6 Portable Pixel Map)\n\
-n FILE : input=FILE.nff, output=FILE.ppm\n\
-t      : show progress in %%\n\
-u      : cheap and fast antialiasing by 2x2 blur\n\
-v VIEW : show picture using viewer VIEW\n\
-j #    : #=0: no antialiasing nor jittering\n\
          #>0, # not square: random jittering # samples/pixel\n\
          #>0, #=N^2: antialiasing by N*N supersamples/pixel\n\
          #<0, #=-N^2: supersample only if contrast>c [default=-j-9]\n\
          #<0, #=-N^2 -d#: if error>c, supersample N-->2N+1 (even) or 3N (odd)\n\
-d #    : area light size (Gaussian jittering), use -j+-LARGE\n\
-c #    : c, threshold for smart supersampling (-j-#) [default=0.02]\n\
-x #    : x size in pixels [default=command `resolution' in NFF file]\n\
-y #    : y size\n\
-r #    : x and y size (resolution)\n\
-s #    : scale view angle (zoom) by # [1]\n\
-S #    : scale x and y size by # [1]\n\
-a #    : pixel aspect ratio (y/x) [1]\n\
-f #    : fog from z-coordinate (no fog in front of this z) [0]\n\
-F #    : fog thickness for attenuation to 1/e (to background color) [0=off]\n\
-b FILE : background image (PPM file), instead of color (command `b') (no fog!)\n\
-X #    : scale background image # times horizontally [1]\n\
-Y #    : scale background image # times vertically [1]\n\
-U #    : move background image up by #*height [0]\n\
-R #    : move background image right by #*width [0]\n\
-B #    : #>0: tile background image\n\
          #<0: chessboard mirror tile (make nonperiodic images continuous)\n\
          +-1: use command `b' for front background\n\
          +-2: calculate background color for front as image average [2]\n\
-l #    : light scaling factor [1] (effective brightness adjusted to -I -N)\n\
-I #    : ambient isotropic light [0.1]\n\
-N #    : ambient light proportional to cos angle(normal,ray) [0.2]\n\
-G #    : gamma for output (does not invert background!) [1]\n\
", Progname) ;
                        exit(-1);
                }
        }

        InitSlabs() ;                   /* Normalizes slab normals...   */
        ReadSceneFile(infilename) ;
        if (xresolutionflag > 0) Xresolution=xresolutionflag ;
        if (yresolutionflag > 0) Yresolution=yresolutionflag ;
        if (relresolutionflag > 0) {
          Yresolution=(Yresolution*relresolutionflag)+0.5;
          Xresolution=(Xresolution*relresolutionflag)+0.5;
          printf("resolution changed into %dx%d\n",Xresolution,Yresolution); }

        /* JK: -a is pixel aspect ratio */
        if (aspectflag!=0) aspect=pixelaspect=aspectflag;
        printf("pixel aspect (ysize/xsize) = %.5f",aspect);
        aspect *= (Flt)Yresolution/(Flt)Xresolution;
        printf(", picture aspect = %.5f\n",aspect);
        if (scale!=1) printf("scale=%.5f\n",scale);

        BuildBoundingSlabs() ;
        times(&pbuf) ;

        if (bgfilename) bg=ReadPic(bgfilename);
        if (bg) {
          Flt xxx=2*tan(Eye.view_angle);
          bg->xrange=xxx;
          bg->yrange=xxx*bg->y/bg->x;
          bg->xrange*=bgxscale;
          bg->yrange*=bgyscale;
          bg->yshift=bgyshift;
          bg->xshift=bgxshift;
          if (bgxscale==0) bg->xrange=xxx*bg->x/Xresolution;
          if (bgyscale==0) bg->yrange=xxx*bg->y/Xresolution;
          printf("bg range=%f %f\n",bg->xrange,bg->yrange); }

        Screen(&Eye, outfilename, Xresolution, Yresolution) ;
        times(&tbuf) ;
        tbuf.tms_utime -= pbuf.tms_utime ;
        tbuf.tms_stime -= pbuf.tms_stime ;
        PrintStatistics(&pbuf, &tbuf) ;
if (viewer) {
  char *call=(char*)malloc(strlen(viewer)+strlen(outfilename)+5);
  sprintf(call,"%s %s &",viewer,outfilename);
  fprintf(stderr,"viewer call status=%d\n",system(call));
  }
return 0;
}

int
PrintStatistics(pbuf, tbuf)
 struct tms *pbuf, *tbuf ;
{

        printf("preprocess time (user code)     %-6d seconds\n",
                        (int)(pbuf -> tms_utime / 60)) ;
        printf("preprocess time (system code)   %-6d seconds\n",
                        (int)(pbuf -> tms_stime / 60)) ;
        printf("tracing time (user code)        %-6d seconds\n",
                        (int)(tbuf -> tms_utime / 60 )) ;
        printf("tracing time (system code)      %-6d seconds\n",
                        (int)(tbuf -> tms_stime / 60 )) ;

        printf("number of rays cast:       %-6d\n", nRays);
        printf("number of shadow rays:     %-6d\n", nShadows);
        printf("number of reflected rays:  %-6d\n", nReflected);
        printf("number of refracted rays:  %-6d\n", nRefracted);

        printf("number of queue inserts:   %-6d\n", totalQueues) ;
        printf("number of queue resets:    %-6d\n", totalQueueResets) ;
        printf("avg number of queues/ray:  %-6g\n", (Flt) totalQueues /
                                                (Flt) totalQueueResets) ;
        printf("max queue size:            %-6d\n", maxQueueSize) ;
        printf("number of bound checks:    %-6d\n", nChecked) ;
        printf("number of bound queued:    %-6d\n", nEnqueued) ;
        printf("number of shadow hits:     %-6d\n", nShadowCacheHits) ;

return 0;
}
