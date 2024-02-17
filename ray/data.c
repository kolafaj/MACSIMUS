/***********************************************************************
 * $Author: markv $
 * $Revision: 1.3 $
 * $Date: 88/09/17 01:23:10 $
 * $Log:	data.c,v $
 * Revision 1.3  88/09/17  01:23:10  markv
 * Added definitions for antialias related variables
 * 
 * Revision 1.2  88/09/12  12:52:03  markv
 * Added counter for shadow cache hits.  Made the max level of 
 * recursion equal to the define MAXLEVEL in defs.h
 * 
 * Revision 1.1  88/09/11  11:02:13  markv
 * Initial revision
 * 
 ***********************************************************************/
#include "defs.h"

int 		yylinecount ;
Viewpoint 	Eye ;
int 		Xresolution = 512 ;
int 		Yresolution = 512 ;
Light		Lights[MAXLIGHTS] ;
int		nLights = 0 ;
Vec		BackgroundColor ;
Surface		* CurrentSurface ;
Object		* Prims[MAXPRIMS] ;
int		nPrims = 0 ;
Flt		rayeps = 1e-6 ;
char *		Progname ;
Object		* Root ;

Flt		minweight = 0.01 ;
int		maxlevel = MAXLEVEL ;
int		nRays = 0 ;
int		nShadows = 0 ;
int		nReflected = 0 ;
int		nRefracted = 0 ;
int		maxQueueSize = 0 ;
int		totalQueues = 0 ;
int		totalQueueResets = 0 ;
int		tickflag = 0 ;
int		filtflag = 0 ;
int             yresolutionflag = 0 ;
int             xresolutionflag = 0 ;
Flt             relresolutionflag = 0 ;
Flt             Gamma = 1;

int		jitterflag = 1 ; /* antialiasing: speed/quality compromise */
int		maxsamples = -9 ;

int		nChecked = 0 ;
int		nEnqueued = 0 ;
int		nShadowCacheHits = 0 ;

/* JK */
Flt             aspect=1;
Flt             aspectflag=0;
Flt             scale=1;
Flt             bgmode=2;
Flt             lightfactor=1;
Flt             isotropiclight=0.1;
Flt             normaldirlight=0.2;
Flt             diffuselight=0;
Flt             threshold=0.02;
Flt             fog=0;
Flt             fogthick=0;

Vec     Slab[] = {
	{1.0, 0.0, 0.0},
	{0.0, 1.0, 0.0},
	{0.0, 0.0, 1.0},
	{1.0, 1.0, 0.0},
	{1.0, 0.0, 1.0},
	{0.0, 1.0, 1.0}
} ;
