/***********************************************************************
 * $Log:	extern.h,v $
 * Revision 1.3  88/09/17  01:22:05  markv
 * Added definitions for new antialiasing variables, plus
 * function definitions for functions which return the chi-squared
 * values.
 * 
 * Revision 1.2  88/09/12  13:11:13  markv
 * Added extern definition for nShadowCacheHits
 * 
 * Revision 1.1  88/09/11  11:00:49  markv
 * Initial revision
 * 
 ***********************************************************************/
extern  int 		yylinecount ;
extern	Viewpoint 	Eye ;
extern	int 		Xresolution ;
extern	int 		Yresolution ;
extern	Light		Lights[] ;
extern	int		nLights ;
extern	Vec		BackgroundColor ;
extern	Surface		* CurrentSurface ;
extern	Object		* Prims[] ;
extern	int		nPrims ;
extern 	Flt		rayeps ;
extern	char *		Progname ;
extern 	int		maxQueueSize ;
extern 	int		totalQueues ;
extern	int		totalQueueResets ;
extern 	int		tickflag ;
extern 	int		filtflag ;
extern 	int		jitterflag ;
extern  int             xresolutionflag ;
extern  int             yresolutionflag ;
extern  Flt             relresolutionflag ;
extern  Flt             Gamma;
extern  int		nChecked ;
extern 	int		nEnqueued ;
extern  int             nShadowCacheHits ;

/* JK */
extern  Flt             aspect;
extern  Flt             aspectflag;
extern  Flt             scale;
extern  Flt             bgmode;
extern  Flt             lightfactor;
extern  Flt             isotropiclight;
extern  Flt             normaldirlight;
extern  Flt             diffuselight;
extern  Flt             threshold;
extern  Flt             fog,fogthick;
extern  char **mainargv;
extern  int mainargc;

extern 	Flt		minweight ;
extern 	int		maxlevel ;
extern	int		maxsamples ;
extern	Flt		variance ;
extern  Flt		maxerror ;
extern 	int		nRays ;
extern	int		nShadows ;
extern	int		nReflected ;
extern 	int		nRefracted ;

/*.....char *          malloc() ;*/
/*.....char *          calloc() ;*/
char *		rindex() ;

extern	Object *	MakeCone() ;
extern	Object *	MakeSphere() ;
extern	Object *	MakePatch() ;
extern	Object *	MakePoly() ;
extern	Object * 	MakeTri() ;

extern 	Flt		VecNormalize() ;
extern  Flt             rndcos() ;
extern	Flt		critchisq() ;
extern	Flt		pochisq() ;
extern	Vec		Slab[] ;
extern	ObjectProcs	NullProcs ;
extern 	Object *	Root ;

#ifdef DUMB_CPP

extern Flt	VecDot() ;
extern Flt 	VecLen() ;

#endif /* DUMB_CPP */
