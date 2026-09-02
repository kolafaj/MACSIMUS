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

/* missing prototypes added 09/2026 (K&R -> ANSI C): */
int yylex(void);
int get_opt(int nargc,char **nargv,char *ostr);
void InitSlabs(void);
void ReadSceneFile(char *infilename);
void BuildBoundingSlabs(void);
void Screen(Viewpoint *view, char *picfile, int xres, int yres);
int compslabs(Object **a, Object **b);
void Shade(int level, Flt weight, Vec P, Vec N, Vec I, Isect *hit, Color col);
void SpecularDirection(Vec I, Vec N, Vec R);
int TransmissionDirection(Surface *m1,Surface *m2, Vec I,Vec N,Vec T);
int Intersect(Ray *ray, Isect *hit, Flt maxdist);
int LookupColorByName(char * name, Vec color);
Object *MakeCone(Vec basepoint, Flt baseradius, Vec apexpoint, Flt apexradius);
Object *MakeSphere(Vec pos, Flt radius);
Object *MakePoly(int npoints, Vec *points);
Object *MakeTri(Vec *point);
int PriorityQueueInsert(Flt key, Object *obj);
int PriorityQueueNull();
int PriorityQueueEmpty();
int PriorityQueueDelete(Flt * key, Object ** obj);
void Trace(int level, Flt weight, Ray *ray, Color color);


#if defined(__STDC__) || defined(__cplusplus)
int yyparse(void);
#else
int yyparse();
#endif

Flt VecNormalize(Vec vec);
Flt             rndcos(void) ;
extern	Vec		Slab[] ;
extern	ObjectProcs	NullProcs ;
extern 	Object *	Root ;
