/***********************************************************************
 * $Log:	defs.h,v $
 * Revision 1.4  88/10/06  12:54:57  markv
 * I forgot (of course) that gcc runs on non-sun machines.  I have
 * defined the FASTPRIMS so that if you are running a sun, you get them
 * otherwise, you don't.
 * 
 * Revision 1.3  88/10/04  14:34:49  markv
 * Added changes to allow for optimized CSG (courtesy of Glassner and
 * Bronsvoort.
 * 
 * Revision 1.2  88/09/12  13:03:28  markv
 * Added space in the Light structure for caching shadow objects.
 * 
 * Revision 1.1  88/09/11  11:00:49  markv
 * Initial revision
 * 
 * Revision 1.1  88/09/09  11:59:55  markv
 * Initial revision
 * 
 ***********************************************************************/

#include "config.h"

typedef double Flt ;
typedef Flt Vec[3] ;
typedef Vec Point ;
typedef Vec Color ;

/*----------------------------------------------------------------------*/

#ifndef DUMB_CPP 

#define MakeVector(x, y, z, v)		(v)[0]=(x),(v)[1]=(y),(v)[2]=(z)
#define VecScale(S,a)	(a)[0] *= S ; (a)[1] *= S ; (a)[2] *= S
#define VecNegate(a)	(a)[0]=0-(a)[0];\
			(a)[1]=0-(a)[1];\
			(a)[2]=0-(a)[2];
#define VecDot(a,b)	((a)[0]*(b)[0]+(a)[1]*(b)[1]+(a)[2]*(b)[2])
#define VecLen(a)	(sqrt(VecDot(a,a)))
#define VecCopy(a,b)	 (b)[0]=(a)[0];(b)[1]=(a)[1];(b)[2]=(a)[2];
#define VecAdd(a,b,c)	 (c)[0]=(a)[0]+(b)[0];\
			 (c)[1]=(a)[1]+(b)[1];\
			 (c)[2]=(a)[2]+(b)[2]
#define VecSub(a,b,c)	 (c)[0]=(a)[0]-(b)[0];\
			 (c)[1]=(a)[1]-(b)[1];\
			 (c)[2]=(a)[2]-(b)[2]
#define VecComb(A,a,B,b,c)	(c)[0]=(A)*(a)[0]+(B)*(b)[0];\
				(c)[1]=(A)*(a)[1]+(B)*(b)[1];\
			 	(c)[2]=(A)*(a)[2]+(B)*(b)[2]
#define VecAddS(A,a,b,c)	 (c)[0]=(A)*(a)[0]+(b)[0];\
				 (c)[1]=(A)*(a)[1]+(b)[1];\
				 (c)[2]=(A)*(a)[2]+(b)[2]
#define VecCross(a,b,c)	 (c)[0]=(a)[1]*(b)[2]-(a)[2]*(b)[1];\
			 (c)[1]=(a)[2]*(b)[0]-(a)[0]*(b)[2];\
			 (c)[2]=(a)[0]*(b)[1]-(a)[1]*(b)[0]
#define VecPrint(msg,v)		printf("%s %g %g %g\n", msg,\
					(v)[0],(v)[1],(v)[2])
#define VecZero(v)	(v)[0]=0.0;(v)[1]=0.0;v[2]=0.0

#endif /* not DUMB_CPP */

/*----------------------------------------------------------------------*/

typedef struct Ray {
        Point P ; /* point of origin */
        Point D ; /* (normalized) direction */
} Ray ;

#define RayPoint(ray,t,point)	VecAddS(t,(ray)->D,(ray)->P,point)

#define max(a,b) 	((a)>(b)?(a):(b))
#define min(a,b) 	((a)<(b)?(a):(b))

/*----------------------------------------------------------------------*/

typedef struct t_surface {
	Color	surf_color ;
	Flt	surf_kd ;
	Flt	surf_ks ;
	Flt	surf_shine ;
	Flt 	surf_kt ;
	Flt	surf_ior ;
} Surface ;

typedef struct t_light {
	Vec	light_pos ;
	Vec	light_pos0 ;
	Flt	light_brightness ;
	struct t_object * light_obj_cache[MAXLEVEL] ;
} Light ;

typedef struct t_viewpoint {
	Vec	view_from ;
	Vec	view_at ;
	Vec	view_up ;
	Flt	view_angle ;
	Flt	view_dist ;
} Viewpoint ;

typedef struct t_object {
	unsigned short 	o_type ;
	Flt	o_dmin[NSLABS] ;
	Flt	o_dmax[NSLABS] ;
	struct t_objectprocs {
		int 	(*print) () ;
		int 	(*intersect) () ;
		int 	(*normal) () ;
	} * o_procs ;
	struct t_surface 	* o_surf ;
	struct t_object		* o_parent ;
	unsigned short		* o_inside ;
	void	* o_data ;
} Object ;

typedef struct t_compositedata {
	unsigned short 	c_size ;
	Object *	c_object[BUNCHINGFACTOR] ;
} CompositeData ;

typedef struct t_objectprocs ObjectProcs ;

typedef struct t_isect {
	Flt 		isect_t ;
	int 		isect_enter ;
	Vec		isect_normal ;
	Object 		* isect_prim ;
	Surface 	* isect_surf ;
} Isect ;

typedef struct t_pixel {
	unsigned char r, g, b, q ;
} Pixel ;

#ifndef PI
#define PI 		(3.14159265358979323844)
#endif /* PI */

#define degtorad(x)	(((Flt)(x))*PI/180.0)

#define		T_COMPOSITE	(0)
#define		T_SPHERE	(1)
#define		T_POLY		(2)
#define		T_CONE		(3)
#define		T_TRI		(4)

#ifdef __GNUC__
#define INLINE	inline

#if defined(FAST_MATH_PRIMS) && defined(sun)
#define sqrt(x)\
	({ double __value, __arg = (x) ;\
	 asm("fsqrtx %1,%0": "=f" (__value): "f" (__arg)) ;\
	 __value ;})
#endif /* FAST_MATH_PRIMS */

#else /* __GNUC__ */

#define INLINE	/* inline not supported  */

#endif /* __GNUC__ */
