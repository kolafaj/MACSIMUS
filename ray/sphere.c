/***********************************************************************
 * $Author: markv $
 * $Revision: 1.1 $
 * $Date: 88/09/11 11:00:44 $
 * $Log:	sphere.c,v $
 * Revision 1.1  88/09/11  11:00:44  markv
 * Initial revision
 * 
 ***********************************************************************/

#include <stdio.h>
#include <math.h>
#include "defs.h"
#include "extern.h"

typedef struct t_spheredata {
	Vec 		sph_center ;
	Flt 		sph_radius ;
	Flt 		sph_radius2 ;
} SphereData ;

int SpherePrint ();
int SphereIntersect ();
int SphereNormal ();

ObjectProcs SphereProcs = {
	SpherePrint,
	SphereIntersect,
	SphereNormal,
} ;

int
SpherePrint(obj)
 Object *obj ;
{
	SphereData * sp ;

	sp = (SphereData *) obj -> o_data ;

	printf("s %g %g %g %g\n", sp -> sph_center[0], 
				   sp -> sph_center[1],
				   sp -> sph_center[2],
				   sp -> sph_radius) ;
return 0;
}

int
SphereIntersect(obj, ray, hit)
 Object * obj ;
 Ray * ray ;
 Isect * hit ;
{

	Flt b, disc, t;
	Point V ;
	SphereData * sp ;

	sp = (SphereData *) obj -> o_data ;

	VecSub((sp->sph_center), ray -> P, V);

	b = VecDot(V, ray -> D);

	disc = b * b - VecDot(V, V) + (sp -> sph_radius2) ;

	if (disc < 0.0)
		return(0);

	disc = sqrt(disc);

	t = (b - disc < rayeps) ? b + disc : b - disc ;

	if (t < rayeps) {
		return(0);
	}

	hit -> isect_t = t ;
	hit -> isect_enter = VecDot(V, V) > sp -> sph_radius2 + rayeps ? 1 : 0 ;
	hit -> isect_prim = obj ;
	hit -> isect_surf = obj -> o_surf ;
        return (1);
}

int
SphereNormal(obj, hit, P, N)
 Object * obj ;
 Isect * hit ;
 Point P, N ;
{
	SphereData * sp ;
	sp = (SphereData *) obj -> o_data ;

	VecSub(P, sp -> sph_center, N);
	(void) VecNormalize(N);
return 0;
}

Object *
MakeSphere(pos, radius)
 Vec pos ;
 Flt radius ;
{
	Object * tmp ;
	int i ;
	SphereData *sp ;

	tmp = (Object *) malloc (sizeof(Object)) ;
	tmp -> o_type = T_SPHERE ;
	tmp -> o_procs = & SphereProcs ;
	tmp -> o_surf = CurrentSurface ;
	sp = (SphereData *) malloc (sizeof(SphereData)) ;
	VecCopy(pos, sp -> sph_center) ;
	sp -> sph_radius = radius ;
	sp -> sph_radius2 = radius * radius ;
	tmp -> o_data = (void *) sp ;

	/*
	 * figure out dmin and dmax values for 
	 * each of the slabs...
	 */
	
	for (i = 0 ; i < NSLABS; i ++) {
		tmp -> o_dmin[i] = VecDot(sp -> sph_center, Slab[i]) 
			- sp -> sph_radius ;
		tmp -> o_dmax[i] = VecDot(sp -> sph_center, Slab[i]) 
			+ sp -> sph_radius ;
	}
	return tmp ;
}
