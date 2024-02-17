/***********************************************************************

 bug fixed (unnormalized axes) by J.Kolafa 9/97

 * $Author: markv $
 * $Revision: 1.3 $
 * $Date: 88/10/04 14:30:02 $
 * $Log:	cone.c,v $
 * Revision 1.3  88/10/04  14:30:02  markv
 * Fixed bug reported by koblas@mips.
 * 
 * Revision 1.2  88/09/15  09:33:02  markv
 * Fixed bug reported by koblas@mips.  Cones which were specified with
 * axis coincident with -y axis had problems.  Caused by incorrect 
 * handling when finding orthogonal axes for the local coordinate system.
 * 
 * Revision 1.1  88/09/11  11:00:38  markv
 * Initial revision
 * 
 ***********************************************************************/
#include <stdio.h>
#include <math.h>
#include "defs.h"
#include "extern.h"

typedef struct t_conedata {
	Vec 		cone_base ;
	Flt		cone_base_radius ;
	Flt		cone_base_d ;
	Vec		cone_apex ;
	Flt		cone_apex_radius ;
	Vec		cone_u ;
	Vec		cone_v ;
	Vec		cone_w ;
	Flt		cone_height ;
	Flt		cone_slope ;
	Flt		cone_min_d ;
	Flt		cone_max_d ;
} ConeData ;

int ConePrint ();
int ConeIntersect ();
int ConeNormal ();

ObjectProcs ConeProcs = {
        ConePrint,
	ConeIntersect,
	ConeNormal,
} ;

int 
ConePrint(obj)
 Object *obj ;
{
	ConeData * cp ;

	cp = (ConeData *) obj -> o_data ;

	printf("c %g %g %g %g %g %g %g %g\n", cp -> cone_base[0], 
				  cp -> cone_base[1],
				  cp -> cone_base[2],
				  cp -> cone_base_radius,
				  cp -> cone_apex[0],
				  cp -> cone_apex[1],
				  cp -> cone_apex[2],
				  cp -> cone_apex_radius) ;
return 0;
}

int
ConeIntersect(obj, ray, hit)
 Object * obj ;
 Ray * ray ;
 Isect * hit ;
{
	Ray tray ;
	ConeData * cd ;
	Vec V, P ;
	Flt a, b, c, d, disc ;
	Flt t1, t2 ;
	int nroots ;

	cd = (ConeData *) (obj -> o_data) ;

	/*
	 * First, we get the coordinates of the ray origin in 
	 * the objects space....
	 */
	
	VecSub(ray -> P, cd -> cone_base, V) ;

	tray.P[0] = VecDot(V, cd -> cone_u) ;
	tray.P[1] = VecDot(V, cd -> cone_v) ;
	tray.P[2] = VecDot(V, cd -> cone_w) ;

	/*
	VecAdd(ray -> P, ray -> D, V) ;
	VecSub(V, cd -> cone_base, V) ;
	*/

	tray.D[0] = VecDot(ray -> D, cd -> cone_u) ;
	tray.D[1] = VecDot(ray -> D, cd -> cone_v) ;
	tray.D[2] = VecDot(ray -> D, cd -> cone_w) ;

	/*
	VecSub(tray.D, tray.P, tray.D) ;
	*/

	a = tray.D[0] * tray.D[0] 
	    + tray.D[1] * tray.D[1] 
	    - cd -> cone_slope * cd -> cone_slope * tray.D[2] * tray.D[2] ;

	b = 2.0 * (tray.P[0] * tray.D[0] + tray.P[1] * tray.D[1] - 
		cd -> cone_slope * cd -> cone_slope * tray.P[2] * tray.D[2]
		- cd -> cone_base_radius * cd -> cone_slope * tray.D[2]) ;

	c = cd -> cone_slope * tray.P[2] + cd -> cone_base_radius ;
	c = tray.P[0] * tray.P[0] + tray.P[1] * tray.P[1] - (c * c) ;

	disc = b * b - 4.0 * a * c ;

	if (disc < 0.0)
		return (0) ;
	
	disc = (Flt) sqrt(disc) ;
	t1 = (-b - disc) / (2.0 * a) ;
	t2 = (-b + disc) / (2.0 * a) ;

	if (t2 < rayeps)
		return (0) ;
	if (t1 < rayeps) {
		nroots = 1 ;
		t1 = t2 ;
	} else {
		nroots = 2 ;
	}
		
	/*
	 * ensure that the points are between the two bounding planes...
	 */
	
	switch(nroots) {
	case 1:
		RayPoint(ray, t1, P) ;
		d = VecDot(cd -> cone_w, P) ;
		if (d >= cd -> cone_min_d && d <= cd -> cone_max_d) {
			hit -> isect_t = t1 ;	
			hit -> isect_prim = obj ;
			hit -> isect_surf = obj -> o_surf ;
			hit -> isect_enter = 0 ;
			return(1) ;
		} else {
			return(0) ;
		}
		break ;
	case 2:
		RayPoint(ray, t1, P) ;
		d = VecDot(cd -> cone_w, P) ;
		if (d >= cd -> cone_min_d && d <= cd -> cone_max_d) {
			hit -> isect_t = t1 ;	
			hit -> isect_prim = obj ;
			hit -> isect_surf = obj -> o_surf ;
			hit -> isect_enter = 1 ;
			return(1) ;
		} else {
			RayPoint(ray, t2, P) ;
			d = VecDot(cd -> cone_w, P) ;
			if (d >= cd -> cone_min_d && d <= cd -> cone_max_d) {
				hit -> isect_t = t2 ;	
				hit -> isect_prim = obj ;
				hit -> isect_surf = obj -> o_surf ;
				hit -> isect_enter = 0 ;
				return(1) ;
			}
		}
		return(0) ;
	}
        return(0) ;
}

int
ConeNormal(obj, hit, P, N)
 Object * obj ;
 Isect * hit ;
 Point P, N ;
{
	Flt t ;
	Vec V ;
	ConeData * cd ;

	cd = (ConeData *) obj -> o_data ;

	/*
         * fill in the real normal...
	 * Project the point onto the base plane.  The normal is
	 * a vector from the basepoint through this point, plus the slope
	 * times the cone_w vector...
	 */

	t = - (VecDot(P, cd -> cone_w) + cd -> cone_base_d) ;
	VecAddS(t, cd -> cone_w, P, V) ;
	VecSub(V, cd -> cone_base, N) ;
	VecNormalize(N) ;
	VecAddS(- cd -> cone_slope, cd -> cone_w, N, N) ;
	VecNormalize(N) ;
return 0;
}

Object *
MakeCone(basepoint, baseradius, apexpoint, apexradius) 
 Vec basepoint, apexpoint ;
 Flt baseradius, apexradius ;
{
	Object * obj ;
	ConeData * cd ;
	Flt dmin, dmax, d , ftmp;
	Vec tmp ;
	int i ;

	obj = (Object *) malloc (sizeof (Object)) ;
	obj -> o_type = T_CONE ;
	obj -> o_procs = & ConeProcs ;
	obj -> o_surf = CurrentSurface ;

	cd = (ConeData *) malloc (sizeof(ConeData)) ;

	VecCopy(basepoint, cd -> cone_base) ;
	VecCopy(apexpoint, cd -> cone_apex) ;

	cd -> cone_base_radius = baseradius ;
	cd -> cone_apex_radius = apexradius ;


	VecSub(apexpoint, basepoint, cd -> cone_w) ;
	cd -> cone_height = VecNormalize(cd -> cone_w) ;
	cd -> cone_slope =  (cd -> cone_apex_radius - cd -> cone_base_radius) /
				(cd -> cone_height) ;
	cd -> cone_base_d = - VecDot(basepoint, cd -> cone_w) ;

	MakeVector(0.0, 0.0, 1.0, tmp) ;

	if (1.0 - fabs(VecDot(tmp, cd ->  cone_w)) < rayeps) {
		MakeVector(0.0, 1.0, 0.0, tmp) ;
	}

	/* find two axes which are at right angles to cone_w
	 */

	VecCross(cd -> cone_w, tmp, cd -> cone_u) ;
        VecCross(cd -> cone_u, cd -> cone_w, cd -> cone_v) ;

        /* added by JK: */
        VecNormalize(cd -> cone_u) ;
        VecNormalize(cd -> cone_v) ;

	cd -> cone_min_d = VecDot(cd -> cone_w, cd -> cone_base) ;
	cd -> cone_max_d = VecDot(cd -> cone_w, cd -> cone_apex) ;

	if (cd -> cone_max_d < cd -> cone_min_d) {
		ftmp = cd -> cone_max_d ;
		cd -> cone_max_d = cd -> cone_min_d ;
		cd -> cone_min_d = ftmp ;
	}

	obj -> o_data = (void *) cd ;

	for (i = 0 ; i < NSLABS ; i ++) {
		dmin = 3.40282347e+38F ;
		dmax = -3.40282347e+38F ;
		d = VecDot(basepoint, Slab[i]) - baseradius ;
		if (d < dmin) dmin = d ; if (d > dmax) dmax = d ;
		d = VecDot(basepoint, Slab[i]) + baseradius ;
		if (d < dmin) dmin = d ; if (d > dmax) dmax = d ;
		d = VecDot(apexpoint, Slab[i]) - apexradius ;
		if (d < dmin) dmin = d ; if (d > dmax) dmax = d ;
		d = VecDot(apexpoint, Slab[i]) + apexradius ;
		if (d < dmin) dmin = d ; if (d > dmax) dmax = d ;

		obj -> o_dmin[i] = dmin ;
		obj -> o_dmax[i] = dmax ;
	}

	return(obj) ;
}
