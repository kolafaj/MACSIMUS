/***********************************************************************
 * $Author: markv $
 * $Revision: 1.2 $
 * $Date: 88/10/04 14:32:25 $
 * $Log:	poly.c,v $
 * Revision 1.2  88/10/04  14:32:25  markv
 * Renamed p1 and p2 to be poly_u and poly_v, which are better names.
 * Also may have solved problems having to do with floating point roundoff
 * when planes are parallel to bounding slabs.
 * 
 * Revision 1.1  88/09/11  11:00:42  markv
 * Initial revision
 * 
 ***********************************************************************/

#include <stdio.h>
#include <math.h>
#include "defs.h"
#include "extern.h"

typedef struct t_polydata {
	int 	poly_npoints ;
	Vec	* poly_point ;
	Vec	poly_normal ;
	Flt 	poly_d ;
	Flt	poly_u, poly_v ;
} PolyData ;

int PolyPrint ();
int PolyIntersect ();
int PolyNormal ();

ObjectProcs PolyProcs = {
	PolyPrint,
	PolyIntersect,
	PolyNormal,
} ;

int 
PolyPrint(obj)
 Object * obj ;
{
	int i ;
	PolyData * pd ;

	pd = (PolyData *) obj -> o_data ;
	printf("p %d\n", pd -> poly_npoints) ;
	for (i = 0 ; i < pd -> poly_npoints ; i++) {
		printf("%g %g %g\n", pd -> poly_point[i][0],
					pd -> poly_point[i][1],
					pd -> poly_point[i][2]) ;
	}
}

/***********************************************************************
 * PolyIntersect(obj, ray, hit)
 * 
 * returns 1 if we hit the polygon, with the hit information in hit.
 * Uses a version of Jordan's theorem to determine whether the point 
 * is inside the polygon.
 * The variable "crosses" will count the number of times that we
 * cross the boundary of the curve.  If it is odd, we are inside.
 *
 ***********************************************************************/

PolyIntersect(obj, ray, hit)
 Object * obj ;
 Ray * ray ;
 Isect * hit ;
{
	Flt n,d,t,m,b ;
	Point V ;
	int i, j, crosses = 0 ;
	int qi, qj ;
	int ri, rj ;
	int u, v ;
	PolyData * pd ;

	pd = (PolyData *) obj -> o_data ;
	n = VecDot(ray -> P, pd -> poly_normal) + pd -> poly_d ;
	d = VecDot(ray -> D, pd -> poly_normal) ;

	if ((Flt) fabs(d) < rayeps) {
		return (0) ;
	}

	t = - n / d ;
	if (t < rayeps)
		return 0 ;

	RayPoint(ray,t,V);

	u = pd -> poly_u ;
	v = pd -> poly_v ;

	for (i = 0 ; i < pd -> poly_npoints ; i++) {

		j = (i + 1) % pd -> poly_npoints ;

		qi = 0 ; qj = 0 ;
		ri = 0 ; rj = 0 ;

		if (pd -> poly_point[i][v] == pd -> poly_point[j][v])
			continue ;		/*ignore horizontal lines */

		/*
		 * If we are both above, or both below the intersection point,
		 * go onto the next one.
		 */

		if (pd -> poly_point[i][v] < V[v])
			qi = 1 ;
		if (pd -> poly_point[j][v] < V[v])
			qj = 1 ;
		if (qi == qj)
			continue ;

		/*
		 * We know one end point was above, and one was below.
		 * If theyare both to the left, then we crossed the line 
		 * to negative infinity, and we continue.
		 */
		if (pd -> poly_point[i][u] < V[u])
			ri = 1 ;
 		if (pd -> poly_point[j][u] < V[u])
			rj = 1 ;

		if (ri & rj) {
			crosses ++ ;
			continue ;
		}

		/*
		 * Otherwise, if we are both to the right, 
		 * we can continue without a cross.
		 */

		if ((ri|rj) == 0)
			continue ;


		/* 
		 * more difficult acceptance...
		 * We have a line segment which occurs with endpoints
		 * in diagonally opposite quadrants.  We must solve
		 * for the intersection, ie, where v = 0.
		 */

		m = (pd -> poly_point[j][v] - pd -> poly_point[i][v]) /
			(pd -> poly_point[j][u] - pd -> poly_point[i][u]) ;
		
		b = (pd -> poly_point[j][v] - V[v]) - 
			m * (pd -> poly_point[j][u] - V[u]);

		if ((-b/m) < rayeps)
			crosses ++ ;
	}

	if (crosses & 1) {
		hit -> isect_t = t ;
		hit -> isect_surf = obj -> o_surf ;
		hit -> isect_prim = obj ;
		hit -> isect_enter = 0 ;
		return(1);
	} else {
		return(0);
	}
}

PolyNormal(obj, hit, P, N)
 Object * obj ;
 Isect * hit ;
 Point P, N ;
{

	PolyData * pd ;
	pd = (PolyData *) obj -> o_data ;
	VecCopy(pd -> poly_normal, N);
}

Object *
MakePoly(npoints, points)
 int npoints ;
 Vec * points ;
{
	Object * obj ;
	PolyData * pd ;
	Vec P1, P2 ;
	Flt d, dmax, dmin ;
	int i, j ;

	obj = (Object *) malloc (sizeof(Object)) ;
	obj -> o_type = T_POLY ;
	obj -> o_procs = & PolyProcs ;
	obj -> o_surf = CurrentSurface ;

	pd = (PolyData *) malloc (sizeof(PolyData)) ;
	pd -> poly_npoints = npoints ;
	pd -> poly_point = points ;

	/*
	 * calculate the normal by giving various cross products...
	 */
	
	VecSub(pd -> poly_point[0], pd -> poly_point[1], P1) ;
	VecSub(pd -> poly_point[2], pd -> poly_point[1], P2) ;

	VecCross(P1, P2, pd -> poly_normal) ;
	VecNormalize(pd -> poly_normal) ;

	if (fabs(pd -> poly_normal[0]) > fabs(pd -> poly_normal[1])
		&& fabs(pd -> poly_normal[0]) > fabs(pd -> poly_normal[2])) {
		pd -> poly_u = 1 ;
		pd -> poly_v = 2 ;
	} else if (fabs(pd -> poly_normal[1]) > fabs(pd -> poly_normal[0]) 
		&& fabs(pd -> poly_normal[1]) > fabs(pd -> poly_normal[2])) {
		pd -> poly_u = 0 ;
		pd -> poly_v = 2 ;
	} else {
		pd -> poly_u = 0 ;
		pd -> poly_v = 1 ;
	}

	pd -> poly_d = - VecDot(pd -> poly_normal, pd -> poly_point[0]) ;

	obj -> o_data = (void *) pd ;

	/*
	 * now, calculate the values of 
	 * the dmin and dmax 'es for the globally defined slabs...
	 */
	
	for (i = 0 ; i < NSLABS ; i ++) {
		dmin = 3.40282347e+38F ;
		dmax = - 3.40282347e+38F ;

		for (j = 0 ; j < pd -> poly_npoints ; j ++) {
			d = VecDot(Slab[i], pd -> poly_point[j]) ;
			if (d < dmin) dmin = d ;
			if (d > dmax) dmax = d ;
		}
		obj -> o_dmin[i] = dmin - rayeps ;
		obj -> o_dmax[i] = dmax + rayeps ;
	}
	return(obj) ;
}
