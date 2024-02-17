/***********************************************************************
 * $Author: markv $
 * $Revision: 1.1 $
 * $Date: 88/11/12 16:23:13 $
 * $Log:	tri.c,v $
 * Revision 1.1  88/11/12  16:23:13  markv
 * Initial revision
 * 
 * TRIs are triangular patches with normals defined at the vertices.
 * When an intersection is found, it interpolates the normal to the 
 * surface at that point.
 *
 * Algorithm is due to Jeff Arenburg, (arenburg@trwrb.uucp) and was 
 * posted to USENET.
 *
 * Basically, for each triangle we calculate an inverse transformation
 * matrix, and use it to determine the point of intersection in the plane
 * of the triangle relative to the "base point" of the triangle.  We then
 * figure its coordinates relative to that base point.  These base points
 * are used to find the barycentric coordinates, and then in normal 
 * interpolation...
 ***********************************************************************/

#include <stdio.h>
#include <math.h>
#include "defs.h"
#include "extern.h"

typedef struct t_spheredata {
	Vec	tri_P[3] ;
	Vec	tri_N[3] ;
	Vec	tri_bb[3] ;
} TriData ;

int TriPrint ();
int TriIntersect ();
int TriNormal ();

ObjectProcs TriProcs = {
	TriPrint,
	TriIntersect,
	TriNormal,
} ;

int 
TriPrint(obj)
 Object *obj ;
{
	TriData *td ;
	int i ;

	td = (TriData *) obj -> o_data ;
	printf("pp 3\n") ;
	for (i = 0 ; i < 3 ; i ++) {
		VecPrint("\t", td -> tri_P[i]) ;
		VecPrint("\t", td -> tri_N[i]) ;
	}
return 0;
}

int
TriIntersect(obj, ray, hit)
 Object * obj ;
 Ray * ray ;
 Isect * hit ;
{
	TriData *td ;
	Flt n, d, dist ;
	Flt r, s, t ;
	Flt a, b ;
	Vec P, Q ;

	td = (TriData *) obj -> o_data ;

	/*
	 * The matrix td -> tri_bb transforms vectors in the world 
	 * space into a space with the following properties.
	 *
	 * 1.  The sides of the triangle are coincident with the
	 *     x and y axis, and have unit length.
	 * 2.  The normal to the triangle is coincident with the 
	 *     z axis.
	 *
	 */

	/*
	 * d is the slope with respect to the z axis.  If d is zero, then
	 * the ray is parallel to the plane of the polygon, and we count 
	 * it as a miss...
	 */

	d = VecDot(ray -> D, td -> tri_bb[2]) ;
	if (fabs(d) < rayeps)
		return 0 ;

	/*
	 * Q is a vector from the eye to the triangles "origin" vertex.
	 * n is then set to be the distance of the tranformed eyepoint
	 * to the plane in the polygon.
	 * Together, n and d allow you to find the distance to the polygon, 
	 * which is merely n / d.
	 */

	VecSub(td -> tri_P[0], ray -> P, Q) ;

	n = VecDot(Q, td -> tri_bb[2]) ;

	dist = n / d ;

	if (dist < rayeps) 
		return 0 ;
	
	/* 
	 * Q is the point we hit.  Find its position relative to the
	 * origin of the triangle.
	 */

	RayPoint(ray, dist, Q) ;
	VecSub(Q, td -> tri_P[0], Q) ;

	a = VecDot(Q, td -> tri_bb[0]) ;
	b = VecDot(Q, td -> tri_bb[1]) ;

	if (a < 0.0 || b < 0.0 || a + b > 1.0) 
		return 0 ;
	
	r = 1.0 - a - b ;
	s = a ;
	t = b ;

	hit -> isect_t = dist ;
	hit -> isect_enter = 0 ;
	hit -> isect_prim = obj ;
	hit -> isect_surf = obj -> o_surf ;

	VecZero(hit->isect_normal) ;
	VecAddS(r, td -> tri_N[0], hit -> isect_normal, hit -> isect_normal) ;
	VecAddS(s, td -> tri_N[1], hit -> isect_normal, hit -> isect_normal) ;
	VecAddS(t, td -> tri_N[2], hit -> isect_normal, hit -> isect_normal) ;
	VecNormalize(hit -> isect_normal) ;

	return(1) ;
}

int
TriNormal(obj, hit, P, N)
 Object * obj ;
 Isect * hit ;
 Point P, N ;
{
	VecCopy(hit -> isect_normal, N) ;
return 0;
}

Object *
MakeTri(point)
 Vec *point ;
{
	Object * o ;
	TriData * td ;
	Vec 	N, P, Q;
	int i, j ;
	Flt dmin, dmax, d ;
	Vec B[3] ;

	o = (Object *) malloc (sizeof(Object)) ;
	o -> o_type = T_TRI ;
	o -> o_procs = & TriProcs ;
	o -> o_surf = CurrentSurface ;

	td = (TriData *) malloc (sizeof(TriData)) ;

	/* 
	 * copy in the points....
	 */
	VecCopy(point[0], td -> tri_P[0]) ;
	VecCopy(point[2], td -> tri_P[1]) ;
	VecCopy(point[4], td -> tri_P[2]) ;

	/*
	 * and the normals, then normalize them...
	 */
	VecCopy(point[1], td -> tri_N[0]) ;
	VecCopy(point[3], td -> tri_N[1]) ;
	VecCopy(point[5], td -> tri_N[2]) ;
	VecNormalize(td -> tri_N[0]) ;
	VecNormalize(td -> tri_N[1]) ;
	VecNormalize(td -> tri_N[2]) ;

	/*
	 * construct the inverse of the matrix...
	 * | P1 |
	 * | P2 |
	 * | N  |
	 * and store it in td -> tri_bb[]
	 */
	
	VecSub(td -> tri_P[1], td -> tri_P[0], B[0]) ;
	VecSub(td -> tri_P[2], td -> tri_P[0], B[1]) ;
	VecCross(B[0], B[1], B[2]) ;
	VecNormalize(B[2]) ;

	InvertMatrix(B, td -> tri_bb) ;

	for (i = 0 ; i < NSLABS ; i ++) {
		dmin = 3.40282347e+38F ;
		dmax = - 3.40282347e+38F ;

		for (j = 0 ; j < 3 ; j ++) {
			d = VecDot(Slab[i], td -> tri_P[j]) ;
			if (d < dmin) dmin = d ;
			if (d > dmax) dmax = d ;
		}
		o -> o_dmin[i] = dmin - rayeps ;
		o -> o_dmax[i] = dmax + rayeps ;
	}

	o -> o_data = (void *) td ;

	return o ;
}

int
InvertMatrix(in, out)
 Vec in[3] ;
 Vec out[3] ;
{
	int i, j, k ;
	Flt tmp, det, sum ;

	out[0][0] = (in[1][1] * in[2][2] - in[1][2] * in[2][1]) ;
	out[1][0] = -(in[0][1] * in[2][2] - in[0][2] * in[2][1]) ;
	out[2][0] = (in[0][1] * in[1][2] - in[0][2] * in[1][1]) ;

	out[0][1] = -(in[1][0] * in[2][2] - in[1][2] * in[2][0]) ;
	out[1][1] = (in[0][0] * in[2][2] - in[0][2] * in[2][0]) ;
	out[2][1] = -(in[0][0] * in[1][2] - in[0][2] * in[1][0]) ;

	out[0][2] = (in[1][0] * in[2][1] - in[1][1] * in[2][0]) ;
	out[1][2] = -(in[0][0] * in[2][1] - in[0][1] * in[2][0]) ;
	out[2][2] = (in[0][0] * in[1][1] - in[0][1] * in[1][0]) ;
	
	det = 
	in[0][0] * in[1][1] * in[2][2] +
	in[0][1] * in[1][2] * in[2][0] +
	in[0][2] * in[1][0] * in[2][1] -
	in[0][2] * in[1][1] * in[2][0] -
	in[0][0] * in[1][2] * in[2][1] -
	in[0][1] * in[1][0] * in[2][2] ;

        if (det==0) perror("singular matrix");
        det = 1 / det ;

	for (i = 0 ; i < 3 ; i ++) {
		for (j = 0 ; j < 3 ; j++) {
			out[i][j] *= det ;
		}
	}
	
#ifdef DEBUG
	for (i = 0 ; i < 3 ; i++) {
		for (j = 0 ; j < 3 ; j++) {
			sum = 0.0 ;
			for (k = 0 ; k < 3 ; k++) {
				sum += in[i][k] * out[k][j] ;
			}
			if (fabs(sum) < rayeps) {
				sum = 0.0 ;
			}
			printf(" %g") ;
		} 
		printf("\n") ;
	}
	printf("\n") ;
#endif /* DEBUG */
return 0;
}
