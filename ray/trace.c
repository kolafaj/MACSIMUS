/***********************************************************************
 * $Author: markv $
 * $Revision: 1.2 $
 * $Date: 88/09/12 13:01:15 $
 * $Log:	trace.c,v $
 * Revision 1.2  88/09/12  13:01:15  markv
 * Fixed Trace to call Intersect with the max dist argument of 3.40282347e+38F.
 * 
 * Revision 1.1  88/09/11  11:00:45  markv
 * Initial revision
 * 
 ***********************************************************************/

#include <stdio.h>
#include <math.h>
#include "defs.h"
#include "extern.h"
#include "pic.h"

void
Trace(level, weight, ray, color) 
 int level;
 Flt weight;
 Ray *ray ;
 Color color ;
{
	Object *prim ;
	Vec P, N ;
	Isect hit ;

	if (level >= maxlevel) {
		color[0] = color[1] = color[2] = 0.0 ;
		return ;
	}
		
	nRays ++ ;

	if (Intersect(ray, &hit, 3.40282347e+38F)) {
		prim = hit.isect_prim ;
		RayPoint(ray, hit.isect_t, P);
		(*prim -> o_procs -> normal) (prim, &hit, P, N);
		if ((VecDot(ray->D, N)) >= 0.0) {
			VecNegate(N);
		}
		Shade(level, weight, P, N, ray -> D, &hit, color);
	} else {
                if (bg && ray->D[2]<-0.05) /* JK: */
                  GetPic(
                    bg,
                    ray->D[0]/(ray->D[2]*bg->xrange),
                    ray->D[1]/(ray->D[2]*bg->yrange),
                    bgmode, color);
                else /* orig: */
                  { VecCopy(BackgroundColor, color) ; }
	}
}
