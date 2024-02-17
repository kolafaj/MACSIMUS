/***********************************************************************
 * $Author: markv $
 * $Revision: 1.2 $
 * $Date: 88/09/12 12:54:48 $
 * $Log:	intersect.c,v $
 * Revision 1.2  88/09/12  12:54:48  markv
 * Added early cutoffs, as suggested by Haines in the RT-News, and 
 * independantly discovered by myself during our correspondence.
 * Now, Intersect takes a max distance, and will not search for 
 * any intersections beyond the max distance.
 * Also, to enable shadow caching, Shadow now returns the object
 * that it found on its way to the light source.
 * 
 * These optimizations speed up shadow testing considerably, and 
 * normal rays some small amount.
 * 
 * Revision 1.1  88/09/11  11:00:40  markv
 * Initial revision
 * 
 ***********************************************************************/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "defs.h"
#include "extern.h"

/*
 * intersect.c
 * Much nicer now, uses the nifty priority queue search
 * as suggested by Kajiya...
 */

Flt		num[NSLABS] ;
Flt		den[NSLABS] ;

/***********************************************************************
 * CheckAndEnqueue(obj, maxdist)
 * Check the current ray (as paramaterized with the num and den 
 * arrays above) against the bounding volume of obj.
 * If we intersect the bounding volume, then insert it into the 
 * priority queue.
 *
 * Note: should be broken into two separate procedures...
 ***********************************************************************/

INLINE
void
CheckAndEnqueue(obj, maxdist)
 Object * obj ;
 Flt maxdist ;
{
	int i ;
	Flt tmp ;
	Flt tmin, tmax ;
	Flt dmin = -3.40282347e+38F ;
	Flt dmax = maxdist ;

	nChecked ++ ;

	for (i = 0 ; i < NSLABS ; i ++) {

		/* enters the slab here...	*/
		tmin = (obj -> o_dmin[i] - num[i]) * den[i] ;
		/* and exits here...		*/
		tmax = (obj -> o_dmax[i] - num[i]) * den[i] ;

		/* but we may have to swap...	*/
		if (tmin > tmax) {
			tmp = tmin ; tmin = tmax ; tmax = tmp ;
		}

		/* if exited closer than we thought, update 	*/
		if (tmax < dmax)
			dmax = tmax ;
		/* if entered farther than we thought, update 	*/
		if (tmin > dmin)
			dmin = tmin ;

		if (dmin > dmax || dmax < rayeps)
			return ;
	}
	PriorityQueueInsert(dmin, obj) ;
	nEnqueued ++ ;
}

/***********************************************************************
 * Intersect(ray, hit, maxdist)
 * 
 * Returns true if we hit something in the root model closer than maxdist.  
 * Returns the closest hit in the "hit" buffer.
 ***********************************************************************/

int
Intersect(ray, hit, maxdist)
 Ray * ray ;
 Isect * hit ;
 Flt maxdist ;
{
	Isect		nhit ;
	int		i ;
	Flt		min_dist = maxdist ;
	Object *	cobj ;
	Object * 	pobj = NULL ;
	CompositeData 	* cdp ;
	Flt		key ;

	/* If the object is simple, then return the hit that it gives you */

	if (Root -> o_type != T_COMPOSITE)
		return (Root -> o_procs -> intersect) (Root, ray, hit) ;

	for (i = 0 ; i < NSLABS ; i ++) {
		num[i] = VecDot(ray -> P, Slab[i]) ;
		den[i] = 1.0 / (1e-33+VecDot(ray -> D, Slab[i])) ;
	}

	/* start with an empty priority queue */
	PriorityQueueNull() ;

	CheckAndEnqueue(Root, maxdist) ;

	for (;;) {

		if (PriorityQueueEmpty())
			break ;

		PriorityQueueDelete(&key, &cobj) ;

		if (key > min_dist) {

			/*
			 * we have already found a primitive
			 * that was closer, we need look no further...
			 */
			 break ;

		} else if (cobj -> o_type == T_COMPOSITE) {
			/* 
			 * if it is in the queue, it got hit.
			 * check each of its children to see if their
			 * bounding volumes get hit.
			 * if so, then push them into the priority
			 * queue...
			 */
			
			cdp = (CompositeData *) cobj -> o_data ;

			for (i = 0 ; i < cdp -> c_size ; i ++ ) {
				CheckAndEnqueue(cdp -> c_object[i], maxdist) ;
			}

		} else {

			/*
			 * we have a primitive 
			 * intersect with the primitive, and possibly
			 * update the nearest hit if it is indeed closer
			 * than the one we currently have...
			 */
			
			if ((cobj -> o_procs -> intersect) (cobj, ray, &nhit)) {
				if (nhit.isect_t < min_dist) {
					pobj = cobj ;
					*hit = nhit ;
					min_dist = nhit.isect_t ;
				}
			}
		}
	}

	if (pobj)
		return 1 ;
	else
		return 0 ;
}

/***********************************************************************
 * Shadow(ray, hit, tmax)
 * 
 * Returns true if we are unshadowed.  Returns the primitive in 
 * the "hit" buffer.
 *
 * Note: the return value of this procedure is a bit strange, as well 
 * as the name.  Should probably be changed.
 ***********************************************************************/

int 
Shadow(ray, hit, tmax) 
 Ray * ray ;
 Isect * hit ;
 Flt tmax ;
{
	if (Intersect(ray, hit, tmax))
		return 0 ;
	else
		return 1 ;
}
