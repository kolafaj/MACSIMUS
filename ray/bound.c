/***********************************************************************
 * $Author: markv $
 * $Revision: 1.3 $
 * $Date: 88/10/31 14:47:22 $
 * $Log:	bound.c,v $
 * Revision 1.3  88/10/31  14:47:22  markv
 * Removed the noisy printout which gave the bounds of the scene...
 * 
 * Revision 1.2  88/09/14  13:54:46  markv
 * Check for overflow of array Prims[].  Should fix problems
 * with rendering the gears at size factor 4.
 * 
 * Revision 1.1  88/09/11  11:00:36  markv
 * Initial revision
 * 
 ***********************************************************************/
#include <stdio.h>
#include <math.h>
#include "defs.h"
#include "extern.h"

/*
 * This function attempts to use median cut
 * to generate tighter bounding volumes than the old
 * code...
 */

int
BuildBoundingSlabs()
{
	int low = 0 ;
	int high, i ;

	high = nPrims ;
	while (SortAndSplit(low, high) == 0) {
		low = high ;
		high = nPrims ;
	}
	fprintf(stderr, "%s: after adding bounding volumes, %d prims\n",
		Progname, nPrims) ;
	fprintf(stderr, "%s: extent of scene\n", Progname) ;
#ifdef NOISY
	for (i = 0 ; i < NSLABS; i++) {
		fprintf(stderr, "%s: <%g -- %g>\n",
				Progname,
				Root -> o_dmin[i],
				Root -> o_dmax[i]) ;
	}
#endif /* NOISY */
return 0;
}

static int Axis ;

int 
compslabs(a, b)
 Object **a, **b ;
{
	Flt am, bm ;

	am = (*a) -> o_dmin[Axis] + (*a) -> o_dmax[Axis] ;
	bm = (*b) -> o_dmin[Axis] + (*b) -> o_dmax[Axis] ;

	if (am < bm)
		return (-1) ;
	else if (am == bm)
		return (0) ;
	else
		return (1) ;
}

int
FindAxis(first, last)
 int first, last ;
{
	Flt	mins[NSLABS] ;
	Flt	maxs[NSLABS] ;
        int     i, j , which=-1 ;
	Flt 	d = -3.40282347e+38F, e ;
	for (i = 0 ; i < NSLABS ; i++) {
		mins[i] = 3.40282347e+38F ;
		maxs[i] = -3.40282347e+38F ;
	}

	for (i = first ; i < last ; i++) {
		for (j = 0 ; j < NSLABS ; j++) {
			if (Prims[i] -> o_dmin[j] < mins[j])
				mins[j] = Prims[i] -> o_dmin [j] ;
			if (Prims[i] -> o_dmin[j] > maxs[j])
				maxs[j] = Prims[i] -> o_dmax [j] ;
		}
	}

	for (i = 0 ; i < NSLABS ; i++) {
		e = maxs[i] - mins[i] ;
		if (e > d) {
			d = e ;
			which = i ;
		}
	}
        if (which<0) perror("FindAxis");
	return(which) ;
}

SortAndSplit(first, last)
 int first, last ;
{
	Object * cp ;
	CompositeData * cd ;
	int size, i, j ;
	Flt dmin, dmax ;
	int m ;

	Axis = FindAxis(first, last) ;

	size = last - first ;

	/*
	 * actually, we could do this faster in several ways.
	 * we could use a logn algorithm to find the median
	 * along the given axis, and then a linear algorithm to 
	 * partition along the axis. Oh well...
	 */

	qsort((char *) (Prims + first), size, sizeof (Object *), compslabs) ;

	if (size <= BUNCHINGFACTOR) {
		/* build a box to contain them */

		cp = (Object *) malloc (sizeof(Object)) ;
		cp -> o_type = T_COMPOSITE ;
		cp -> o_procs = & NullProcs ; 	/* die if you call any 	*/
		cp -> o_surf = NULL ;		/* no surface...	*/
		cd = (CompositeData *) malloc (sizeof(CompositeData)) ;
		cd -> c_size = size ;

		for(i = 0 ; i < size ; i++) {
			cd -> c_object[i] = Prims[first + i] ;
		}

		for (i = 0 ; i < NSLABS ; i++ ) {
			dmin = 3.40282347e+38F ;
			dmax = -3.40282347e+38F ;
			for (j = 0 ; j < size ; j++) {
				if (cd -> c_object[j] -> o_dmin[i] < dmin)
					dmin = cd -> c_object[j] -> o_dmin[i] ;
				if (cd -> c_object[j] -> o_dmax[i] > dmax)
					dmax = cd -> c_object[j] -> o_dmax[i] ;
			}
			cp -> o_dmin[i] = dmin ;
			cp -> o_dmax[i] = dmax ;
		}
		cp -> o_data = (void *) cd ;
		Root = cp ;
		if (nPrims < MAXPRIMS) {
			Prims[nPrims ++] = cp ;
			return (1) ;
		} else {
			fprintf(stderr, "too many primitives, max is %d\n",
				MAXPRIMS) ;
				exit(0);
		}
	} else {
		m = (first + last) / 2 ;
		SortAndSplit(first, m) ;
		SortAndSplit(m , last) ;
		return (0) ;
	}
}
	

int
InitSlabs()
{
	int i ;

	for (i = 0 ; i < NSLABS ; i ++) {
		VecNormalize(Slab[i]) ;
	}
return 0;
}
