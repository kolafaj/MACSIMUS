/***********************************************************************
 * $Author: markv $
 * $Revision: 1.1 $
 * $Date: 88/09/11 11:00:46 $
 * $Log:	vector.c,v $
 * Revision 1.1  88/09/11  11:00:46  markv
 * Initial revision
 * 
 ***********************************************************************/
#include <stdio.h>
#include <math.h>
#include "defs.h"
#include "extern.h"

Flt 
VecNormalize(vec)
 Vec vec ;
{
	Flt len ;
	len = (Flt) VecLen(vec);
	vec[0]/=len ;
	vec[1]/=len ;
	vec[2]/=len ;
	return(len) ;
}

#ifdef DUMB_CPP
/*
 * Some machines can't handle all the vector operations, so if we define
 * DUMB_CPP, we replace them with equivalent function calls...
 */

MakeVector(x, y, z, v)
 Flt x, y, z ;
 Vec v ;
{
	v[0] = x ; v[1] = y ; v[2] = z ;
}

VecNegate(v)
 Vec v ;
{
	v[0] = -v[0] ;
	v[1] = -v[1] ;
	v[2] = -v[2] ;
}

Flt 
VecDot(a, b)
 Vec a, b ;
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] ;
}

Flt
VecLen(a)
 Vec a;
{
	return sqrt(VecDot(a, a)) ;
}

VecCopy(a, b) 
 Vec a, b ;
{
	b[0] = a[0] ;
	b[1] = a[1] ;
	b[2] = a[2] ;
}

VecAdd(a, b, c)
 Vec a, b, c ;
{
	c[0] = a[0] + b[0] ;
	c[1] = a[1] + b[1] ;
	c[2] = a[2] + b[2] ;
}

VecSub(a, b, c)
 Vec a, b, c ;
{
	c[0] = a[0] - b[0] ;
	c[1] = a[1] - b[1] ;
	c[2] = a[2] - b[2] ;
}

VecComb(A, a, B, b, c)
 Flt A, B ;
 Vec a, b, c ;
{
	c[0] = A * a[0] + B * b[0] ;	
	c[1] = A * a[1] + B * b[1] ;	
	c[2] = A * a[2] + B * b[2] ;	
}

VecAddS(A, a, b, c)
 Flt A ;
 Vec a, b, c ;
{
	c[0] = A * a[0] + b[0] ;	
	c[1] = A * a[1] + b[1] ;	
	c[2] = A * a[2] + b[2] ;	
}

VecCross(a, b, c)
 Vec a, b, c ;
{
	c[0] = a[1] * b[2] - a[2] * b[1] ;
	c[1] = a[2] * b[0] - a[0] * b[2] ;
	c[2] = a[0] * b[1] - a[1] * b[0] ;
}

#endif /* DUMB_CPP */
