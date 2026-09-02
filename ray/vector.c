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

Flt VecNormalize(Vec vec) /************************************ VecNormalize */
{
	Flt len ;
	len = (Flt) VecLen(vec);
	vec[0]/=len ;
	vec[1]/=len ;
	vec[2]/=len ;
	return(len) ;
}
