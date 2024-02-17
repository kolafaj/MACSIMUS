/***********************************************************************
 * $Author: markv $
 * $Revision: 1.1 $
 * $Date: 88/09/11 11:00:39 $
 * $Log:	error.c,v $
 * Revision 1.1  88/09/11  11:00:39  markv
 * Initial revision
 * 
 ***********************************************************************/
#include <stdio.h>
#include <math.h>
#include "defs.h"
#include "extern.h"

/*
 * various routines to print error messages and die...
 */

int NullPrint() ;
int NullIntersect() ;
int NullNormal() ;

ObjectProcs NullProcs = {
	NullPrint,
	NullIntersect,
	NullNormal
} ;

NullPrint()
{
	fprintf(stderr, "%s: called (* print)(...), dying...\n", Progname) ;
	abort() ;
}

NullIntersect()
{
	fprintf(stderr, "%s: called (* intersect)(...), dying...\n", Progname) ;
	abort() ;
}

NullNormal()
{
	fprintf(stderr, "%s: called (* normal)(...), dying...\n", Progname) ;
	abort() ;
}
