/***********************************************************************
 * $Author: markv $
 * $Revision: 1.1 $
 * $Date: 88/09/11 11:00:42 $
 * $Log:	pqueue.c,v $
 * Revision 1.1  88/09/11  11:00:42  markv
 * Initial revision
 * 
 ***********************************************************************/

#include <stdio.h>
#include <math.h>
#include "defs.h"
#include "extern.h"

typedef struct t_qelem {
	Flt	q_key ;
	Object	* q_obj ;
} Qelem ;

static int 	Qsize ;
static Qelem	Q[PQSIZE] ;

int
PriorityQueueNull()
{
	Qsize = 0 ;
	totalQueueResets ++ ;
#ifdef DEBUG	
	printf("resetting\n") ;
#endif /* DEBUG */
return 0;
}

PriorityQueueEmpty()
{
        return (Qsize == 0) ;
}

int
PriorityQueueInsert(key, obj)
 Flt key ;
 Object * obj ;
{
	int i ; 
	Qelem tmp ;

	totalQueues ++ ;
#ifdef DEBUG
	printf("inserting element, key = %g\n", key) ;
#endif
 	Qsize ++ ;
	if (Qsize > maxQueueSize)
		maxQueueSize = Qsize ;
	if (Qsize >= PQSIZE) {
		fprintf(stderr, "%s: exhausted priority queue space\n", 
			Progname) ;
		exit(1) ;
	}
	Q[Qsize].q_key = key ;
	Q[Qsize].q_obj = obj ;

	i = Qsize ;
	while (i > 1 && Q[i].q_key < Q[i/2].q_key) {
		tmp = Q[i] ;
		Q[i] = Q[i/2] ;
		Q[i/2] = tmp ;
		i = i / 2 ;
	}
return 0;
}

int
PriorityQueueDelete(key, obj)
 Flt * key ;
 Object ** obj ;
{
	Qelem tmp ;
	int i, j ;

	if (Qsize == 0) {
		fprintf(stderr, "%s: priority queue is empty\n",
			Progname) ;
		exit(1) ;
	}

	*key = Q[1].q_key ;
	*obj = Q[1].q_obj ;

#ifdef DEBUG
	printf("deleting element, key = %g\n", *key) ;
#endif

	Q[1] = Q[Qsize] ;
	Qsize -- ;

	i = 1 ;

	while (2 * i <= Qsize) {

		if (2 * i == Qsize) {
			j = 2 * i ;
		} else if (Q[2*i].q_key < Q[2*i+1].q_key) {
			j = 2 * i ;
		} else {
			j = 2 * i + 1 ;
		}

		if (Q[i].q_key > Q[j].q_key) {
			tmp = Q[i] ;
			Q[i] = Q[j] ;
			Q[j] = tmp ;
			i = j ;
		} else {
			break ;
		}
	}
return 0;
}
