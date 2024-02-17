/***********************************************************************
 * $Log:	config.h,v $
 * Revision 1.2  88/09/12  13:02:06  markv
 * Changed BUNCHINGFACTOR to something larger.
 * Added defines for compilers which don't support void.
 * Added a conditional define for SHADOW_CACHING
 * 
 * Revision 1.1  88/09/11  11:00:48  markv
 * Initial revision
 * 
 ***********************************************************************/

#ifdef SYS_V
#define index 	strchr
#define rindex 	strrchr
#endif /* SYS_V */

#define		NSLABS		(3)
#define		BUNCHINGFACTOR	(4)
/*.....#define		PQSIZE		(100)*/
#define		PQSIZE		(200) /* changed JK */

#define		XMAX		(1024)
#define		YMAX		(1024)
#define		MAXLIGHTS	(48) /* changed JK */
#define         MAXPRIMS        (4194304) /* changed JK */
#define		MAXLEVEL 	(5)

/***********************************************************************
 * If your compiler doesn't grok the void type, then define NO_VOID
 * here...
 ***********************************************************************/

/* #define NO_VOID */
#ifdef		NO_VOID
#define		void		char
#endif 		/* NO_VOID */

/***********************************************************************
 * Shadow caching is an interesting optimization.  I'd leave it in, 
 * I think it can really only help.
 ***********************************************************************/

#define SHADOW_CACHING
