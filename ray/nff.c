/*.....# line 2 "nff.y"*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"
#include "extern.h"
#include <string.h>

/* added by JK */
char *get_txt(char *a, char *b)
{
static char str[256];
strcpy(str,a);
strcat(str," ");
strcat(str,b);
return str;
}

extern char yytext[] ;
extern FILE * yyin ;
Vec * pl, * plist ;

/*.....# line 15 "nff.y"*/
typedef union
#ifdef __cplusplus
	YYSTYPE
#endif
 {
	Vec	vec ;
	Vec *	vecl ;
	double	flt ;
	Object * obj ;
} YYSTYPE;
# define VIEWPOINT 257
# define FROM 258
# define AT 259
# define UP 260
# define ANGLE 261
# define HITHER 262
# define RESOLUTION 263
# define LIGHT 264
# define BACKGROUND 265
# define SURFACE 266
# define CONE 267
# define SPHERE 268
# define POLYGON 269
# define PATCH 270
# define NUM 271
# define TOKEN 272

#include <malloc.h>
#include <memory.h>
#include <unistd.h>
/*.....#include <values.h>*/

#ifdef __cplusplus

#ifndef yyerror
	void yyerror(const char *);
#endif
#ifndef yylex
	extern "C" int yylex(void);
#endif
	int yyparse(void);

#endif
#define yyclearin yychar = -1
#define yyerrok yyerrflag = 0
extern int yychar;
extern int yyerrflag;
YYSTYPE yylval;
YYSTYPE yyval;
typedef int yytabelem;
#ifndef YYMAXDEPTH
#define YYMAXDEPTH 150
#endif
#if YYMAXDEPTH > 0
int yy_yys[YYMAXDEPTH], *yys = yy_yys;
YYSTYPE yy_yyv[YYMAXDEPTH], *yyv = yy_yyv;
#else	/* user does initial allocation */
int *yys;
YYSTYPE *yyv;
#endif
static int yymaxdepth = YYMAXDEPTH;
# define YYERRCODE 256

/*.....# line 193 "nff.y"*/


yyerror(str)
 char * str ;
{
	fprintf(stderr, "%s: error at line %d\n", 
		Progname, yylinecount) ;
	fprintf(stderr, "%s: %s\n", Progname, str) ;
	exit(-1) ;
}

int
ReadSceneFile(str)
 char *str ;
{
	if (str == NULL) 
		yyin = stdin ;
	else {
		if ((yyin = fopen(str, "r")) == NULL) {
			fprintf(stderr, "%s: cannot open %s\n", Progname,
						str) ;
			exit(-1) ;
		}
	}
	if (yyparse() == 1) {
		fprintf(stderr, "%s: invalid input specification\n", Progname);
		exit(-1) ;
	}
	fprintf(stderr, "%s: %d prims, %d lights\n", 
			Progname, nPrims, nLights) ;
	fprintf(stderr, "%s: inputfile = \"%s\"\n", Progname, str) ;
	fprintf(stderr, "%s: resolution %dx%d\n", Progname, Xresolution,
		Yresolution) ;
return 0;
}

yytabelem yyexca[] ={
-1, 1,
	0, -1,
	-2, 0,
	};
# define YYNPROD 28
# define YYLAST 97
yytabelem yyact[]={

    11,    12,    13,    18,    19,    20,    21,    60,    24,    28,
    24,    58,    55,    49,    34,     5,     3,    23,    47,    26,
    41,    40,    10,     9,     8,     7,     6,     4,     2,     1,
    27,    27,    17,    29,    16,    15,    14,     0,    32,    33,
     0,    35,     0,     0,     0,    36,    52,    37,    38,    39,
     0,     0,    22,    43,    44,    45,     0,     0,    25,     0,
    48,     0,     0,    50,    51,    30,    31,     0,    54,     0,
     0,     0,    56,    57,     0,     0,    59,     0,    61,    62,
     0,    42,     0,     0,     0,    46,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,    53 };
yytabelem yypact[]={

  -241,-10000000,-10000000,  -243,  -264,  -261,-10000000,-10000000,-10000000,-10000000,
-10000000,  -261,  -263,  -263,-10000000,-10000000,-10000000,-10000000,  -261,  -261,
  -261,  -261,  -245,  -261,-10000000,-10000000,-10000000,  -261,-10000000,  -261,
  -261,  -261,-10000000,-10000000,  -261,  -261,  -261,  -261,  -261,-10000000,
-10000000,-10000000,  -247,-10000000,-10000000,  -261,  -261,  -261,  -261,  -261,
  -261,-10000000,-10000000,  -249,  -261,  -261,-10000000,  -251,  -261,  -256,
  -261,  -261,-10000000 };
yytabelem yypgo[]={

     0,    46,    19,    36,    35,    34,    32,    17,    29,    28,
    27,    26,    25,    24,    23,    22,    21,    18,    20 };
yytabelem yyr1[]={

     0,     8,    10,    10,    11,    11,    11,    11,    15,    15,
    15,    15,     9,    12,    13,    14,     3,     4,    16,     5,
    18,     6,     2,     2,     1,    17,    17,     7 };
yytabelem yyr2[]={

     0,     5,     4,     0,     2,     2,     2,     3,     2,     2,
     2,     2,    29,     5,     5,    15,    11,     7,     1,     9,
     1,     9,     7,     3,     7,     5,     0,     3 };
yytabelem yychk[]={

-10000000,    -8,    -9,   257,   -10,   258,   -11,   -12,   -13,   -14,
   -15,   264,   265,   266,    -3,    -4,    -5,    -6,   267,   268,
   269,   270,    -1,    -7,   271,    -1,    -2,    -7,   272,    -2,
    -1,    -1,    -7,    -7,   259,    -7,    -7,    -7,    -7,    -7,
   -16,   -18,    -1,    -7,    -7,    -7,    -1,   -17,   -17,   260,
    -7,    -7,    -1,    -1,    -7,   261,    -7,    -7,   262,    -7,
   263,    -7,    -7 };
yytabelem yydef[]={

     0,    -2,     3,     0,     1,     0,     2,     4,     5,     6,
     7,     0,     0,     0,     8,     9,    10,    11,     0,     0,
     0,     0,     0,     0,    27,    13,    14,     0,    23,     0,
     0,     0,    18,    20,     0,     0,     0,     0,     0,    17,
    26,    26,     0,    24,    22,     0,     0,    19,    21,     0,
     0,    16,    25,     0,     0,     0,    15,     0,     0,     0,
     0,     0,    12 };
typedef struct
#ifdef __cplusplus
	yytoktype
#endif
{ char *t_name; int t_val; } yytoktype;
#ifndef YYDEBUG
#	define YYDEBUG	0	/* don't allow debugging */
#endif

#if YYDEBUG

yytoktype yytoks[] =
{
	"VIEWPOINT",	257,
	"FROM",	258,
	"AT",	259,
	"UP",	260,
	"ANGLE",	261,
	"HITHER",	262,
	"RESOLUTION",	263,
	"LIGHT",	264,
	"BACKGROUND",	265,
	"SURFACE",	266,
	"CONE",	267,
	"SPHERE",	268,
	"POLYGON",	269,
	"PATCH",	270,
	"NUM",	271,
	"TOKEN",	272,
	"-unknown-",	-1	/* ends search */
};

char * yyreds[] =
{
	"-no such reduction-",
	"scene : camera elementlist",
	"elementlist : elementlist element",
	"elementlist : /* empty */",
	"element : light",
	"element : background",
	"element : surface",
	"element : object",
	"object : cone",
	"object : sphere",
	"object : polygon",
	"object : ppatch",
	"camera : VIEWPOINT FROM point AT point UP point ANGLE num HITHER num RESOLUTION num num",
	"light : LIGHT point",
	"background : BACKGROUND primcolor",
	"surface : SURFACE primcolor num num num num num",
	"cone : CONE point num point num",
	"sphere : SPHERE point num",
	"polygon : POLYGON num",
	"polygon : POLYGON num pointlist",
	"ppatch : PATCH num",
	"ppatch : PATCH num pointlist",
	"primcolor : num num num",
	"primcolor : TOKEN",
	"point : num num num",
	"pointlist : pointlist point",
	"pointlist : /* empty */",
	"num : NUM",
};
#endif /* YYDEBUG */
/* 
 *	Copyright 1987 Silicon Graphics, Inc. - All Rights Reserved
 */

/* #ident	"@(#)yacc:yaccpar	1.10" */
#ident	"$Revision: 1.10 $"

/*
** Skeleton parser driver for yacc output
*/
#include <stddef.h>

/*
** yacc user known macros and defines
*/
#define YYERROR		goto yyerrlab
#define YYACCEPT	return(0)
#define YYABORT		return(1)
#ifdef __cplusplus
#define YYBACKUP( newtoken, newvalue )\
{\
	if ( yychar >= 0 || ( yyr2[ yytmp ] >> 1 ) != 1 )\
	{\
                yyerror( get_txt("uxlibc:78", "syntax error - cannot backup") );\
		goto yyerrlab;\
	}\
	yychar = newtoken;\
	yystate = *yyps;\
	yylval = newvalue;\
	goto yynewstate;\
}
#else
#define YYBACKUP( newtoken, newvalue )\
{\
	if ( yychar >= 0 || ( yyr2[ yytmp ] >> 1 ) != 1 )\
	{\
                yyerror( get_txt("uxlibc:78", "Syntax error - cannot backup") );\
		goto yyerrlab;\
	}\
	yychar = newtoken;\
	yystate = *yyps;\
	yylval = newvalue;\
	goto yynewstate;\
}
#endif
#define YYRECOVERING()	(!!yyerrflag)
#define YYNEW(type)	malloc(sizeof(type) * yynewmax)
#define YYCOPY(to, from, type) \
	(type *) memcpy(to, (char *) from, yynewmax * sizeof(type))
#define YYENLARGE( from, type) \
	(type *) realloc((char *) from, yynewmax * sizeof(type))
#ifndef YYDEBUG
#	define YYDEBUG	1	/* make debugging available */
#endif

/*
** user known globals
*/
int yydebug;			/* set to 1 to get debugging */

/*
** driver internal defines
*/
#define YYFLAG		(-10000000)

/*
** global variables used by the parser
*/
YYSTYPE *yypv;			/* top of value stack */
int *yyps;			/* top of state stack */

int yystate;			/* current state */
int yytmp;			/* extra var (lasts between blocks) */

int yynerrs;			/* number of errors */
int yyerrflag;			/* error recovery flag */
int yychar;			/* current input token number */



/*
** yyparse - return 0 if worked, 1 if syntax error not recovered from
*/
#if defined(__STDC__) || defined(__cplusplus)
int yyparse(void)
#else
int yyparse()
#endif
{
	register YYSTYPE *yypvt;	/* top of value stack for $vars */

	/*
	** Initialize externals - yyparse may be called more than once
	*/
	yypv = &yyv[-1];
	yyps = &yys[-1];
	yystate = 0;
	yytmp = 0;
	yynerrs = 0;
	yyerrflag = 0;
	yychar = -1;

#if YYMAXDEPTH <= 0
	if (yymaxdepth <= 0)
	{
		if ((yymaxdepth = YYEXPAND(0)) <= 0)
		{
#ifdef __cplusplus
                        yyerror(get_txt("uxlibc:79", "yacc initialization error"));
#else
                        yyerror(get_txt("uxlibc:79", "Yacc initialization error"));
#endif
			YYABORT;
		}
	}
#endif

	goto yystack;
	{
		register YYSTYPE *yy_pv;	/* top of value stack */
		register int *yy_ps;		/* top of state stack */
		register int yy_state;		/* current state */
		register int  yy_n;		/* internal state number info */

		/*
		** get globals into registers.
		** branch to here only if YYBACKUP was called.
		*/
	yynewstate:
		yy_pv = yypv;
		yy_ps = yyps;
		yy_state = yystate;
		goto yy_newstate;

		/*
		** get globals into registers.
		** either we just started, or we just finished a reduction
		*/
	yystack:
		yy_pv = yypv;
		yy_ps = yyps;
		yy_state = yystate;

		/*
		** top of for (;;) loop while no reductions done
		*/
	yy_stack:
		/*
		** put a state and value onto the stacks
		*/
#if YYDEBUG
		/*
		** if debugging, look up token value in list of value vs.
		** name pairs.  0 and negative (-1) are special values.
		** Note: linear search is used since time is not a real
		** consideration while debugging.
		*/
		if ( yydebug )
		{
			register int yy_i;

			printf( "State %d, token ", yy_state );
			if ( yychar == 0 )
				printf( "end-of-file\n" );
			else if ( yychar < 0 )
				printf( "-none-\n" );
			else
			{
				for ( yy_i = 0; yytoks[yy_i].t_val >= 0;
					yy_i++ )
				{
					if ( yytoks[yy_i].t_val == yychar )
						break;
				}
				printf( "%s\n", yytoks[yy_i].t_name );
			}
		}
#endif /* YYDEBUG */
		if ( ++yy_ps >= &yys[ yymaxdepth ] )	/* room on stack? */
		{
			int yynewmax;
			ptrdiff_t yys_off;

			/* The following pointer-differences are safe, since
			 * yypvt, yy_pv, and yypv all are a multiple of
			 * sizeof(YYSTYPE) bytes from yyv.
			 */
			ptrdiff_t yypvt_off = yypvt - yyv;
			ptrdiff_t yy_pv_off = yy_pv - yyv;
			ptrdiff_t yypv_off = yypv - yyv;

			int *yys_base = yys;
#ifdef YYEXPAND
			yynewmax = YYEXPAND(yymaxdepth);
#else
			yynewmax = 2 * yymaxdepth;	/* double table size */
			if (yymaxdepth == YYMAXDEPTH)	/* first time growth */
			{
				void *newyys = YYNEW(int);
				void *newyyv = YYNEW(YYSTYPE);
				if (newyys != 0 && newyyv != 0)
				{
					yys = YYCOPY(newyys, yys, int);
					yyv = YYCOPY(newyyv, yyv, YYSTYPE);
				}
				else
					yynewmax = 0;	/* failed */
			}
			else				/* not first time */
			{
				yys = YYENLARGE(yys, int);
				yyv = YYENLARGE(yyv, YYSTYPE);
				if (yys == 0 || yyv == 0)
					yynewmax = 0;	/* failed */
			}
#endif
			if (yynewmax <= yymaxdepth)	/* tables not expanded */
			{
#ifdef __cplusplus
                                yyerror( get_txt("uxlibc:80", "yacc stack overflow") );
#else
                                yyerror( get_txt("uxlibc:80", "Yacc stack overflow") );
#endif
				YYABORT;
			}
			yymaxdepth = yynewmax;

			/* reset pointers into yys */
			yys_off = yys - yys_base;
			yy_ps = yy_ps + yys_off;
			yyps = yyps + yys_off;

			/* reset pointers into yyv */
			yypvt = yyv + yypvt_off;
			yy_pv = yyv + yy_pv_off;
			yypv = yyv + yypv_off;
		}
		*yy_ps = yy_state;
		*++yy_pv = yyval;

		/*
		** we have a new state - find out what to do
		*/
	yy_newstate:
		if ( ( yy_n = yypact[ yy_state ] ) <= YYFLAG )
			goto yydefault;		/* simple state */
#if YYDEBUG
		/*
		** if debugging, need to mark whether new token grabbed
		*/
		yytmp = yychar < 0;
#endif
		if ( ( yychar < 0 ) && ( ( yychar = yylex() ) < 0 ) )
			yychar = 0;		/* reached EOF */
#if YYDEBUG
		if ( yydebug && yytmp )
		{
			register int yy_i;

			printf( "Received token " );
			if ( yychar == 0 )
				printf( "end-of-file\n" );
			else if ( yychar < 0 )
				printf( "-none-\n" );
			else
			{
				for ( yy_i = 0; yytoks[yy_i].t_val >= 0;
					yy_i++ )
				{
					if ( yytoks[yy_i].t_val == yychar )
						break;
				}
				printf( "%s\n", yytoks[yy_i].t_name );
			}
		}
#endif /* YYDEBUG */
		if ( ( ( yy_n += yychar ) < 0 ) || ( yy_n >= YYLAST ) )
			goto yydefault;
		if ( yychk[ yy_n = yyact[ yy_n ] ] == yychar )	/*valid shift*/
		{
			yychar = -1;
			yyval = yylval;
			yy_state = yy_n;
			if ( yyerrflag > 0 )
				yyerrflag--;
			goto yy_stack;
		}

	yydefault:
		if ( ( yy_n = yydef[ yy_state ] ) == -2 )
		{
#if YYDEBUG
			yytmp = yychar < 0;
#endif
			if ( ( yychar < 0 ) && ( ( yychar = yylex() ) < 0 ) )
				yychar = 0;		/* reached EOF */
#if YYDEBUG
			if ( yydebug && yytmp )
			{
				register int yy_i;

				printf( "Received token " );
				if ( yychar == 0 )
					printf( "end-of-file\n" );
				else if ( yychar < 0 )
					printf( "-none-\n" );
				else
				{
					for ( yy_i = 0;
						yytoks[yy_i].t_val >= 0;
						yy_i++ )
					{
						if ( yytoks[yy_i].t_val
							== yychar )
						{
							break;
						}
					}
					printf( "%s\n", yytoks[yy_i].t_name );
				}
			}
#endif /* YYDEBUG */
			/*
			** look through exception table
			*/
			{
				register int *yyxi = yyexca;

				while ( ( *yyxi != -1 ) ||
					( yyxi[1] != yy_state ) )
				{
					yyxi += 2;
				}
				while ( ( *(yyxi += 2) >= 0 ) &&
					( *yyxi != yychar ) )
					;
				if ( ( yy_n = yyxi[1] ) < 0 )
					YYACCEPT;
			}
		}

		/*
		** check for syntax error
		*/
		if ( yy_n == 0 )	/* have an error */
		{
			/* no worry about speed here! */
			switch ( yyerrflag )
			{
			case 0:		/* new error */
#ifdef __cplusplus
                                yyerror( get_txt("uxlibc:81", "syntax error") );
#else
                                yyerror( get_txt("uxlibc:81", "Syntax error") );
#endif
				goto skip_init;
			yyerrlab:
				/*
				** get globals into registers.
				** we have a user generated syntax type error
				*/
				yy_pv = yypv;
				yy_ps = yyps;
				yy_state = yystate;
				yynerrs++;
				/* FALLTHRU */
			skip_init:
			case 1:
			case 2:		/* incompletely recovered error */
					/* try again... */
				yyerrflag = 3;
				/*
				** find state where "error" is a legal
				** shift action
				*/
				while ( yy_ps >= yys )
				{
					yy_n = yypact[ *yy_ps ] + YYERRCODE;
					if ( yy_n >= 0 && yy_n < YYLAST &&
						yychk[yyact[yy_n]] == YYERRCODE)					{
						/*
						** simulate shift of "error"
						*/
						yy_state = yyact[ yy_n ];
						goto yy_stack;
					}
					/*
					** current state has no shift on
					** "error", pop stack
					*/
#if YYDEBUG
#	define _POP_ "Error recovery pops state %d, uncovers state %d\n"
					if ( yydebug )
						printf( _POP_, *yy_ps,
							yy_ps[-1] );
#	undef _POP_
#endif
					yy_ps--;
					yy_pv--;
				}
				/*
				** there is no state on stack with "error" as
				** a valid shift.  give up.
				*/
				YYABORT;
			case 3:		/* no shift yet; eat a token */
#if YYDEBUG
				/*
				** if debugging, look up token in list of
				** pairs.  0 and negative shouldn't occur,
				** but since timing doesn't matter when
				** debugging, it doesn't hurt to leave the
				** tests here.
				*/
				if ( yydebug )
				{
					register int yy_i;

					printf( "Error recovery discards " );
					if ( yychar == 0 )
						printf( "token end-of-file\n" );
					else if ( yychar < 0 )
						printf( "token -none-\n" );
					else
					{
						for ( yy_i = 0;
							yytoks[yy_i].t_val >= 0;
							yy_i++ )
						{
							if ( yytoks[yy_i].t_val
								== yychar )
							{
								break;
							}
						}
						printf( "token %s\n",
							yytoks[yy_i].t_name );
					}
				}
#endif /* YYDEBUG */
				if ( yychar == 0 )	/* reached EOF. quit */
					YYABORT;
				yychar = -1;
				goto yy_newstate;
			}
		}/* end if ( yy_n == 0 ) */
		/*
		** reduction by production yy_n
		** put stack tops, etc. so things right after switch
		*/
#if YYDEBUG
		/*
		** if debugging, print the string that is the user's
		** specification of the reduction which is just about
		** to be done.
		*/
		if ( yydebug )
			printf( "Reduce by (%d) \"%s\"\n",
				yy_n, yyreds[ yy_n ] );
#endif
		yytmp = yy_n;			/* value to switch over */
		yypvt = yy_pv;			/* $vars top of value stack */
		/*
		** Look in goto table for next state
		** Sorry about using yy_state here as temporary
		** register variable, but why not, if it works...
		** If yyr2[ yy_n ] doesn't have the low order bit
		** set, then there is no action to be done for
		** this reduction.  So, no saving & unsaving of
		** registers done.  The only difference between the
		** code just after the if and the body of the if is
		** the goto yy_stack in the body.  This way the test
		** can be made before the choice of what to do is needed.
		*/
		{
			/* length of production doubled with extra bit */
			register int yy_len = yyr2[ yy_n ];

			if ( !( yy_len & 01 ) )
			{
				yy_len >>= 1;
				yyval = ( yy_pv -= yy_len )[1];	/* $$ = $1 */
				yy_state = yypgo[ yy_n = yyr1[ yy_n ] ] +
					*( yy_ps -= yy_len ) + 1;
				if ( yy_state >= YYLAST ||
					yychk[ yy_state =
					yyact[ yy_state ] ] != -yy_n )
				{
					yy_state = yyact[ yypgo[ yy_n ] ];
				}
				goto yy_stack;
			}
			yy_len >>= 1;
			yyval = ( yy_pv -= yy_len )[1];	/* $$ = $1 */
			yy_state = yypgo[ yy_n = yyr1[ yy_n ] ] +
				*( yy_ps -= yy_len ) + 1;
			if ( yy_state >= YYLAST ||
				yychk[ yy_state = yyact[ yy_state ] ] != -yy_n )
			{
				yy_state = yyact[ yypgo[ yy_n ] ];
			}
		}
					/* save until reenter driver code */
		yystate = yy_state;
		yyps = yy_ps;
		yypv = yy_pv;
	}
	/*
	** code supplied by user is placed in this switch
	*/
	switch( yytmp )
	{
		
case 1:
/*.....# line 30 "nff.y"*/
{
		int i, l ;
		for (l = 0 ; l < nLights ; l++) {
		        VecCopy(Lights[l].light_pos,Lights[l].light_pos0);
			Lights[l].light_brightness = lightfactor*(1-isotropiclight-normaldirlight)/nLights ;
			/* JK: was 1/sqrt((Flt)nLights) */
			for (i = 0 ; i < MAXLEVEL ; i++) {
				Lights[l].light_obj_cache[i] = NULL ;
			}
		}
	} break;
case 7:
/*.....# line 50 "nff.y"*/
{
		char buf[80] ;
		if (nPrims >= MAXPRIMS) {
			sprintf(buf, "max objects = %d", MAXPRIMS) ;
			yyerror(buf) ;
		}
	} break;
case 12:
/*.....# line 71 "nff.y"*/
{
		VecCopy(yypvt[-11].vec, Eye.view_from) ;
		VecCopy(yypvt[-9].vec, Eye.view_at) ;
		VecCopy(yypvt[-7].vec, Eye.view_up) ;
		Eye.view_angle = degtorad(yypvt[-5].flt/2.0 ) ;
		Eye.view_dist = yypvt[-3].flt ;
        if (xresolutionflag > 0) Xresolution=xresolutionflag ;
        else Xresolution = (int) yypvt[-1].flt ;
        if (yresolutionflag > 0) Yresolution=yresolutionflag ;
        else Yresolution = (int) yypvt[-0].flt ;
	} break;
case 13:
/*.....# line 88 "nff.y"*/
{
		if (nLights>=MAXLIGHTS) {
		  fprintf(stderr,"too many lights\n"); exit(-1); }
		VecCopy(yypvt[-0].vec, Lights[nLights].light_pos) ;
		/* fill in brightness of the light, after we 
		 * know the number of lights sources in the scene
		 */
		nLights ++ ;
	} break;
case 14:
/*.....# line 98 "nff.y"*/
{
		VecCopy(yypvt[-0].vec, BackgroundColor) ;
	} break;
case 15:
/*.....# line 104 "nff.y"*/
{
		CurrentSurface = (Surface *) malloc (sizeof(Surface)) ;
		VecCopy(yypvt[-5].vec, CurrentSurface -> surf_color) ;
		CurrentSurface -> surf_kd = yypvt[-4].flt ;
		CurrentSurface -> surf_ks = yypvt[-3].flt ;
		CurrentSurface -> surf_shine = yypvt[-2].flt ;
		CurrentSurface -> surf_kt = yypvt[-1].flt ;
		CurrentSurface -> surf_ior = yypvt[-0].flt ;
	} break;
case 16:
/*.....# line 116 "nff.y"*/
{
		yyval.obj = MakeCone(yypvt[-3].vec, yypvt[-2].flt, yypvt[-1].vec, yypvt[-0].flt) ;
		Prims[nPrims++] = yyval.obj ;

	} break;
case 17:
/*.....# line 123 "nff.y"*/
{
		yyval.obj = MakeSphere(yypvt[-1].vec, yypvt[-0].flt) ;
		Prims[nPrims++] = yyval.obj ;
	} break;
case 18:
/*.....# line 130 "nff.y"*/
{
		plist = (Vec *) calloc((int) yypvt[-0].flt, sizeof(Vec)) ;
		pl = plist ;
	} break;
case 19:
/*.....# line 135 "nff.y"*/
{
		yyval.obj = MakePoly((int) yypvt[-2].flt, plist) ;
		Prims[nPrims++] = yyval.obj ;
	} break;
case 20:
/*.....# line 141 "nff.y"*/
{
		if ((int) yypvt[-0].flt != 3)
			fprintf(stderr, "patches must have 3 vertices...\n") ;
		plist = (Vec *) calloc(2 * (int) yypvt[-0].flt, sizeof(Vec)) ;
		pl = plist ;
	} break;
case 21:
/*.....# line 148 "nff.y"*/
{
		yyval.obj = MakeTri(plist) ;
		Prims[nPrims++] = yyval.obj ;
	} break;
case 22:
/*.....# line 155 "nff.y"*/
{
		yyval.vec[0] = yypvt[-2].flt ;
		yyval.vec[1] = yypvt[-1].flt ;
		yyval.vec[2] = yypvt[-0].flt ;
	} break;
case 23:
/*.....# line 161 "nff.y"*/
{
		char buf[80] ;

		if (LookupColorByName(yytext, yyval.vec) == 0) {
			sprintf(buf, "cannot find color \"%s\"\n",
				yytext) ;
			yyerror(buf) ;
		}
	} break;
case 24:
/*.....# line 173 "nff.y"*/
{
		yyval.vec[0] = yypvt[-2].flt ;
		yyval.vec[1] = yypvt[-1].flt ;
		yyval.vec[2] = yypvt[-0].flt ;
	} break;
case 25:
/*.....# line 181 "nff.y"*/
{
		VecCopy(yypvt[-0].vec, (*pl)) ;
		pl ++ ;
	} break;
case 27:
/*.....# line 189 "nff.y"*/
{
			yyval.flt = atof(yytext) ;
		} break;
	}
	goto yystack;		/* reset registers in driver code */
}
