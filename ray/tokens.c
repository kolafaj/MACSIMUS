#include <stdio.h>
# define U(x) ((x)&0377)
# define NLSTATE yyprevious=YYNEWLINE
# define BEGIN yybgin = yysvec + 1 +
# define INITIAL 0
# define YYLERR yysvec
# define YYSTATE (yyestate-yysvec-1)
# define YYOPTIM 1
# define YYLMAX 200
# define output(c) putc(c,yyout)
# define input() (((yytchar=yysptr>yysbuf?U(*--yysptr):getc(yyin))==10?(yylineno++,yytchar):yytchar)==EOF?0:yytchar)
# define unput(c) {yytchar= (c);if(yytchar=='\n')yylineno--;*yysptr++=yytchar;}
# define yymore() (yymorfg=1)
# define ECHO fprintf(yyout, "%s",yytext)
# define REJECT { nstr = yyreject(); goto yyfussy;}
int yyleng; extern unsigned char yytext[];
int yymorfg;
extern unsigned char *yysptr, yysbuf[];
int yytchar;
FILE *yyin, *yyout ;
extern int yylineno;
struct yysvf {
	struct yywork *yystoff;
	struct yysvf *yyother;
	int *yystops;};
struct yysvf *yyestate;
extern struct yysvf yysvec[], *yybgin;
#include <stdio.h>
#include "defs.h"
/*.....#include "y.tab.h"*/
#include "y_tab.h"
#include "extern.h"
# define YYNEWLINE 10
yylex(){
int nstr; extern int yyprevious;
while((nstr = yylook()) >= 0)
yyfussy: switch(nstr){
case 0:
if(yywrap()) return(0); break;
case 1:
		;
break;
case 2:
		;
break;
case 3:
		yylinecount ++ ;
break;
case 4:
		return VIEWPOINT ;
break;
case 5:
	return VIEWPOINT ;
break;
case 6:
		return FROM ;
break;
case 7:
		return AT ;
break;
case 8:
		return UP ;
break;
case 9:
		return ANGLE ;
break;
case 10:
		return HITHER ;
break;
case 11:
	return RESOLUTION ;
break;
case 12:
		return LIGHT ;
break;
case 13:
		return LIGHT ;
break;
case 14:
		return BACKGROUND ;
break;
case 15:
	return BACKGROUND ;
break;
case 16:
		return SURFACE ;
break;
case 17:
		return SURFACE ;
break;
case 18:
		return CONE ;
break;
case 19:
		return CONE ;
break;
case 20:
		return SPHERE ;
break;
case 21:
		return SPHERE ;
break;
case 22:
		return POLYGON ;
break;
case 23:
		return POLYGON ;
break;
case 24:
		return PATCH ;
break;
case 25:
		return PATCH ;
break;
case 26:
	return NUM ;
break;
case 27:
	return TOKEN ;
break;
case 28:
		return yytext[0] ;
break;
case -1:
break;
default:
fprintf(yyout,"bad switch yylook %d",nstr);
} return(0); }
/* end of yylex */

yywrap()
{
	return 1 ;
}
int yyvstop[] = {
0,

26,
0,

26,
0,

28,
0,

1,
28,
0,

3,
0,

28,
-2,
0,

26,
28,
0,

26,
28,
0,

26,
27,
28,
0,

27,
28,
0,

27,
28,
0,

14,
27,
28,
0,

18,
27,
28,
0,

16,
27,
28,
0,

27,
28,
0,

12,
27,
28,
0,

22,
27,
28,
0,

27,
28,
0,

20,
27,
28,
0,

27,
28,
0,

4,
27,
28,
0,

-2,
0,

2,
0,

26,
0,

26,
0,

26,
27,
0,

27,
0,

27,
0,

7,
27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

24,
27,
0,

27,
0,

27,
0,

27,
0,

8,
27,
0,

27,
0,

26,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

19,
27,
0,

6,
27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

9,
27,
0,

27,
0,

27,
0,

13,
27,
0,

25,
27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

10,
27,
0,

27,
0,

27,
0,

21,
27,
0,

27,
0,

27,
0,

27,
0,

23,
27,
0,

27,
0,

17,
27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

27,
0,

5,
27,
0,

15,
27,
0,

11,
27,
0,
0};
# define YYTYPE unsigned char
struct yywork { YYTYPE verify, advance; } yycrank[] = {
	{0,0},	{0,0},	{1,3},	{0,0},
	{0,0},	{0,0},	{6,22},	{0,0},
	{0,0},	{0,0},	{1,4},	{1,5},
	{0,0},	{0,0},	{6,22},	{6,23},
	{0,0},	{0,0},	{0,0},	{0,0},
	{0,0},	{0,0},	{0,0},	{0,0},
	{0,0},	{0,0},	{0,0},	{0,0},
	{0,0},	{0,0},	{0,0},	{0,0},
	{0,0},	{0,0},	{0,0},	{0,0},
	{1,6},	{0,0},	{0,0},	{0,0},
	{0,0},	{0,0},	{0,0},	{0,0},
	{0,0},	{0,0},	{1,7},	{1,8},
	{0,0},	{1,9},	{0,0},	{0,0},
	{0,0},	{6,22},	{0,0},	{0,0},
	{0,0},	{0,0},	{2,6},	{0,0},
	{0,0},	{0,0},	{0,0},	{0,0},
	{0,0},	{0,0},	{1,10},	{0,0},
	{2,7},	{2,8},	{6,22},	{7,24},
	{0,0},	{7,25},	{7,25},	{7,25},
	{7,25},	{7,25},	{7,25},	{7,25},
	{7,25},	{7,25},	{7,25},	{8,24},
	{8,24},	{8,24},	{8,24},	{8,24},
	{8,24},	{8,24},	{8,24},	{8,24},
	{8,24},	{0,0},	{0,0},	{0,0},
	{0,0},	{0,0},	{1,11},	{1,12},
	{1,13},	{0,0},	{0,0},	{1,14},
	{12,31},	{1,15},	{31,47},	{52,64},
	{18,39},	{1,16},	{29,46},	{35,51},
	{15,34},	{1,17},	{16,35},	{1,18},
	{1,19},	{11,29},	{1,20},	{1,21},
	{2,11},	{2,12},	{2,13},	{11,30},
	{13,32},	{2,14},	{14,33},	{2,15},
	{19,40},	{20,42},	{21,43},	{2,16},
	{32,48},	{19,41},	{33,49},	{2,17},
	{8,26},	{2,18},	{2,19},	{34,50},
	{2,20},	{2,21},	{9,24},	{17,36},
	{9,27},	{9,27},	{9,27},	{9,27},
	{9,27},	{9,27},	{9,27},	{9,27},
	{9,27},	{9,27},	{36,52},	{37,53},
	{39,54},	{17,37},	{17,38},	{40,55},
	{41,56},	{9,28},	{9,28},	{9,28},
	{9,28},	{9,28},	{9,28},	{9,28},
	{9,28},	{9,28},	{9,28},	{9,28},
	{9,28},	{9,28},	{9,28},	{9,28},
	{9,28},	{9,28},	{9,28},	{9,28},
	{9,28},	{9,28},	{9,28},	{9,28},
	{9,28},	{9,28},	{9,28},	{43,57},
	{46,58},	{47,59},	{48,60},	{9,28},
	{49,61},	{9,28},	{9,28},	{9,28},
	{9,28},	{9,28},	{9,28},	{9,28},
	{9,28},	{9,28},	{9,28},	{9,28},
	{9,28},	{9,28},	{9,28},	{9,28},
	{9,28},	{9,28},	{9,28},	{9,28},
	{9,28},	{9,28},	{9,28},	{9,28},
	{9,28},	{9,28},	{9,28},	{10,28},
	{10,28},	{10,28},	{10,28},	{10,28},
	{10,28},	{10,28},	{10,28},	{10,28},
	{10,28},	{50,62},	{51,63},	{53,65},
	{54,66},	{55,67},	{56,68},	{57,69},
	{10,28},	{10,28},	{10,28},	{10,28},
	{10,28},	{10,28},	{10,28},	{10,28},
	{10,28},	{10,28},	{10,28},	{10,28},
	{10,28},	{10,28},	{10,28},	{10,28},
	{10,28},	{10,28},	{10,28},	{10,28},
	{10,28},	{10,28},	{10,28},	{10,28},
	{10,28},	{10,28},	{58,70},	{59,71},
	{62,72},	{63,73},	{10,28},	{64,74},
	{10,28},	{10,28},	{10,28},	{10,28},
	{10,28},	{10,28},	{10,28},	{10,28},
	{10,28},	{10,28},	{10,28},	{10,28},
	{10,28},	{10,28},	{10,28},	{10,28},
	{10,28},	{10,28},	{10,28},	{10,28},
	{10,28},	{10,28},	{10,28},	{10,28},
	{10,28},	{10,28},	{26,44},	{65,75},
	{66,76},	{26,45},	{26,45},	{26,45},
	{26,45},	{26,45},	{26,45},	{26,45},
	{26,45},	{26,45},	{26,45},	{44,45},
	{44,45},	{44,45},	{44,45},	{44,45},
	{44,45},	{44,45},	{44,45},	{44,45},
	{44,45},	{67,77},	{68,78},	{69,79},
	{71,80},	{72,81},	{75,82},	{76,83},
	{77,84},	{78,85},	{79,86},	{80,87},
	{82,88},	{83,89},	{85,90},	{86,91},
	{87,92},	{89,93},	{91,94},	{92,95},
	{93,96},	{94,97},	{95,98},	{96,99},
	{0,0}};
struct yysvf yysvec[] = {
	{0,	0,	0},
	{yycrank+-1,	0,		yyvstop+1},
	{yycrank+-23,	yysvec+1,	yyvstop+3},
	{yycrank+0,	0,		yyvstop+5},
	{yycrank+0,	0,		yyvstop+7},
	{yycrank+0,	0,		yyvstop+10},
	{yycrank+-5,	0,		yyvstop+12},
	{yycrank+25,	0,		yyvstop+15},
	{yycrank+35,	0,		yyvstop+18},
	{yycrank+96,	0,		yyvstop+21},
	{yycrank+171,	0,		yyvstop+25},
	{yycrank+7,	yysvec+10,	yyvstop+28},
	{yycrank+7,	yysvec+10,	yyvstop+31},
	{yycrank+13,	yysvec+10,	yyvstop+35},
	{yycrank+12,	yysvec+10,	yyvstop+39},
	{yycrank+7,	yysvec+10,	yyvstop+43},
	{yycrank+9,	yysvec+10,	yyvstop+46},
	{yycrank+46,	yysvec+10,	yyvstop+50},
	{yycrank+7,	yysvec+10,	yyvstop+54},
	{yycrank+16,	yysvec+10,	yyvstop+57},
	{yycrank+17,	yysvec+10,	yyvstop+61},
	{yycrank+25,	yysvec+10,	yyvstop+64},
	{yycrank+0,	yysvec+6,	yyvstop+68},
	{yycrank+0,	0,		yyvstop+70},
	{yycrank+0,	yysvec+8,	yyvstop+72},
	{yycrank+0,	yysvec+7,	yyvstop+74},
	{yycrank+249,	0,		0},
	{yycrank+0,	yysvec+9,	yyvstop+76},
	{yycrank+0,	yysvec+10,	yyvstop+79},
	{yycrank+7,	yysvec+10,	yyvstop+81},
	{yycrank+0,	yysvec+10,	yyvstop+83},
	{yycrank+7,	yysvec+10,	yyvstop+86},
	{yycrank+22,	yysvec+10,	yyvstop+88},
	{yycrank+23,	yysvec+10,	yyvstop+90},
	{yycrank+23,	yysvec+10,	yyvstop+92},
	{yycrank+8,	yysvec+10,	yyvstop+94},
	{yycrank+38,	yysvec+10,	yyvstop+96},
	{yycrank+47,	yysvec+10,	yyvstop+98},
	{yycrank+0,	yysvec+10,	yyvstop+100},
	{yycrank+41,	yysvec+10,	yyvstop+103},
	{yycrank+55,	yysvec+10,	yyvstop+105},
	{yycrank+46,	yysvec+10,	yyvstop+107},
	{yycrank+0,	yysvec+10,	yyvstop+109},
	{yycrank+86,	yysvec+10,	yyvstop+112},
	{yycrank+259,	0,		0},
	{yycrank+0,	yysvec+44,	yyvstop+114},
	{yycrank+80,	yysvec+10,	yyvstop+116},
	{yycrank+82,	yysvec+10,	yyvstop+118},
	{yycrank+89,	yysvec+10,	yyvstop+120},
	{yycrank+83,	yysvec+10,	yyvstop+122},
	{yycrank+125,	yysvec+10,	yyvstop+124},
	{yycrank+126,	yysvec+10,	yyvstop+126},
	{yycrank+8,	yysvec+10,	yyvstop+128},
	{yycrank+110,	yysvec+10,	yyvstop+130},
	{yycrank+121,	yysvec+10,	yyvstop+132},
	{yycrank+132,	yysvec+10,	yyvstop+134},
	{yycrank+132,	yysvec+10,	yyvstop+136},
	{yycrank+116,	yysvec+10,	yyvstop+138},
	{yycrank+161,	yysvec+10,	yyvstop+140},
	{yycrank+160,	yysvec+10,	yyvstop+142},
	{yycrank+0,	yysvec+10,	yyvstop+144},
	{yycrank+0,	yysvec+10,	yyvstop+147},
	{yycrank+163,	yysvec+10,	yyvstop+150},
	{yycrank+149,	yysvec+10,	yyvstop+152},
	{yycrank+163,	yysvec+10,	yyvstop+154},
	{yycrank+192,	yysvec+10,	yyvstop+156},
	{yycrank+188,	yysvec+10,	yyvstop+158},
	{yycrank+203,	yysvec+10,	yyvstop+160},
	{yycrank+221,	yysvec+10,	yyvstop+162},
	{yycrank+207,	yysvec+10,	yyvstop+164},
	{yycrank+0,	yysvec+10,	yyvstop+166},
	{yycrank+206,	yysvec+10,	yyvstop+169},
	{yycrank+207,	yysvec+10,	yyvstop+171},
	{yycrank+0,	yysvec+10,	yyvstop+173},
	{yycrank+0,	yysvec+10,	yyvstop+176},
	{yycrank+211,	yysvec+10,	yyvstop+179},
	{yycrank+206,	yysvec+10,	yyvstop+181},
	{yycrank+223,	yysvec+10,	yyvstop+183},
	{yycrank+226,	yysvec+10,	yyvstop+185},
	{yycrank+215,	yysvec+10,	yyvstop+187},
	{yycrank+216,	yysvec+10,	yyvstop+189},
	{yycrank+0,	yysvec+10,	yyvstop+191},
	{yycrank+218,	yysvec+10,	yyvstop+194},
	{yycrank+213,	yysvec+10,	yyvstop+196},
	{yycrank+0,	yysvec+10,	yyvstop+198},
	{yycrank+229,	yysvec+10,	yyvstop+201},
	{yycrank+226,	yysvec+10,	yyvstop+203},
	{yycrank+215,	yysvec+10,	yyvstop+205},
	{yycrank+0,	yysvec+10,	yyvstop+207},
	{yycrank+228,	yysvec+10,	yyvstop+210},
	{yycrank+0,	yysvec+10,	yyvstop+212},
	{yycrank+224,	yysvec+10,	yyvstop+215},
	{yycrank+225,	yysvec+10,	yyvstop+217},
	{yycrank+225,	yysvec+10,	yyvstop+219},
	{yycrank+221,	yysvec+10,	yyvstop+221},
	{yycrank+238,	yysvec+10,	yyvstop+223},
	{yycrank+229,	yysvec+10,	yyvstop+225},
	{yycrank+0,	yysvec+10,	yyvstop+227},
	{yycrank+0,	yysvec+10,	yyvstop+230},
	{yycrank+0,	yysvec+10,	yyvstop+233},
	{0,	0,	0}};
struct yywork *yytop = yycrank+339;
struct yysvf *yybgin = yysvec+1;
unsigned char yymatch[] = {
00  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,011 ,012 ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
011 ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
'0' ,'0' ,'0' ,'0' ,'0' ,'0' ,'0' ,'0' ,
'0' ,'0' ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,
'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,
'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,
'A' ,'A' ,'A' ,01  ,01  ,01  ,01  ,'A' ,
01  ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,
'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,
'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,
'A' ,'A' ,'A' ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
0};
unsigned char yyextra[] = {
0,0,1,0,0,0,0,0,
0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,
0};
/*
 * (c) Copyright 1990, OPEN SOFTWARE FOUNDATION, INC.
 * ALL RIGHTS RESERVED
 */
/*
 * OSF/1 Release 1.0
*/
/*
#
# IBM CONFIDENTIAL
# Copyright International Business Machines Corp. 1989
# Unpublished Work
# All Rights Reserved
# Licensed Material - Property of IBM
#
#
# US Government Users Restricted Rights - Use, duplication or
# disclosure restricted by GSA ADP Schedule Contract with IBM Corp.
#
*/
/* @(#)ncform	1.3  com/lib/l,3.1,8951 9/7/89 18:48:47 */
int yylineno =1;
# define YYU(x) x
# define NLSTATE yyprevious=YYNEWLINE
unsigned char yytext[YYLMAX+1];
struct yysvf *yylstate [YYLMAX], **yylsp, **yyolsp;
unsigned char yysbuf[YYLMAX];
unsigned char *yysptr = yysbuf;
int *yyfnd;
extern struct yysvf *yyestate;
int yyprevious = YYNEWLINE;
yylook(){
	register struct yysvf *yystate, **lsp;
	register struct yywork *yyt;
	struct yysvf *yyz;
	int yych, yyfirst;
	struct yywork *yyr;
# ifdef LEXDEBUG
	int debug;
# endif
	unsigned char *yylastch;
	/* start off machines */
# ifdef LEXDEBUG
	debug = 0;
# endif
	yyfirst=1;
	if (!yymorfg)
		yylastch = yytext;
	else {
		yymorfg=0;
		yylastch = yytext+yyleng;
		}
	for(;;){
		lsp = yylstate;
		yyestate = yystate = yybgin;
		if (yyprevious==YYNEWLINE) yystate++;
		for (;;){
# ifdef LEXDEBUG
			if(debug)fprintf(yyout,"state %d\n",yystate-yysvec-1);
# endif
			yyt = yystate->yystoff;
			if(yyt == yycrank && !yyfirst){  /* may not be any transitions */
				yyz = yystate->yyother;
				if(yyz == 0)break;
				if(yyz->yystoff == yycrank)break;
				}
			*yylastch++ = yych = input();
            if (yylastch >= yytext + YYLMAX) {
                fprintf(yyout, "Maximum token length exceeded\n");
                yytext[YYLMAX] = 0;
                return 0;
            }
			yyfirst=0;
		tryagain:
# ifdef LEXDEBUG
			if(debug){
				fprintf(yyout,"char ");
				allprint(yych);
				putchar('\n');
				}
# endif
			yyr = yyt;
			if ( yyt > yycrank){
				yyt = yyr + yych;
				if (yyt <= yytop && yyt->verify+yysvec == yystate){
					if(yyt->advance+yysvec == YYLERR)	/* error transitions */
						{unput(*--yylastch);break;}
					*lsp++ = yystate = yyt->advance+yysvec;
					goto contin;
					}
				}
# ifdef YYOPTIM
			else if(yyt < yycrank) {		/* r < yycrank */
				yyt = yyr = yycrank+(yycrank-yyt);
# ifdef LEXDEBUG
				if(debug)fprintf(yyout,"compressed state\n");
# endif
				yyt = yyt + yych;
				if(yyt <= yytop && yyt->verify+yysvec == yystate){
					if(yyt->advance+yysvec == YYLERR)	/* error transitions */
						{unput(*--yylastch);break;}
					*lsp++ = yystate = yyt->advance+yysvec;
					goto contin;
					}
				yyt = yyr + YYU(yymatch[yych]);
# ifdef LEXDEBUG
				if(debug){
					fprintf(yyout,"try fall back character ");
					allprint(YYU(yymatch[yych]));
					putchar('\n');
					}
# endif
				if(yyt <= yytop && yyt->verify+yysvec == yystate){
					if(yyt->advance+yysvec == YYLERR)	/* error transition */
						{unput(*--yylastch);break;}
					*lsp++ = yystate = yyt->advance+yysvec;
					goto contin;
					}
				}
			if ((yystate = yystate->yyother) && (yyt= yystate->yystoff) != yycrank){
# ifdef LEXDEBUG
				if(debug)fprintf(yyout,"fall back to state %d\n",yystate-yysvec-1);
# endif
				goto tryagain;
				}
# endif
			else
				{unput(*--yylastch);break;}
		contin:
# ifdef LEXDEBUG
			if(debug){
				fprintf(yyout,"state %d char ",yystate-yysvec-1);
				allprint(yych);
				putchar('\n');
				}
# endif
			;
			}
# ifdef LEXDEBUG
		if(debug){
			fprintf(yyout,"stopped at %d with ",*(lsp-1)-yysvec-1);
			allprint(yych);
			putchar('\n');
			}
# endif
		while (lsp-- > yylstate){
			*yylastch-- = 0;
			if (*lsp != 0 && (yyfnd= (*lsp)->yystops) && *yyfnd > 0){
				yyolsp = lsp;
				if(yyextra[*yyfnd]){		/* must backup */
					while(yyback((*lsp)->yystops,-*yyfnd) != 1 && lsp > yylstate){
						lsp--;
						unput(*yylastch--);
						}
					}
				yyprevious = YYU(*yylastch);
				yylsp = lsp;
				yyleng = yylastch-yytext+1;
                if (yyleng >= YYLMAX) {
                    fprintf(yyout, "Maximum token length exceeded\n");
                    yytext[YYLMAX] = 0;
                    return 0;
                }
				yytext[yyleng] = 0;
# ifdef LEXDEBUG
				if(debug){
					fprintf(yyout,"\nmatch ");
					sprint(yytext);
					fprintf(yyout," action %d\n",*yyfnd);
					}
# endif
				return(*yyfnd++);
				}
			unput(*yylastch);
			}
		if (yytext[0] == 0  /* && feof(yyin) */)
			{
			yysptr=yysbuf;
			return(0);
			}
		yyprevious = yytext[0] = input();
		if (yyprevious>0)
			output(yyprevious);
		yylastch=yytext;
# ifdef LEXDEBUG
		if(debug)putchar('\n');
# endif
		}
	}
yyback(p, m)
	int *p;
{
if (p==0) return(0);
while (*p)
	{
	if (*p++ == m)
		return(1);
	}
return(0);
}
	/* the following are only used in the lex library */
yyinput(){
	return(input());
	}

int yyoutput(c)
  int c; {
	output(c);
        return 0;
	}

int yyunput(c)
   int c; {
	unput(c);
        return 0;
        }
