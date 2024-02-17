/***********************************************************************
 * $Author: markv $
 * $Revision: 1.2 $
 * $Date: 88/09/12 12:53:47 $
 * $Log:	color.c,v $
 * Revision 1.2  88/09/12  12:53:47  markv
 * Fixed problem in LookupColorbyName, had return ; and return(0).
 * [ Thank you lint! ]
 * 
 * Revision 1.1  88/09/11  11:00:37  markv
 * Initial revision
 * 
 ***********************************************************************/

#include <stdio.h>
#include "defs.h"

#define		NCOLORS		(142)

typedef struct t_color_entry {
	char *	ce_name ;
	Vec 	ce_color ;
} ColorEntry ;

#define LESS_THAN -1
#define GREATER_THAN 1
#define EQUAL_TO 0

/*
 * Note: These colors must be in sorted order, because we binary search
 * for them.
 *
 * They were swiped from the X-11 distribution.  Sorry....
 */

ColorEntry Colors[] = {
	{"Aquamarine", {.439216, .858824, .576471}},
	{"Black", {0, 0, 0}},
	{"Blue", {0, 0, 1}},
	{"BlueViolet", {.623529, .372549, .623529}},
	{"Brown", {.647059, .164706, .164706}},
	{"CadetBlue", {.372549, .623529, .623529}},
	{"Coral", {1, .498039, 0}},
	{"CornflowerBlue", {.258824, .258824, .435294}},
	{"Cyan", {0, 1, 1}},
	{"DarkGreen", {.184314, .309804, .184314}},
	{"DarkOliveGreen", {.309804, .309804, .184314}},
	{"DarkOrchid", {.6, .196078, .8}},
	{"DarkSlateBlue", {.419608, .137255, .556863}},
	{"DarkSlateGray", {.184314, .309804, .309804}},
	{"DarkSlateGrey", {.184314, .309804, .309804}},
	{"DarkTurquoise", {.439216, .576471, .858824}},
	{"DimGray", {.329412, .329412, .329412}},
	{"DimGrey", {.329412, .329412, .329412}},
	{"Firebrick", {.556863, .137255, .137255}},
	{"ForestGreen", {.137255, .556863, .137255}},
	{"Gold", {.8, .498039, .196078}},
	{"Goldenrod", {.858824, .858824, .439216}},
	{"Gray", {.752941, .752941, .752941}},
	{"Green", {0, 1, 0}},
	{"GreenYellow", {.576471, .858824, .439216}},
	{"Grey", {.752941, .752941, .752941}},
	{"IndianRed", {.309804, .184314, .184314}},
	{"Khaki", {.623529, .623529, .372549}},
	{"LightBlue", {.74902, .847059, .847059}},
	{"LightGray", {.658824, .658824, .658824}},
	{"LightGrey", {.658824, .658824, .658824}},
	{"LightSteelBlue", {.560784, .560784, .737255}},
	{"LimeGreen", {.196078, .8, .196078}},
	{"Magenta", {1, 0, 1}},
	{"Maroon", {.556863, .137255, .419608}},
	{"MediumAquamarine", {.196078, .8, .6}},
	{"MediumBlue", {.196078, .196078, .8}},
	{"MediumForestGreen", {.419608, .556863, .137255}},
	{"MediumGoldenrod", {.917647, .917647, .678431}},
	{"MediumOrchid", {.576471, .439216, .858824}},
	{"MediumSeaGreen", {.258824, .435294, .258824}},
	{"MediumSlateBlue", {.498039, 0, 1}},
	{"MediumSpringGreen", {.498039, 1, 0}},
	{"MediumTurquoise", {.439216, .858824, .858824}},
	{"MediumVioletRed", {.858824, .439216, .576471}},
	{"MidnightBlue", {.184314, .184314, .309804}},
	{"Navy", {.137255, .137255, .556863}},
	{"NavyBlue", {.137255, .137255, .556863}},
	{"Orange", {.8, .196078, .196078}},
	{"OrangeRed", {1, 0, .498039}},
	{"Orchid", {.858824, .439216, .858824}},
	{"PaleGreen", {.560784, .737255, .560784}},
	{"Pink", {.737255, .560784, .560784}},
	{"Plum", {.917647, .678431, .917647}},
	{"Red", {1, 0, 0}},
	{"Salmon", {.435294, .258824, .258824}},
	{"SeaGreen", {.137255, .556863, .419608}},
	{"Sienna", {.556863, .419608, .137255}},
	{"SkyBlue", {.196078, .6, .8}},
	{"SlateBlue", {0, .498039, 1}},
	{"SpringGreen", {0, 1, .498039}},
	{"SteelBlue", {.137255, .419608, .556863}},
	{"Tan", {.858824, .576471, .439216}},
	{"Thistle", {.847059, .74902, .847059}},
	{"Turquoise", {.678431, .917647, .917647}},
	{"Violet", {.309804, .184314, .309804}},
	{"VioletRed", {.8, .196078, .6}},
	{"Wheat", {.847059, .847059, .74902}},
	{"White", {.988235, .988235, .988235}},
	{"Yellow", {1, 1, 0}},
	{"YellowGreen", {.6, .8, .196078}},
	{"aquamarine", {.439216, .858824, .576471}},
	{"black", {0, 0, 0}},
	{"blue", {0, 0, 1}},
	{"blue_violet", {.623529, .372549, .623529}},
	{"brown", {.647059, .164706, .164706}},
	{"cadet_blue", {.372549, .623529, .623529}},
	{"coral", {1, .498039, 0}},
	{"cornflower_blue", {.258824, .258824, .435294}},
	{"cyan", {0, 1, 1}},
	{"dark_green", {.184314, .309804, .184314}},
	{"dark_olive_green", {.309804, .309804, .184314}},
	{"dark_orchid", {.6, .196078, .8}},
	{"dark_slate_blue", {.419608, .137255, .556863}},
	{"dark_slate_gray", {.184314, .309804, .309804}},
	{"dark_slate_grey", {.184314, .309804, .309804}},
	{"dark_turquoise", {.439216, .576471, .858824}},
	{"dim_gray", {.329412, .329412, .329412}},
	{"dim_grey", {.329412, .329412, .329412}},
	{"firebrick", {.556863, .137255, .137255}},
	{"forest_green", {.137255, .556863, .137255}},
	{"gold", {.8, .498039, .196078}},
	{"goldenrod", {.858824, .858824, .439216}},
	{"gray", {.752941, .752941, .752941}},
	{"green", {0, 1, 0}},
	{"green_yellow", {.576471, .858824, .439216}},
	{"grey", {.752941, .752941, .752941}},
	{"indian_red", {.309804, .184314, .184314}},
	{"khaki", {.623529, .623529, .372549}},
	{"light_blue", {.74902, .847059, .847059}},
	{"light_gray", {.658824, .658824, .658824}},
	{"light_grey", {.658824, .658824, .658824}},
	{"light_steel_blue", {.560784, .560784, .737255}},
	{"lime_green", {.196078, .8, .196078}},
	{"magenta", {1, 0, 1}},
	{"maroon", {.556863, .137255, .419608}},
	{"medium_aquamarine", {.196078, .8, .6}},
	{"medium_blue", {.196078, .196078, .8}},
	{"medium_forest_green", {.419608, .556863, .137255}},
	{"medium_goldenrod", {.917647, .917647, .678431}},
	{"medium_orchid", {.576471, .439216, .858824}},
	{"medium_sea_green", {.258824, .435294, .258824}},
	{"medium_slate_blue", {.498039, 0, 1}},
	{"medium_spring_green", {.498039, 1, 0}},
	{"medium_turquoise", {.439216, .858824, .858824}},
	{"medium_violet_red", {.858824, .439216, .576471}},
	{"midnight_blue", {.184314, .184314, .309804}},
	{"navy", {.137255, .137255, .556863}},
	{"navy_blue", {.137255, .137255, .556863}},
	{"orange", {.8, .196078, .196078}},
	{"orange_red", {1, 0, .498039}},
	{"orchid", {.858824, .439216, .858824}},
	{"pale_green", {.560784, .737255, .560784}},
	{"pink", {.737255, .560784, .560784}},
	{"plum", {.917647, .678431, .917647}},
	{"red", {1, 0, 0}},
	{"salmon", {.435294, .258824, .258824}},
	{"sea_green", {.137255, .556863, .419608}},
	{"sienna", {.556863, .419608, .137255}},
	{"sky_blue", {.196078, .6, .8}},
	{"slate_blue", {0, .498039, 1}},
	{"spring_green", {0, 1, .498039}},
	{"steel_blue", {.137255, .419608, .556863}},
	{"tan", {.858824, .576471, .439216}},
	{"thistle", {.847059, .74902, .847059}},
	{"turquoise", {.678431, .917647, .917647}},
	{"violet", {.309804, .184314, .309804}},
	{"violet_red", {.8, .196078, .6}},
	{"wheat", {.847059, .847059, .74902}},
	{"white", {.988235, .988235, .988235}},
	{"yellow", {1, 1, 0}},
	{"yellow_green", {.6, .8, .196078}}
} ;

int
LookupColorByName(name, color)
 char * name ;
 Vec color ;
{
	int rc ;
	rc = BinarySearch(name, 0, NCOLORS - 1 , Colors) ;
	if (rc < 0) {
		return(0) ;
	}

	VecCopy(Colors[rc].ce_color, color) ;
	return 1 ;
}


int 
BinarySearch(name, l, h, array)
 char * name ;
 int l, h ;
 ColorEntry array[] ;
{
	int m, rc ;
	if (l > h)
		return(-1) ;
	
	m = (l + h) / 2 ;

	rc = strcmp(name, array[m].ce_name) ;
	if (rc == 0)
		return m ;
	else if (rc < 0)
		return BinarySearch(name, l, m-1, array) ;
	else
		return BinarySearch(name, m + 1, h, array) ;
}
