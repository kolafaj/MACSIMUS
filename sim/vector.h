#ifndef _H_VECTOR_
#define _H_VECTOR_

/*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                              vectors                              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vector macros mnemonics: V = <vector> | <operator vector>, O = <operator>
example: VO(v,*=2) multiplies vector v by 2
*/

#include "simopt.h"

#ifdef FLOAT
typedef float real;
#else
typedef double real;
#endif

#ifdef TWODIM
#include "vector2d.h"
#else
#include "vector3d.h"
#endif

#endif
