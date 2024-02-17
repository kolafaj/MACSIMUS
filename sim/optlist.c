/* optionlist for cook 
   and utilities working with .cfg files (plb2cfg, blendmed) 
   to be directly #included after #include "options.h"
*/   

#define x badoption
int optionlist[32] =
/*` a   b c d  e f g h  i j  k l m n o    p q r s t u v w x y          z { | } ~ */
{-9,0,-15,9,0,-1,0,x,0,-9,0,-16,0,2,0,1,229,0,2,0,0,0,3,1,7,0x7fffffff,0,x,x,x,0,100};
/*@ A   B C D  E F G H  I J  K  L M N O   P Q R S T U V W X Y          Z [ \ ] ^ _*/
/*  options -ajq are set to option -_ in the first step */
#undef x

