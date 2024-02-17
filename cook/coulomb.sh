
MODULE="$MODULE elst ewald"

echo
echo "===== Select algorithm for ELECTROSTATIC calculations:"
echo " -3 = Gaussian charges (on spring or nonpolarizable) + Ewald"
echo ">-2 = The same as -1 with more accurate erfc for 1-2 and 1-3 interactions"
echo "      For models with partial charges close together (DEFAULT)"
echo " -1 = Ewald summation"
echo "      All r-space terms (erfc) are calculated by hyperbolic splines"
echo "  0 = Cut-and-shift approximation (MACSIMUS style)"
echo "  2 = Deprecated: COULOMB=0 implemented by splines"
echo "  3 = Cut-and-shift approximation by Fennell, Gezelter, JCP 124, 234104 (2006)"
echo "      Not suitable for models with large partial charges very close together"
echo "NOTE: Consider also switches QQTAB and SPLINES, see ../cook/generic/simopt.h"
read

COULOMB=$REPLY
if [ "$COULOMB" = "" ] ; then
  COULOMB="-2"
fi
echo "#define COULOMB" ${COULOMB} >> simopt.h
case "$COULOMB" in
-3 )
   PROJECT=$PROJECT"-Ewald+GC"
   NAME=$NAME"gc"
   echo "#define SPLINE 3" >> simopt.h
   echo "#define GAUSSIANCHARGES" >> simopt.h

   echo
   echo "===== Select Gaussian charges model (charge-charge table) type:"
   echo "  0 = Central atoms have zero charge (all charges on a spring, as BK3 water)"
   echo "  1 = There are charges on the atom and on a spring (e.g., KB ions)"
   echo "> 2 = Optimized for both cases (e.g., ions in water) or nonpolarizable (DEFAULT)"
   read
   QQTAB=$REPLY
   if [ "$QQTAB" = "" ] ; then
     QQTAB="2"
   fi
   echo "#define QQTAB" ${QQTAB} >> simopt.h

   ;;
-1 | -2 )
   PROJECT=$PROJECT"-Ewald"
   NAME=$NAME"ew"
   echo "#define SPLINE -2" >> simopt.h
   ;;
 0 )
   PROJECT=$PROJECT"-cutelst"$COULOMB
   NAME=$NAME"ce"
   ;;
 2 )
   PROJECT=$PROJECT"-cutelst"$COULOMB
   NAME=$NAME"ce"
   echo "#define SPLINE 2" >> simopt.h
   ;;
 3 ) NAME=$NAME"fg"
   echo "#define SPLINE 2" >> simopt.h
   ;;
 * ) echo "ERROR: unsupported algorithm"; exit ;;
esac

if [ "$COULOMB" != "0" ] ; then
  echo "WARNING:"
  echo "Real-space charge-charge interactions are calculated by splines."
  echo "You may need to check the automatically generated simopt.h"
fi
