
MODULE="$MODULE elst ewald"

cat <<EOF
===== Select algorithm for ELECTROSTATIC calculations:
 -3 = Gaussian charges (on spring or nonpolarizable) + Ewald
>-2 = The same as -1 with more accurate erfc for 1-2 and 1-3 interactions
      For models with partial charges close together (DEFAULT)
 -1 = Ewald summation
      All r-space terms (erfc) are calculated by hyperbolic splines
  0 = Cut-and-shift approximation (MACSIMUS style)
  2 = As above via quadratic splines (slightly faster but less accurate)
  3 = Cut-and-shift approximation by Fennell, Gezelter, JCP 124, 234104 (2006)
      Not suitable for models with large partial charges very close together
NOTE: Consider also switches QQTAB and SPLINES, see ../cook/generic/simopt.h
EOF
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
