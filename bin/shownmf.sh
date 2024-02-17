#!/bin/bash
if [ $# -le  0 ]; then
  echo "Call by:"
  echo "  shownmf.sh MOLNAME[.nmf]"
  echo "MOLNAME.nm####.plb files must be present, too"
  echo "See also:"
  echo "  nmf.sh nmf4guest.sh"
  exit
fi

N=${1/%.nmf/}
cat $N.nmf
# header of 3 lines assumed in $N.nmf, n=L-3
L=`wc -l $N.nmf | cut -d' ' -f1`
if (( $L < 9 )) ; then
  echo "$N.nmf is not long enough (incl. assumed 3 header lines)"
  exit
fi

# find the first nonzero vibration
# general molecule - skip 3 translations + 3 rotations
H=-h6
# linear molecule
(( $L == 9 )) && H=-h5
# PATCH for CO2, which is linear
[ "$N" == "CO2" ] && H=-h5

[ -e kill.sh ] && rm kill.sh
export GSHOWKILL=1
export SHOWNMF=$N.nmf
echo show $H -\\ -d30 $N $N.nm%04d.plb -I\$wwir
show $H -\\ -d30 $N $N.nm%04d.plb -I\$wwir
[ -e kill.sh ] && sh kill.sh
