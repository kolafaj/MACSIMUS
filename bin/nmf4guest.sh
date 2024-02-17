#!/bin/bash
if [ $# -le  0 ]; then
  echo "Normal mode frequences for students (and Midnight Commander). Call by:"
  echo "  nmf4guest.sh MOLNAME[.che] [BLEND-OPTIONS]"
  echo "If che.sh is bound to .che in mc.ext, then use as:"
  echo "  CHE=2 mc"
  echo "Clicking a che-file will optimize the molecule, then clicking [Go]"
  echo "or hotkey . it will calculate and show the vibrations."
  echo "If available, plb/MOLNAME.plb and plb/MOLNAME.mol will be used."
  echo "Works in a temporary directory which is erased after exit."
  echo "BUG: Too slow for more than 100 atoms!"
  exit
fi

D=`mktemp -d`
echo "Temporary directory = $D"
N=`basename $1`
N=${N/%.che/}

if [ $N == PEPTIDE ] ; then
  echo "-----------------------------------------------------------------------------"
  echo "Enter sequence of aminoacids using 3-letter codes, lowercase, space-separated"
  echo "  For protonated acids use: asph gluh"
  echo "  For neutral amines use: hisn lysn argn"
  echo "  Can be abbreviated using 1-letter code (ala ala lys -> A2K)"
  echo "  N-teminus will be acetylated"
  echo "  C-teminus will be methylated"
  echo "Example:"
  echo "  A3PT hisn"
  echo "will make peptide Acetyl-ALA-ALA-ALA-PRO-THR-HIS(neutral)-Methyl"
  read
  PEPT=$REPLY
  echo "Enter force field:"
  echo "  1=charmm21 (default)"
  echo "  2=charmm22"
  echo "  3=gromos96"
  read
  FF=charmm21
  [ $REPLY == 2 ] && FF=charmm22
  [ $REPLY == 3 ] && FF=gromos96
  E=-e40
  makepept $FF ace $PEPT ct1 > PEPTIDE.che
fi

if [ ! -e $N.che ] ; then
  echo "cannot find $N.che"
  rmdir $D
  exit
else
  cp $N.che $D
fi

if [ -e plb/$N.plb ] ; then
  cp plb/$N.plb $D
  cp plb/$N.mol $D
fi

if ! cd $D ; then
  echo cannot cd $D
  exit
fi

if [ "${BLENDPATH}" == "" ] ; then
  export BLENDPATH=/home/jiri/macsimus/blend/data
fi

blend $E -i16 -N -M-1 -P999 -m99999 -p -g10 $2 $N.che
nice shownmf.sh $N.nmf
