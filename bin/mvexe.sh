#!/bin/bash
# move either $2.exe [$3.exe ...] or $2 [$3 ...] (whichever exists) to $1
# part of MACSIMUS installation script compile.sh

# WARNING: because of stupid and every year changing way how .exe
# files are treated under CygWin (unclear meaning of "if [ -e FILE ]"
# and "if [ -e FILE.exe ]"), this fool-proof version just tries to
# move both FILE and FILE.exe.  If everything goes fine, just one of
# them exists, and the other is ignored with confusing but harmless
# error message.

DIR=$1

echo "CygWin hack:"
echo "  Under Unix, ignore messages as XXX.exe does not exist"
echo "  Under CygWin, ignore messages as XXX does not exist"

while [ $# -ge 2 ] ; do
  shift

#  if [ -e $1 ] ; then
#    # UNIX
#    echo mv $1 ${DIR}
#    mv $1 ${DIR}
#  else
#    if [ -e $1.exe ]; then
#    # CygWin
#      echo mv $1.exe ${DIR}
#      mv $1.exe ${DIR}
#    else
#      echo no $1 nor $1.exe
#    fi
#  fi

# "if [ -e $1 ]"  does not work under CygWin, doing it the hard way
# (ignore all warnings)

  echo mv $1 $1.exe ${DIR}
  mv $1 $1.exe ${DIR}

done
