#!/bin/bash
if [ $# -lt 1 ]; then 
  echo "Call by:"
  echo "  ./yes.sh text [default]"
  echo "[aAyYzZ] accepted as yes"
  echo "Example:"
  echo "  yes \"erase xxx.* (A/n)\" y"
  exit
fi
echo $1
read
R=$REPLY
[ x$R == x ] && R=$2
[ "x${R:0:1}" == xa ] && exit 0
[ "x${R:0:1}" == xA ] && exit 0
[ "x${R:0:1}" == xy ] && exit 0
[ "x${R:0:1}" == xY ] && exit 0
[ "x${R:0:1}" == xz ] && exit 0
[ "x${R:0:1}" == xZ ] && exit 0
exit 1
