#!/bin/bash
if [ $# -lt  1 ]; then
  echo "Make NAME and move it to macsimus/bin, to be used by install.sh"
  echo "Compatible with both Unix and CygWin (NAME.exe)"
  echo "Call by: makeone.sh NAME"
  exit
fi
# Note that the meaning of command "if [ -e FILE.exe ]" is unclear in CygWin,
# therefore we try to move both FILE and FILE.exe
make $1 || exit 1
[ -e $1.exe ] && mv $1.exe $BIN
[ -e $1 ] && mv $1 $BIN
exit 0
