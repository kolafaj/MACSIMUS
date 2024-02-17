#!/bin/bash
if [ $# -lt  1 ]; then
  echo "Make NAME and move it to macsimus/bin, to be used by install.sh"
  echo "Compatible with both Unix and CywWin (NAME.exe)"
  echo "Call by: makeone.sh NAME"
  exit
fi
# Note that the meaning of command "if [ -e FILE.exe ]" is unclear in CygWin,
# therefore we try to move both FILE and FILE.exe
make $1
mv $1 $1.exe ../bin 2> /dev/null
