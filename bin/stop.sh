#!/bin/bash
if [ $# -lt  1 ]; then
  echo "Stop cook* simulation(s) by creating the .stp file"
  echo "  stop.sh {SIMNAME|LOC-FILE|PRT-FILE} [...]"
  echo "Example (stop all):"
  echo "  stop.sh *.loc"
  echo "See also:"
  echo "  stp.sh stpall.sh"
  exit
fi

for f in $* ; do
  N=${f/%.loc/} 
  N=${N/%.prt/} 
  if [ -e $N.loc ] ; then
    echo $N.loc found - touch $N.stp
    touch $N.stp
  else
    echo $N.loc not found - ignored
  fi
done
