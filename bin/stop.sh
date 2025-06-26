#!/bin/bash
if [ $# -lt  1 ]; then
  echo "Stop cook* simulation(s) by creating the SIMNAME.stp file(s)"
  echo "  stop.sh [-w] {SIMNAME|LOC-FILE|PRT-FILE} [...]"
  echo "-w = wait until all loc-files disappear"
  echo "Example (stop all):"
  echo "  stop.sh *.loc"
  echo "See also:"
  echo "  stp.sh stpall.sh"
  exit
fi

W=$1

[ "$W" == "-w" ] && shift

for f in "$@" ; do
  N="${f/%.loc/}"
  N="${N/%.prt/}"
  if [ -e "$N.loc" ] ; then
    echo "$N.loc" found - touch "$N.stp"
    touch "$N.stp"
  else
    echo "$N.loc" not found - ignored
  fi
done

[ "$W" == "-w" ] || exit

echo "waiting until all loc-files disappear and erasing the stp-files"

while true ; do
  sleep 3
  S="STOP"
  for f in "$@" ; do
    N="${f/%.loc/}"
    N="${N/%.prt/}"
    [ -e "$N.loc" ] && S="CONT"
  done
  [ $S == "STOP" ] && exit
done

for f in "$@" ; do                                                                                                      
  N="${f/%.loc/}"
  N="${N/%.prt/}"
  rm "$N.stp"
done
