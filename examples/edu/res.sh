#!/bin/bash

while [ x$R == x ] ; do

  echo
  echo "======================= Analysis of results ========================="
  echo "Select one of:"
  echo "  1=show the trajectory"
  echo "  2=convegence profiles"
  echo "  3=radial distribution functions"
  echo "    Hints: use right mouse click to show info"
  echo "           use left mouse to zoom to a rectangle"
  echo "  4=running coordination number"
  echo "  0=quit"
  read

  I='-I%i'
  [ $1 = solution ] && I='-I@i'
  
  case $REPLY in
    1 )
      nice show $1 $I ;;
    2 )
      nice showcp -p5 $1 Tkin Epot $2 ;;
    3 ) 
      export PLOTINIT=W
      [ -e $1.rdf ] && rdfg $1.rdf -g -p ;;
    4 ) 
      export PLOTINIT=W
      [ -e $1.rdf ] && rdfg $1.rdf -c -p ;;
    0 ) exit ;;
    * ) echo "wrong number"
  esac

done
