#!/bin/bash
if [ $# -lt 1 ]; then
  cat <<EOF
Helper script for .che files, to be used by mc or start. Call by:
  CHE= che.sh NAME.che [BLEND-OPTIONS]
CHE value:
(none) = ask
     0 = blend -g5 -r2 BLEND-OPTIONS NAME.che # optimize and show/edit
     1 = nmf.sh NAME.che BLEND-OPTIONS        # normal modes of vibration, in place
     2 = nmf4guest.sh NAME.che BLEND-OPTIONS  # normal modes of vibration
         # in temporary directory and optional plb/*.{plb,mol}
Example:
  CHE=1 che.sh fulerene.che
EOF
  exit
fi

M=$1
shift

if [ x$CHE == x ] ; then
  cat <<EOF
Environment variable CHE not specified, select one of:
  0 = (or just Enter) optimize and show/edit molecule, command:
      blend -g5 -r2 $* $M
  1 : normal mode vibrations (in place), command:
      nmf.sh $M $*
  2 : as above, for guest account (in temporary directory), command:
      nmf4guest.sh $M $*
EOF
  read
  CHE=$REPLY
fi

case $CHE in
  1 )
    echo "CHE=1: starting nmf.sh (in current directory)"
    nmf.sh $M $* ;;
  2 )
    echo "CHE=2: starting nmf4guest.sh (in /tmp, erased after run)"
    nmf4guest.sh $M $* ;;
  * )
    echo "default or CHE=0: starting blend only"
    blend -g5 -r2 $* $M ;;
esac
