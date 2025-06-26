#!/bin/bash
if [ $# -lt 1 ]; then
  cat <<EOF
Install MACSIMUS. Call by:
  ./install.sh COMPILER
where COMPILER is the C compiler name:
  gcc = GCC compiler - recommended
        optimization options = -O2 -ffast-math
        -O3 may work, but has not been extensively tested
  cc = generic compiler without options
  icc = Intel compiler - not tested recently, old info:
        faster code than gcc, less reliable on AMD
        optimization option = -O2 (flag -fast is not reliable)
        locale must be set, LANG=C is recommended
  pgcc = Portland Group, Inc. (PGI) compiler - not tested recently, old info:
        bug: random numbers do not work
To implement a new compiler, edit file "sys/compile.mmk".

In linux, the following packages are needed:
  gcc make g++ libncurses-dev libx11-dev

Example (Ubuntu, Debian):
  sudo apt install gcc make g++ libncurses-dev libx11-dev
  ./install.sh gcc
EOF
  exit
fi

# makemake option: now unix and linux are equivalent
UNIX=linux
CC="$1"
export BIN=$PWD/bin/

# OLD, now spaces are allowed
# # no spaces in the path allowed
# if [ "${PWD/ /z}" != "${PWD}" ] ; then
#   echo "Sorry, I cannot handle spaces in PWD=\"$PWD\"!"
#   echo "Please install MACSIMUS in a path not containing spaces"
#   exit
# fi

# check if the selected compiler exists
if ! which "$CC" ; then
  echo "ERROR: compiler $1 not found"
  exit
fi

# checking diff and ed
if ! which diff ; then
  echo "utility diff not found"
  exit
fi
if ! which ed ; then
  echo "utility ed not found"
  exit
fi

# check bin/makeone.sh
if [ ! -e bin/makeone.sh ]; then
  echo "ERROR bin/makeone.sh not found (MACSIMUS incorrectly unzipped?)"
  exit
fi
chmod +x bin/makeone.sh

# using lib64 if exists (probably no longer needed)
# if [ -d /usr/X11R6/lib64 ] ; then
# cat <<EOF | ed graph.mmk
# ,s:/usr/X11R6/lib:/usr/X11R6/lib64:g
# w
# q
# EOF
# echo /usr/X11R6/lib64 used
# fi

# environment for cook/configure.sh
cat > makemake.sh <<EOF
export CC="$CC"
# this file is used by script cook/configure.sh
EOF

cp sys/*.mmk .

# export MACSIMUS home directory for makemake
cat > home.mmk <<EOF
# This setup file is !included by all project metamake.mmk files.
# Non-! statements will become a part of the generated makefile files.

# Read the compiler setup first:
!include sys/compile.mmk

# define MACSIMUS home directory:
EOF
# NEW: spaces etc. allowed in PWD
echo "!home=\"${PWD}\"" >> home.mmk
echo "home.mmk generated"

echo ">>>>>>>>>>>>>>"
echo
echo ">>> int4 check"
cd gen/
$1 -o int4 int4.c
if ! ./int4 ; then
  echo "This is file $PWD/int4.h:"
  echo "----------------------------------------------------"
  cat "$PWD/int4.h"
  echo "----------------------------------------------------"

  echo "Continue anyway (y/N)?"
  read
  if [ "$REPLY" != "y" ] ; then
    exit
  fi
fi
cd ..

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo
echo ">>> Compiling the makemake utility"
cd c/
rm *.o
"$CC" -o makemake makemake.c
if ! ./makemake $UNIX "$CC" ; then
  echo "makemake problem, exit"
  exit
fi
rm makemake # will be recompiled

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo
echo ">>> Compiling general utilities from subdirectory c/"
# do not forget to include these modules to metamake.mmk
for f in \
  makemake \
  tabproc extract comment field \
  mergetab tab start sortcite filttab lemon eqfield sum \
  lattice sumetc blocktab eddata ice naclcryst histogr runsum exp2pow \
  runint binrepl liat lc repl ren line oneline sparsetab
do
  ../bin/makeone.sh $f
done
  echo "$CC" -o ev -O2 -Wall -DCALC=0 ev.c -I../gen -lm -lncurses
  "$CC" -o ev -O2 -Wall -DCALC=0 ev.c -I../gen -lm -lncurses
  echo "$CC" -o evu -O2 -Wall -DCALC=1 ev.c -I../gen -lm -lncurses
  "$CC" -o evu -O2 -Wall -DCALC=1 ev.c -I../gen -lm -lncurses
  mv ev ev.exe evu evu.exe ../bin 2> /dev/null
cd ..

# support of fitting by plot
cd c/fit4plot
touch fitinc.c
rm *.o *.ol *.ox
../../bin/makemake $UNIX "$CC" long
mv makefile makefile.long
../../bin/makemake $UNIX "$CC" high
mv makefile makefile.high
../../bin/makemake $UNIX "$CC"
mv makefile makefile.double
cd ../..

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo
echo ">>> Compiling blend (lj version), pdb, and similar from blend/"
cd blend/
if ! ../bin/makemake $UNIX "$CC" lj ; then
  echo "makemake problem, exit"
  exit
fi

rm *.o
# do not forget to include these modules to metamake.mmk
for f in \
  blend pdb makepept blefilt ramachan pdb2pdb bonds
do
  ../bin/makeone.sh $f
done

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo
echo ">>> Installing bin versions of common force fields from blend/data/"
cd data
rm gromos96.bin charmm19.bin charmm21.bin charmm22.bin
export BLENDPATH=.
../../bin/blend -b gromos96.par
../../bin/blend -b charmm19.par
../../bin/blend -b charmm21.par
../../bin/blend -b charmm22.par
cd ../..

# compiling graphical utilities
cd show/
if ! ../bin/makemake $UNIX "$CC" ; then
  echo "makemake problem, exit"
  exit
fi
rm *.o
for f in \
  show plot stereo plbmatch
do
  ../bin/makeone.sh $f
done
cd ..

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo
echo ">>> Compiling different simulation utilities from util/"
cd util/
if ! ../bin/makemake $UNIX "$CC" ; then
  echo "makemake problem, exit"
  exit
fi
rm *.o
# do not forget to include these modules to metamake.mmk
for f in \
  autocorr cppak molcfg plb2diff plbmsd plbpak rdfg sfourier showcp \
  atomdist spectrum staprt \
  coordn cp2cp cutprt plbcut plbcenter hbonds density densprof plbfilt \
  plbmerge plbinfo plb2nbr plb2cryst shownear plbstack plbreplicate plbsites \
  plb2asc plbcut plb2plb plbconv asc2plb cfg2atm plbsmooth \
  plb2cfg cfg2plb cfg2atm plbframe plbscale plbbox plbtran \
  m2m mol2mol molren cpacol
do
  ../bin/makeone.sh $f
done
cd ..

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo
echo ">>> Compiling the raytracer + related software from ray/"
echo "(more utilities available upon request)"
cd ray/
if ! ../bin/makemake $UNIX "$CC" ; then
  echo "makemake problem, exit"
  exit
fi
rm *.o
for f in \
  ray ppm2ps ppminfo ppmscale ppm2ppm
do
  ../bin/makeone.sh $f
done
cd ..

echo
echo "Comparing (using diff):"
echo "  the contents of \"${PWD}/\" (> instbin.lst)"
echo "  with the expected listing \"sys/bin.lst\""
echo "----------------------------------------------------------------------------"
cd bin
# LC_COLLATE=C added by T.Trnka (../bin.lst is of C locale order)
LC_COLLATE=C \ls -1 > ../instbin.lst
cd ..
# this is because of CygWin (removing .exe)
cat <<EOF | ed instbin.lst
,s:\.exe::g
w
q
EOF
echo
if diff sys/bin.lst instbin.lst > diff.lst ; then
  echo "=== MACSIMUS UTILITIES HAVE BEEN SUCCESSFULLY COMPILED ==="
else
  cat diff.lst
  LT=`fgrep '<' diff.lst | wc -c`
  GT=`fgrep '>' diff.lst | wc -c`
  if [ "$LT" == "0" ] ; then
    echo "=== All registered MACSIMUS have been successfully compiled ==="
  else
    echo "ERRORS DETECTED, symbol '<' above denotes not compiled executables"
  fi
  if [ "$GT" != "0" ] ; then
    echo "=== but there are additional executables marked by '>' (check bin.lst) ==="
  else
    echo "=== there is a problem with diff ==="
  fi
fi
cat <<EOF

>>> Precompile selected cook version (or versions: use string; e.g.. "eE"):
  f = cookfree: free (vacuum, Cartesian) boundary conditions,
                pair sum of Lennard-Jones and point charges
  e = cookewslc: periodic boundary conditions, slab support, linked-cell list,
                 Lennard-Jones, Ewald (point charges), serial version
                 [This is the default needed to run the installation test]
  E = cookewslcP1: as above, thread-parallel version
  c = cookceslc: periodic boundary conditions, slab support, linked-cell list,
                 Lennard-Jones, smoothed cutoff electrostatics, serial version
  C = cookceslcP1: as above, thread-parallel version
  a = fec (all serial versions)
  A = fecEC (all serial and parallel)
  n = nothing (cannot run the installation test)
EOF

read
KEY=$REPLY
[ "$KEY" == "a" ] && KEY=fec
[ "$KEY" == "A" ] && KEY=fecEC
[ "$KEY" == "" ] && KEY=e

FAILED=""

while [ "$KEY" != "" ] ; do

  echo ">>>>>>>>>>>>>>>>>>>"
  echo
  echo ">>> compiling key=${KEY:0:1}"

  case ${KEY:0:1} in
    f )
      pushd cook/exe/free
      rm *.o
      if ! ../../../bin/makemake $UNIX "$CC" ; then
        echo "makemake problem, exit"
        exit
      fi
      ../../../bin/makeone.sh cookfree || FAILED=$FAILED' f=cookfree'
      popd
      ;;
    e )
      pushd cook/exe/ewald
      rm *.o
      if ! ../../../bin/makemake $UNIX "$CC" ; then
        echo "makemake problem, exit"
        exit
      fi
      ../../../bin/makeone.sh cookewslc || FAILED=$FAILED' e=cookewslc'
      popd
      RUNTEST=1
      ;;
    E )
      pushd cook/exe/ewaldP1
      rm *.o
      if ! ../../../bin/makemake $UNIX "$CC" pthread ; then
        echo "makemake problem, exit"
        exit
      fi
      ../../../bin/makeone.sh cookewslcP1 || FAILED=$FAILED' E=cookewslcP1'
      popd
      ;;
    c )
      pushd cook/exe/cutelst
      rm *.o
      if ! ../../../bin/makemake $UNIX "$CC" ; then
        echo "makemake problem, exit"
        exit
      fi
      ../../../bin/makeone.sh cookceslc || FAILED=$FAILED' c=cookceslc'
      popd
      ;;
    C )
      pushd cook/exe/cutelstP1
      rm *.o
      if ! ../../../bin/makemake $UNIX "$CC" pthread ; then
        echo "makemake problem, exit"
        exit
      fi
      ../../../bin/makeone.sh cookceslcP1 || FAILED=$FAILED' C=cookceslc'
      popd
      ;;
    n )
      ;;
    * )
      echo "unknown character in the control string $KEY"
  esac
  KEY=${KEY:1:9}
done

if [ "$FAILED" == "" ] ; then
   echo
   echo "=== All cook* versions have been successfully compiled ==="
   echo
else
  echo "=== THE FOLLOWING COMPILATIONS HAVE FAILED:"
  echo "===$FAILED"
fi

echo "Type Enter to continue installation."
read

cat <<EOF

----------------------------------------------------------------------------
MACSIMUS offers the following commands:
  ev    - command-line-oriented calculator, e.g.:
          ev "exp(-1.1*kcal/298/R)"
  evu   - as above with units, e.g.:
          evu "exp(-1.1[kcal/mol]/298[K]/R)"
  start - extension-based start (MACSIMUS files and more), e.g.:
          start simul.plb
NB: render option in show requires start installed
Get help by running them without parameters.
* Installation will copy the executables to ~/bin/
  and initialization .files to ~/.
* Not-installation will remove the temporary files
  $PWD/bin/{ev,evu,start}.
Install ev,evu,start (y/N)?
EOF

read
if [ "$REPLY" == "y" ] ; then
  echo "copying customizable initialization files .evdata .evudata .startdata to ~/"
  pushd sys
  cp -ip .evdata .evudata .startdata ~/
  popd
else
  echo "removing executables $PWD/bin/{ev,evu,start}"
  rm bin/ev bin/evu bin/start
fi

if [ -d ~/.config/mc ] ; then
  echo
  echo "INFO: ~/.config/mc found, it looks you are using Midnight Commander"
  if [ -e ~/.config/mc/mc.ext ] ; then
    echo "INFO: extension association file ~/.config/mc/mc.ext found!"
    echo "To activate the MACSIMUS extensions (incl. common graphic files via"
    echo "display/ImageMagick), you should install macsimus/mc.ext:"
    echo "  a = Append macsimus/mc.ext to ~/.config/mc/mc.ext"
    echo "  r = Replace ~/.config/mc/mc.ext by macsimus/mc.ext (with backup~)"
    echo "> n = Nothing (default)"
  else
    echo "INFO: no extension association file ~/.config/mc/mc.ext found!"
    echo "To activate the MACSIMUS extensions (incl. common graphic files via"
    echo "display/ImageMagick), you should install macsimus/mc.ext:"
    echo "  c = Copy macsimus/mc.ext to ~/.config/mc/mc.ext"
    echo "> n = Nothing (default)"
  fi

  read
  case $REPLY in
    a | A ) cat mc.ext >>  ~/.config/mc/mc.ext ;;
    c | C | r | R )
      [ -e  ~/.config/mc/mc.ext ] && mv  ~/.config/mc/mc.ext ~/.config/mc/mc.ext~
      cp mc.ext ~/.config/mc/mc.ext
      ;;
    * ) ;;
  esac
fi

if [ -d ~/.mc ] ; then
  echo
  echo "INFO: ~/.mc found, it looks you are using an OBSOLETE Midnight Commander"
  if [ -e ~/.mc/bindings ] ; then
    echo "INFO: extension association file ~/.mc/bindings found!"
    echo "To activate the MACSIMUS extensions (incl. common graphic files via"
    echo "display/ImageMagick), you should install macsimus/mc.ext:"
    echo "  a = Append macsimus/mc.ext to ~/.mc/bindings"
    echo "  r = Replace ~/.mc/bindings by macsimus/mc.ext (with backup~)"
    echo "> n = Nothing (default)"
  else
    echo "INFO: no extension association file ~/.mc/bindings found!"
    echo "To activate the MACSIMUS extensions (incl. common graphic files via"
    echo "display/ImageMagick), you should install macsimus/mc.ext:"
    echo "  c = Copy macsimus/mc.ext to ~/.mc/bindings"
    echo "> n = Nothing (default)"
  fi

  read
  case $REPLY in
    a | A ) cat sys/mc.ext >>  ~/.mc/bindings ;;
    c | C | r | R )
      [ -e  ~/.mc/bindings ] && mv ~/.mc/bindings ~/.mc/bindings~
      cp sys/mc.ext ~/.mc/bindings
      ;;
    * ) ;;
  esac
fi

echo "MACSIMUS directory: ${PWD}"
cat sys/FINISHINST.txt

if [ "$RUNTEST" == "1" ]; then
  X=' later'
  echo "Would you like to start a simple test of MACSIMUS now (y/N)?"
  read
  if [ "$REPLY" == "y" ] ; then
    export BLENDPATH="${PWD}/blend/data"
    export PATH="${PWD}/bin:$PATH"
    pushd examples/inst
    echo "NB: the environment has been set within the script."
    ./test.sh
    X=" again"
    popd
  fi

  echo
  echo "You can run the test$X as ./test.sh from directory"
  echo "${PWD}/examples/inst/"
  echo

  echo "MORE: To show normal mode vibrations of selected molecules:"
  pushd blend/che
  cat <<EOF
  cd "$PWD"
Then EITHER from the command prompt as:
  CHE=2 che.sh MOLECULE.che
OR using the Midnight Commander (if the MACSIMUS associations are installed)
  CHE=2 mc
and double click a molecule.
If minimized, press button [finish] or hotkey .
EOF
  popd
  echo "MACSIMUS directory: ${PWD}"
  cat sys/FINISHINST.txt
fi
