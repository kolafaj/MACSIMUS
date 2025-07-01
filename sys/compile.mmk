!! *** part of MACSIMUS - basic configuration file for makemake ***
!! See c/makemake.c for more info.
!!
!! This file is !included (by command !include) from files `metamake'
!! in subdirectories of the MACSIMUS home directory during processing
!! by utility `makemake'.
!! - If graphics is needed, switch x11 should be !defined.
!! - Preprocessor switches (-D) should be passed in makefile variable $O
!!
!! For a particular purpose, this file can be as simple as:
!!   CC = <C compiler name>
!!   LIBOPT = <linker options (incl. the X11 library)>
!!   OPT = -DX11 $O <more C compiler options>

!! use OPT = ... -DC99 in case of problems with doubly defined fmin,fmax

!! *** X11 graphics support ***
!if x11
!! This preprocessor switch is needed by MACSIMUS code to use X11 graphics:
OPTX11 := -DX11
!! This is the X library:
LIBOPTX11 := -lX11
!if gui
!! OBSOLETE: gtk-based GUI for show, plot, blend
OPTX11 += -DGUI
!endif
!! In case of X11 linking problems, one of these may work instead:
!!   LIBOPTX11 -L/usr/X11R6/lib/
!!   LIBOPTX11 -L/usr/X11R6/lib64/
!endif

!! default compiling options
!if shm
!! OBSOLETE
EXTOPT := -DSHM
!endif
OPT = $O $(OPTX11) $(EXTOPT)

!! default linking options
!if pthread
EXTLIBOPT := -lpthread
!endif
!if static
EXTLIBOPT := -static
!endif
!if pthread & static
EXTLIBOPT := -lpthread -static
!endif
LIBOPT := $(LIBOPTX11) -lm ${EXTLIBOPT}

!! cancel the default for CC to force a user to select a compiler
CC = undefined-compiler

!! *** generic C compiler ***
!if cc
CC := cc
!endif

!! *** GNU C compiler (recommended) ***
!if gcc
CC := gcc
OPT = -O2 -ffast-math -Wall $O $(OPTX11) $(EXTOPT)
!endif

!! *** GNU C compiler (obsolete - now gcc colorized by default) ***
!if colorgcc
CC := colorgcc
OPT = -O2 -ffast-math -Wall $O $(OPTX11) $(EXTOPT)
!endif

!! *** Intel C compiler ***
!if icc
CC := icc
OPT = -fast $O $(OPTX11) $(EXTOPT)
!! OPT = -O2 $O $(OPTX11) $(EXTOPT)
!! flag -fast sometimes causes "Trace/breakpoint trap" just after a run starts
!endif

!! *** Portland Group C compiler ***
! problems!
!if pgcc
CC := pgcc
OPT = -O2 -fast $O $(OPTX11) $(EXTOPT)
!endif

!! *** SGI (obsolete) ***
!if sgi
CC := cc
OPT = -fullwarn -O2 $O $(OPTX11) $(EXTOPT)
!endif

!! *** Digital (obsolete) ***
!if digital
OPT = -O2 $O -std1 -warnprotos -w0 $(OPTX11) $(EXTOPT)
# digital cc bug: some prototypes incorrectly reported as missing
!endif

# do not export these modules to the list of dependencies for make
!ignore = shmalloc.c xxxdef.c xxxplus.c xxxdip.c vircdiat.c seasr.c shakev30.c shakev4.c shakev9.c spctcf.c
