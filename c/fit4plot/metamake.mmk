!! 
!! to prepare all makefiles required by plot/fit, do (cf. makemake.sh)
!!   makemake linux gcc ; mv makefile makefile.double    # mapped to ctrl-t
!!   makemake linux gcc long ; mv makefile makefile.long # mapped to t
!!   makemake linux gcc high ; mv makefile makefile.high # mapped to T

!include "../../home.mmk"
!include "../../compile.mmk"

!dir = gen c/fit4plot

!if high
!objext=.xo
CC = g++
CPLUS = g++
OPT = -O1 -DPRECISION=8 -Wall
!fit.high : fit ground minimize gjlineq qmin highprec
!else
!if long
!objext=.lo
OPT = -O2 -DPRECISION=2 -Wall
!fit.long : fit ground minimize gjlineq qmin
!else
OPT = -O2 -DPRECISION=1 -Wall
!fit.double : fit ground minimize gjlineq qmin
!myfit.double : myfit ground minimize gjlineq qmin
!endif
!endif
