!! makemake gcc

!include "../home.mmk"

! required by cygwin:
OPT += -DSCR -I/usr/include/ncurses

!dir = gen c

!! several modules
!tabproc: tabproc ground
!lr : lr ground linregr
!ice : ice ground statics varfile

!! one module (can be also compiled directly by the 1st line of the SOURCE.c)
!makemake : LIBOPT= : makemake
!mergetab : LIBOPT= -lm : mergetab
!tab : LIBOPT= : tab
!start : LIBOPT= : start
!sortcite : LIBOPT= : sortcite
!lattice : lattice
!filttab : LIBOPT= -lm : OPT=-O0 : filttab
!sumetc : sumetc
!lemon : LIBOPT= : lemon
!eqfield : LIBOPT= : eqfield
!sum : LIBOPT= : sum
!blocktab : blocktab
!eddata : eddata
!histogr : histogr
!runsum : runsum
!runint : runint
!binrepl : binrepl
!grouptab : grouptab
!lc : LIBOPT= : lc
!liat : LIBOPT= : liat
!ren : ren
!naclcryst : naclcryst

!private "private.mmk"
