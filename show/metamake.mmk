!! makemake gcc

!include "../home.mmk"

!define x11
# O = -DSCR
!include "../compile.mmk"

!dir = gen sim show

#############################################################################

!plot : plot ploterr xdraw mydraw ground
!show : show ground xdraw
!stereo : stereo ground
!plbmatch : plbmatch ground

!private "private.mmk"
