!! makemake gcc

!define x11
!include "../home.mmk"

!dir = gen sim show

#############################################################################

!plot : plot ploterr xdraw mydraw ground
!show : show ground xdraw
!stereo : stereo ground
!plbmatch : plbmatch ground

!private "private.mmk"
