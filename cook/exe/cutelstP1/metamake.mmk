!! makemake gcc
!! This file was generated by configure.sh
!! See MACSIMUS/generic/metamake.mmk for more details

O = -DSCR
!include "../../../home.mmk"
!! project directories, relative to !home (see "home.mmk")
!dir = gen sim sim/lj cook cook/exe/cutelstP1
!! project dependencies, will be expanded to makefile by "makemake"
!cookceslcP1 : ground statics varfile bitfile cputime conjgrad gjlineq asksig \
  simcg simgear interpot setss intrapot sitesite norm xforces linregr \
  constrd rhs simglob pakcp simils siminit main simpot simmeas simmeasx \
  forces simdef \
  elst ewald linklist

