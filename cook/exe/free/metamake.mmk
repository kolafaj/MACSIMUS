!! makemake gcc
!! This file was generated by configure.sh
!! See MACSIMUS/generic/metamake for more details

O = -DSCR
!include "../../../home.mmk"
!include "../../../compile.mmk"
!! !dir is relative to !home (see "home.mmk")
!dir = gen sim sim/lj cook cook/exe/free
!cookfree : ground statics varfile bitfile cputime conjgrad gjlineq asksig \
  simcg simgear interpot setss intrapot sitesite norm xforces linregr \
  constrd rhs simglob pakcp simils siminit main simpot simmeas simmeasx \
  forces simdef \
  water

