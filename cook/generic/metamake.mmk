!! COMMENTED EXAMPLE
!! use script "macsimus/COOK/configure.sh" to create "metamake"

!! additional option passed to compilation
O = -DSCR

!! define project home
!include "../../home.mmk"

!! define the compiler and compiler-specific options
!include "../../../compile.mmk"

!! project directories, relative to !home (see "home.mmk")
!dir = gen sim sim/busing cook "cook/gcpol/P0a"

!! project dependencies, will be expanded to makefile by "makemake"
!cookgcp0slc+a : ground statics varfile bitfile cputime conjgrad gjlineq asksig \
  simcg simgear interpot setss intrapot sitesite norm xforces linregr \
  constrd rhs simglob pakcp simils siminit main simpot simmeas simmeasx \
  forces simdef \
  elst ewald linklist

