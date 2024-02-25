!! makemake linux

!include "../home.mmk"

!dir = gen ray

!! several modules
!ray : data main sphere vector shade trace intersec \
        screen pic poly bound error pqueue cone \
        color tri getopt nff tokens

!! one module (can be also compiled directly by the 1st line of the SOURCE.c)
!ppm2ps : ppm2ps
!ppminfo : ppminfo
!ppmscale : ppmscale
!ppm2ppm : LIBOPT=: ppm2ppm
