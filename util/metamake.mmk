!! makemake linux gcc
!include "../home.mmk"

# -DReverse causes some utilities to support reversed endian, now rarely needed
# -DSCR extra code for scrolling screens of output (internal `less')
# O = -DSCR -DReverse

!dir = gen sim util

#############################################################################

!! several modules
!autocorr : autocorr ground varfile statics
!cppak : cppak pakcp bitfile ground
!molcfg : molcfg ground
!plb2diff : plb2diff ground linregr sdfit varfile statics
!plbmsd : plbmsd ground linregr
!plbpak : plbpak bitfile ground
!rdfg : rdfg ground varfile
!rdfgv1 : rdfgv1 ground varfile
!rdfconv : rdfconv ground varfile
!sfourier : sfourier ground fft
!showcp : showcp pakcp bitfile ground varfile statics linregr
!tcov : tcov pakcp bitfile ground
!spectrum : spectrum fft ground solve2
!staprt : staprt ground varfile statics
!cfg2plb : cfg2plb ground varfile
!cfg2asc : cfg2asc ground varfile
!asc2cfg : asc2cfg ground varfile
!cfg2atm : cfg2atm ground varfile
!plb2cfg : plb2cfg ground varfile
!cfgconv : cfgconv ground varfile

!! one module (can be also compiled directly by the 1st line of the SOURCE.c)
!atomdist : atomdist
!coordn: coordn
!cp2cp: cp2cp
!cutprt: cutprt
!hbonds: hbonds
!density: density
!densprof: densprof
!plbfilt: plbfilt
!plbmerge: plbmerge
!plbinfo: plbinfo
!plb2nbr: plb2nbr
!plb2cryst: plb2cryst
!shownear: shownear
!plbtran : plbtran
!plbstack: plbstack
!plbreplicate: plbreplicate
!plbsites: plbsites
!plbscale : plbscale

!plb2asc: LIBOPT= : plb2asc
!plbcut: LIBOPT= : plbcut
!plb2plb: LIBOPT= : plb2plb
!plbconv: LIBOPT= : plbconv
!asc2plb: LIBOPT= : asc2plb
!plbframe: LIBOPT= : plbframe
!plbbox : LIBOPT= : plbbox

!private private.mmk
