!! makemake gcc lj ; rm *.o blend; make blend
!! makemake gcc busing ; rm *.o blendbus; make blendbus
!! makemake gcc lj polar ; rm *.o blendpol; make blendpol
!! makemake gcc busing polar ; rm *.o blendpolbus; make blendpolbus
!! makemake gcc metal ; rm *.o blendmetal; make blendmetal
!! See below for more names! Adding rep8 or rep12 to exp6,busing,buck
!! adds -D/r^8 or -D/r^12 terms
!!
!! Synopsis:
!!   makemake COMPILER SITESITE MODIFIER [MODIFIER ..]
!! COMPILER:
!!   cc = general linux/unix
!!   gcc = GNU C
!!   colorgcc = gcc with a colorizer (OLD - now in gcc)
!!   icc = Intel
!! SITESITE (site-site potential and combining rule):
!!   lj = Lennard-Jones, see MACSIMUS/sim/lj
!!   buck = Buckingham version of exp-6, see MACSIMUS/sim/buck
!!   exp6 = LJ-like version of exp-6 (UNFINISHED), see MACSIMUS/sim/exp6
!!   busing = Busing version of exp-6, see MACSIMUS/sim/busing
!!   goldErkoc = Erkoc gold efective potential
!!   metal = tight-binding metal potentials (export to ble-files only)
!! MODIFIERs:
!!   polar = polar (linear or shell-core), exe=blendpol
!!   rep12 = add D/r^12 terms to exp-6/buck/busing potentials
!!   gui = GTK-based GUI (OLD VERSIONS ONLY)

!define x11
!include "../home.mmk"

!dir = gen sim blend

OREP = -DPOW=0
!if rep12
OREP = -DPOW=12
!endif
!if rep8
OREP = -DPOW=8
!endif

# preprocessor conditionals (use in statement O = -DXXX before !dir):
#   -DPOLAR : force field with polarizability
#   -DPOW={0,8,12} : [do not] add term D/r^POW to the exp6/Buckingham/Busing FF
#   -DSCR : extra code for scrolling screens of output (internal `less')
#   -DX11 : X windows graphics, needed, added later

# SS_TABLE is #defined in MACSIMUS/sim/FIELD/sitesite.h and refers to the name
# of the table containing atom potential parameters: RvdW, EvdW, parm[SS_PARMS]

!if lj
# see MACSIMUS/sim/lj
!dir = sim/lj
!if polar
O = -DPOLAR
!blendpol : blend blendpar blendmed blendedt blendgen blendimp blendmin \
  blenddep intrapot sitesite ground varfile xdraw
!else
!blend    : blend blendpar blendmed blendedt blendgen blendimp blendmin \
  blenddep intrapot sitesite ground varfile xdraw
!endif
!endif

!if buck
# see MACSIMUS/sim/buck
!dir = sim/buck
!if polar
O = -DPOLAR $(OREP)
!blendpolbuck   : blend blendpar blendmed blendedt blendgen blendimp blendmin \
  blenddep intrapot sitesite ground varfile xdraw
!else
O = $(OREP)
!blendbuck      : blend blendpar blendmed blendedt blendgen blendimp blendmin \
  blenddep intrapot sitesite ground varfile xdraw
!endif
!endif

!if exp6
# see MACSIMUS/sim/exp6
!dir = sim/exp6
!if polar
O = -DPOLAR $(OREP)
!blendpolexp6   : blend blendpar blendmed blendedt blendgen blendimp blendmin \
  blenddep intrapot sitesite ground varfile xdraw
!else
O = $(OREP)
!blendexp6      : blend blendpar blendmed blendedt blendgen blendimp blendmin \
  blenddep intrapot sitesite ground varfile xdraw
!endif
!endif

!if busing
!dir = sim/busing
!if polar
O = -DPOLAR $(OREP)
!blendpolbus : blend blendpar blendmed blendedt blendgen blendimp blendmin \
  blenddep intrapot sitesite ground varfile xdraw
!else
O = $(OREP)
!blendbus : blend blendpar blendmed blendedt blendgen blendimp blendmin \
  blenddep intrapot sitesite ground varfile xdraw
!endif
!endif

!if goldErkoc
# see MACSIMUS/sim/goldErkoc
!dir = gen sim sim/goldErkoc blend
!blendgoldErkoc : blend blendpar blendmed blendedt blendgen blendimp blendmin \
  blenddep intrapot sitesite ground varfile xdraw
!endif

!if metal
# see MACSIMUS/sim/metal
!dir = gen sim sim/metal blend
!blendmetal : blend blendpar blendmed blendedt blendgen blendimp blendmin \
  blenddep intrapot sitesite ground varfile xdraw
!endif

!! ------------------------------------------------------------------------

!pdb: LIBOPT=-lm: ground forscanf pdb pdbbasic pdbpdb pdbinfo pdbmol pdbconv pdbrep
!showpro: LIBOPT=-lm: showpro
!makepept: LIBOPT=: makepept
!blefilt: LIBOPT=: blefilt
!ramachan: LIBOPT=-lm: ramachan
!pdb2pdb: LIBOPT=-lm: pdb2pdb
!bonds: LIBOPT=-lm: bonds
!waterdep: LIBOPT=-lm: waterdep

!private "private.mmk"
