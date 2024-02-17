/* support for the Widom insertion particle method */

#ifdef POLAR
#  warning "WIDOM+POLAR not a good idea (unless uncharged molecule is inserted)"
#endif /*# POLAR */

#ifdef LINKCELL

#  undef DOM
#  define DOM(X) { \
  En.el=-En.el; LJqqm(X,ls,ssi,qi,&Ureal); En.el=-En.el; \
  ls->sp=sp; \
  LJqqm(X,ls,ssivirt,qivirt,&Uvirt); \
  ls->sp=mol->sp; }

void Eidchange(int sp,molecule_t *mol,vector *rpf,vector *rp) /*** Eidchange */
/*
  Energy difference caused by particle identity change from mol->sp to sp 
  Molecules mol->sp and sp must have the same structure - they should
    differ in charges and LJ terms only
  See also Eidchange in sim/forces.c
*/
{
  linklist_t *ls,*ls0,*l,*last;
  real DX,dy,DY,dyq,dz;
  double Ureal=0,Uvirt=0,qi,qivirt;
  vector cellr;
  vector r0;
  int ix,ixm,ix0,ix1;
  int iy,iym,iy0,iy1;
  int iz,izm,iz0,iz1;
  vector cell,rL; /* cell size and rL=1/cell */
  int n,i;
  sitesite_t *ssi,*ssivirt;
  siteinfo_t *si=spec[mol->sp]->si,*sivirt=spec[sp]->si;
  int (*is)[3];

  if (mol->ns != spec[sp]->ns)
    ERROR(("Eidchange: nsvirt=%d != nsreal=%d",spec[sp]->ns,mol->ns))

  VVV(cell,=box.L,/No.cell)
  VVV(rL,=No.cell,/box.L)
  En.el=En.pot=0;

  n=mol-cfg; /* of the molecule grown */
  allocarray(is,mol->ns);

  /* TRICK:
    cells are in the same order in array site_storage as the sites in rp 
    then, ls=the cell of the 1st site in mol
  */
  ls0 = site_storage + mol->ir/sizeof(vector);

  /* REMOVE temporarily the molecule from the list and prepare is[] */
  ls=ls0;
  loop (i,0,mol->ns) {
    if (ls->n!=n) ERROR(("Eidchange: ls->n=%d != n=%d",ls->n,n))

    VVV(cellr,=ls->r,*rL) /* in cell units */
    VV(is[i],=cellr)

    last=NULL;
    looplist (l,lastlist[is[i][0]][is[i][1]][is[i][2]]) {
    if (l==ls) {
      if (last) last->next=l->next;
      else lastlist[is[i][0]][is[i][1]][is[i][2]]=lastlist[is[i][0]][is[i][1]][is[i][2]]->next; 
      goto lsremoved; }
    last=l; }
  ERROR(("mol.site=%d.%d not in list[%d][%d][%d]\n\
*** l->r=%g %g %g",n,i,is[i][0],is[i][1],is[i][2],VARG(l->r)))
   lsremoved:; 
    ls++; }

  ls=ls0;
  loop (i,0,mol->ns) {
    VVV(cellr,=ls->r,*rL) /* in cell units */

    VV(r0,=ls->r)

    ssi=sstab[si[i].st];
    ssivirt=sstab[sivirt[i].st];
    qi=si[i].charge;
    qivirt=sivirt[i].charge;

    ixm = ix0 = (int)(cellr[0]+No.cell[0]-box.cutoff*rL[0])-No.cell[0];
    if (ixm<0) { ls->r[0]+=box.L[0]; ixm+=No.cell[0]; }

    if (ixm<0 || ixm>No.cell[0])
      ERROR(("Einsert: cutoff too large or internal error"))
    ix1=cellr[0]+box.cutoff*rL[0];

    loopto (ix,ix0,ix1) {
      if (ix>is[i][0]) DX=ix-cellr[0];
      else if (ix<is[i][0]) DX=ix+1-cellr[0];
      else DX=0;

      if (ixm==No.cell[0]) { ls->r[0]-=box.L[0]; ixm-=No.cell[0]; }
      dyq=box.cq-Sqr(DX*cell[0]);

      if (dyq<0) {
        if (dyq/box.cq<-1e-13) ERROR(("dyq=%g<0",dyq))
        dyq=0; }

      dy=sqrt(dyq)*rL[1];

      iym = iy0 = (int)(cellr[1]+No.cell[1]-dy)-No.cell[1];
      if (iym<0) { ls->r[1]+=box.L[1]; iym+=No.cell[1]; }
      if (iym<0 || iym>No.cell[1])
        ERROR(("lcforce: cutoff too large or internal error"))
      iy1 = (int)(cellr[1]+dy);

      loopto (iy,iy0,iy1) {
        if (iym==No.cell[1]) {
          ls->r[1]-=box.L[1]; iym-=No.cell[1]; }
        if (iy>is[i][1]) DY=iy-cellr[1];
        else if (iy<is[i][1]) DY=iy+1-cellr[1];
        else DY=0;
        dz=dyq-Sqr(DY*cell[1]);

        if (dz<0) {
          if (dz/box.cq<-1e-13) ERROR(("dz=%g<0",dz))
          dz=0; }

        dz=sqrt(dz)*rL[2];

#  ifdef SLIT
        /* wall or slit pore: not periodic in z */
        iz0 = (int)(cellr[2]-dz);
        if (iz0<0) iz0=0;
        iz1 = (int)(cellr[2]+dz);
        if (iz1>=No.cell[2]) iz1=No.cell[2]-1;

        loopto (izm,iz0,iz1) DOM(lastlist[ixm][iym][izm])
#  else /*# SLIT */
        /* periodic b.c. in z */
        izm = iz0 = (int)(cellr[2]+No.cell[2]-dz)-No.cell[2];
	if (izm<0) { ls->r[2]+=box.L[2]; izm+=No.cell[2]; }
	if (izm<0 || izm>No.cell[2])
	  ERROR(("lcforce: cutoff too large or internal error"))
        iz1 = (int)(cellr[2]+dz);

        loopto (iz,iz0,iz1) {
          if (izm==No.cell[2]) {
            ls->r[2]-=box.L[2]; izm-=No.cell[2]; }
          DOM(lastlist[ixm][iym][izm]) 
	  izm++; }
#  endif /*#!SLIT */
        ls->r[2]=r0[2];
 
        iym++; } /*iy*/

      ls->r[1]=r0[1];

      ixm++; } /*ix*/
    ls->r[0]=r0[0];

    ls++; }

  ls=ls0;
  loop (i,0,mol->ns) {
    /* RETURN ls to the list (at head) */
    ls->next=lastlist[is[i][0]][is[i][1]][is[i][2]];
    lastlist[is[i][0]][is[i][1]][is[i][2]]=ls;
    ls++; }

  //  put3(Uvirt,Ureal,En.el)

  free(is);
  En.pot=Uvirt-Ureal+En.el;
} /* Eidchange */

#  undef DOM
#  define DOM(X) LJqqm(X,&lsloc,ssi,qi,&U);

void Einsert(molecule_t *mol,vector *rpf,vector *rp) /************** Einsert */
{
  linklist_t lsloc;
  real DX,dy,DY,dyq,dz;
  double U=0,qi;
  vector cellr;
  real *r0;
  int ix,ixm,ixs,ix0,ix1;
  int iy,iym,iys,iy0,iy1;
  int iz,izm,izs,iz0,iz1;
  vector cell,rL; /* cell size and rL=1/cell */
  int i;
  sitesite_t *ssi;
  siteinfo_t *si=spec[mol->sp]->si;

  VVV(cell,=box.L,/No.cell)
  VVV(rL,=No.cell,/box.L)

  En.el=En.pot=0;
  lsloc.sp=mol->sp;
  lsloc.next=NULL;
  lsloc.n=No.s; /* impossible value - never equal (don't use negative because of FROM) */

  loop (i,0,mol->ns) {
    r0=rof(mol,rp)[i];
    VV(lsloc.r,=r0)
    lsloc.f=rof(mol,rpf)[i];
    ssi=sstab[si[i].st];
    qi=si[i].charge;
    
    lsloc.si=&si[i];

    VVV(cellr,=lsloc.r,*rL) /* in cell units */
    ixs=cellr[0];
    iys=cellr[1];
    izs=cellr[2];

    ixm = ix0 = (int)(cellr[0]+No.cell[0]-box.cutoff*rL[0])-No.cell[0];
    if (ixm<0) { lsloc.r[0]+=box.L[0]; ixm+=No.cell[0]; }

    if (ixm<0 || ixm>No.cell[0])
      ERROR(("Einsert: cutoff too large or internal error"))
    ix1=cellr[0]+box.cutoff*rL[0];

    loopto (ix,ix0,ix1) {
      if (ix>ixs) DX=ix-cellr[0];
      else if (ix<ixs) DX=ix+1-cellr[0];
      else DX=0;

      if (ixm==No.cell[0]) { lsloc.r[0]-=box.L[0]; ixm-=No.cell[0]; }
      dyq=box.cq-Sqr(DX*cell[0]);

      if (dyq<0) {
        if (dyq/box.cq<-1e-13) ERROR(("dyq=%g<0",dyq))
        dyq=0; }

      dy=sqrt(dyq)*rL[1];

      iym = iy0 = (int)(cellr[1]+No.cell[1]-dy)-No.cell[1];
      if (iym<0) { lsloc.r[1]+=box.L[1]; iym+=No.cell[1]; }
      if (iym<0 || iym>No.cell[1])
        ERROR(("lcforce: cutoff too large or internal error"))
      iy1 = (int)(cellr[1]+dy);

      loopto (iy,iy0,iy1) {
        if (iym==No.cell[1]) {
          lsloc.r[1]-=box.L[1]; iym-=No.cell[1]; }
        if (iy>iys) DY=iy-cellr[1];
        else if (iy<iys) DY=iy+1-cellr[1];
        else DY=0;
        dz=dyq-Sqr(DY*cell[1]);

        if (dz<0) {
          if (dz/box.cq<-1e-13) ERROR(("dz=%g<0",dz))
          dz=0; }

        dz=sqrt(dz)*rL[2];

#  ifdef SLIT
        /* wall or slit pore: not periodic in z */
        iz0 = (int)(cellr[2]-dz);
        if (iz0<0) iz0=0;
        iz1 = (int)(cellr[2]+dz);
        if (iz1>=No.cell[2]) iz1=No.cell[2]-1;

        loopto (izm,iz0,iz1) DOM(lastlist[ixm][iym][izm])
#  else /*# SLIT */
        /* periodic b.c. in z */
        izm = iz0 = (int)(cellr[2]+No.cell[2]-dz)-No.cell[2];
	if (izm<0) { lsloc.r[2]+=box.L[2]; izm+=No.cell[2]; }
	if (izm<0 || izm>No.cell[2])
	  ERROR(("lcforce: cutoff too large or internal error"))
        iz1 = (int)(cellr[2]+dz);

        loopto (iz,iz0,iz1) {
          if (izm==No.cell[2]) {
            lsloc.r[2]-=box.L[2]; izm-=No.cell[2]; }
          DOM(lastlist[ixm][iym][izm]) 
	  izm++; }
#  endif /*#!SLIT */
        lsloc.r[2]=r0[2];
 
        iym++; } /*iy*/

      lsloc.r[1]=r0[1];

      ixm++; } /*ix*/
    //    lsloc.r[0]=r0[0]; /* restore shifted - last not needed */
  }

  En.pot=U+En.el;
} /* Einsert */

#else /*# LINKCELL */

void Einsert(molecule_t *mol,vector *rpf,vector *rp) /************** Einsert */
/* energy of a (virtual) particle inserted; Ewald not supported */
{
  int m;

  En.el=En.pot=0;
  loop (m,0,No.N)
    En.pot+=(*pot[mol->sp][cfg[m].sp]) (
	       rof(mol,rpf),rof(cfg+m,rpf),
	       mol, cfg+m, rp );
  En.pot+=En.el;
}


void Eidchange(int sp,molecule_t *mol,vector *rpf,vector *rp) /*** Eidchange */
/*
  energy difference caused by particle identity change from mol->sp to sp
*/
{
  int m;
  int oldsp=mol->sp;

  /* old energy */
  En.el=En.pot=0;
  loop (m,0,No.N) if (mol!=&cfg[m])
    En.pot-=(*pot[mol->sp][cfg[m].sp]) (
	       rof(mol,rpf),rof(cfg+m,rpf),
	       mol, cfg+m, rp );
  En.pot-=En.el;

  En.el=0;
  /* new energy */
  mol->sp=sp;
  loop (m,0,No.N) if (mol!=&cfg[m])
    En.pot+=(*pot[mol->sp][cfg[m].sp]) (
	       rof(mol,rpf),rof(cfg+m,rpf),
	       mol, cfg+m, rp );
  En.pot+=En.el;

  mol->sp=oldsp;
}
#endif /*#!LINKCELL */
