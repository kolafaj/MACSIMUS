/* make show
 */
int NS; /* actual size based (selection/ full), plb only */
int movesel,render=0;
struct animated_s {
  int gif;
  char *options;
} animated={0,CONVERT};

void preparedump(int sel) /************************************ preparedump */
{
  if (dumpmode==NONE) ERROR(("internal"))
  if (!dumpname) {
    if (strchr(plbname,'%'))
/*.....      ERROR(("output not supported for a series of plb files"))*/
      prt("! WARNING: dump from a series of plb files may cause problems");
    alloc(dumpname,strlen(plbname)+21); }
  if (!percentfn[0]) strcpy(percentfn,plbname);
  if (option('_'))
    sprintf(dumpname,"%s.%d%c.%s",percentfn,pos,'a'+dumpindex++,dumpext[dumpmode]);
  else
    sprintf(dumpname,"%s.%04d.%s",percentfn,dumpindex++,dumpext[dumpmode]);
  dumpprepared=1;

  NS=ns;
  if (sel) {
    /* applies to plb + selection only */
    int i;

    NS=0;
    loop (i,0,ns) if (site[i].mark&4) NS++;
    if (!NS) {
      prt("WARNING: empty selection, using whole configuration");
      NS=ns; }
    prt("%d sites to be dumped",NS); }
}

void opendump(void) /********************************************** opendump */
{
  if (dumpprepared) {
    dump=fopen(dumpname,dumpmode>=NFF ? "wt": "wb");
    if (!dump) {
      ERROR(("cannot write to %s",dumpname))
      return; }

    switch (dumpmode) {
      case BONDEPS:
      case EPS: {
        time_t t;
        time(&t);

        fprintf(dump,"\
%%!PS-Adobe-3.0 EPSF-3.0\n\
%%%%BoundingBox: 0 0 %d %d\n\
%%%%Creator: MACSIMUS/show\n\
%%%%Title: %s\n\
%%%%CreationDate: %s\
%%%%Pages: 1\n\
%%%%LanguageLevel: 1\n\
0 setgray 1 setlinejoin 1 setlinewidth\n\
/DS %d string def\n",maxxn,maxyn,molname,ctime(&t),maxxn*3); }
        break;

      case PLB:
        hdr[0]=NS;
        fwrite(hdr,sizeof(float),2,dump);
        if (hdr[1]==-3) fwrite(L,sizeof(float),3,dump);
        break;

      case ATM:
        fprintf(dump,"%d\n",NS);
        break;

      case POV:
      case NFF: {
        fvector l[2],up={0,1,0};
        int i,ibg;
        vector bg={0,0,0}; // default if !ballbackground

        if (ballbackground)

# if 0
        if (xwindowhints.background && xwindowhints.background[0]) {
          sscanf(xwindowhints.background+1,"%x",&ibg);
          bg[0]=(ibg&0xff0000)/(double)0xff0000;
          bg[1]=(ibg&0xff00)/(double)0xff00;
          bg[2]=(ibg&0xff)/(double)0xff; }
# else
        if (ballbackground) {
          sscanf(ballbackground,"%x",&ibg);
          bg[0]=(ibg&0xff0000)/(double)0xff0000;
          bg[1]=(ibg&0xff00)/(double)0xff00;
          bg[2]=(ibg&0xff)/(double)0xff; }
# endif

        if (option('w')) bg[0]=bg[1]=bg[2]=1;

        NFFfloor=3e33;
        switch (option('t')) {
          default: if (dumpstat==ONE) break;
          case 0: NFFfloor=-NFFfloor;
          case 1: break; }

        if (dumpstat!=CONT) {
          /* two lights */
          loop (i,0,2) {
            l[i][0]=eye[0]+eye[2]*(-0.15+0.45*i);
            l[i][1]=eye[1]+eye[2];
            l[i][2]=eye[2]*1.3; }

          /* zview=2*yrg gives angle=28 */
          if (dumpstat!=CONT) {
            if (dumpmode==NFF) fprintf(dump,"\
# %s, frame %d dumped by show, %f %f\n\
v\n\
from %f %f %f\n\
at %f %f 0\n\
up %f %f %f\n\
angle %f\n\
hither 1\n\
resolution %d %d\n\
l %f %f %f\n\
l %f %f %f\n\
b %f %f %f\n\
",
                                       molname,pos,
                                       Scale[0],Scale[1],
                                       VEC(eye),
                                       eye[0],eye[1],
                                       VEC(up),
                                       (360/PI)*atan(maxxn/Scale[0]/eye[2]/2),
                                       maxxn,maxyn,
                                       VEC(l[0]),
                                       VEC(l[1]),
                                       VEC(bg));
            else fprintf(dump,"\
//%s, frame %d dumped by show, %f %f\n\
#include \"show.inc\"\n\
camera {location <%f,%f,%f> look_at <%f,%f,0> angle %f}\n\
light_source {<%f,%f,%f> color rgb 0.3}\n\
light_source {<%f,%f,%f> color rgb .9}\n\
light_source {<%f,%f,%f> color rgb .9}\n\
background {color rgb <%f,%f,%f>}\n\
",
                         molname,pos,
                         Scale[0],Scale[1],
                         iVEC(eye),
                         eye[0],eye[1],
                         /*.....iVEC(up), must be perpendicular!!! */
                         (360/PI)*atan(maxxn/Scale[0]/eye[2]/2),
                         iVEC(eye) /* =camera (attempt to improve poor ambient light) */,
                         iVEC(l[0]),
                         iVEC(l[1]),
                         VEC(bg)); }
        } } /* NFF */
    default:; }
    if (dumpstat==START) dumpstat=CONT;
    dumpprepared=0;

  if (NS!=ns) {
    char *dot=NULL,*c;
    int i,j,*ren;
    FILE *mol;

    for (c=dumpname; *c; c++) if (*c=='.') dot=c;
    if (!dot) ERROR(("no . in %s",dumpname))

    strcpy(dot,".mol");
    mol=fopen(dumpname,"wt");
    fprintf(mol,"selection\n\
parameter_set = %s\n\
number_of_atoms = %d\n\
atoms\n\
! i  atom-id  a-type charge chir nbonds bound_atoms\n",
            parameter_set,NS);

    allocarrayzero(ren,ns);
    j=0;
    loop (i,0,ns) if (site[i].mark&4) ren[i]=j++;
    if (j!=NS) ERROR(("number of marked atoms unexpectedly changed"))

    j=0;
    loop (i,0,ns) if (site[i].mark&4) {
      int k,h,nnbr=0;

      loop (k,0,nc) loop (h,0,2) if (bond[k][h]==i)
        if (site[bond[k][1-h]].mark&4) nnbr++;

      // BUG: chirality not copied, 0 used
      fprintf(mol,"%4d %s %s %9.6f 0 %d",j,
              site[i].id,site[i].type,site[i].charge,nnbr);
      loop (k,0,nc) loop (h,0,2) if (bond[k][h]==i)
        if (site[bond[k][1-h]].mark&4) fprintf(mol," %d",ren[bond[k][1-h]]);

      fprintf(mol,"\n");
      j++; }

    free(ren);

    /* table of exported sites, may be useful */
    fprintf(mol,"\n! original numbering of sites:\n");

    loop (i,0,ns) if (site[i].mark&4)
      fprintf(mol,"! %d SITE\n",i);

    fclose(mol);

    strcpy(dot,".gol");
    mol=fopen(dumpname,"wt");
    fprintf(mol,"!sphere#1\n! selection\n%d\n",NS);
    loop (i,0,ns) if (site[i].mark&4)
      fprintf(mol,"%s %.2f\n",colortab[site[i].c/COLLEN].col,site[i].r*0.01);
    fclose(mol); }
  }
}

void closedump(int closeseries) /********************************* closedump */
{
  if (dump) {
    if (dumpmode>=EPS) fputs("\nshowpage\n%%EOF\n",dump);
    if (dumpmode==POV) if (NFFfloor>-2.999e33) fprintf(dump,"plane {y,%f texture{FLOOR}}\n",NFFfloor);
    if (dumpmode==NFF) {
      int i;
      fvector p[4];

      loop (i,0,4) {
        p[i][0]=20*eye[2]*(1-2*(i>1));
        p[i][1]=NFFfloor;
        p[i][2]=p[i][0]*(1-2*(i&1)); }

      if (NFFfloor>-1e9) fprintf(dump,"f 0.8 0.8 0.6 0.9 0.1 1 0 0\n\
p 4\n\
%f %f %f\n\
%f %f %f\n\
%f %f %f\n\
%f %f %f\n",
VEC(p[0]),
VEC(p[1]),
VEC(p[2]),
VEC(p[3]));
    }
    fclose(dump);
    dump=NULL;
    fprintf(stderr,"%s dumped%s\n",dumpname,dumpstat==CONT?" (merged)":"");

    /* original code:
    fprintf(stderr,"%s dumped%s\n",dumpname,
            closeseries && dumpstat==SERIES?", series closed":"");
    */
    if (closeseries && dumpstat==SERIES) WARNING(("series closed in unexpected way"))

    if (dumpstat!=SERIES && render) {
      if (system(string("start %s &",dumpname)))
        fprintf(stderr,"start %s failed\n",dumpname); }

    if (!closeseries && dumpstat==SERIES) preparedump(ns!=NS);
    else dumpmode=NONE; }
  else {
    fprintf(stderr,"(merged) series closed\n");
    if (closeseries && animated.gif) {
      if (dumpstat==SERIES && dumpmode==PIC) {
        char *s=string("convert %s %s.????.%s %s.gif",
                       animated.options,
                       percentfn,dumpext[dumpmode],
                       percentfn);
        int retcode;

        fprintf(stderr,"Making animated gif by command:\n%s\nplease wait...",s);
        if (option('_')) ERROR(("option -_ not allowed for making animated GIFS"))
        retcode=system(string("convert %s %s.????.%s %s.gif",
                           animated.options,
                           percentfn,dumpext[dumpmode],percentfn));
        fprintf(stderr,"done retcode=%d\n",retcode);
      }
    }

    dumpmode=NONE;
    dumpstat=ONE;
    dumpindex--;
  }
}

void calcNFFfloor(double r) /********************************** calcNFFfloor */
{
  if (NFFfloor>r) NFFfloor=r;
}

void dumpNFFcol(int c) /***************************************** dumpNFFcol */
{
  int i;
  unsigned char *rgb=colortab[c/=COLLEN].rgb;

  fprintf(dump,"f");
  loop (i,0,3) fprintf(dump," %.3f",NFFq.opaque*rgb[i]/255.);

  fprintf(dump," %s\n",colortab[c].NFFcol);
}

void dumponeplb(int ifrot) /************************************* dumponeplb */
{
  int locns=0;

  if (dump) {
    ERROR(("dump opened - closing"))
    closedump(1); }
  dumpmode=PLB;
  preparedump(movesel);

  opendump();
  if (ifrot) {
    int i;
    fvector cr,R;

    loop (i,0,ns) if (ns==NS || site[i].mark&4) {
      locns++;

      VVV(cr,=cfg[i],-center)

      R[0]=rot[0][0]*cr[0]+rot[0][1]*cr[1]+rot[0][2]*cr[2]+eye[0];
      R[1]=rot[1][0]*cr[0]+rot[1][1]*cr[1]+rot[1][2]*cr[2]+eye[1];
      R[2]=rot[2][0]*cr[0]+rot[2][1]*cr[1]+rot[2][2]*cr[2];

      fwrite(R,sizeof(fvector),1,dump); }
  }
  else {
    int i;

    loop (i,0,ns) if (ns==NS || site[i].mark&4) {
      locns++;
      fwrite(cfg+i,sizeof(fvector),1,dump); }
  }
  closedump(0);

  if (locns!=NS) WARNING(("selection changed during write"))
}
