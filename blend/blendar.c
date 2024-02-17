/* msd + angle list */

/******************************* independent angles for entropy estimation */

static int minsite(torsion_t *t)
{
  int *i=t->indx;

  if (i[0]<i[1] && i[0]<i[2] && i[0]<i[3]) return i[0];
  if (i[1]<i[2] && i[1]<i[3]) return i[1];
  if (i[2]<i[3]) return i[2];

  return i[3];
}

static unsigned incl;
static FILE *info;
static int list[10];

static void makecycle(int c,torsion_t *t,int a,int b)
{
  int i;

  loop (i,0,4) list[i]=t->indx[i];
  list[4]=a; list[5]=b;
  loop (i,c,10) list[i]=list[i-c];
  loop (i,0,c) site[list[i]].nest++;
}

static void removedihedral(int i)
{
  torsion_t *t;

  for (t=d0; t; t=t->next) if (t->parm.K[0])
    if ((t->indx[0]==list[i] && t->indx[1]==list[i+1] && t->indx[2]==list[i+2] && t->indx[3]==list[i+3])
     || (t->indx[3]==list[i] && t->indx[2]==list[i+1] && t->indx[1]==list[i+2] && t->indx[0]==list[i+3])) {
      t->parm.K[0]=0;
      fprintf(info,"D %d %d %d %d  %s %s %s %s\n",
              t->indx[0],t->indx[1],t->indx[2],t->indx[3],
              site[t->indx[0]].id,
               site[t->indx[1]].id,
              site[t->indx[2]].id,
              site[t->indx[3]].id);
      incl++; }
}

static void removeangle(int i)
{
  angle_t *a;

  for (a=a0; a; a=a->next) if (a->parm.K) if (a->indx[1]==list[i+1])
    if ((a->indx[0]==list[i] && a->indx[2]==list[i+2])
        || (a->indx[2]==list[i] && a->indx[0]==list[i+2])) {
       a->parm.K=a->parm.K2=0;
      fprintf(info,"A %d %d %d  %s %s %s\n",
              a->indx[0],a->indx[1],a->indx[2],
              site[a->indx[0]].id,site[a->indx[1]].id,site[a->indx[2]].id);
      incl++; }
}

static void removeimproper(int i)
/* here i is site, not index in list[] */
{
  torsion_t *t;

  for (t=i0; t; t=t->next) if (t->indx[0]==i && t->parm.K[0]) {
    t->parm.K[0]=0;
    fprintf(info,"I %d %d %d %d  %s %s %s %s\n",
            t->indx[0],t->indx[1],t->indx[2],t->indx[3],
            site[t->indx[0]].id,
             site[t->indx[1]].id,
            site[t->indx[2]].id,
            site[t->indx[3]].id);
    incl++; }
}

void anglemsd(species_t *spec) /********************************** anglemsd */
{
  int ns=spec->ns,size=ns*sizeof(vector),i,j,k,l;
  unsigned it,cyc5=0,cyc6=0,sw;
  vector *r,*r0;
  FILE *msd=NULL;
  angle_t *a;
  torsion_t *t,*tt;
  FILE *ang;
  static vector dummy;

  site=spec->site;

  if (spec->Xopt.A) {
    strcpy(spec->ext,".ang");
    info=fopen(spec->fn,"wt");
    fprintf(info,"*** list of removed Angles Dihedrals Impropers ***\n");
    incl=0;

    /* concatenating with aromatics */
    for (t=d0; t; t=t->next) if (t->next==NULL) { t->next=ar0; break; }

    fprintf(info,"removing angles around 3- and 4-bonded atoms\n");
    /* removing redundant angles for 3 and 4 bonded atoms */
    for (a=a0; a; a=a->next) {
      i=a->indx[0], j=a->indx[1], k=a->indx[2];
      if (site[j].nnbr>=3) {
        if ((site[j].nbr[0]==i && site[j].nbr[1]==k)
         || (site[j].nbr[0]==k && site[j].nbr[1]==i)) {
          a->parm.K=a->parm.K2=0;
          fprintf(info,"A %d %d %d  %s %s %s\n",
                  i,j,k,
                  site[i].id,site[j].id,site[k].id);
          incl++; }
      } }

    /* removing in 6- and 5-cycles */
    for (t=d0; t; t=t->next) {
      for (tt=t->next; tt; tt=tt->next) {
        if ((i=t->indx[0]==tt->indx[0] && t->indx[3]==tt->indx[3])
         ||   (t->indx[0]==tt->indx[3] && t->indx[3]==tt->indx[0])) {
          int it=minsite(t),itt=minsite(tt);
          Min(it,itt)
          if (t->indx[0]==it || t->indx[3]==it) {
            fprintf(info,"removing in 6-cycle [%d]\n",cyc6++);
            if (i) makecycle(6,t,tt->indx[2],tt->indx[1]);
            else makecycle(6,t,tt->indx[1],tt->indx[2]);
            i=incl;
            removedihedral(1);
            removedihedral(2);
            removedihedral(4);
            removedihedral(5);
            removeangle(2);
            removeangle(5);
            loop (j,0,6) {
              if (incl==i+6) break;
              removedihedral(j);
              if (incl==i+6) break;
              removeangle(j); }
            if (incl!=i+6) fprintf(info,"COULD NOT REMOVE %d CONSTRAINTS\n",i+6-incl);
/*.....         loop (i,0,6) removeimproper(i);*/
            }
          } }
      for (a=a0; a; a=a->next) {
        if ((t->indx[0]==a->indx[0] && t->indx[3]==a->indx[2])
         || (t->indx[0]==a->indx[2] && t->indx[3]==a->indx[0])) {
          int it=minsite(t);

          Min(it,a->indx[1])
          if (a->indx[1]==it) {
            fprintf(info,"removing in 5-cycle [%d]\n",cyc5++);
            makecycle(5,t,a->indx[1],0);
            i=incl;
            removedihedral(1);
            removedihedral(2);
            removedihedral(4);
            removedihedral(0);
            removeangle(2);
            removeangle(0);
            loop (j,0,5) {
              if (incl==i+6) break;
              removedihedral(j);
              if (incl==i+6) break;
              removeangle(j); }
            if (incl!=i+6) fprintf(info,"COULD NOT REMOVE %d CONSTRAINTS\n",i+6-incl);
/*.....         loop (i,0,5) removeimproper(i);*/
            }
          } } }

    fprintf(info,"impropers in merged rings:\n");
    loop (i,0,ns) if (site[i].nest>1) {
      loop (k,0,site[i].nnbr) {
        l=site[i].nbr[k];
        if (site[l].nest>1) {
          j=incl;
          removeimproper(i);
          if (j==incl) removeimproper(l);
          if (j!=incl) site[i].nest--,site[l].nest--; } } }

    fprintf(info,"\n%d 5-cycles and %d 6-cycles detected\n",cyc5,cyc6);
    fprintf(info,"%d constraints removed\n\n",incl);
    fprintf(info,"*** list of included Angles Dihedrals Impropers ***\n");
    incl=0; }

  ralloc(r0,size);
  ralloc(r,size);
  measure=2;

  if (spec->Xopt.D) {
    strcpy(spec->ext,".msd");
    msd=fopen(spec->fn,"wt"); }

  spec->opt_w=0; /* unify! */

  for (it=0;;it++) {
    if (spec->frame>spec->Xopt.toframe) break;
    if (abs(option('r'))!=4) if (it>=1) break;
    if (!read3D(spec,-1)) break;
    spec->frame+=spec->Xopt.byframe;
    loop (i,0,ns) VV(r[i],=site[i].r)
    if (it) {
      vector dr;
      double summsd=0,summ=0;

      if (msd) {
        loop (i,0,ns) {
          double m=atom[site[i].type].mass;

          VVV(dr,=r[i],-r0[i])
          summsd += SQR(dr)*m;
          summ += m; }
        fprintf(msd,"%10.4f\n",summsd/summ); } }
    else
      copy(r0,r,size);

    if (spec->Xopt.A) {
      sprintf(spec->ext,".%d",it+1);
      ang=fopen(spec->fn,"wt");

      for (a=a0; a; a=a->next) if (a->parm.K) {
        i=a->indx[0]; j=a->indx[1]; k=a->indx[2];
        anglepot(r[i],r[j],r[k], dummy,dummy,dummy, &a->parm);
        fprintf(ang,"%12.8f\n",phi*(180/PI));
        if (info) {
          fprintf(info,"a %d %d %d  %s %s %s\n",
                  i,j,k,
                  site[i].id,site[j].id,site[k].id);
          incl++; } }

      /*** dihedrals and impropers ***/
      loop (sw,0,2)
        for (t=sw?i0:d0; t; t=t->next) if (t->parm.K[0]) {
          i=t->indx[0]; j=t->indx[1]; k=t->indx[2]; l=t->indx[3];
          (sw?improperpot:dihedralpot)(r[i],r[j],r[k],r[l],
                                       dummy,dummy,dummy,dummy, &t->parm);
          fprintf(ang,"%12.8f\n",phi*(180/PI));
          if (info) {
            fprintf(info,"%c %d %d %d %d  %s %s %s %s\n","di"[sw],
                    i,j,k,l,
                    site[i].id,site[j].id,site[k].id,site[l].id);
            incl++; } }

      fclose(ang);
      ang=NULL;
      if (info) {
        fprintf(info,"... total %d included, theory=%d\n\n",
                incl,2*ns-cyc5-cyc6-5*spec->nclust);
        fclose(info); } }

     info=NULL; }

  if (msd) fclose(msd);

  release(r0);
}
