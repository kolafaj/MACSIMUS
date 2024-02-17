
/*** bonds ***/
if (option('b')) for (b=b0; b; b=b->next) {
  i=b->indx[0]; j=b->indx[1];
#ifdef OMITKEPT
  if (!(keepmask & site[i].keep & site[j].keep))
#endif
    Ubond += bondpot(r[i],r[j], f[i],f[j], &b->parm); }

if (option('i'))
/*** angles - check pyramidal error of tetrahedral angles ***/
for (a=a0; a; a=a->next) {
  i=a->indx[0]; j=a->indx[1]; k=a->indx[2];
#ifdef OMITKEPT
  if (!(keepmask & site[i].keep & site[j].keep & site[k].keep))
#endif
    {
      double U1=anglepot(r[i],r[j],r[k], f[i],f[j],f[k], &a->parm);
      Uangle += U1; 
      if (U1>Uanglemax) { 
        Uanglemax=U1; ianglemax=j; } } }
else     
/*** angles - do not check pyramidal error ***/
for (a=a0; a; a=a->next) {
  i=a->indx[0]; j=a->indx[1]; k=a->indx[2];
#ifdef OMITKEPT
  if (!(keepmask & site[i].keep & site[j].keep & site[k].keep))
#endif
    Uangle += anglepot(r[i],r[j],r[k], f[i],f[j],f[k], &a->parm); }

/*** dihedrals ***/
for (t=d0; t; t=t->next) {
  i=t->indx[0]; j=t->indx[1]; k=t->indx[2]; l=t->indx[3];
#ifdef OMITKEPT
  if (!(keepmask & site[i].keep & site[j].keep & site[k].keep & site[l].keep))
#endif
    Udih += dihedralpot(r[i],r[j],r[k],r[l], f[i],f[j],f[k],f[l], &t->parm); }

/*** cisdihedrals ***/
for (t=cis0; t; t=t->next) {
  i=t->indx[0]; j=t->indx[1]; k=t->indx[2]; l=t->indx[3];
#ifdef OMITKEPT
  if (!(keepmask & site[i].keep & site[j].keep & site[k].keep & site[l].keep))
#endif
    Udih += dihedralpot(r[i],r[j],r[k],r[l], f[i],f[j],f[k],f[l], &t->parm); }

/*** impropers ***/
for (t=i0; t; t=t->next) {
  i=t->indx[0]; j=t->indx[1]; k=t->indx[2]; l=t->indx[3];
#ifdef OMITKEPT
  if (!(keepmask & site[i].keep & site[j].keep & site[k].keep & site[l].keep))
#endif
    Uimp += improperpot(r[i],r[j],r[k],r[l], f[i],f[j],f[k],f[l], &t->parm); }

/*** aromatic dihedrals ***/
for (t=ar0; t; t=t->next) {
  i=t->indx[0]; j=t->indx[1]; k=t->indx[2]; l=t->indx[3];
#ifdef OMITKEPT
  if (!(keepmask & site[i].keep & site[j].keep & site[k].keep & site[l].keep))
#endif
    Uaro += dihedralpot(r[i],r[j],r[k],r[l], f[i],f[j],f[k],f[l], &t->parm); }
