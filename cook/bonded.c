sp=spec[m1->sp];

if (sp->intra) {

  /*** bonds ***/
  for (b=sp->bond; b; b=b->next) {
    i=b->indx[0]; j=b->indx[1];
    U+=bondpot(r1[i],r1[j], f1[i],f1[j], &b->parm); }
  DEBUGU(bonds)

  /*** angles ***/
  for (an=sp->angle; an; an=an->next) {
    i=an->indx[0]; j=an->indx[1]; k=an->indx[2];
    U+=anglepot(r1[i],r1[j],r1[k], f1[i],f1[j],f1[k], &an->parm); }
  DEBUGU(angles)

  /*** dihedrals ***/
  for (t=sp->dihedral; t; t=t->next) {
    i=t->indx[0]; j=t->indx[1]; k=t->indx[2]; l=t->indx[3];
    U+=dihedralpot(r1[i],r1[j],r1[k],r1[l], f1[i],f1[j],f1[k],f1[l], &t->parm); }
  DEBUGU(dihedrals)

  /*** impropers ***/
  for (t=sp->improper; t; t=t->next) {
    i=t->indx[0]; j=t->indx[1]; k=t->indx[2]; l=t->indx[3];
    U+=improperpot(r1[i],r1[j],r1[k],r1[l], f1[i],f1[j],f1[k],f1[l], &t->parm); }
  DEBUGU(impropers)

  /*** aromatic dihedrals ***/
  for (t=sp->aromatics; t; t=t->next) {
    i=t->indx[0]; j=t->indx[1]; k=t->indx[2]; l=t->indx[3];
    U+=dihedralpot(r1[i],r1[j],r1[k],r1[l], f1[i],f1[j],f1[k],f1[l], &t->parm); }
  DEBUGU(aromatics)
  }
