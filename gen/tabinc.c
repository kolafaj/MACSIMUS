/*
  05/2020 return value of to_colon() decreased by 1:
    now points to the position of ':'
    (which is replaced by 0, so strlen(to_colon()) returns always 0)
  04/2016 bug fixed: PI now (again) recognized
*/

int max_column=0; /* max ID found in :EXPR */

char *to_colon(char *a) /****************************************** to_colon */
/*
  translate #1->A, #0->n, x->A, etc., stop at `:'
  NEW: if over Z, then continues by c27, c28, ... (just #->c)
  returns the position of : or NULL
*/
{
  char *e;

  if (!a) return NULL;

  /* patch: PI problem */
  for (e=a; *e; e++) if (*e=='P' && e[1]=='I') *e='p',e[1]='i';

  for (e=a; ; a++, e++) {
    *e=*a;
    if (*a==0) return NULL;
    if (*a=='@') *e='n';
    if (strchr("xyz",*a) && (!isalnum(a[1]))) *e=*a-'x'+'A'; /* x->A ... */
    if (isupper(*e)) Max(max_column,*e-'@')
    if (strchr("#c",*a) && isdigit(a[1])) {
      int n=atoi(a+1);

      Max(max_column,n);
      if (n<27) {
        a++;
        while (isdigit(a[1])) a++; /* shift */
        if (n) *e='@'+n; else *e='n'; }
      else if (n>N)
        ERROR(("column #%d out of range (max. %d)",n,N))
      else
        *e='c'; }
    if (*a==':') {
      *a=*e=0;
      return a; }
  }
}

void columnidlist(void) /************************************** columnidlist */
/*
  making a list of id's for Calc(); n,A, etc. are first (=fastest)
*/
{
  int i;
  struct _Idlist_s *id;

  for (i=N-1; i>=0; i--) {
    alloczero(id,sizeof(struct _Idlist_s)+7); /* bit wasteful - but padded anyway*/
    id->next=_Id.head;
    _Id.head=id;
    id->val=0;
    id->id[0]='@'+i;
    if (i>26) sprintf(id->id,"c%d",i);
    else id->id[0]='@'+i; /* A..Z; n (i=0) see below: */ }

  /* column 0=line count, denoted by n */
  count=&id->val;
  id->id[0]='n';
}
