/* cc -O2 -o exp2pow exp2pow.c
 */
#include "../gen/include.h"

#define xDEBUG

struct items_s {
  struct items_s *next,*prev;
  char *name;
  enum type_t { UNKNOWN, NUM, ID, OP, LPAR, RPAR } type;
  int prec; /* OP precedence: -1: lower than ^, 0=^, +1: higher */
  int nest; /* +1 inside parentheses */
  int mark;
} *head,*last;

int prec;
char *expop="^";

char *additem(char *name,char *end,enum type_t type)
{
  struct items_s *new;

  allocone(new);
  new->next=NULL;
  if (last) {
    last->next=new;
    new->prev=last;
    last=new; }
  else {
    head=last=new;
    head->prev=NULL; }
  new->prec=prec;
  new->type=type;
  alloczero(new->name,end-name+1);
  memcpy(new->name,name,end-name);

  return end;
}

char *isid(char *a)
/* returns a pointer after identifier end */
{
  if (isalpha(*a) || *a=='_') {
    while (isalnum(*a) ||  *a=='_') a++;
    return a; }

  return NULL;
}

char *isnum(char *a)
/* returns a pointer after number end */
/* hexadecimal numbers not accepted */
{
  if (isdigit(*a) || (*a=='.' && isdigit(a[1]))) {
    while (isdigit(*a) || *a=='.') a++;
    if (*a && strchr("lLuU",*a)) a++; /* trailing integer modifiers */
    else if (*a && strchr("eE",*a)) {
      /* exp format */
      a++;
      if (*a && strchr("+-",*a)) a++;
      while (isdigit(*a)) a++; }
    return a; }
  return NULL;
}

int isint(struct items_s *item)
{
  if (item->type!=NUM) return 0;
  if (strchr(item->name,'e') || strchr(item->name,'E') || strchr(item->name,'.') ) return 0;

  return 1;
}

char *isop(char *a)
/* returns a pointer after operator end */
{
  char *ops[]={ /* 1st char denotes the precedence wrt ^ */
    " ^",
    "+++", "+--", /* should not happen */
    "-+=","-+=","-*=","-/=","-%=","->=","-<=","-!=","-==",
    "->>","-<<", /* bitwise | & ^ not considered */
    "-+", "--", "-*", "-/", "-%", /* binary */
    "+.",
    "-,",
    NULL };
  char **o;

  for (o=ops; *o; o++)
    if (!memcmp(*o+1,a,strlen(*o+1))) {
      switch ((*o)[0]) {
        case '+': prec= 1; break;
        case '-': prec=-1; break;
        default: prec=0; }
      return a+strlen(*o+1); }      

  return NULL;
}      

void line2items(char *line)
/* Splits a line into items. Unary +/- are detached from numbers */
{
  char *c=line,*e;
  
  head=last=NULL;

  for (;;) {

    if (*c==0) return;

    while (*c<=' ') {
      if (*c==0) return;
      c++; }

    prec=0;

    if (e=isid(c))            c=additem(c,e,ID);
    else if (e=isop(c))       c=additem(c,e,OP);
    else if (e=isnum(c))      c=additem(c,e,NUM);
    else if (strchr("([",*c)) c=additem(c,c+1,LPAR);
    else if (strchr(")]",*c)) c=additem(c,c+1,RPAR);
    else                      c=additem(c,c+1,UNKNOWN); }
}

struct items_s *addafter(struct items_s *l,char *name,int type)
{
  struct items_s *new;

  if (!l) Error("internal: add after empty");

  allocone(new);
  new->prev=l;
  new->next=l->next;
  if (l->next) l->next->prev=new;
  l->next=new;
  new->name=strdup(name);
  new->type=type;

  return new;
}

struct items_s *addbefore(struct items_s *l,char *name,int type)
{
  struct items_s *new;

  if (!l) Error("internal: add before empty");

  allocone(new);
  new->next=l;
  new->prev=l->prev;
  if (l->prev) l->prev->next=new;
  l->prev=new;
  new->name=strdup(name);
  new->type=type;

  if (l==head) head=new;

  return new;
}

struct items_s *delete(struct items_s *l)
{
  struct items_s *n=l->next;

  if (n) n->prev=l->prev;
  if (l->prev) l->prev->next=n; else head=n;

  return n;
}

void log2log (struct items_s *head)
{
  struct items_s *l,*pow,*beg,*end;

  looplist (l,head) {
    if (!strcmp(l->name,"log")) l->name="log10";
    if (!strcmp(l->name,"ln")) l->name="log";
  }
}


void exp2pow(char *expr)
{
  struct items_s *l,*pow,*beg,*end;

  head=last=NULL;
  line2items(expr);
  log2log(head);

  for (;;) {
    int nest=0;

    looplist (l,head) {
      l->mark=0;
      l->nest=nest;
      if (l->type==LPAR) { nest++; l->nest=nest; }
      if (l->type==RPAR) nest--; }

    looplist (pow,head) if (pow->type==OP && pow->prec==0) break;
    if (!pow) break;

    if (strcmp(pow->name,expop))
      Error("internal: name does not match pow operator marked");

    /* AFTER ^ */
    end=pow->next;
    /* one unary + allowed => removed because C does not allow it */
    if (!strcmp("+",end->name)) end=delete(end);
    /* unary - */
    while (end && !strcmp("-",end->name)) {
      end->mark++; end=end->next; }
    /* unary * */
    while (end && !strcmp("*",end->name)) {
      end->mark++; end=end->next; }

    for (;end;end=end->next) {
      if (end->type==ID || end->type==NUM) end->mark++;
      if (end->type==OP && end->prec>0) end->mark++;
      if (end->type==LPAR) end->mark++;
      if (end->nest>pow->nest) end->mark++; }
	
    for (beg=pow->prev; beg; beg=beg->prev) {
      if (beg->type==ID || beg->type==NUM) beg->mark++;
      if (beg->type==OP && beg->prec>0) beg->mark++;
      if (beg->type==RPAR) beg->mark++;
      if (beg->nest>pow->nest) beg->mark++; }
	
    pow->mark++;

    for (l=pow; l; l=l->prev) {
      if (l->mark) beg=l;
      else break; }

    for (l=pow; l; l=l->next) {
      if (l->mark) end=l;
      else break; }

#ifdef DEBUG
    looplist (l,head)  {
      printf("%-12s nest=%d mark=%d type=%c prec=%d %s\n",l->name,l->nest,l->mark,"?NIO()"[l->type],l->prec,l==beg?"BEG":l==end?"END":"");
    }
    printf("===========\n");
#endif

    if (isint(pow->next)) {
      /* integer power */
      /* negative should not happen because - is detached */
      int i=atoi(pow->next->name),j;
      struct items_s *xbeg=beg,*xend;
      xend=pow->prev;
      delete(end);
      delete(pow);
      end=xend;

      if (i==0) {
	for (; xbeg!=end; ) xbeg=delete(xbeg); 
	xbeg->name="1"; }
      else {

	for (j=1; j<=i; j<<=1);
	j>>=1;
	while ((j>>=1)) { 
	  struct items_s *x,*e;

	  e=addafter(end,"*",OP);	e->prec=-1;
	  for (x=beg; ;x=x->next) {
	    e=addafter(e,x->name,x->type); e->prec=x->prec;
	    if (x==end) break; }
	  end=addafter(e,")",RPAR);
	  beg=addbefore(beg,"(",LPAR);
	  
	  if (j & i) {
	    beg=addbefore(beg,"(",LPAR);
	    end=addafter(end,"*",OP); end->prec=-1;
	    for (x=xbeg; ;x=x->next) {
	      end=addafter(end,x->name,x->type); end->prec=x->prec;
	      if (x==xend) break; } 
	    end=addafter(end,")",RPAR); } } }
    } else {
      pow->name=","; 
      pow->prec=1;
      end=addafter(end,")",RPAR);
      beg=addbefore(beg,"(",LPAR);
      beg=addbefore(beg,"pow",ID); }
  }

}

char line[16000];

int main(int narg,char **arg)
{
  int iarg;

  getsbufsize=16000;
  
  if (narg<2) {
    fprintf(stderr,"Converts expression with power written by ^ into C-syntax. Call by:\n\
  exp2pow \"EXPR\" [\"EXPR\"..]\n\
  exp2pow - < FILE-WITH-EXPR\n\
EXPR is an arithmetic expression in the usual syntax with power written by ^\n\
bit operators as & | ! etc. do not have the correct precedence and must be ()\n\
^ (bit) can be written as ^^\n\
a^#, where # is a positive decimal number, is translated into multiplications\n\
other expressions call pow(,)\n\
In addition, log is translated into log10 and ln to log\n\
Example:\n\
  exp2pow \"a+b*A^2+c*A^d\"\n\
gives\n\
  a+b*(A*A)+c*pow(A,d)\n\
Bugs:\n\
  the parser is simple, so use parentheses in case of doubts\n\
  A^-3 is translated to pow(A,-3), use 1/A^3 to emit a more efficient code\n\
  the output code uses more () than necessary\n");
    exit(0); }

  if (strcmp(arg[1],"-")) loop (iarg,1,narg) {
    struct items_s *l;

    exp2pow(arg[iarg]);

    looplist (l,head) printf("%s",l->name); 
    printf("\n"); }

  else while (gets(line)) {
    struct items_s *l;

    exp2pow(line);

    looplist (l,head) printf("%s",l->name); 
    printf("\n"); }

  return 0;
}
