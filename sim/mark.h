extern struct mark_s {
  double LJ,el/*,bond*/; /* bonded not included -- later?! */
  char id1[32],id2[32];
  char *set; } *mark_tab; 

/* 
   mark_tab[mark_n+2] : 
   mark_tab[0] = reference molecule
   mark_tab[1..mark_n] = molecules to count all LJ+el pairs with the ref 
   mark_tab[mark_n+1].set=NULL is the sentinel 
*/

extern double mark_el,U; /* U needed by MARK_PATCH in functions not
                            having internal U*/

void mark_init(void);
void mark_do(int no,int endian);

/* to be included to all xxxM() in interpot.c (incl. those that cannot
   give any contribution because of  mark_el=En.el; statement) */
#define MARK_PATCH { \
  struct mark_s *mt=mark_tab; \
  int i=(vector*)r1-a[0]->rp,j=(vector*)r2-a[0]->rp,n; \
  if (!mt->set[i]) n=j,j=i,i=n; \
  if (mt->set[i]) for (mt++; mt->set; mt++)\
    if (mt->set[j]) mt->LJ+=U,mt->el+=En.el-mark_el; \
  mark_el=En.el; }
