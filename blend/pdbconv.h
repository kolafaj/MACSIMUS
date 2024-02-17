void connect(site_t *s1,site_t *s2); /* creates bond between sites s1 and s2 */

site_t *patch(site_t *pat,site_t *site,int merge); /* patch pat to site */

/*.....void mergeres(residue_t *r1,residue_t *r2);*/

int renumber(site_t *site,int offset); /* renumber sites by adding `offset' */

site_t *findid(residue_t *res,char *id1,char *id2); /* finds 1st `id1' [id2] from res->from... */

int delete(site_t *del); /* deletes site *del and all its references */

int altlocations(void); /* solve alternate locations */

void findCYSCYS(void); /* find CYS-CYS bonds by distances */

int sitebyoption(char *id,residue_t *res,int opt); /* option -d,-i support */

typedef int remove_f(site_t *);

int freeH(site_t *s); /* returns 1 if s is free H */

int removesites(remove_f *condition,char *msg); /* removes sites by condition */

