#ifndef loop

/* C-style loop (last not included) */
#define loop(I,FROM,TO) for ((I)=(FROM); (I)<(TO); (I)++)

/* Pascal etc. style loop (last included) */
#define loopto(I,FROM,TO) for ((I)=(FROM); (I)<=(TO); (I)++)

/* loop over linked list.  Example:

   struct mylist_s {
     struct mylist_s *next; // next item; NULL terminates the list
     int DATA;
   } *head=NULL,*l;

   int i;

   // building a list in reverse order
   loop (i,0,10) {
     allocone(l);
     l->DATA=i;
     l->next=head;
     head=l; }

   // building a list in normal order
   loop (i,0,10) {
     if (head) {
       allocone(l->next);
       l=l->next; }
     else {
       allocone(head);
       l=head; }
     l->next=NULL;
     l->DATA=i; }

   // loop over the list
   looplist (l,head) printf("%d\n",l->DATA);

   // remove the list (no ; at end)
   freelist(l,head)

  // SEE ALSO: macsimus/c/alg/linklist.c
*/
#define looplist(PTR,HEAD) for ((PTR)=(HEAD); (PTR); (PTR)=(PTR)->next)
#define freelist(PTR,HEAD) for ((PTR)=(HEAD); (PTR); (PTR)=(HEAD)) { \
    (HEAD)=(PTR)->next; \
    free(PTR); (PTR)=NULL; } \
  (HEAD)=NULL;

#endif
