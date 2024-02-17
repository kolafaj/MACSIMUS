/* empty module userforces.c, #included by sim/xforces.c

   If a user provides MACSIMUS/cook/PROJECT/userforces.c, it will be
   used instead of this default empty one.

   Do NOT place userforces.c to MACSIMUS/sim, because it is read by xforces.c.
*/

void userforces(ToIntPtr B, ToIntPtr A) /************************ userforces */
{
  double U=0;

  En.pot += En.usr = U;

  //  if (measure) En.vir += f*r;
}
