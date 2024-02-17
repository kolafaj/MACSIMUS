/* #include code for Shake
   (called before Verlet; note that there are FORCES in p[*])

1/ calculate and print force on the center-of-mass

2/ calculate and print the torque (force momentum) (with respect to
   the center-of-mass)

3/ print position of center-of-mass and velocity

4/ clear constraint forces for groups
*/

/* momentum (torque) */
if (mn->anchor & (ANCHOR_r | ANCHOR_v | ANCHOR_x | ANCHOR_f | ANCHOR_m)) {
  vector cm={0,0,0},molf={0,0,0},molv={0,0,0},dr,lever,moltorque;
  double m;

  loop (i,0,ns) {
    VV(cm,+=si[i].mass*r[i])
    VV(molf,+=p[i]) }
  m=spec[mn->sp]->mass;
  VO(cm,/=m)

  if (mn->anchor & ANCHOR_r) printanchor(1,cm,mn->xyz);

  if (mn->anchor & ANCHOR_v) {
    loop (i,0,ns) VV(molv,+=si[i].mass*r1[i])
    VO(molv,/=m)

    printanchor(1,molv,mn->xyz); }

  /* acceleration */
  if (mn->anchor & ANCHOR_x) printanchor(1./m,molf,mn->xyz);

  /* forces */
  if (mn->anchor & ANCHOR_f) printanchor(FORCEUNIT,molf,mn->xyz);

  /* torque */
  if (mn->anchor & ANCHOR_m) {
    VO(moltorque,=0)
    loop (i,0,ns) {
      VVV(dr,=r[i],-cm)
      VECT(lever,dr,p[i])

      VV(moltorque,+=lever) }

    printanchor(TORQUEUNIT,moltorque,mn->xyz); }
}

/* clearing CF of groups */
if (mn->anchor & ANCHOR_g) {
  struct anchorgroup_s *group;

  looplist (group,mn->group) VO(group->cf,=0)
}
