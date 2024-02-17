/* order of sites = HHO */
{
  QQX(&waters.ss.HH,r1[0],r2[0],f1[0],f2[0],waters.qq.HH GLOB);
  QQX(&waters.ss.HH,r1[0],r2[1],f1[0],f2[1],waters.qq.HH GLOB);
  QQX(&waters.ss.HH,r1[1],r2[0],f1[1],f2[0],waters.qq.HH GLOB);
  QQX(&waters.ss.HH,r1[1],r2[1],f1[1],f2[1],waters.qq.HH GLOB);
  QQX(&waters.ss.OH,r1[2],r2[0],f1[2],f2[0],waters.qq.MH GLOB);
  QQX(&waters.ss.OH,r1[0],r2[2],f1[0],f2[2],waters.qq.MH GLOB);
  QQX(&waters.ss.OH,r1[2],r2[1],f1[2],f2[1],waters.qq.MH GLOB);
  QQX(&waters.ss.OH,r1[1],r2[2],f1[1],f2[2],waters.qq.MH GLOB);
LJQQX(&waters.ss.OO,r1[2],r2[2],f1[2],f2[2],waters.qq.MM GLOB);
}
