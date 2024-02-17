/* order of sites = HOMH */
{
  QQX(&waters.ss.HH,r1[0],r2[0],f1[0],f2[0],waters.qq.HH GLOB);
  QQX(&waters.ss.HH,r1[0],r2[3],f1[0],f2[3],waters.qq.HH GLOB);
  QQX(&waters.ss.HH,r1[3],r2[0],f1[3],f2[0],waters.qq.HH GLOB);
  QQX(&waters.ss.HH,r1[3],r2[3],f1[3],f2[3],waters.qq.HH GLOB);
  QQX(&waters.ss.dummy,r1[2],r2[0],f1[2],f2[0],waters.qq.MH GLOB);
  QQX(&waters.ss.dummy,r1[0],r2[2],f1[0],f2[2],waters.qq.MH GLOB);
  QQX(&waters.ss.dummy,r1[2],r2[3],f1[2],f2[3],waters.qq.MH GLOB);
  QQX(&waters.ss.dummy,r1[3],r2[2],f1[3],f2[2],waters.qq.MH GLOB);
  QQX(&waters.ss.dummy,r1[2],r2[2],f1[2],f2[2],waters.qq.MM GLOB);
  LJX(&waters.ss.OO,r1[1],r2[1],f1[1],f2[1] GLOB);
  /* measuring OH correlation functions is in waters.c */
}
