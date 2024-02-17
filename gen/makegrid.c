/* used by nsk and plot */

void makegrid(int lbl,REAL density,                            /*** makegrid */
              int xaxis,REAL x0,REAL x1,REAL y0,REAL y1,int at)
/*
   warning: if !xaxis then x and y are swapped: x0,x1 is always the axis
   to be labeled 
   lbl=0 : no grid
   lbl=1 : light
   lbl=2 : more dense grid
   density : normally 3--4 for better graphics, 1--2 for worse
*/
{
  if (lbl) {
    REAL dd=1,x1x0=x1-x0;
    int ifmt=0, j=0,ix,at8,at9;
    static char fmt[8];
    int l,lfrom,lto;
    REAL ladd;

    setlinestyle(4,0x4444,1); /* user pattern: . . . . . */

    if (at<100) at8=at+8,at9=at8+1;
    else at8=at-2,at9=at8-1;

    if (x1x0<=0) return;

    while (x1x0 > 10*density*dd) { dd*=10; ifmt--; }
    while (x1x0 <    density*dd) { dd/=10; ifmt++; }

    if (ifmt>16) ifmt=16;
    sprintf(fmt,"%%.%df",ifmt);
    if (ifmt<0 || (fabs(x0)<1e-5 && fabs(x1)<1e-5) ) strcpy(fmt,"%g");

    ladd=fint(x0/(dd*10))*10;
    lfrom=Int(x0/dd-1e-4-ladd);
    lto=Int(x1/dd+1e-4-ladd);

    loopto (l,lfrom,lto) {
      REAL xi=(l+ladd)*dd;
      setcolor(lbl==2?LIGHTGRAY:DARKGRAY);
      if (xi+dd/2<x1+0.0001*x1x0)
	if (xaxis) lline(xi+dd/2,y0, xi+dd/2,y1);
	else lline(y0,xi+dd/2, y1,xi+dd/2);
      setcolor(l%10==0 || lbl==2 ? WHITE: LIGHTGRAY);
      if (l%5==0 || (lto-lfrom)<11) {
	char nr[24];
	int ii;

	sprintf(nr,fmt,Val(xi));
	if (xaxis) {
	  ii=SX(xi);
	  ix=ii-strlen(nr)*4+1;
	  if (ix>=j) {
	    outtextxy(ix,at,nr);
	    j=ix+strlen(nr)*8+8; }
	  putpixel(ii,at8,WHITE);
	  putpixel(ii-1,at8,WHITE);
	  putpixel(ii+1,at8,WHITE);
	  putpixel(ii,at9,WHITE); }
	else
	  outtextxy(at,SY(xi)-4,nr);
        }
      if (xi>x0-0.0001*x1x0)
	if (xaxis) lline(xi,y0,xi,y1);
	else lline(y0,xi,y1,xi); }
    setlinestyle(0,0,1);
    }
}
