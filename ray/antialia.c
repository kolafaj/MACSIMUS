{
float (**rgb)[3];
int nsmart=0;

rgb=(float(**)[3])malloc((yres+2)*sizeof(rgb[0]));
if (!rgb) exit(-1);
rgb++;
for (y = -1 ; y < yres+1 ; y++) {
  rgb[y]=(float(*)[3])malloc((xres+2)*sizeof(rgb[0][0]));
  if (!rgb[y]) { fprintf(stderr,"no heap\n"); exit(-1); }
  rgb[y]++; }
printf("picture allocated, normal scan:\n");

for (y = -1 ; y <= yres ; y++) {
  ylen = ( ((Flt) (2 * y) / (Flt) yres) - 1.0 ) * (aspect/scale) ;
  for (x = -1 ; x <= xres ; x++) {
    xlen = ( ((Flt) (2 * x) / (Flt) xres) - 1.0 ) / scale ;
    
    avg[0] = avg[1] = avg[2] = 0.0 ;
    
    VecComb(xlen , leftvec,
	    ylen , upvec, ray.D) ;
    VecAdd(ray.D, viewvec, ray.D) ;
    VecNormalize(ray.D);
    Trace(0, 1.0, &ray, color);
    VecAdd(color, avg, avg) ;
    avg[0] = min(1.0, avg[0]) ;
    avg[1] = min(1.0, avg[1]) ;
    avg[2] = min(1.0, avg[2]) ;
    rgb[y][x][0]=avg[0]; rgb[y][x][1]=avg[1]; rgb[y][x][2]=avg[2];
  }
  if (tickflag) fprintf(stderr, "%4.1f%%\r",100.0*y/yres) ;
}
if (tickflag) fprintf(stderr, "\n") ;
printf("antialiasing scan:\n") ;

for (y = 0 ; y < yres ; y++) {
  ylen = ( ((Flt) (2 * y) / (Flt) yres) - 1.0 ) * (aspect/scale) ;
  for (x = 0 ; x < xres ; x++) {
#define DIF(DX,DY) \
 (fabs(rgb[y-DY][x-DX][0]-2*rgb[y][x][0]+rgb[y+DY][x+DX][0])*0.3 \
 +fabs(rgb[y-DY][x-DX][1]-2*rgb[y][x][1]+rgb[y+DY][x+DX][1])*0.59\
 +fabs(rgb[y-DY][x-DX][2]-2*rgb[y][x][2]+rgb[y+DY][x+DX][2])*0.11)

    double contrast=(DIF(0,1)+DIF(1,0))*.125+(DIF(1,1)+DIF(1,(-1)))*0.1;
    
    xlen = ( ((Flt) (2 * x) / (Flt) xres) - 1.0 ) / scale ;

    if (contrast>threshold) {
      /* square jitter pattern */
      int xx,yy;
      if (square&1) {
	avg[0]=rgb[y][x][0]; avg[1]=rgb[y][x][1]; avg[2]=rgb[y][x][2]; }
      else
	avg[0] = avg[1] = avg[2] = 0.0 ;
      for (xx=0; xx<square; xx++)
	for (yy=0; yy<square; yy++) {
	  if ((square&1) && xx==square/2 && yy==square/2) continue;
	  VecComb(xlen + (Flt)(2*xx+1-square)/(Flt)(square) * xwidth, leftvec,
		  ylen + (Flt)(2*yy+1-square)/(Flt)(square) * ywidth, upvec, ray.D) ;
	  VecAdd(ray.D, viewvec, ray.D) ;
	  VecNormalize(ray.D);
	  Trace(0, 1.0, &ray, color);
	  VecAdd(color, avg, avg) ; }
      nsmart++;
      avg[0]*=(1./(square*square)); 
      avg[1]*=(1./(square*square)); 
      avg[2]*=(1./(square*square)); }
    else {
      avg[0]=rgb[y][x][0]; avg[1]=rgb[y][x][1]; avg[2]=rgb[y][x][2]; }
    avg[0] = min(1.0, avg[0]) ;
    avg[1] = min(1.0, avg[1]) ;
    avg[2] = min(1.0, avg[2]) ;
    PicWritePixel(pic, avg) ;
    
                }
  if (tickflag) fprintf(stderr, "%4.1f%%\r",100.0*y/yres) ;
}
if (tickflag) fprintf(stderr, "\a\r\n") ;
fprintf(stderr, "%.1f%% pixels re-sampled\n",100.*nsmart/(xres*yres)) ;
}
