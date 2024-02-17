/***********************************************************************
 * $Log:	pic.h,v $
 * Revision 1.1  88/09/11  11:00:50  markv
 * Initial revision
 * 
 * Revision 1.1  88/09/09  11:59:56  markv
 * Initial revision
 * 
 ***********************************************************************/
typedef struct Pic {
	char	* filename ;
	FILE	* filep ;
	int	x, y ;
} Pic ;

/* JK: picture from pic (ppm) file (for background) */
typedef struct pic_s {
  int x,y;
  Flt xrange,yrange;
  Flt xshift,yshift;
  unsigned char (*c)[3]; } pic_t;

extern pic_t *bg; /* background picture */
pic_t *ReadPic(char *fn);
void GetPic(pic_t *pic,Flt x,Flt y,int mode,Vec c);

