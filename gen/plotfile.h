/* -1 is returned if the file does not exist, 0 on success */

int drawfile(char *fn,int colx,int coly,int color,int style);

int minmaxfile(char *fn,int colx,int coly,int cont,
	       double *X0,double *X1,double *Y0,double *Y1);

