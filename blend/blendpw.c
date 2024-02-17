#define Rcos 0.587728588833071
#define Rsin 0.759061990794090

const vector relHH[12][2]={
  { {0, Rcos, Rsin},{0, Rcos,-Rsin} },
  { {0,-Rcos, Rsin},{0,-Rcos,-Rsin} },
  { {0, Rsin, Rcos},{0,-Rsin, Rcos} },
  { {0, Rsin,-Rcos},{0,-Rsin,-Rcos} },
  { { Rcos, Rsin,0},{ Rcos,-Rsin,0} },
  { {-Rcos, Rsin,0},{-Rcos,-Rsin,0} },
  { { Rsin, Rcos,0},{-Rsin, Rcos,0} },
  { { Rsin,-Rcos,0},{-Rsin,-Rcos,0} },
  { { Rcos,0, Rsin},{ Rcos,0,-Rsin} },
  { {-Rcos,0, Rsin},{-Rcos,0,-Rsin} },
  { { Rsin,0, Rcos},{-Rsin,0, Rcos} },
  { { Rsin,0,-Rcos},{-Rsin,0,-Rcos} } };

static void makeHH(vector *rO,int n)
{
int i;
loop (i,0,2)
  VVV(rO[i+1],=rO[0],+relHH[n][i])
}
