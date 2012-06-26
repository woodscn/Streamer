#include <stdio.h>
int pnpoly(int npol, float *xp, float *yp, float x, float y)
{
  int i, j, c = 0;
  for (i = 0, j = npol-1; i < npol; j = i++) {
    printf("i = %d, j = %d\n",i,j);
    if ((((yp[i] <= y) && (y < yp[j])) ||
	 ((yp[j] <= y) && (y < yp[i]))) &&
	(x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
      c = !c;
  }
  return c;
}

int main()
{
  printf( "I am alive!  Beware.\n" );
  float yp[4];
  float xp[4];
  yp[0] = 0.0;
  yp[1] = 0.0;
  yp[2] = 1.0;
  yp[3] = 1.0;
  xp[0] = 0.0;
  xp[1] = 1.0;
  xp[2] = 1.0;
  xp[3] = 0.0;
  printf("%f\n",yp[0]);
  printf("%i\n",pnpoly(4,xp,yp,.25,.25));
  return 0;

}
