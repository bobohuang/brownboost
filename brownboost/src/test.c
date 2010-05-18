#include "R.h"	
#include "Rmath.h"



void solvede(double* r, double* s, 
	     double* h, double* y, 
	     double* c, double* output)
{
  // do some math. //
  printf("in the c code\n");
  static double stest[2];
  stest[0] = r[0]*s[0];
  stest[1] = r[1]*s[1];

  // Use a normal curve function. //
  double normTest = 0;
  double x = 0.5;
  double mu = 0;
  double sigma = 1;
  normTest = pnorm(x, mu, sigma, 0, 0);
  printf("normTest: %lf\n", normTest);

  // Set the output. //
  output[0] = stest[0];
  output[1] = stest[1];

  return;
}

