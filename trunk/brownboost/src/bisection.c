#include "R.h"	
//#include "R.internals"
#include "Rmath.h"


// vector_plus_const //
// takes a vector and a constant 
// and returns new memory with result.
double* vpc (double* v, double c, int n) 
{
  int i = 0;
  double* result = malloc(n * sizeof(double));
  for(i = 0; i < n; i++) {
    result[i] = v[i] + c;
  }
  return(result);
}

// vector_plus_vector //
// takes a vector and a vector
// and returns new memory with result.
double* vpv (double* v1, double* v2, int n) 
{
  int i = 0;
  double* result = malloc(n * sizeof(double));
  for(i = 0; i < n; i++) {
    result[i] = v1[i] + v2[i];
  }
  return(result);
}


// vector_minus_vector //
// takes a vector and a vector
// and returns new memory with result.
double* vmv (double* v1, double* v2, int n) 
{
  int i = 0;
  double* result = malloc(n * sizeof(double));
  for(i = 0; i < n; i++) {
    result[i] = v1[i] - v2[i];
  }
  return(result);
}


// vector_times_vector //
// takes two vectors and returns a new one //
double* vtv (double* v1, double* v2, int n) 
{
  int i = 0;
  double* result = malloc(n * sizeof(double));
  for(i = 0; i < n; i++) {
    result[i] = v1[i] * v2[i];
  }
  return(result);
}


// vector_times_const //
// takes two vectors and returns a new one //
double* vtc (double* v, double c, int n) 
{
  int i = 0;
  double* result = malloc(n * sizeof(double));
  for(i = 0; i < n; i++) {
    result[i] = v[i] * c;
  }
  return(result);
}


// vector_divided_by_const //
// takes two vectors and returns a new one //
double* vdc (double* v, double c, int n) 
{
  int i = 0;
  double* result = malloc(n * sizeof(double));
  for(i = 0; i < n; i++) {
    result[i] = v[i] / c;
  }
  return(result);
}


// vector_sum //
// takes one vector and returns the sum of elements //
double vsum (double* v, int n) 
{
  int i = 0;
  double result = 0.0;
  for(i = 0; i < n; i++) {
    result += v[i];
  }
  return(result);
}



// special dot product function ... c == -1 //
double* dotl (double* x, double* b, double v, int n) {
  // R Code //
  // (v[1]*l[[1]] + v[2]*l[[2]]) //
  // Where v = c(alpha, tee) and l == [[b, -1]]

  double* inner = vtc(b, x[0], n);
  double c = (x[1] * v);
  double* result = vpc(inner, c, n);
  //printf("dotl inner: %lf\t%lf\t%lf\n", inner[0], inner[1], inner[2]);
  //printf("dotl c: %lf", c);
  //printf("dotl result: %lf\t%lf\t%lf\n", result[0], result[1], result[2]);
  free(inner);
  return(result);
}


// square each element of a vector //
double* vpow (double* v, int e, int n) {
  int i = 0;
  double* result = malloc(n * sizeof(double));
  for(i = 0; i < n; i++) {
    result[i] = pow(v[i], e);
  }
  return(result);
}


// exp each element of a vector //
double* expv (double* v, int n) {
  int i = 0;
  double* result = malloc(n * sizeof(double));
  for(i = 0; i < n; i++) {
    result[i] = exp(v[i]);
  }
  return(result);
}


// returns the sign of a double //
int fastsign (double x) {
  return((x > 0) - (x < 0));
  //  if (x > 0) return 1;
  //if (x < 0) return -1;
  //return 0;
}


// The error function
double errfun (double* x) {
  // R code (2*pnorm(a*sqrt(2)) - 1) //  Wrong tail here???
  return( (double) (2.0 * pnorm((*x) * M_SQRT2, 0, 1, 0, 0) - 1.0) * -1 );
}


// The error function for vectors
double* verrfun (double* v, int n) {
  // R code (2*pnorm(a*sqrt(2)) - 1) //
  int i = 0;
  double* result = malloc(n * sizeof(double));
  for(i = 0; i < n; i++) {
    result[i] = errfun(&v[i]);
  }
  return(result);
}


// The first equation we're trying to zero... //
// a = r *c and b = h(x) * y and v = -1s and x is the point (alpha, tee) //
double zero1 (double* a, double* b, 
	      double* v, double* x,
	      double* c, int n) 
{
  // R code
  //  f1 <- sum(b * exp( -(1/c) * (a + dotl(z, v))^2))
  double result = 0;
  double cd = -(1.00/(*c));
  double* inner1 = dotl(x, b, *v, n);    // dot product
  double* inner2 = vpv(a, inner1, n);   // add vectors
  double* inner3 = vpow(inner2, 2, n);         // square vector
  double* inner4 = vtc(inner3, cd, n); // const multiplier
  double* inner5 = expv(inner4, n);        // exp vector
  double* inner6 = vtv(inner5, b, n);
  double final  = vsum(inner6, n);

  /*
  printf("\nzero1\n");
  printf("x: %lf\t%lf\n", x[0], x[1]);
  printf("b: %lf\t%lf\t%lf\n", b[0], b[1], b[2]);
  printf("v: %lf\n", v[0]);
  printf("%lf\t%lf\t%lf\n", inner1[0], inner1[1], inner1[2]);
  printf("%lf\t%lf\t%lf\n", inner2[0], inner2[1], inner2[2]);
  printf("%lf\t%lf\t%lf\n", inner3[0], inner3[1], inner3[2]);
  printf("%lf\t%lf\t%lf\n", inner4[0], inner4[1], inner4[2]);
  printf("%lf\t%lf\t%lf\n", inner5[0], inner5[1], inner5[2]);
  printf("%lf\t%lf\t%lf\n", inner6[0], inner6[1], inner6[2]);
  printf("final: %lf\n", final);
  printf("done\n");
  */

  free(inner1);
  free(inner2);
  free(inner3);
  free(inner4);
  free(inner5);
  free(inner6);
  return(final);
}



// The second equation we're trying to zero... //
// a = r *c and b = h(x) * y and v = -1s and x is the point (alpha, tee) //
double zero2 (double* a, double* b, 
	      double* v, double* x,
	      double* c, int n) 
{
  // R code
  // f2 <- sum(erf((a + dotl(z, v))/sqrt(c)) - erf(a/sqrt(c)))
  double result = 0;
  double sqrtc = sqrt(*c);
  double* inner1 = vdc(a, sqrtc, n); 
  double* inner2 = verrfun(inner1, n);
  double* inner3 = dotl(x, b, *v, n);
  double* inner4 = vpv(inner3, a, n);
  double* inner5 = vdc(inner4, sqrtc, n);
  double* inner6 = verrfun(inner5, n);
  double* inner7 = vmv(inner6, inner2, n);
  result = vsum(inner7, n);

  /*
  printf("\nzero2\n");
  printf("x: %lf\t%lf\n", x[0], x[1]);
  printf("b: %lf\t%lf\t%lf\n", b[0], b[1], b[2]);
  printf("v: %lf\n", v[0]);
  printf("sqrtc: %lf\n", sqrtc);
  printf("%lf\t%lf\t%lf\n", inner1[0], inner1[1], inner1[2]);  // a/sqrt(c)
  printf("%lf\t%lf\t%lf\n", inner2[0], inner2[1], inner2[2]);  // erf(a/sqrt(c))
  printf("%lf\t%lf\t%lf\n", inner3[0], inner3[1], inner3[2]);  
  printf("%lf\t%lf\t%lf\n", inner4[0], inner4[1], inner4[2]);
  printf("%lf\t%lf\t%lf\n", inner5[0], inner5[1], inner5[2]);
  printf("%lf\t%lf\t%lf\n", inner6[0], inner6[1], inner6[2]);
  printf("%lf\t%lf\t%lf\n", inner7[0], inner7[1], inner7[2]);
  printf("final: %lf\n", result);
  printf("done\n");
  */

  free(inner1);
  free(inner2);
  free(inner3);
  free(inner4);
  free(inner5);
  free(inner6);
  free(inner7);
  return(result);
}


// a = r *c and b = h(x) * y and v = -1s and x is the point //
void bigfun (double* a, double* b, 
	     double* v, double* x,
	     double* c,
	     int* sign, int* n)
{
  double result = 0;  
  sign[0] = fastsign(zero1(a, b, v, x, c, *n));
  sign[1] = fastsign(zero2(a, b, v, x, c, *n));
  return;
}


double* getNewPoints (double* a, double* b, 
		      double* v, double* c, int n)
{
  GetRNGstate();
  double* result = malloc(4 * sizeof(double));
  int* signx1 = malloc(2 * sizeof(int));  
  int* signx2 = malloc(2 * sizeof(int));  
  double x1[2];  x1[0] = 0;  x1[1] = 0;
  double x2[2];  x2[0] = 0;  x2[1] = 0;
  signx1[0] = -1;  signx1[1] = 1;
  signx2[0] = 1;  signx2[1] = -1;
  while ((signx1[0] == signx2[0] 
	  || signx1[1] == signx2[1]
	  || signx1[0] != signx1[1]
	  || signx2[0] != signx2[1])) {
    x1[0] = unif_rand() * 2;  x1[1] = unif_rand() * 2;
    x2[0] = -unif_rand() * 2;  x2[1] = -unif_rand() * 2;
    bigfun(a, b, v, x1, c, signx1, &n);
    bigfun(a, b, v, x2, c, signx2, &n);
  }
  free(signx1);
  free(signx2);
  result[0] = x1[0]; result[1] = x1[1]; result[2] = x2[0]; result[3] = x2[1];
  PutRNGstate();
  return(result); // remember to free this memory!
}


double* decideOnNewPoint(double* pts, double* a, double* b, 
			  double* v, double* c, int n)
{
  int i = 0;
  int signv[2];                            // the sign of the new point set
  double newerr[2];
  int oldsign[2];                          // sign of incoming points
  double olderror[4];
  double newpts[14];                          // set of new points
  double* midpt1 = vpv(pts, &pts[2], 2);      // the midpoint
  double* mid    = vdc(midpt1, 2, 2);       
  olderror[0] = zero1(a, b, v, pts, c, n);
  olderror[1] = zero2(a, b, v, pts, c, n);
  olderror[2] = zero1(a, b, v, &pts[2], c, n);
  olderror[3] = zero2(a, b, v, &pts[2], c, n);
  newpts[0] = mid[0];  newpts[1] = mid[1];
  newpts[2] = mid[0];  newpts[3] = pts[3];
  newpts[4] = pts[2];  newpts[5] = mid[1];
  newpts[6] = mid[0];  newpts[7] = pts[1];
  newpts[8] = pts[0];  newpts[9] = mid[1];
  newpts[10] = pts[2];  newpts[11] = pts[1];
  newpts[12] = pts[0];  newpts[13] = pts[3];
                                                      // *For Sure* // 
  bigfun(a, b, v, &newpts[i], c, oldsign, &n);  // other sign is opposite //
  free(midpt1);
  free(mid);

  for (i = 0; i < 14; i += 2) {
   
    // get new signs for this point .. save in signv
    bigfun(a, b, v, &newpts[i], c, signv, &n);
    newerr[0] = zero1(a, b, v, &newpts[i], c, n);
    newerr[1] = zero2(a, b, v, &newpts[i], c, n);
    // if the signs are the same and oppposite of the old sign //
    // then let's take it. //
    if ((signv[0] == signv[1]) && (signv[0] != oldsign[0])
	&& pow(newerr[0], 2) + pow(newerr[1],2) < pow(olderror[2],2) + pow(olderror[3],2)) {
      pts[2] = newpts[i];  pts[3] = newpts[i+1];
      return(pts);
    // Otherwise, if they are the same, then opposite of other point //
    } else if (signv[0] == signv[1]
	       &&  pow(newerr[0], 2) + pow(newerr[1],2) < pow(olderror[0],2) + pow(olderror[1],2)) {
      pts[0] = newpts[i];  pts[1] = newpts[i+1];
      return(pts);
    }
  }
  printf("Couldn't find any new points!\n");
  exit(1);
  return NULL;
}



void solvede(double* r, double* s, 
	     double* h, double* y, 
	     double* c,    int* n, 
	     double* output)
{
  // set up a = r *c and b = h(x) * y. //
  double* a = vpc(r, *s, *n);
  double* b = vtv(h, y, *n);
  
  // v is the list containing b, and -1 //
  double* v1 = malloc(1 * sizeof(double));  v1[0] = -1.0;
  int tries = 0;
  int loops = 0;
  double* points = NULL;

  while (tries < 100) {

    // get starting points //
    points = getNewPoints (a, b, v1, c, *n); // needs to be freed.

    while (loops < 111) {
      points = decideOnNewPoint(points, a, b, v1, c, *n);

      // Are the new points sufficently close to the previous? //
      if (fabs(points[0] - points[2]) < 0.0000000001
	  && fabs(points[1] - points[3]) < 0.0000000001) {
	output[0] = points[0];
	output[1] = points[1];
	printf("points: %lf\t%lf\t%lf\t%lf\n", points[0], points[1], points[2], points[3]);
	printf("zero1: %lf\n", zero1(a, b, v1, output, c, *n));
	printf("zero2: %lf\n", zero2(a, b, v1, output, c, *n));
	return;
      }

      loops += 1;
    }

    free(points);
    loops = 0;
    tries += 1;
  }

  // Set the output. //
  output[0] = -1;
  output[1] = -1;

  free(a);
  free(b);
  free(v1);

  return;
}



  



