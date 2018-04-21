#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gmath.h"
#include "matrix.h"

double *calculate_normal(double x0, double y0, double z0,
			 double x1, double y1, double z1,
			 double x2, double y2, double z2) {

  double *N = (double *)malloc(3 * sizeof(double));
  
  N[0] = (y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0);
  N[1] = (z1 - z0) * (x2 - x0) - (x1 - x0) * (z2 - z0);
  N[2] = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);

  return N;
}

double dot_product(double *v0, double *v1) {
  return v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];
}
