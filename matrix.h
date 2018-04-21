#ifndef MATRIX_H
#define MATRIX_H
#include <assert.h>
#include <math.h>
#define HERMITE 0
#define BEZIER 1
#define PI  3.1415
#define TAO 6.2832
#define CONV PI / 180


struct matrix {
  double **m;
  int rows, cols;
  int lastcol;
} matrix;

//Basic matrix manipulation routines
struct matrix *new_matrix(int rows, int cols);
void free_matrix(struct matrix *m);
void grow_matrix(struct matrix *m, int newcols);
void copy_matrix(struct matrix *a, struct matrix *b);

//functions students need to complete
void print_matrix(char *name, struct matrix *m);
void ident(struct matrix *m);
//a*b -> new
struct matrix *matrix_mult(struct matrix *a, struct matrix *b);

//transformation routines
struct matrix * make_translate(double x, double y, double z);
struct matrix * make_scale(double x, double y, double z);
struct matrix * make_rotX(double theta);
struct matrix * make_rotY(double theta);
struct matrix * make_rotZ(double theta);

struct matrix *make_bezier();
struct matrix *make_hermite();

#endif
