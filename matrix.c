#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"


/*-------------- void print_matrix() --------------
Inputs:  struct matrix *m 
Returns: 

print the matrix
*/
void print_matrix(char *name, struct matrix *m) {
  printf("%s: (%d x %d), %d cols allocated;\n", name, m->rows, m->lastcol+1, m->cols);
  
  int r, c;
  for (r = 0; r < m->rows; r++) {
    for (c = 0; c <= m->lastcol; c++) {
      printf("%.2f ", m->m[r][c]);
    }
    printf("\n");
  }
  printf("\n");
}

/*-------------- void ident() --------------
Inputs:  struct matrix *m <-- assumes m is a square matrix
Returns: 

turns m in to an identity matrix
*/
void ident(struct matrix *m) {
  m->lastcol = m->rows-1;

  int r, c;
  
  for (r = 0; r < m->rows; r++) {
    for (c = 0; c < m->rows; c++) {
      if (r == c)
	m->m[r][c] = 1;
      else
	m->m[r][c] = 0;
    }
  }
}


  

/*-------------- void matrix_mult() --------------
Inputs:  struct matrix *a
         struct matrix *b 
Returns: 

a*b -> new
*/
struct matrix *matrix_mult(struct matrix *a, struct matrix *b) {
  /*  printf("a: ");
  print_matrix(a);
  printf("b: ");
  print_matrix(b);*/
  
  assert(a->lastcol == b->rows-1);
  struct matrix *out = new_matrix(a->rows, b->lastcol + 1);
  out->lastcol = b->lastcol;
  
  int r, r2, c;
  for (r = 0; r < a->rows; r++) {
    for (c = 0; c <= b->lastcol; c++) {
      double sum = 0;
      
      for (r2 = 0; r2 < b->rows; r2++) {
	sum += a->m[r][r2] * b->m[r2][c];
      }

      out->m[r][c] = sum;
    }
  }

  /*  printf("result: ");
  print_matrix(out);
  printf("\n");*/
  
  return out;
}



/*===============================================
  These Functions do not need to be modified
  ===============================================*/

/*-------------- struct matrix *new_matrix() --------------
Inputs:  int rows
         int cols 
Returns: 

Once allocated, access the matrix as follows:
m->m[r][c]=something;
if (m->lastcol)... 
*/
struct matrix *new_matrix(int rows, int cols) {
  double **tmp;
  int i, j;
  struct matrix *m;

  tmp = (double **)malloc(rows * sizeof(double *));
  for (i=0;i<rows;i++) {
      tmp[i]=(double *)malloc(cols * sizeof(double));
      for (j = 0; j < cols; j++) {
	tmp[i][j] = 0;
      }
  }

  m=(struct matrix *)malloc(sizeof(struct matrix));
  m->m=tmp;
  m->rows = rows;
  m->cols = cols;
  m->lastcol = -1;

  return m;
}


/*-------------- void free_matrix() --------------
Inputs:  struct matrix *m 
Returns: 

1. free individual rows
2. free array holding row pointers
3. free actual matrix
*/
void free_matrix(struct matrix *m) {

  int i;
  for (i=0;i<m->rows;i++) {
      free(m->m[i]);
    }
  free(m->m);
  free(m);
}


/*======== void grow_matrix() ==========
Inputs:  struct matrix *m
         int newcols 
Returns: 

Reallocates the memory for m->m such that it now has
newcols number of collumns
====================*/
void grow_matrix(struct matrix *m, int newcols) {
  
  int i;
  for (i=0;i<m->rows;i++) {
      m->m[i] = realloc(m->m[i],newcols*sizeof(double));
  }
  m->cols = newcols;
}


/*-------------- void copy_matrix() --------------
Inputs:  struct matrix *a
         struct matrix *b 
Returns: 

copy matrix a to matrix b
*/
void copy_matrix(struct matrix *a, struct matrix *b) {
  assert(a->lastcol < b->cols);

  int r, c;

  for (r=0; r < a->rows; r++) 
    for (c=0; c <= a->lastcol; c++)  
      b->m[r][c] = a->m[r][c];  

  b->lastcol = a->lastcol;

  
}

struct matrix *make_translate(double x, double y, double z) {
  struct matrix *trans = new_matrix(4, 4);
  ident(trans);
  trans->m[0][3] = x;
  trans->m[1][3] = y;
  trans->m[2][3] = z;
  return trans;
}

struct matrix *make_scale(double x, double y, double z) {
  struct matrix *trans = new_matrix(4, 4);
  ident(trans);
  trans->m[0][0] = x;
  trans->m[1][1] = y;
  trans->m[2][2] = z;
  return trans;
}

struct matrix *make_rotX(double theta) {
  struct matrix *trans = new_matrix(4, 4);
  ident(trans);
  trans->m[1][1] = cos(theta * CONV);
  trans->m[1][2] = -1 * sin(theta * CONV);
  trans->m[2][1] = sin(theta * CONV);
  trans->m[2][2] = cos(theta * CONV);
  return trans;
}

struct matrix *make_rotY(double theta) {
  struct matrix *trans = new_matrix(4, 4);
  ident(trans);
  trans->m[2][2] = cos(theta * CONV);
  trans->m[2][0] = -1 * sin(theta * CONV);
  trans->m[0][2] = sin(theta * CONV);
  trans->m[0][0] = cos(theta * CONV);
  return trans;
}

struct matrix *make_rotZ(double theta) {
  struct matrix *trans = new_matrix(4, 4);
  ident(trans);
  trans->m[0][0] = cos(theta * CONV);
  trans->m[0][1] = -1 * sin(theta * CONV);
  trans->m[1][0] = sin(theta * CONV);
  trans->m[1][1] = cos(theta * CONV);
  return trans;
}

/*======== struct matrix * make_bezier() ==========
  Inputs:
  Returns: The correct 4x4 matrix that can be used
  to generate the coefiecients for a bezier curve
  ====================*/
struct matrix * make_bezier() {
  struct matrix *bezier = new_matrix(4, 4);
  ident(bezier);
  bezier->m[1][0] = -3;
  bezier->m[1][1] = 3;
  bezier->m[2][0] = 3;
  bezier->m[2][1] = -6;
  bezier->m[2][2] = 3;
  bezier->m[3][0] = -1;
  bezier->m[3][1] = 3;
  bezier->m[3][2] = -3;
  
  return bezier;
}

/*======== struct matrix * make_hermite() ==========
  Inputs:
  Returns:
>
  The correct 4x4 matrix that can be used to generate
  the coefiecients for a hermite curve
  ====================*/
struct matrix * make_hermite() {
  struct matrix *hermite = new_matrix(4, 4);
  hermite->m[0][0] = 1;
  hermite->m[1][2] = 1;
  hermite->m[2][0] = -3;
  hermite->m[2][1] = 3;
  hermite->m[2][2] = -2;
  hermite->m[2][3] = -1;
  hermite->m[3][0] = 2;
  hermite->m[3][1] = -2;
  hermite->m[3][2] = 1;
  hermite->m[3][3] = 1;
  hermite->lastcol = 3;

  return hermite;

  
}
