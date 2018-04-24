#include <stdio.h>
#include <stdlib.h>

#include "ml6.h"
#include "gmath.h"
#include "display.h"
#include "draw.h"
#include "matrix.h"

//function of step2:
#define OFFSET 0 //(double) rand() / (double) RAND_MAX
#define COLOR_VAR 80 //variance of colors in solid

/*======== void add_point() ==========
Inputs:   struct matrix * points
         int x
         int y
         int z 
Returns: 
adds point (x, y, z) to points and increment points.lastcol
if points is full, should call grow on points
====================*/
void add_point(struct matrix * points, double x, double y, double z) {
  assert(points->rows == 4);

  //printf("point: (%.1f, %.1f, %.1f)\n", x, y, z);
  
  int col = points->lastcol + 1;

  if (col >= points->cols) {
    grow_matrix(points, col * 2);
  }

  points->m[0][col] = x;
  points->m[1][col] = y;
  points->m[2][col] = z;
  points->m[3][col] = 1;

  points->lastcol++;
}

/*======== void add_edge() ==========
Inputs:   struct matrix * points
          int x0, int y0, int z0, int x1, int y1, int z1
Returns: 
add the line connecting (x0, y0, z0) to (x1, y1, z1) to points
should use add_point
====================*/

int check_even(struct matrix *points) {
  if (points->lastcol % 2 == 0) {

    add_point(points, 0, 0, 0);
    return 1;
  }
  return 0;
}

void add_edge(struct matrix * points, 
	       double x0, double y0, double z0, 
	       double x1, double y1, double z1) {
  check_even(points);
  
  add_point(points, x0, y0, z0);
  add_point(points, x1, y1, z1);
}

void add_polygon(struct matrix * points,
		 double x0, double y0, double z0,
		 double x1, double y1, double z1,
		 double x2, double y2, double z2) {

  //  printf("Polygon:\n");
  add_point(points, x0, y0, z0);
  add_point(points, x1, y1, z1);
  add_point(points, x2, y2, z2);
  //  printf("End polygon\n");
}

void add_matrix(struct matrix *points, struct matrix *mo_points) {

  assert(points->rows == 4 && mo_points->rows == 4);
  check_even(points);

  int col;
  for (col = 0; col <= mo_points->lastcol; col++) {

    add_point(points, mo_points->m[0][col],
	      mo_points->m[1][col], mo_points->m[2][col]);
  }

}	   

/*======== void draw_lines() ==========
Inputs:   struct matrix * points
         screen s
         color c 
Returns: 
Go through points 2 at a time and call draw_line to add that line
to the screen
====================*/
void draw_lines(struct matrix * points, screen s, color c, depthmap d) {
  assert(points->rows == 4);
  
  int i;

  for (i = 0; i < (points->lastcol + 1) / 2; i++) {
    draw_line(points, 2*i, s, c, d);
  }
}


void draw_shallow(int x0, int y0, double z0, int x1, int y1, double z1, screen s, color c, depthmap d) {
  //ternary operators account for both positive and negative slopes

  int dx = abs(x1 - x0);
  int ddx = (x1 >= x0)? 1: -1;

  int dy = abs(y0 - y1);
  int ddy =  (y1 >= y0)? 1: -1;

  double dzdx = (z1 - z0) / dx;
  double z = z0 - dzdx;
  
  int line_sum = 2 * dy - dx;
  int y = y0;
  int x = x0 - ddx;

  while ((x += ddx) != x1) {
    plot(s, c, d, x, y, (z += dzdx));

    if (line_sum > 0) {
      y += ddy;
      line_sum -= 2 * dx;
    }

    line_sum += 2 * dy;
  }
}

void draw_steep(int x0, int y0, double z0, int x1, int y1, double z1, screen s, color c, depthmap d) {
  
  //x and y switched from abovex
  int dx = abs(x1 - x0);
  int ddx = (x1 >= x0)? 1: -1;

  int dy = abs(y0 - y1);
  int ddy =  (y1 >= y0)? 1: -1;

  double dzdy = (z1 - z0) / dy;
  double z = z0 - dzdy;
  
  int line_sum = 2 * dx - dy;
  int x = x0;

  int y = y0 - ddy;

  while ((y += ddy) != y1) {
    plot(s, c, d, x, y, (z += dzdy));
    
    if (line_sum > 0) {
      x += ddx;
      line_sum -= 2 * dy;
    }
    
    line_sum += 2 * dx;
  }
}

void draw_line(struct matrix *points, int pos, screen s, color c, depthmap d) {
  int x0 = points->m[0][pos];
  int y0 = points->m[1][pos];
  double z0 = points->m[2][pos];
  
  int x1 = points->m[0][pos+1];
  int y1 = points->m[1][pos+1];
  double z1 = points->m[2][pos+1];


  if (x0 == x1) {
    if (y0 == y1)
      return;
    int x = x0;
    int dy = y1 > y0? 1: -1;
    int y = y0;
    double dz = (z1 - z0) / abs(y1 - y0);
    double z = z0;
    while (y != y1) {
      plot(s, c, d, x, y, z);
      y += dy;
      z += dz;
    }
    plot(s, c, d, x, y, z);
  } else if (abs(y1 - y0) < abs(x1 - x0)) {
    draw_shallow(x0, y0, z0, x1, y1, z1, s, c, d);
  } else {
    draw_steep(x0, y0, z0, x1, y1, z1, s, c, d);
  }
  /*color cr;
  cr.red = 255;
  cr.green = 0;
  cr.blue = 0;
  plot(s, cr, x0, y0);
  plot(s, cr, x1, y1);*/
}
/*======== int draw_polygon() ========
  Inputs:  double x0, y0, z0
           double x1, y1, z1
	   double x2, y2, z2
	   screen s, color c
  Returns: 1 if front-facing (drawn), 0 if back-facing (ignored)
  
  draws triangle with three given vertices
  =========*/
int draw_polygon(struct matrix *points, int pos, screen s, color c, depthmap d) {
  //  printf("Poly: (%.1f, %.1f, %.1f)\n(%.1f, %.1f, %.1f)\n(%.1f, %.1f, %.1f)\n", x0, y0, z0, x1, y1, z1, x2, y2, z2);
  /*back-face culling assumes that all 3d shapes are solids, and thus only the front faces must be drawn.
     This is computed by taking the cross product of v1 - v0 and v2 - v0
     If the z component is positive, it faces forwards and should be displayed
     if negative, it faces away and will be hidden
     if zero, it is perpendicular to the plane of viewing (also invisible)

     (v1 - v0) x (v2 - v0) = <smth, smth, (x1 - x0)(y2 - y0) - (y1 - y0)(x2 - x0)>*/


  /*  double view_angle[3] = {0, 0, 1};
  
  if (dot_product(view_angle, calculate_normal(points, pos)) > 0) {
  return 0;*/
  struct matrix *edges = new_matrix(4, 6);
  add_edge(edges, points->m[0][pos+0], points->m[1][pos+0], points->m[2][pos+0],
	   points->m[0][pos+1], points->m[1][pos+1], points->m[2][pos+1]);
  add_edge(edges, points->m[0][pos+1], points->m[1][pos+1], points->m[2][pos+1],
	   points->m[0][pos+2], points->m[1][pos+2], points->m[2][pos+2]);
  add_edge(edges, points->m[0][pos+2], points->m[1][pos+2], points->m[2][pos+2],
	   points->m[0][pos+0], points->m[1][pos+0], points->m[2][pos+0]);
  draw_lines(edges, s, c, d);
  free(edges);
  
  return 1;
}

int draw_filled_triangle(struct matrix *points, int pos,
			 screen s, color c, depthmap d) {
  double view_angle[3] = {0, 0, 1};
  
  if (dot_product(view_angle, calculate_normal(points, pos)) > 0) {
    return 0;
  }

  
  int *left, *middle, *right;
  int p0[3] = {points->m[0][pos], points->m[1][pos], points->m[2][pos]};
  int p1[3] = {points->m[0][pos+1], points->m[1][pos+1], points->m[2][pos+1]};
  int p2[3] = {points->m[0][pos+2], points->m[1][pos+2], points->m[2][pos+2]};
  if (p0[0] >= p1[0] && p0[0] >= p2[0]) {
    right = p0;
    if (p1[0] >= p2[0]) {
      middle = p1;
      left = p2;
    } else {
      middle = p2;
      left = p1;
    }
  } else if (p1[0] >= p2[0]) {
    right = p1;
    if (p0[0] >= p2[0]) {
      middle = p0;
      left = p2;
    } else {
      middle = p2;
      left = p0;
    }
  } else {
    right = p2;
    if (p0[0] >= p1[0]) {
      middle = p0;
      left = p1;
    } else {
      middle = p1;
      left = p0;
    }
  }
  
  int x; //iterator
  struct matrix *edges = new_matrix(4, 50);
  
  if (right[0] == left[0])
    return 0;
  double y_red = left[1];
  double dydx_red = 1.0 * (right[1] - left[1]) / (right[0] - left[0]);
  double z_red = left[2];
  double dzdx_red = 1.0 * (right[2] - left[2]) / (right[0] - left[0]);

  if (middle[0] != left[0]) {
    double y_green = left[1];
    double dydx_green = 1.0 * (middle[1] - left[1]) / (middle[0] - left[0]);
    double z_green = left[2];
    double dzdx_green = 1.0 * (middle[2] - left[2]) / (middle[0] - left[0]);

    for (x = left[0]; x < floor(middle[0]); x++) {
      y_red += dydx_red;
      z_red += dzdx_red;
      
      y_green += dydx_green;
      z_green += dzdx_green;
      
      add_edge(edges, x, y_red, z_red, x, y_green, z_green);
    }
  }

  if (right[0] != middle[0]) {
    double y_blue = middle[1];
    double dydx_blue = 1.0 * (right[1] - middle[1]) / (right[0] - middle[0]);
    double z_blue = middle[2];
    double dzdx_blue = 1.0 * (right[2] - middle[2]) / (right[0] - middle[0]);

    for (x = middle[0]; x < floor(right[0]); x++) {
      y_red += dydx_red;
      z_red += dzdx_red;
      
      y_blue += dydx_blue;
      z_blue += dzdx_blue;
      
      add_edge(edges, x, y_red, z_red, x, y_blue, z_blue);
    }
  }

  draw_lines(edges, s, c, d);
  free(edges);
    
  return 1;
}
  
int min(int a, int b) {
  return a>b? b: a;
}


/*======== int draw_polygons() =======
  Inputs: struct matrix * points
          screen s, color c
  Returns: number of polygons drawn

  go thru points 3 each and draw corresponding polygon
  ==================*/
int draw_polygons(struct matrix *points, screen s, color c, depthmap d) {
  int r0 = min(COLOR_VAR, c.red);
  int r1 = min(COLOR_VAR, 256 - c.red);
  int g0 = min(COLOR_VAR, c.green);
  int g1 = min(COLOR_VAR, 256 - c.green);
  int b0 = min(COLOR_VAR, c.blue);
  int b1 = min(COLOR_VAR, 256 - c.blue);
  int j = rand() % 32;
  
  int i;
  int count = 0;
  color c2;
  for (i = 0; i < (points->lastcol + 1) / 3; i++) {
    c2 = get_color(c.red + ((263 * i + j) % (r0 + r1)) - r0,
		   c.green + ((149 * i + j) % (g0 + g1)) - g0,
		   c.blue + ((257 * i + j) % (b0 + b1)) - b0);
    count += draw_filled_triangle(points, 3 * i, s, c2, d);
  }

  return count;
}

		 

/*======== void add_circle() ==========
  Inputs:   struct matrix * points
            double cx
            double cy
            double r
            double step
  Returns:
>
  Adds the circle at (cx, cy) with radius r to points
  ====================*/
void add_circle(struct matrix * points,
                 double cx, double cy, double cz,
                 double r, double step ) {
  double steps = -1 * step;
  while ((steps += step) <= 1) {
    add_edge(points,
	     cx + r * cos(TAO * steps), cy + r * sin(TAO * steps), cz,
	     cx + r * cos(TAO * (steps + step)), cy + r * sin(TAO * (steps + step)), cz);
  }
}

/*
Adds the curve bounded by the 4 points passsed as parameters
of type specified in type (see matrix.h for curve type constants)
to the matrix points
====================*/
void add_curve(struct matrix *points,
                double x0, double y0,
                double x1, double y1,
                double x2, double y2,
                double x3, double y3,
                double step, int type ) {
  struct matrix *coeff_trans;
  if (type == BEZIER) {
    coeff_trans = make_bezier();
  } else if (type == HERMITE) {
    coeff_trans = make_hermite();
  } else {
    printf("Error invalid curve type %d, expected 0 (HERMITE) or 1 (BEZIER)\n", type);
    return;
  }

  struct matrix *coeff = new_matrix(4, 2);
  coeff->m[0][0] = x0;
  coeff->m[0][1] = y0;
  coeff->m[1][0] = x1;
  coeff->m[1][1] = y1;
  coeff->m[2][0] = x2;
  coeff->m[2][1] = y2;
  coeff->m[3][0] = x3;
  coeff->m[3][1] = y3;
  coeff->lastcol = 1;

  coeff = matrix_mult(coeff_trans, coeff);

  double steps = -1 * step;
  struct matrix *time = new_matrix(1, 4);
  struct matrix *prev_point = NULL;
  struct matrix *point = NULL;
  while ((steps += step) <= 1) {
    time->m[0][0] = 1;
    time->m[0][1] = steps;
    time->m[0][2] = steps * steps;
    time->m[0][3] = steps * steps * steps;
    time->lastcol = 3;
    point = matrix_mult(time, coeff);

    if (prev_point)
      add_edge(points, prev_point->m[0][0], prev_point->m[0][1], 0,
	       point->m[0][0], point->m[0][1], 0);

    prev_point = point;
  }
}
	     
void add_box(struct matrix *points, double x0, double y0, double z0, double x_depth, double y_depth, double z_depth) {
  
  double vals[3] = {x0, y0, z0};
  double depths[3] = {x_depth, -1 * y_depth, -1 * z_depth};
  int i, j, a, b, c;
  for (i = 0; i < 3; i++) {
    for (a = 0; a < 2; a++) {
      
      double vertices[4][3];
      for (b = 0; b < 2; b++) {
  	for (c = 0; c < 2; c++) {
	  j = a? b + 2 * c: c + 2 * b;
  	  vertices[j][i] = vals[i] + a * depths[i];
  	  vertices[j][(i + 1) % 3] = vals[(i + 1) % 3] + b * depths[(i + 1) % 3];
  	  vertices[j][(i + 2) % 3] = vals[(i + 2) % 3] + c * depths[(i + 2) % 3];
  	}
      }

      add_polygon(points,
  		  vertices[0][0], vertices[0][1], vertices[0][2],
  		  vertices[2][0], vertices[2][1], vertices[2][2],
  		  vertices[1][0], vertices[1][1], vertices[1][2]);
      add_polygon(points,
  		  vertices[3][0], vertices[3][1], vertices[3][2],
  		  vertices[1][0], vertices[1][1], vertices[1][2],
  		  vertices[2][0], vertices[2][1], vertices[2][2]);
    }
  }
}

//spherical to cylindrical coordinates. theta, phi in [0, 1]
double *sphere_cartesian(double x, double y, double z, double r, double phi, double theta) {
  double *xyz = malloc(3 * sizeof(double));
  
  xyz[0] = x + r * cos(TAO * theta) * sin(PI * phi);
  xyz[1] = y + r * sin(TAO * theta) * sin(PI * phi);
  xyz[2] = z + r * cos(PI * phi);

  return xyz;
}

//tracks coordinates along the surface of a torus, theta phi in [0, 1]
double *torus_cartesian(double x, double y, double z, double r, double R, double phi, double theta) {
  double *xyz = malloc(3 * sizeof(double));

  xyz[0] = x + cos(TAO * phi) * (R + r * cos(TAO * theta));
  xyz[1] = y + sin(TAO * phi) * (R + r * cos(TAO * theta));
  xyz[2] = z + r * sin(TAO * theta);

  return xyz;
}
 
void add_sphere(struct matrix *points, double cx, double cy, double cz, double r, double step) {

  double phi = -1 * step;
  double theta = -1 * step;

  while ((phi += step) <= 1) {
    theta = -1 * step;
    while ((theta += step) <= 1) {
      double *point1 = sphere_cartesian(cx, cy, cz, r, phi, theta);
      double *point2 = sphere_cartesian(cx, cy, cz, r, phi, theta + step);
      double *point3 = sphere_cartesian(cx, cy, cz, r, phi + step, theta);
      double *point4 = sphere_cartesian(cx, cy, cz, r, phi + step, theta + step);
      add_polygon(points,
		  point1[0], point1[1], point1[2],
		  point2[0], point2[1], point2[2],
		  point3[0], point3[1], point3[2]);
      add_polygon(points,
		  point3[0], point3[1], point3[2],
		  point2[0], point2[1], point2[2],
		  point4[0], point4[1], point4[2]);
    }
  }
}

void add_torus(struct matrix *points, double cx, double cy, double cz, double r, double R, double step) {

  double phi = -1 * step;
  double theta = -1 * step;

  while ((phi += step) <= 1) {
    theta = -1 * step;
    while ((theta += step) <= 1) {
      double *point1 = torus_cartesian(cx, cy, cz, r, R, phi, theta);
      double *point2 = torus_cartesian(cx, cy, cz, r, R, phi, theta + step);
      double *point3 = torus_cartesian(cx, cy, cz, r, R, phi + step, theta);
      double *point4 = torus_cartesian(cx, cy, cz, r, R, phi + step, theta + step);
      add_polygon(points,
		  point1[0], point1[1], point1[2],
		  point2[0], point2[1], point2[2],
		  point3[0], point3[1], point3[2]);
      add_polygon(points,
		  point3[0], point3[1], point3[2],
		  point2[0], point2[1], point2[2],
		  point4[0], point4[1], point4[2]);
    }
  }
}


