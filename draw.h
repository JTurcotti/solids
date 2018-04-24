#ifndef DRAW_H
#define DRAW_H

#include "matrix.h"
#include "ml6.h"
#include <time.h>
#include <stdlib.h>

void add_point(struct matrix *points,
	       double x, double y, double z);
void add_edge(struct matrix *points, 
	       double x0, double y0, double z0, 
	       double x1, double y1, double z1);
void add_polygon(struct matrix *points,
		double x0, double y0, double z0,
		double x1, double y1, double z1,
		double x2, double y2, double z2);

		
void draw_lines(struct matrix *points, screen s, color c, depthmap d);
void draw_line(struct matrix *points, int pos, screen s, color c, depthmap d);


int draw_polygons(struct matrix *points, screen s, color c, depthmap d);
int draw_polygon(struct matrix *points, int pos, screen s, color c, depthmap d);
int draw_filled_triangle(struct matrix *points, int pos, screen s, color c, depthmap d);


void add_circle(struct matrix *points,
		double cx, double cy, double cz,
		double r, double step);

void add_curve(struct matrix *points,
	       double x0, double y0,
	       double x1, double y1,
	       double x2, double y2,
	       double x3, double y3,
	       double step, int type );

void add_box(struct matrix *points,
	     double x, double y, double z,
	     double x_depth, double y_depth, double z_depth);

void add_sphere(struct matrix *points,
		double cx, double cy, double cz,
		double r, double step);

void add_torus(struct matrix *points,
		double cx, double cy, double cz,
	       double slicer, double bigr, double step);


#endif
