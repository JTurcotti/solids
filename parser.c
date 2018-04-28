#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "stack.h"
#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "matrix.h"
#include "parser.h"
#inclyde "preparser.h"

/*======== void parse_file () ==========
Inputs:   char * filename 

Returns: void

Goes through the file named filename and performs all of the actions listed in that file.
The file follows the following format:
     Every command is a single character that takes up a line
     Any command that requires arguments must have those arguments in the second line.
     The commands are as follows:
         line: draw a line to the screen
	    takes 6 arguemnts (x_0, y_0, z_0) 
	                      (x_1, y_1, z_1)
	 polygon: draw a triangle to the screen
	    takes 9 arguments (x_0, y_0, z_0) 
	                      (x_1, y_1, z_1) 
			      (x_2, y_2, z_2) 
	 circle: draw a circle to the screen
	    takes 4 arguments (x_center, y_center, z_center)
	                       radius
	 bezier: draw a bezier curve to the screen from endpoints and rates of change -
	    takes 8 arguments (x_0, y_0) (x_1, y_1)
	                      (dx/dt_0, dy/dt_0) (dx/dt_1, dy/dt_1)
	 hermite: draw hermite curve to screen from endpoints and control points -
	    takes 8 arguments (x_0, y_0)
                              (ctrlx_0, ctrly_0) (ctrlx_1, ctrly_1)
			      (x_1, y_1)
	 box: draw box to screen from minimal coordinate vertex and widths -
	    takes 6 arguments (x_min, y_min, z_min)
	                      (x_width, y_width, z_width)
	 sphere: draw sphere from center and radius -
	    takes 4 arguments (x_center, y_center, z_center)
	                       radius
	 torus: draw torus from center, radius of tube (minor radius), and the distance from the center to the center of the tube (major radius) -
	    takes 5 arguments (x_center, y_center, z_center)
	                       radius_minor
			       radius_major
	 push: copy the top coordinate system and push the copy to the stack
	 pop: remove the top coordinate system from the stack
	 scale: create a scale matrix, 
	    then multiply the top coordinate system by it
	    takes 3 arguments (sx, sy, sz)
	 translate: create a translation matrix, 
	    then multiply the top coordinate system by it
	    takes 3 arguments (tx, ty, tz)
	 rotate: create an rotation matrix,
	    then multiply the top coordinate system by it
	    takes 2 arguments (axis, theta) axis should be x y or z
	 color: set the drawing color to desired color
	    takes 3 arguments (red, green, blue)
	 background: set the background to desired color
	    takes 3 arguments (red, green, blue)
	 clear: clears the screen
	 save: draws and saves the screen to a file -
	    takes 1 argument (file name)
	 quit: end parsing

====================*/
void parse_file ( char * input) {
  int nargs = 0;
  int ndisplay = 0;

  struct matrix *edges = new_matrix(1, 1);
  struct matrix *polygons = new_matrix(1, 1);

  struct stack *stack = new_stack();

  screen s;
  depthmap d;

  FILE *f;
  char line[256], argline[256];

  color c = get_color(255, 255, 255);
  color back = get_color(0, 0, 0);
  
  clear_screen(s, back);
  clear_depthmap(d);

  if ( strcmp(input, "stdin") == 0 ) 
    f = stdin;
  else
    f = fopen(input, "r");
  
  while (fgets(line, 255, f)) {
    if (strchr(line, '\r')) {
      fprintf(stderr, "Please use a script file with Unix-style (\\n) line endings\n");
      break;
    }
    
    line[strlen(line)-1]='\0'; //remove new line
    if (strlen(line) == 0 ||
	line[0] == '#' ||
	(line[0] == '/' && line[1] == '/')) {
      continue;
    }

    printf(":%s\n", line);
    if (!strcmp(line, "quit")) {
      break;
    } else if (!strcmp(line, "line")) {
      double *args = malloc(6 * sizeof(double));
      edges = new_matrix(4, 1);
      
      if (!fgets(argline, 255, f) ||
	  ((nargs = sscanf(argline, "%lf %lf %lf %lf %lf %lf",
			   args, args+1, args+2, args+3, args+4, args+5)) != 6)) {
	printf("Error: 'line' requires 6 arguments of type double, found %d\n", nargs);
      } else {
	add_edge(edges, args[0], args[1], args[2], args[3], args[4], args[5]);
      }
      draw_lines(matrix_mult(peek(stack), edges), s, c, d);
      free(args);
      free(edges);
    } else if (!strcmp(line, "polygon")) {
      double *args = malloc(9 * sizeof(double));
      polygons = new_matrix(4, 3);

      if (!fgets(argline, 255, f) ||
	  ((nargs = sscanf(argline, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
			   args, args+1, args+2, args+3, args+4, args+5,
			   args+6, args+7, args+8)) != 9)) {
	printf("Error: 'polygon' requires 9 arguments of type double, found %d\n", nargs);
      } else {
	add_polygon(polygons, args[0], args[1], args[2], args[3], args[4],
		    args[5], args[6], args[7], args[8]);
      }
      draw_polygons(matrix_mult(peek(stack), polygons), s, c, d);
      free(args);
      free(polygons);
    } else if (!strcmp(line, "box")) {
      double *args = malloc(6 * sizeof(double));
      polygons = new_matrix(4, 36);

      if (!fgets(argline, 255, f) ||
	  ((nargs = sscanf(argline, "%lf %lf %lf %lf %lf %lf", args, args+1, args+2, args+3, args+4, args+5)) != 6)) {
	printf("Error: 'box' requires 6 arguments of type double, found %d\n", nargs);
      } else {
	add_box(polygons, args[0], args[1], args[2], args[3], args[4], args[5]);
      }
      draw_polygons(matrix_mult(peek(stack), polygons), s, c, d);
      free(args);
    } else if (!strcmp(line, "sphere")) {
      double *args = malloc(4 * sizeof(double));
      polygons = new_matrix(4, 128);

      if (!fgets(argline, 255, f) ||
	  ((nargs = sscanf(argline, "%lf %lf %lf %lf", args, args+1, args+2, args+3)) != 4)) {
	printf("Error: 'sphere' requires 4 arguments of type double, found %d\n", nargs);
      } else {
	add_sphere(polygons, args[0], args[1], args[2], args[3], STEP_SIZE * 5);
      }
      draw_polygons(matrix_mult(peek(stack), polygons), s, c, d);
      free(args);
    } else if (!strcmp(line, "torus")) {
      double *args = malloc(5 * sizeof(double));
      polygons = new_matrix(4, 128);
      
      if (!fgets(argline, 255, f) ||
	  ((nargs = sscanf(argline, "%lf %lf %lf %lf %lf", args, args+1, args+2, args+3, args+4)) != 5)) {
	printf("Error: 'torus' requires 5 arguments of type double, found %d\n", nargs);
      } else {
	add_torus(polygons, args[0], args[1], args[2], args[3], args[4], STEP_SIZE * 3);
      }
      draw_polygons(matrix_mult(peek(stack), polygons), s, c, d);
      free(args);
    } else if (!strcmp(line, "circle")) {
      double *args = malloc(4 * sizeof(double));
      edges = new_matrix(4, 64);

      if (!fgets(argline, 255, f) ||
	  ((nargs = sscanf(argline, "%lf %lf %lf %lf", args, args+1, args+2, args+3)) != 4)) {
	printf("Error: 'circle' requires 4 arguments of type double, found %d\n", nargs);
      } else {
	add_circle(edges, args[0], args[1], args[2], args[3], STEP_SIZE);
      }
      draw_lines(matrix_mult(peek(stack), edges), s, c, d);
      free(args);
    } else if (!strcmp(line, "bezier")) {
      double *args = malloc(8 * sizeof(double));
      edges = new_matrix(4, 64);

      if (!fgets(argline, 255, f) ||
	  ((nargs = sscanf(argline, "%lf %lf %lf %lf %lf %lf %lf %lf",
			   args, args+1, args+2, args+3, args+4, args+5,
			   args+6, args+7)) != 8)) {
	printf("Error: 'bezier' requires 8 arguments of type double, found %d\n", nargs);
      } else {
	add_curve(edges, args[0], args[1], args[2], args[3], args[4],
		  args[5], args[6], args[7], STEP_SIZE, BEZIER);
      }
      draw_lines(matrix_mult(peek(stack), edges), s, c, d);
      free(args);
    } else if (!strcmp(line, "hermite")) {
      double *args = malloc(8 * sizeof(double));
      edges = new_matrix(4, 64);

      if (!fgets(argline, 255, f) ||
	  ((nargs = sscanf(argline, "%lf %lf %lf %lf %lf %lf %lf %lf",
			   args, args+1, args+2, args+3, args+4, args+5,
			   args+6, args+7)) != 8)) {
	printf("Error: 'hermite' requires 8 arguments of type double, found %d\n", nargs);
      } else {
	add_curve(edges, args[0], args[1], args[2], args[3], args[4],
		  args[5], args[6], args[7], STEP_SIZE, HERMITE);
      }      
      draw_lines(matrix_mult(peek(stack), edges), s, c, d);
      free(args);
    } else if (!strcmp(line, "push")) {
      push(stack);
    } else if (!strcmp(line, "pop")) {
      pop(stack);
    } else if (!strcmp(line, "scale")) {
      double *args = malloc(3 * sizeof(double));
      
      if (!fgets(argline, 255, f) ||
	  ((nargs = sscanf(argline, "%lf %lf %lf", args, args+1, args+2)) != 3)) {
	printf("Error: 'scale' requires 3 arguments of type double, found %d\n", nargs);
      } else {
	apply_to_stack(stack, make_scale(args[0], args[1], args[2]));
      }
      
      free(args);
      
    } else if (!strcmp(line, "translate") || !strcmp(line, "move")) {
      double *args = malloc(3 * sizeof(double));
      
      if (!fgets(argline, 255, f) ||
	  ((nargs = sscanf(argline, "%lf %lf %lf", args, args+1, args+2)) != 3)) {
	printf("Error: 'translate' requires 3 arguments of type double, found %d\n", nargs);
      } else {
	apply_to_stack(stack, make_translate(args[0], args[1], args[2]));
      }

      free(args);
    } else if (!strcmp(line, "rotate")) {
      char *axis = malloc(1);
      double *theta = malloc(sizeof(double));
      if (!fgets(argline, 255, f) || (sscanf(argline, "%c %lf", axis, theta) != 2)) {
	printf("Error: 'rotate' requires both an axis and an angle\n");
      } else {
	if (*axis == 'x' || *axis == 'X') {
	  apply_to_stack(stack, make_rotX(*theta));
	} else if (*axis == 'y' || *axis == 'Y') {
	  apply_to_stack(stack, make_rotY(*theta));
	} else if (*axis == 'z' || *axis == 'Z') {
	  apply_to_stack(stack, make_rotZ(*theta));
	} else {
	  printf("Error: %c is not a valid axis\n", *axis);
	}
      }
    } else if (!strcmp(line, "color")) {
      int rgb[3];
      if (!fgets(argline, 255, f) || (sscanf(argline, "%d %d %d", rgb, rgb+1, rgb+2) != 3)) {
	printf("Error: 'color' requires three integeres\n");
      } else {
	c = get_color(rgb[0], rgb[1], rgb[2]);
      }
    } else if (!strcmp(line, "background")) {
      int rgb[3];
      if (!fgets(argline, 255, f) || (sscanf(argline, "%d %d %d", rgb, rgb+1, rgb+2) != 3)) {
	printf("Error: 'background' requires three integeres\n");
      } else {
	back = get_color(rgb[0], rgb[1], rgb[2]);
      }
    } else if (!strcmp(line, "display")) {
      char *filename = malloc(16);
      sprintf(filename, "display-%d.png", ndisplay++);
      save_extension(s, filename);
      free(filename);
    } else if (!strcmp(line, "save")) {
      char *filename = malloc(32);
      if (!fgets(argline, 255, f) || (sscanf(argline, "%s", filename) == 0)) {
	printf("Error: 'save' requires a filename, none given\n");
      } else {
	printf("saving as %s\n", filename);
	save_extension(s, filename);
      }
      free(filename);
    } else if (!strcmp(line, "quit")) {
      printf("Parse terminated\n");
      exit(0);
    } else {
      printf("Error: unrecognized command '%s'\n", line);
    }
			  
  }
}
  
//*/
