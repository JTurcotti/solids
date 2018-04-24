#ifndef DISPLAY_H
#define DISPLAY_H

void plot(screen s, color c, depthmap d, int x, int y, double z);
void clear_screen(screen s, color c);
void clear_depthmap(depthmap d);
void save_ppm(screen s, char *file);
void save_extension(screen s, char *file);
void display(screen s);
color get_color(int r, int g, int b);

#endif
