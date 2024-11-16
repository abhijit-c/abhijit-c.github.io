#include <complex.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define BOX_SIZE 2.0
#define NUM_FRAMES 240

int main(int argc, char **argv)
{
  if (argc != 2)
  {
    perror("USAGE: ./julia.out resolution \n");
    return -1;
  }
  int resolution = atoi(argv[1]);
  int w = 4*resolution; double w_grid[w];
  for (int i = 0; i < w; i++) { w_grid[i] = -BOX_SIZE + 2*i*BOX_SIZE/w; }
  int h = 4*resolution; double h_grid[h];
  for (int i = 0; i < h; i++) { h_grid[i] = -BOX_SIZE + 2*i*BOX_SIZE/h; }

  #pragma omp parallel for
  for (int frame = 1; frame <= NUM_FRAMES; frame++)
  {
    char filename[50]; sprintf(filename, "%05d.pgm", frame);
    FILE* current_frame = fopen(filename, "w");
    fprintf(current_frame, "P2\n\n%d %d\n255\n", w, h);

    double theta = frame*(2*M_PI/NUM_FRAMES);
    for (int i = 0; i < w; i++)
    {
      for (int j = 0; j < h; j++)
      {
        complex double z = w_grid[i] + I*h_grid[j];
        complex double c = 0.7885*cexp(I*theta);
        int n = 255;
        while (cabs(z) < 10 && n > 0)
        {
          z = cpow(z,2) + c;
          n = n - 5;
        }
        fprintf(current_frame, "%d ", n);
      }
      fprintf(current_frame, "\n");
    }

    fclose(current_frame);

    char cmd_buffer[50];
    sprintf(cmd_buffer, "convert %05d.pgm %05d.png", frame, frame);
    system(cmd_buffer);

    printf("%s\n", filename);
  }

  system("convert -delay 3.33 -loop 0 *.png animation.gif");
  system("make clean");

  return 0;
}
