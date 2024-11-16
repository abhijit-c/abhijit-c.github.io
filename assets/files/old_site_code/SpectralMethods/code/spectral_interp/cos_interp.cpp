/*
 * FILE: dff_cos.cpp
 * AUTHOR: Abhijit Chowdhary 2019/06/07.
 * ----------------------------------
 * An example showing how we can use the Fourier transform to numerically
 * approximate the derivative, with spectral accuracy.
 */
#include <stdlib.h>
#include <math.h>

#include <fftw3.h>

#include "interpft.h"

// Numerically differentiate cos(x) on [0, 2*pi] via spectral differentiation.
int main(int argc, char **argv)
{
  if (argc != 3)
  {
    printf("USAGE: ./spectral_differentiation.out num_grid num_interp\n\
            SAMPLE: ./spectral_differentiation.out 16 4096");
    exit(1);
  }
  int N = atoi(argv[1]), M = atoi(argv[2]);

  printf("# Interpolated cosine on [0, 2pi] to %d grid points\n", M);
  printf("# num_pts\t||phi(x)-cos(x)||_\\infty\n");

  fftw_complex *fxx, *out;
  out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * M);
  while (2*N < M)
  {
    fxx = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);

    for (int k = 0; k < N; k++)
    { //Build [0,2*pi] and map cos() onto it.
      fxx[k][0] = cos(k*2*M_PI / N); fxx[k][1] = 0;
    }
    for (int k = 0; k < M; k++) { out[k][0] = 0; out[k][1] = 0; }
    interpft(fxx, N, out, M);

    double err = 0.0;
    for (int k = 0; k < M; k++)
    {
      err = fmax(err, fabs( out[k][0] - cos(k*2*M_PI/M) ));
    }
    printf("%d\t%e\n", N, err);
    N = N + 1;

    fftw_free(fxx); 
  }

  fftw_free(out);
}
