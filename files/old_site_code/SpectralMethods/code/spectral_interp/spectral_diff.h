/*
 * FILE: spectral_diff.h
 * AUTHOR: Abhijit Chowdhary 2019/06/07.
 * -------------------------------------
 * Provides a function to differentiate a function applied to the grid [0,2*pi].
 */

#ifndef INTERPFT_H
#define INTERPFT_H

#include <fftw3.h>

/* 
 * FUNCTION: spectral_diff
 * -----------------------
 * Differentiates the function f onto a finer grid via spectral interpolation.
 *
 * INPUT:
 * fftw_complex &fxx: An array of complex numbers describing the function which
 *                    the interpolation is to be taken on applied to the grid
 *                    [0,2*pi].
 * int N: Spacing of the above grid, i.e. fxx[k] = f(2*pi*k/N).
 * int M: Spacing of output grid, i.e. out[k] = f'(2*pi*k/M). M >= N must be
 *        true, if not M is set to N and out is realloced to size N. If this
 *        fails, function returns having done nothing.
 * OUTPUT:
 * fftw_complex &out: Output array of size M where out[k] = f'(2*pi*k/M).
 */
inline void 
spectral_diff(fftw_complex *fxx, int N, fftw_complex *out, int M)
{
  if (M < N)
  {
    fftw_complex *new_out = (fftw_complex *)realloc(out, N);
    if (new_out == NULL) { return; }
    else { out = new_out; M = N; }
  }
  fftw_plan p = fftw_plan_dft_1d(N, fxx, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p); fftw_destroy_plan(p);

  int nyqst = N/2;
  for (int k = 0; k < N; ++k)
  { // Build [y(1:nyqst), zeros(M-N), y(nyqst+1:N)]
    // Multiply by M/N for oversampling, then by 1/M for normalized ifft.
    int w = (k < nyqst) ? k : k-N;
    double temp = out[k][0]; out[k][0] = -w*out[k][1]/N; out[k][1] = w*temp/N;
    if (k >= nyqst)
    {
      out[k+(M-N)][0] = out[k][0]; out[k+(M-N)][1] = out[k][1]; 
      out[k][0] = 0; out[k][1] = 0;
    }
  }
  if (N % 2 == 0)
  { // If even divide energy to symmetrize frequencies.
    out[M-nyqst][0] /= 2; out[M-nyqst][1] /= 2;
    out[nyqst][0] = out[M-nyqst][0]; out[nyqst][1] = out[M-nyqst][1];
  }

  fftw_plan p = fftw_plan_dft_1d(M, out, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p); fftw_destroy_plan(p);
}

#endif
