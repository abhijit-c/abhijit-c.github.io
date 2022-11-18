#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <fftw3.h>
#include <omp.h>

double alpha;

double u_exact(double t, double x) { return t*sin(x); }

// Manufactured so that u(t,x) = u_exact(t,x)
double f(double t, double x) { return (alpha*t + 1)*sin(x); }

double I(double x) { return u_exact(0,x); }

int main(int argc, char **argv)
{
  if (argc != 5)
  {
    printf("USAGE: ./binary alpha T number_of_spacialpoints CFL\n");
    exit(1);
  }
  /*
   * Solve u_t = (alpha)u_xx + f(t,x) for x in [0, 2pi], t in [0, T], with
   * initial condition I(x) and dirichet zero boundary conditions.
   */
  alpha = atof(argv[1]);
  double T = atof(argv[2]);
  int N = atoi(argv[3]);
  double CFL = atof(argv[4]);

  //fftw_init_threads();
  //fftw_plan_with_nthreads(omp_get_max_threads());

  fftw_complex *u      = (fftw_complex *)malloc(sizeof(fftw_complex)*N);
  fftw_complex *u_diff = (fftw_complex *)malloc(sizeof(fftw_complex)*N);

  fftw_plan fw, bk;
  fw = fftw_plan_dft_1d(N, u, u_diff, FFTW_FORWARD, FFTW_ESTIMATE);
  bk = fftw_plan_dft_1d(N, u_diff, u_diff, FFTW_BACKWARD, FFTW_ESTIMATE);

  for (int k = 0; k < N; k++) { u[k][0] = I(k*2*M_PI/N); u[k][1] = 0.0; }

  struct timespec start, finish;
  clock_gettime(CLOCK_MONOTONIC, &start);

  double dx = 2*M_PI/N, dt = CFL*dx*dx/alpha, t = 0.0;
  while (t + dt/2 <= T)
  {
    fftw_execute(fw);
    //#pragma omp parallel for
    for (int k = 0; k < N; k++)
    {
      int w = (k < N/2) ? k : N-k;
      u_diff[k][0] *= -(double)w/N; u_diff[k][1] *= -(double)w/N;
    }
    fftw_execute(bk);
    //#pragma omp parallel for
    for (int k = 1; k < N-1; k++) 
    {
      u[k][0] = u[k][0] + dt*( alpha*u_diff[k][0] + f(t, k*2*M_PI/N) );
      u[k][1] = u[k][1] + dt*( alpha*u_diff[k][1] + f(t, k*2*M_PI/N) );
    }
    u[0][0] = u[N-1][0] = 0; u[0][1] = u[N-1][1] = 0;
    t += dt;
  }
  clock_gettime(CLOCK_MONOTONIC, &finish);
  double t_elap = (finish.tv_sec - start.tv_sec);
  t_elap += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

  double err = 0.0;
  for (int k = 0; k < N; k++)
  {
    err = fmax( err, fabs( u[k][0] - u_exact(T, k*2*M_PI/N) ) );
  }

  printf("Solution computed in %fs with absolute error %f.\n", t_elap, err);
  fftw_destroy_plan(fw); fftw_destroy_plan(bk); 
  fftw_free(u), fftw_free(u_diff);
  return 0;
}
