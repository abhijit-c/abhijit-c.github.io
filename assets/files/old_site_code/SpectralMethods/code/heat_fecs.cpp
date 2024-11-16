#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

  double *u      = (double *)malloc(sizeof(double)*N);
  double *u_temp = (double *)malloc(sizeof(double)*N);
  for (int k = 0; k < N; k++) { u[k] = I(k*2*M_PI/N); }

  struct timespec start, finish;
  clock_gettime(CLOCK_MONOTONIC, &start);

  double dx = 2*M_PI/N, dt = CFL*dx*dx/alpha, t = 0.0;
  while (t + dt/2 <= T)
  {
    #pragma omp parallel for
    for (int k = 1; k < N-1; k++) 
    {
      u_temp[k] = u[k] + CFL*(u[k+1] -2*u[k] + u[k-1]) + dt*f(t, k*2*M_PI/N);
    }
    u_temp[0] = u_temp[N-1] = 0;
    #pragma omp parallel for
    for (int k = 0; k < N; k++) { u[k] = u_temp[k]; }
    t += dt;
  }
  clock_gettime(CLOCK_MONOTONIC, &finish);
  double t_elap = (finish.tv_sec - start.tv_sec);
  t_elap += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

  double err = 0.0;
  for (int k = 0; k < N; k++)
  {
    err = fmax( err, fabs( u[k] - u_exact(T, k*2*M_PI/N) ) );
  }

  printf("Solution computed in %fs with absolute error %f\n", t_elap, err);
  free(u); free(u_temp);
  return 0;
}
