#include <stdio.h>
#include <stdlib.h>
#include <math.h>

enum
{
    N         = 100,
    MAX_ITER  = 10000,
};

#define EPS 1e-6

// D = DX = DY = 1
double h = 1.0 / N; // step

double
exact_res(double x, double y)
{
    double u = (sinh(M_PI*y)/sinh(M_PI))*sin(M_PI*x);
    return u;
}

int
solve_slae_via_jacobi(double kx, double ky)
{
    double u[N+1][N+1] = {{0.0}};
    double u_new[N+1][N+1] = {{0.0}};

    for (int i = 0; i <= N; ++i){
        double x = i*h;
        u_new[i][N] = u[i][N] = sin(M_PI*x);
    }

    int iter = 0;
    double max_dif;
    do {
        max_dif = 0.0;
        for (int i = 1; i < N; ++i){
           for (int j = 1; j < N; ++j){
                u_new[i][j] = (kx*(u[i+1][j] + u[i-1][j]) + ky*(u[i][j+1] + u[i][j-1])) / (2*(kx+ky));
                double dif = fabs(u[i][j] - u_new[i][j]);
                if (max_dif < dif){
                    max_dif = dif;
                }
            }
        }

        for (int i = 0; i <= N; ++i){
           for (int j = 0; j <= N; ++j){
                u[i][j] = u_new[i][j];
            }
        }
    } while ((++iter < MAX_ITER) && (max_dif > EPS));

    if (iter == MAX_ITER){
        printf("Doesn't converge\n");
        return 1;
    }

    printf("Iters: %d\n", iter);
    for (int i = 0; i < N; ++i){
       for (int j = 0; j < N; ++j){
            double x = i*h; double y = j*h;
            printf("u(%8.6lf,%8.6lf) | %8.6lf | %8.6lf\n", x, y, u_new[i][j], exact_res(x, y));
        }
    }

    return 0;
}

int
main(void)
{
    printf("JACOBI\n");

    double kx; double ky;
    if (scanf("%lf%lf", &kx, &ky) != 2){
        return 1;
    }

    solve_slae_via_jacobi(kx, ky);

    return 0;
}
