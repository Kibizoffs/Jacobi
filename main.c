#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/wait.h>

enum
{
    N         = 100,
    MAX_ITER  = 66666,
};

#define EPS 1e-6

// D = DX = DY = 1
double h = 1.0 / N; // step

double
exact_result(double x, double y)
{
    double u = (sinh(M_PI*y)/sinh(M_PI))*sin(M_PI*x);
    return u;
}

int
solve_slae_via_jacobi(double u[N+1][N+1], double kx, double ky)
{
    if ((kx / ky < 0.1) || (kx / ky > 10)){
        printf(
            "Caution! There is a change of significant anisotropy\n" \
            "Would you like t"
        ");
    }
    
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

        for (int i = 1; i < N; ++i){
           for (int j = 1; j < N; ++j){
                u[i][j] = u_new[i][j];
            }
        }

    } while ((++iter < MAX_ITER) && (max_dif > EPS));

    if (iter == MAX_ITER){
        printf("Doesn't converge\n");
        return -1;
    }

    return iter;
}

void
output_result(double u[N+1][N+1], int iters)
{
    printf("--------\n");

    char csv_name[16] = "result.csv";
    FILE *csv_fd = fopen(csv_name, "w");

    for (int i = 0; i < N; ++i){
       for (int j = 0; j < N; ++j){
            double x = i*h; double y = j*h;
            
            if ((i % 10 == 0) && (j % 10 == 0)){
                printf("u(%8.6lf,%8.6lf) | %8.6lf | %8.6lf\n", x, y, u[i][j], exact_result(x, y));
            }

            if (csv_fd){
                fprintf(csv_fd, "%8.6lf %8.6lf %8.6lf\n", x, y, u[i][j]);
            }
        }
    }

    printf(
        "--------\n" \
        "Iters: %d\n" \
        "\n",
        iters
    );

    if (csv_fd) {
        fclose(csv_fd);
    }

    if (fork() == 0){
        execlp("gnuplot", "gnuplot", "plot.gnu", NULL);
        exit(1);
    }
}

int
main(void)
{
    printf(
        "JACOBI\n" \
        "Please, input kx and ky: "
    );

    double kx; double ky;
    if (scanf("%lf%lf", &kx, &ky) != 2){
        return 1;
    }

    double u[N+1][N+1] = {{0.0}};
    int iters = solve_slae_via_jacobi(u, kx, ky);

    output_result(u, iters);

    return 0;
}
