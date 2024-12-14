#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

enum
{
    MAX_TESTS_AMOUNT = 10, // the max amount of tests, however we need only 2 in the task
    N                = 100, // the amount of steps by x and y
    MAX_ITERS        = 66666, // the max amout of iterations
};

#define PLOT 1 // "1" to enable gnuplot; "0" to disable gnuplot
#define EPS 1e-6
#define RESULT_CSV "result.csv"
#define OMEGA_CSV "omega.csv"
#define RESULT_GNU "result.gnu"
#define OMEGA_GNU "omega.gnu"

double h = 1.0 / N; // delta-step by x and y
double *glob_tests = NULL;
pid_t glob_gnuplot_pid = 0;

// Cleanup resources related to iteration
void cleanup_iter(void) {
    if (glob_gnuplot_pid > 0) {
        kill(glob_gnuplot_pid, SIGTERM);
        waitpid(glob_gnuplot_pid, NULL, 0);
        glob_gnuplot_pid = 0;
    }
    //remove(RESULT_CSV);
    //remove(OMEGA_CSV);
}

// Cleanup all resources
void cleanup_all(void) {
    if (glob_tests) {
        free(glob_tests);
        glob_tests = NULL;
    }
    cleanup_iter();
}

// Signal handler for SIGINT
void sigint_handler(int sig) {
    cleanup_all();
    printf("\nProgram terminated by user\n");
    exit(EXIT_FAILURE);
}

// Input test cases
void input_tests(unsigned int *tests_c) {
    printf("Please, input amount of tests (max %d): ", MAX_TESTS_AMOUNT);
    if ((scanf("%u", tests_c) != 1) || *tests_c > MAX_TESTS_AMOUNT) {
        printf("Err: invalid number of tests\n");
        cleanup_all();
        exit(EXIT_FAILURE);
    }

    glob_tests = malloc((*tests_c) * 3 * sizeof(double));
    if (!glob_tests) {
        printf("Err: can't allocate memory\n");
        cleanup_all();
        exit(EXIT_FAILURE);
    }

    printf("Leave 3rd argument as 0 (for Jacobi) or as w (w = 1 for Gauss-Seidel; 1 < w < 2 for SOR; w = -1 for optimal w)\n");
    for (unsigned int i = 0; i < *tests_c; ++i) {
        printf("Input kx and ky and w; for test #%u: ", i + 1);
        if (scanf("%lf%lf%lf", &glob_tests[3 * i], &glob_tests[3 * i + 1], &glob_tests[3 * i + 2]) != 3) {
            printf("Err: invalid input for test #%u\n", i + 1);
            cleanup_all();
            exit(EXIT_FAILURE);
        }
    }
}

// Analytical solution
double solve_analyt(double x, double y, double ky) {
    return (sinh(M_PI * y / sqrt(ky)) / sinh(M_PI) / sqrt(ky)) * sin(M_PI * x);
}

// Solve SLAE using Jacobi method
int solve_slae_via_jacobi(double u[N + 1][N + 1], double kx, double ky) {
    double u_new[N + 1][N + 1] = {{0.0}};
    for (int i = 0; i <= N; ++i) {
        double x = i * h;
        u_new[i][N] = u[i][N] = sin(M_PI * x);
    }

    int iter = 0;
    double max_dif;
    do {
        max_dif = 0.0;
        for (int i = 1; i < N; ++i) {
            for (int j = 1; j < N; ++j) {
                u_new[i][j] = (kx*(u[i+1][j] + u[i-1][j]) + ky*(u[i][j+1] + u[i][j-1])) / (2*(kx+ky));
                double dif = fabs(u[i][j] - u_new[i][j]);
                if (max_dif < dif) {
                    max_dif = dif;
                }
            }
        }

        for (int i = 1; i < N; ++i) {
            for (int j = 1; j < N; ++j) {
                u[i][j] = u_new[i][j];
            }
        }

    } while ((++iter < MAX_ITERS) && (max_dif > EPS));

    if (iter == MAX_ITERS) {
        return -1;
    }

    return iter;
}

// Solve SLAE using method with w
int solve_slae_via_w(double u[N + 1][N + 1], double kx, double ky, double w) {
    for (int i = 0; i <= N; ++i) {
        double x = i * h;
        u[i][N] = sin(M_PI * x);
    }

    int iter = 0;
    double max_dif;
    do {
        max_dif = 0.0;
        for (int i = 1; i < N; ++i) {
            for (int j = 1; j < N; ++j) {
                double old_val = u[i][j];
                double tmp = (kx*(u[i+1][j] + u[i-1][j]) + ky*(u[i][j+1] + u[i][j-1])) / (2*(kx+ky));
                u[i][j] = (1 - w)*old_val + w*tmp;

                double dif = fabs(u[i][j] - old_val);
                if (max_dif < dif) {
                    max_dif = dif;
                }
            }
        }

    } while ((++iter < MAX_ITERS) && (max_dif > EPS));

    if (iter == MAX_ITERS) {
        return -1;
    }

    return iter;
}

// Output results and plot using gnuplot
void output_result(double u[N + 1][N + 1], int iters, int cur_test) {
    FILE *csv_fd = fopen(RESULT_CSV, "w");
    if (!csv_fd) {
        perror("Err: can't open CSV file\n");
        cleanup_all();
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            double x = i * h;
            double y = j * h;
            if ((i % 25 == 0) && (j % 25 == 0)) {
                printf("u(%8.6lf, %8.6lf) | %8.6lf | %8.6lf\n", x, y, u[i][j], solve_analyt(x, y, glob_tests[3 * cur_test + 1]));
            }
            fprintf(csv_fd, "%8.6lf %8.6lf %8.6lf\n", x, y, u[i][j]);
        }
    }

    fclose(csv_fd);

    printf(
        "--------\n"
        "Test #%d:\n"
        "Iterations = %d\n\n"
        "Type enter to continue...\n",
        cur_test + 1, iters);

    if (PLOT) {
        glob_gnuplot_pid = fork();
        if (glob_gnuplot_pid < 0) {
            perror("Err: can't fork process for gnuplot\n");
            cleanup_all();
            exit(EXIT_FAILURE);
        } else if (glob_gnuplot_pid == 0) {
            execlp("gnuplot", "gnuplot", RESULT_GNU, NULL);
            perror("Err: can't execute gnuplot\n");
            exit(EXIT_FAILURE);
        }
        wait(NULL);
    }
}

int main(void) {
    signal(SIGINT, sigint_handler);

    unsigned int tests_c;
    input_tests(&tests_c);

    for (unsigned int cur_test = 0; cur_test < tests_c; ++cur_test) {
        double kx = glob_tests[3 * cur_test];
        double ky = glob_tests[3 * cur_test + 1];
        double w = glob_tests[3 * cur_test + 2];

        double u[N + 1][N + 1] = {{0.0}};
        int iters = -1;
        if (w >= 0){
            if (w == 0.0){
                iters = solve_slae_via_jacobi(u, kx, ky);
            }
            else if ((1 <= w) && (w < 2)){
                iters = solve_slae_via_w(u, kx, ky, w);
            }

            if (iters == -1) {
                printf("Test #%d doesn't converge. Try increasing MAX_ITERS\n", cur_test + 1);
                continue;
            }

            output_result(u, iters, cur_test);
        }
        else if (w == -1.0){
            FILE *csv_fd = fopen(OMEGA_CSV, "w");
            if (!csv_fd) {
                perror("Err: can't open CSV file\n");
                cleanup_all();
                exit(EXIT_FAILURE);
            }

            double omega = 1.0;
            int optimal_iters = MAX_ITERS;
            double optimal_omega = -1.0;

            printf("This may take time...\n");
            while (omega < 2.0){
                iters = solve_slae_via_w(u, kx, ky, omega);
                fprintf(csv_fd, "%.2lf %d\n", omega, iters);
                if ((optimal_iters > iters) && (iters != -1)){
                    optimal_iters = iters;
                    optimal_omega = omega;
                }
                omega += 0.01;
                memset(u, 0, sizeof u);
            }

            fclose(csv_fd);

            printf("Optimal omega: %lf\nIters: %d\n", optimal_omega, optimal_iters);

            if (PLOT) {
                glob_gnuplot_pid = fork();
                if (glob_gnuplot_pid < 0) {
                    perror("Err: can't fork process for gnuplot\n");
                    cleanup_all();
                    exit(EXIT_FAILURE);
                } else if (glob_gnuplot_pid == 0) {
                    execlp("gnuplot", "gnuplot", OMEGA_GNU, NULL);
                    perror("Err: can't execute gnuplot\n");
                    exit(EXIT_FAILURE);
                }
                wait(NULL);
            }
        }

        if (cur_test != tests_c - 1){
            while (getchar() != '\n');
        }
    }

    cleanup_all();
    return 0;
}
