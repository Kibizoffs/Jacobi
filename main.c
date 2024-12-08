#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

enum
{
    MAX_TESTS_AMOUNT = 10,
    N                = 100,
    MAX_ITER         = 66666,
};

#define PLOT 1 // "1" to enable gnuplot; "0" to disable gnuplot
#define EPS 1e-6
#define CSV_FILE "result.csv"
#define PLOT_SCRIPT "plot.gnu"

double s = 1; // length
double h = 1.0 / N; // step by x and y

double *glob_tests = NULL;
pid_t glob_gnuplot_pid = 0;

// Cleanup resources related to iteration
void cleanup_iter(void) {
    if (glob_gnuplot_pid > 0) {
        kill(glob_gnuplot_pid, SIGTERM);
        waitpid(glob_gnuplot_pid, NULL, 0);
        glob_gnuplot_pid = 0;
    }
    if (remove(CSV_FILE) != 0) {
        perror("Error removing CSV file");
    }
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

    glob_tests = malloc((*tests_c) * 2 * sizeof(double));
    if (!glob_tests) {
        printf("Err: can't allocate memory");
        cleanup_all();
        exit(EXIT_FAILURE);
    }

    for (unsigned int i = 0; i < *tests_c; ++i) {
        printf("Input kx and ky for test #%u: ", i + 1);
        if (scanf("%lf%lf", &glob_tests[2 * i], &glob_tests[2 * i + 1]) != 2) {
            printf("Err: invalid input for test #%u.\n", i + 1);
            cleanup_all();
            exit(EXIT_FAILURE);
        }
    }
}

// Analytical solution
double solve_analyt(double x, double y) {
    return (sinh(M_PI * y) / sinh(M_PI)) * sin(M_PI * x);
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
                u_new[i][j] = (kx * (u[i + 1][j] + u[i - 1][j]) +
                               ky * (u[i][j + 1] + u[i][j - 1])) /
                              (2 * (kx + ky));
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

    } while ((++iter < MAX_ITER) && (max_dif > EPS));

    if (iter == MAX_ITER) {
        printf("Doesn't converge\n");
        return -1;
    }

    return iter;
}

// Output results and plot using gnuplot
void output_result(double u[N + 1][N + 1], int iters, int cur_test) {
    FILE *csv_fd = fopen(CSV_FILE, "w");
    if (!csv_fd) {
        perror("Err: can't open CSV file");
        cleanup_all();
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            double x = i * h;
            double y = j * h;
            if ((i % 10 == 0) && (j % 10 == 0)) {
                printf("u(%8.6lf, %8.6lf) | %8.6lf | %8.6lf\n", x, y, u[i][j], solve_analyt(x, y));
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
            perror("Err: can't fork process for gnuplot");
            cleanup_all();
            exit(EXIT_FAILURE);
        } else if (glob_gnuplot_pid == 0) {
            execlp("gnuplot", "gnuplot", PLOT_SCRIPT, NULL);
            perror("Err: can't execute gnuplot");
            exit(EXIT_FAILURE);
        }
        wait(NULL);
    }
}

// Main function
int main(void) {
    signal(SIGINT, sigint_handler);

    unsigned int tests_c;
    input_tests(&tests_c);

    for (unsigned int cur_test = 0; cur_test < tests_c; ++cur_test) {
        double kx = glob_tests[2 * cur_test];
        double ky = glob_tests[2 * cur_test + 1];

        double u[N + 1][N + 1] = {{0.0}};
        int iters = solve_slae_via_jacobi(u, kx, ky);

        if (iters == -1) {
            printf("Test #%d doesn't converge\n", cur_test + 1);
            continue;
        }

        output_result(u, iters, cur_test);

        while (getchar() != '\n');
    }

    cleanup_all();
    return 0;
}
