/*
 * =============================================================
 * Brusselator Reaction-Diffusion System
 * Turing Pattern Formation on a 2D Periodic Grid
 * =============================================================
 * Computational Physics PEC - Exercise 1
 *
 * Model:
 *   du/dt = Du * Lap(u) + a - u + u^2 * v
 *   dv/dt = Dv * Lap(v) + b - u^2 * v
 *
 * Method: Explicit Euler, periodic boundary conditions,
 *         Gaussian initial noise (Box-Muller transform).
 *
 * Compile:  gcc main.c -o simulacion -lm
 * Run:      ./simulacion
 * =============================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Portable definition of PI (not guaranteed by all C standards) */
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

/* ==================== SIMULATION PARAMETERS ==================== */

#define L               100     /* Grid size: L x L lattice points      */
#define DU              1.0     /* Diffusion coefficient for u (activator)*/
#define DV              10.0    /* Diffusion coefficient for v (inhibitor) */
#define A_PARAM         0.1     /* Brusselator feed parameter a           */
#define B_PARAM         0.9     /* Brusselator drain parameter b          */
#define DT              0.01    /* Time step (Euler)                      */
#define T_MAX           200.0   /* Total simulation time                  */
#define N_STEPS         20000   /* Number of steps: T_MAX / DT            */
#define SIGMA           0.01    /* Std. deviation of initial Gaussian noise*/
#define N_FRAMES        10      /* Total field snapshots to save          */
#define FRAME_INTERVAL  2000    /* Save one frame every FRAME_INTERVAL steps*/

/* ================================================================ */

/*
 * Global field arrays.
 * We use a double-buffer scheme: u/v hold the current state,
 * u_new/v_new hold the updated state before we swap them.
 *
 * Declared static to avoid stack overflow (100*100*4 doubles ~ 3.2 MB).
 */
static double u    [L][L];
static double v    [L][L];
static double u_new[L][L];
static double v_new[L][L];


/* ==============================================================
 * gaussian_noise(sigma)
 * --------------------------------------------------------------
 * Returns a random sample from the normal distribution N(0, sigma)
 * using the Box-Muller transform:
 *
 *   Given U1, U2 ~ Uniform(0,1):
 *   Z = sqrt(-2 * ln(U1)) * cos(2*pi*U2)   ->  N(0, 1)
 *
 * Multiplying by sigma gives a sample from N(0, sigma).
 * ============================================================== */
double gaussian_noise(double sigma)
{
    double u1, u2, z;

    /* Draw U1 > 0 to avoid log(0) */
    do {
        u1 = (double)rand() / ((double)RAND_MAX + 1.0);
    } while (u1 <= 0.0);

    u2 = (double)rand() / ((double)RAND_MAX + 1.0);

    z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return sigma * z;
}


/* ==============================================================
 * initialize_fields(u0, v0)
 * --------------------------------------------------------------
 * Sets the initial condition as the homogeneous steady state
 * (u0, v0) plus a small Gaussian perturbation on each grid point:
 *
 *   u[i][j] = u0 + noise(0, sigma)
 *   v[i][j] = v0 + noise(0, sigma)
 *
 * The noise seeds the spatial Fourier modes that the Turing
 * instability will selectively amplify.
 * ============================================================== */
void initialize_fields(double u0, double v0)
{
    int i, j;
    for (i = 0; i < L; i++) {
        for (j = 0; j < L; j++) {
            u[i][j] = u0 + gaussian_noise(SIGMA);
            v[i][j] = v0 + gaussian_noise(SIGMA);
        }
    }
}


/* ==============================================================
 * compute_laplacian(field, i, j)
 * --------------------------------------------------------------
 * Computes the discrete 2D Laplacian at grid point (i, j)
 * using a 5-point stencil with dx = 1:
 *
 *   Lap(f)[i][j] = f[i+1][j] + f[i-1][j]
 *                + f[i][j+1] + f[i][j-1]
 *                - 4 * f[i][j]
 *
 * Periodic (toroidal) boundary conditions are applied via
 * modular indexing so the grid wraps around on all sides.
 * ============================================================== */
double compute_laplacian(double field[L][L], int i, int j)
{
    int ip = (i + 1) % L;          /* right  neighbour (wraps) */
    int im = (i - 1 + L) % L;      /* left   neighbour (wraps) */
    int jp = (j + 1) % L;          /* upper  neighbour (wraps) */
    int jm = (j - 1 + L) % L;      /* lower  neighbour (wraps) */

    return field[ip][j] + field[im][j]
         + field[i][jp] + field[i][jm]
         - 4.0 * field[i][j];
}


/* ==============================================================
 * compute_amplitude(u0)
 * --------------------------------------------------------------
 * Measures the spatial RMS deviation of u from its steady state:
 *
 *   A(t) = sqrt( (1/N) * SUM_{i,j} (u[i][j] - u0)^2 )
 *
 * where N = L * L.  A(t) starts near sigma, grows exponentially
 * in the Turing unstable regime, and saturates as nonlinear
 * effects saturate the pattern.
 * ============================================================== */
double compute_amplitude(double u0)
{
    int i, j;
    double sum = 0.0;
    double diff;

    for (i = 0; i < L; i++) {
        for (j = 0; j < L; j++) {
            diff  = u[i][j] - u0;
            sum  += diff * diff;
        }
    }
    return sqrt(sum / (double)(L * L));
}


/* ==============================================================
 * save_field(frame)
 * --------------------------------------------------------------
 * Writes the current field u[i][j] to "u_<frame>.dat".
 * Format: L rows, each with L space-separated floating-point
 * values.  Suitable for reading with NumPy or gnuplot.
 * ============================================================== */
void save_field(int frame)
{
    char filename[32];
    int i, j;

    sprintf(filename, "u_%d.dat", frame);
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        fprintf(stderr, "ERROR: cannot open %s for writing.\n", filename);
        return;
    }

    for (i = 0; i < L; i++) {
        for (j = 0; j < L; j++) {
            fprintf(fp, "%.8f", u[i][j]);
            if (j < L - 1) fprintf(fp, " ");
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    printf("    [Frame %2d saved] %s\n", frame, filename);
}


/* ==============================================================
 * MAIN
 * ============================================================== */
int main(void)
{
    int    i, j, step, frame;
    double t, amplitude;
    double u0, v0;

    /* ---- Seed the random number generator with the current time ---- */
    srand((unsigned int)time(NULL));

    /* ---- Compute homogeneous steady state --------------------------------
     * The Brusselator has a unique fixed point:
     *   u0 = a + b
     *   v0 = b / (a + b)^2
     * For a=0.1, b=0.9: u0 = 1.0, v0 = 0.9
     * --------------------------------------------------------------------- */
    u0 = A_PARAM + B_PARAM;
    v0 = B_PARAM / (u0 * u0);

    /* ---- Print simulation parameters ---- */
    printf("================================================\n");
    printf("  Brusselator Reaction-Diffusion Simulation\n");
    printf("  Turing Pattern Formation  |  Grid %d x %d\n", L, L);
    printf("================================================\n");
    printf("Parameters:\n");
    printf("  Du          = %.4f\n",  DU);
    printf("  Dv          = %.4f\n",  DV);
    printf("  a           = %.4f\n",  A_PARAM);
    printf("  b           = %.4f\n",  B_PARAM);
    printf("  dt          = %.4f\n",  DT);
    printf("  T_max       = %.1f\n",  T_MAX);
    printf("  N_steps     = %d\n",    N_STEPS);
    printf("  sigma       = %.4f\n",  SIGMA);
    printf("  Frame every = %d steps\n", FRAME_INTERVAL);
    printf("\nHomogeneous steady state:\n");
    printf("  u0 = %.6f\n", u0);
    printf("  v0 = %.6f\n", v0);
    printf("================================================\n\n");

    /* ---- Initialize fields ---- */
    printf("Initializing fields with Gaussian noise (sigma = %.4f)...\n\n",
           SIGMA);
    initialize_fields(u0, v0);

    /* ---- Open output files ---- */
    FILE *fp_amp  = fopen("amplitud.dat", "w");
    FILE *fp_time = fopen("time.dat", "w");
    if (!fp_amp || !fp_time) {
        fprintf(stderr, "ERROR: cannot open output files. "
                        "Check write permissions.\n");
        return 1;
    }

    /* ---- Time evolution ---- */
    frame = 0;
    printf("Starting time evolution (%d steps)...\n\n", N_STEPS);
    printf("  %-8s  %-10s  %-15s\n",
           "Step", "Time", "Amplitude A(t)");
    printf("  %-8s  %-10s  %-15s\n",
           "--------", "----------", "---------------");

    for (step = 0; step <= N_STEPS; step++) {

        t = (double)step * DT;

        /* -- Compute spatial amplitude of perturbation -- */
        amplitude = compute_amplitude(u0);

        /* -- Write amplitude and time to files -- */
        fprintf(fp_time, "%.6f\n", t);
        fprintf(fp_amp,  "%.10f\n", amplitude);

        /* -- Save field snapshot every FRAME_INTERVAL steps -- */
        if (step % FRAME_INTERVAL == 0 && frame < N_FRAMES) {
            save_field(frame);
            frame++;
        }

        /* -- Print progress every 2000 steps -- */
        if (step % 2000 == 0) {
            printf("  %-8d  %-10.2f  %-15.8f\n", step, t, amplitude);
        }

        /* -- On the last step, only record; do not update fields -- */
        if (step == N_STEPS) break;

        /* ============================================================
         * EXPLICIT EULER UPDATE
         * For each grid point compute the discrete reaction-diffusion
         * equations and advance u and v by one time step dt.
         * ============================================================ */
        for (i = 0; i < L; i++) {
            for (j = 0; j < L; j++) {

                double uij = u[i][j];
                double vij = v[i][j];
                double u2v = uij * uij * vij;   /* u^2 * v (shared term) */

                double lap_u = compute_laplacian(u, i, j);
                double lap_v = compute_laplacian(v, i, j);

                /*  du/dt = Du * Lap(u) + a - u + u^2*v  */
                u_new[i][j] = uij + DT * (DU * lap_u + A_PARAM - uij + u2v);

                /*  dv/dt = Dv * Lap(v) + b - u^2*v      */
                v_new[i][j] = vij + DT * (DV * lap_v + B_PARAM - u2v);
            }
        }

        /* -- Swap buffers: new state becomes current state -- */
        for (i = 0; i < L; i++) {
            for (j = 0; j < L; j++) {
                u[i][j] = u_new[i][j];
                v[i][j] = v_new[i][j];
            }
        }

    } /* end time loop */

    fclose(fp_amp);
    fclose(fp_time);

    /* ---- Final summary ---- */
    printf("\n================================================\n");
    printf("Simulation complete!\n");
    printf("Generated output files:\n");
    printf("  u_0.dat ... u_9.dat    field snapshots of u(x,y)\n");
    printf("  amplitud.dat           A(t) at every time step\n");
    printf("  time.dat               t   at every time step\n");
    printf("================================================\n");

    return 0;
}
