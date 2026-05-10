/*
 * Brusselator 2D — Simulación de Patrones de Turing
 * PEC Física Computacional — Ejercicio 1
 *
 * du/dt = Du·Lap(u) + a - u + u²v
 * dv/dt = Dv·Lap(v) + b - u²v
 *
 * Compile:  gcc main.c -o simulacion -lm
 * Run:      ./simulacion
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

/* Double-buffer scheme; static para evitar desbordamiento de pila con L=100 */
static double u    [L][L];
static double v    [L][L];
static double u_new[L][L];
static double v_new[L][L];


/* Transformada Box-Muller: devuelve una muestra de N(0, sigma) */
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


/* Condición inicial: estado estacionario (u0, v0) + ruido gaussiano pequeño */
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


/* Laplaciano discreto 5 puntos con condiciones de contorno periódicas (dx=1) */
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


/* A(t) = RMS espacial de (u - u0); mide el crecimiento del patrón */
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


/* Escribe u[i][j] en u_<frame>.dat como matriz L×L separada por espacios */
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
    printf("  [Frame %d] %s\n", frame, filename);
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

    /* Estado estacionario: u0 = a+b, v0 = b/(a+b)^2 */
    u0 = A_PARAM + B_PARAM;
    v0 = B_PARAM / (u0 * u0);

    printf("Brusselator 2D | red %dx%d\n", L, L);
    printf("  Du=%.1f  Dv=%.1f  a=%.2f  b=%.2f\n", DU, DV, A_PARAM, B_PARAM);
    printf("  dt=%.3f  T=%.0f  N=%d  sigma=%.3f\n", DT, T_MAX, N_STEPS, SIGMA);
    printf("  u0=%.4f  v0=%.4f\n\n", u0, v0);

    initialize_fields(u0, v0);

    FILE *fp_amp  = fopen("amplitud.dat", "w");
    FILE *fp_time = fopen("time.dat", "w");
    if (!fp_amp || !fp_time) {
        fprintf(stderr, "ERROR: no se pueden abrir los archivos de salida.\n");
        return 1;
    }

    frame = 0;
    printf("Ejecutando %d pasos...\n", N_STEPS);

    for (step = 0; step <= N_STEPS; step++) {

        t = (double)step * DT;

        amplitude = compute_amplitude(u0);
        fprintf(fp_time, "%.6f\n", t);
        fprintf(fp_amp,  "%.10f\n", amplitude);

        if (step % FRAME_INTERVAL == 0 && frame < N_FRAMES) {
            save_field(frame);
            frame++;
        }

        if (step % 2000 == 0)
            printf("  paso %5d  t=%6.1f  A=%.6f\n", step, t, amplitude);

        if (step == N_STEPS) break;

        /* Actualización por Euler explícito */
        for (i = 0; i < L; i++) {
            for (j = 0; j < L; j++) {

                double uij = u[i][j];
                double vij = v[i][j];
                double u2v = uij * uij * vij;   /* término compartido u²v */

                double lap_u = compute_laplacian(u, i, j);
                double lap_v = compute_laplacian(v, i, j);

                u_new[i][j] = uij + DT * (DU * lap_u + A_PARAM - uij + u2v);
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

    printf("\nListo. Archivos generados: u_0.dat...u_9.dat, amplitud.dat, time.dat\n");
    return 0;
}
