#include "uvp.h"
#include <stdlib.h>
#include <math.h>

void calculate_fg(double Re, double GX, double GY, double alpha, double dt,
		double dx, double dy, int imax, int jmax, double **U, double **V,
		double **F, double **G) {

}

void calculate_rs(double dt, double dx, double dy, int imax, int jmax,
		double **F, double **G, double **RS) {
	int i, j;
	for (int i = 1; i <= imax; i++) {
		for (int j = 1; j <= jmax; j++) {
			RS[i][j] = ((F[i][j] - F[i - 1][j]) / dx
					+ (G[i][j] - G[i][j - 1]) / dy) / dt;
		}
	}
}

void calculate_dt(double Re, double tau, double *dt, double dx, double dy,
		int imax, int jmax, double **U, double **V) {

	int i, j;
	double max_U, max_V, c1, c2, c3, minc;
	/* in case of negative tau the function is doing nothing because dt is already read from the input file */
	if (tau < 0) {
		return;
	}
	/* fine maximum of U and V */
	max_U = U[0][0];
	max_V = V[0][0];
	for (i = 0; i <= imax + 1; i++) {
		for (j = 0; j <= jmax + 1; j++) {
			if (abs(U[i][j]) > max_U) {
				max_U = abs(U[i][j]);
			}
			if (abs(V[i][j]) > max_V) {
				max_V = abs(V[i][j]);
			}
		}
	}

	c2 = dx / max_U;
	c3 = dy / max_V;
	c1 = Re / 2 / (pow(dx, -2) + pow(dy, -2));

	/* take the min out of the c values */
	minc = c1;
	if (c2 < minc) {
		minc = c2;
	}
	if (c3 < minc) {
		minc = c3;
	}
	*dt = tau * minc;

	return;
}

void calculate_uv(double dt, double dx, double dy, int imax, int jmax,
		double **U, double **V, double **F, double **G, double **P) {

	int i, j;

	for (i = 1; i <= imax - 1; i++) {
		for (j = 1; j <= jmax; j++) {
			U[i][j] = F[i][j] - dt / dx * (P[i + 1][j] - P[i][j]);
		}
	}
	for (i = 1; i <= imax; i++) {
		for (j = 1; j <= jmax - 1; j++) {
			V[i][j] = G[i][j] - dt / dy * (P[i][j + 1] - P[i][j]);
		}
	}
}
