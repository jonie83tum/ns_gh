#include "boundary_val.h"

void boundaryvalues(int imax, int jmax, double **U, double **V) {
	int k;
	for (k = 1; k <= jmax; k++) {
		U[0][k] = 0;
		U[imax][k] = 0;
		V[0][k] = -V[1][k];
		V[imax + 1][k] = -V[imax][k];
	}
	for (k = 1; k <= imax; k++) {
		V[k][0] = 0;
		V[k][jmax] = 0;
		U[k][0] = -U[k][1];
		U[k][jmax + 1] = -U[k][jmax];
	}
}
