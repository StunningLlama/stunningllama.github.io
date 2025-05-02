#include <cmath>
#include <algorithm>

extern "C" {


double logmean(double x, double y)
{
	if (x <= 0 || y <= 0)
		return 0;

	//My approximation
	if (std::abs((x-y)/(x+y)) <  1e-3)
		return (2/3.0)*std::sqrt(x*y) + (1/6.0)*(x+y);

	return (x-y)/std::log(x/y);
}

void iterateSimulation(
		int nx, int ny,
		double dt, double ds,
		double Hz_dissipation, double absorbing_coeff,
		double mu_electron, double mu_hole,
		double D_electron, double D_hole,
		double q_n, double q_p, double E_sat,

		double* Hz, double* Hz_laplacian,
		double* Ex, double* Ey,
		double* mu_z, double* absorptivity,
		double* epsx, double* epsy,
		double* absorptivity_x, double* absorptivity_y,
		double* emfx, double* emfy,
		double* cmfx_n, double* cmfx_p,
		double* cmfy_n, double* cmfy_p,
		double* mobility_factor,

		double* rho_n, double* rho_p, double* rho_abs, double* rho_back, double* rho_free,
		double* R, double* K, double* conducting,

		double* Jx_n, double* Jx_p, double* Jx_abs,
		double* Jy_n, double* Jy_p, double* Jy_abs,
		double* conducting_x, double* conducting_y)
{
	// Laplacian of Hz
	for (int i = 1; i < nx - 2; i++) {
		for (int j = 1; j < ny - 2; j++) {
			Hz_laplacian[i * ny + j] = (
					Hz[(i + 1) * ny + j] + Hz[i * ny + j + 1] +
					Hz[(i - 1) * ny + j] + Hz[i * ny + j - 1] -
					4 * Hz[i * ny + j]) / (ds * ds);
		}
	}

	// Update Hz
	for (int i = 0; i < nx - 1; i++) {
		for (int j = 0; j < ny - 1; j++) {
			int ij = i * ny + j;
			double sigma = absorptivity[ij] * mu_z[ij] * absorbing_coeff;

			double denom = 1 + 0.5 * dt * sigma / mu_z[ij];
			double numer = Hz[ij] * (1 - 0.5 * dt * sigma / mu_z[ij]) +
					((- (Ey[(i + 1) * ny + j] - Ey[ij]) + (Ex[ij + 1] - Ex[ij])) * dt / (ds * mu_z[ij])) +
					Hz_dissipation * dt * Hz_laplacian[ij];

			Hz[ij] = numer / denom;
		}
	}

	// Update charges and mobilities
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {
			int ij = i * ny + j;
			int ijL = (i - 1) * ny + j;
			int ijB = i * ny + j - 1;

			double gen = conducting[ij] * R[ij] * (K[ij] - rho_n[ij] * rho_p[ij] / (q_n * q_p));

			rho_n[ij] -= (Jx_n[ij] - Jx_n[ijL] + Jy_n[ij] - Jy_n[ijB]) * dt / ds - dt * q_n * gen;
			rho_p[ij] -= (Jx_p[ij] - Jx_p[ijL] + Jy_p[ij] - Jy_p[ijB]) * dt / ds - dt * q_p * gen;
			rho_abs[ij] -= (Jx_abs[ij] - Jx_abs[ijL] + Jy_abs[ij] - Jy_abs[ijB]) * dt / ds;

			rho_free[ij] = rho_abs[ij] + rho_n[ij] + rho_p[ij] + rho_back[ij];

			double E_avg = std::sqrt(0.5 * (Ex[ij] * Ex[ij] + Ex[ijL] * Ex[ijL] + Ey[ij] * Ey[ij] + Ey[ijB] * Ey[ijB]));
			mobility_factor[ij] = std::min(1.0, E_sat / (E_avg > 0 ? E_avg : 1));
		}
	}

	// Update Ex
	for (int i = 0; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {
			int ij = i * ny + j;
			double ex_prev = Ex[ij];

			double mf = std::min(mobility_factor[(i + 1) * ny + j], mobility_factor[ij]);

			double sigma_n = conducting_x[ij] * mf * mu_electron * logmean(-rho_n[(i + 1) * ny + j], -rho_n[ij]);
			double sigma_p = conducting_x[ij] * mf * mu_hole * logmean(rho_p[(i + 1) * ny + j], rho_p[ij]);

			Jx_abs[ij] = 0;

			Jx_n[ij] = conducting_x[ij] * (-mf * D_electron * (rho_n[(i + 1) * ny + j] - rho_n[ij]) / ds +
					sigma_n * (emfx[ij] + cmfx_n[ij] / q_n));

			Jx_p[ij] = conducting_x[ij] * (-mf * D_hole * (rho_p[(i + 1) * ny + j] - rho_p[ij]) / ds +
					sigma_p * (emfx[ij] + cmfx_p[ij] / q_p));

			double sigma = sigma_n + sigma_p + absorptivity_x[ij] * epsx[ij] * absorbing_coeff;
			double jx = Jx_abs[ij] + Jx_n[ij] + Jx_p[ij];

			Ex[ij] = (Ex[ij] * (1 - 0.5 * dt * sigma / epsx[ij]) + ((Hz[ij] - Hz[ij - 1]) / ds - jx) * dt / epsx[ij]) /
					(1 + 0.5 * dt * sigma / epsx[ij]);

			Jx_abs[ij] += 0.5 * absorptivity_x[ij] * epsx[ij] * absorbing_coeff * (ex_prev + Ex[ij]);
			Jx_n[ij]  += 0.5 * sigma_n * (ex_prev + Ex[ij]);
			Jx_p[ij]  += 0.5 * sigma_p * (ex_prev + Ex[ij]);
		}
	}

	// Update Ey
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 0; j < ny - 1; j++) {
			int ij = i * ny + j;
			double ey_prev = Ey[ij];

			double mf = std::min(mobility_factor[i * ny + j + 1], mobility_factor[ij]);

			double sigma_n = conducting_y[ij] * mf * mu_electron * logmean(-rho_n[i * ny + j + 1], -rho_n[ij]);
			double sigma_p = conducting_y[ij] * mf * mu_hole * logmean(rho_p[i * ny + j + 1], rho_p[ij]);

			Jy_abs[ij] = 0;

			Jy_n[ij] = conducting_y[ij] * (-mf * D_electron * (rho_n[i * ny + j + 1] - rho_n[ij]) / ds +
					sigma_n * (emfy[ij] + cmfy_n[ij] / q_n));

			Jy_p[ij] = conducting_y[ij] * (-mf * D_hole * (rho_p[i * ny + j + 1] - rho_p[ij]) / ds +
					sigma_p * (emfy[ij] + cmfy_p[ij] / q_p));

			double sigma = sigma_n + sigma_p + absorptivity_y[ij] * epsy[ij] * absorbing_coeff;
			double jy = Jy_abs[ij] + Jy_n[ij] + Jy_p[ij];

			Ey[ij] = (Ey[ij] * (1 - 0.5 * dt * sigma / epsy[ij]) - ((Hz[ij] - Hz[(i - 1) * ny + j]) / ds + jy) * dt / epsy[ij]) /
					(1 + 0.5 * dt * sigma / epsy[ij]);

			Jy_abs[ij] += 0.5 * absorptivity_y[ij] * epsy[ij] * absorbing_coeff * (ey_prev + Ey[ij]);
			Jy_n[ij]  += 0.5 * sigma_n * (ey_prev + Ey[ij]);
			Jy_p[ij]  += 0.5 * sigma_p * (ey_prev + Ey[ij]);
		}
	}

	// Enforce boundary conditions
	for (int j = 0; j < ny - 1; j++) {
		Ey[j] = 0;
		Ey[(nx - 1) * ny + j] = 0;
	}

	for (int i = 0; i < nx - 1; i++) {
		Ex[i * ny + 0] = 0;
		Ex[i * ny + (ny - 1)] = 0;
	}
}

void downscale(double* source, double* dest, int steps, int nx, int ny) {
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			dest[(steps * nx + i) * ny + j] = source[i * ny + j];
		}
	}
	int nx_d = nx;
	int ny_d = ny;
	for (int k = steps-1; k >= 0; k--) {
		nx_d = nx_d/2;
		ny_d = ny_d/2;

		for (int i = 0; i < nx_d; i++) {
			for (int j = 0; j < ny_d; j++) {
				dest[(k * nx + i) * ny + j] = 0.25*(dest[((k+1) * nx + 2*i) * ny + 2*j]+dest[((k+1) * nx + (2*i+1)) * ny + 2*j]+dest[((k+1) * nx + 2*i) * ny + 2*j+1]+dest[((k+1) * nx + (2*i+1)) * ny + 2*j+1]);
			}
		}
	}
}

void JacobiIteration(int nx, int ny, int steps, int xbound, int ybound, double alpha, int fineness, bool calcPhi,
	double* MG_phi1,
	double* MG_phi2,
	double* MG_rho,
	double* MG_epsx,
	double* MG_epsy,
	double* MG_eps_avg

) {
if (calcPhi) {
	for (int p = 0; p < steps; p++) {
		for (int i = 1; i < xbound - 1; i++) {
			for (int j = 1; j < ybound - 1; j++) {
				int ij = i * ny + j;
				MG_phi2[ij] = (0.1*MG_phi1[i * ny + j] + (MG_phi1[(i - 1) * ny + j] + MG_phi1[(i + 1) * ny + j] +
						MG_phi1[i * ny + (j - 1)] + MG_phi1[i * ny + (j + 1)] +
						MG_rho[(fineness * nx + i) * ny + j] * alpha) / 4.0) / 1.1;
			}
		}
		for (int i = 1; i < xbound - 1; i++) {
			for (int j = 1; j < ybound - 1; j++) {
				int ij = i * ny + j;
				MG_phi1[ij] = (0.1*MG_phi2[i * ny + j] + (MG_phi2[(i - 1) * ny + j] + MG_phi2[(i + 1) * ny + j] +
						MG_phi2[i * ny + (j - 1)] + MG_phi2[i * ny + (j + 1)] +
						MG_rho[(fineness * nx + i) * ny + j] * alpha) / 4.0) / 1.1;
			}
		}
	}
} else {
	for (int i = 1; i < xbound - 1; i++) {
		for (int j = 1; j < ybound - 1; j++) {
			int ij = i * ny + j;
			MG_eps_avg[ij] = MG_epsx[(fineness * nx + i - 1) * ny + j] + MG_epsx[(fineness * nx + i) * ny + j] +
					MG_epsy[(fineness * nx + i) * ny + j - 1] + MG_epsy[(fineness * nx + i) * ny + j];
		}
	}

	for (int p = 0; p < steps; p++) {
		for (int i = 1; i < xbound - 1; i++) {
			for (int j = 1; j < ybound - 1; j++) {
				int ij = i * ny + j;
				MG_phi2[ij] = (0.1*MG_phi1[i * ny + j] + (MG_phi1[(i - 1) * ny + j] * MG_epsx[(fineness * nx + i - 1) * ny + j] +
						MG_phi1[(i + 1) * ny + j] * MG_epsx[(fineness * nx + i) * ny + j] +
						MG_phi1[i * ny + (j - 1)] * MG_epsy[(fineness * nx + i) * ny + j - 1] +
						MG_phi1[i * ny + (j + 1)] * MG_epsy[(fineness * nx + i) * ny + j] +
						MG_rho[(fineness * nx + i) * ny + j] * alpha) / MG_eps_avg[ij]) / 1.1;
			}
		}

		for (int i = 1; i < xbound - 1; i++) {
			for (int j = 1; j < ybound - 1; j++) {
				int ij = i * ny + j;
				MG_phi1[ij] = (0.1*MG_phi2[i * ny + j] + (MG_phi2[(i - 1) * ny + j] * MG_epsx[(fineness * nx + i - 1) * ny + j] +
						MG_phi2[(i + 1) * ny + j] * MG_epsx[(fineness * nx + i) * ny + j] +
						MG_phi2[i * ny + (j - 1)] * MG_epsy[(fineness * nx + i) * ny + j - 1] +
						MG_phi2[i * ny + (j + 1)] * MG_epsy[(fineness * nx + i) * ny + j] +
						MG_rho[(fineness * nx + i) * ny + j] * alpha) / MG_eps_avg[ij]) / 1.1;
			}
		}
	}
}
}

void multigridSolve(int nx, int ny, int log2_resolution, bool correctEfield, bool computePhi,
		double* MG_phi1,
		double* MG_phi2,
		double* MG_rho,
		double* MG_rho0,
		double* MG_epsx,
		double* MG_epsy,
		double* MG_eps_avg,
		double* epsx,
		double* epsy,
		double* Ex,
		double* Ey,
		double* phi,
		double* rho_free,
		double ds,
		double width
) {
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			MG_phi1[i * ny + j] = 0;
			for (int k = 0; k <= log2_resolution; k++) {
				MG_rho[(k * nx + i) * ny + j] = 0;
			}
		}
	}

	if (correctEfield) {
		for (int i = 1; i < nx - 1; i++) {
			for (int j = 1; j < ny - 1; j++) {
				MG_rho0[i * ny + j] = ((Ex[i * ny + j] * epsx[i * ny + j] - Ex[(i - 1) * ny + j] * epsx[(i - 1) * ny + j] +
						Ey[i * ny + j] * epsy[i * ny + j] - Ey[i * ny + j - 1] * epsy[i * ny + j - 1]) / ds) -
						rho_free[i * ny + j];
			}
		}
	}

	if (computePhi) {
		for (int i = 1; i < nx - 1; i++) {
			for (int j = 1; j < ny - 1; j++) {
				MG_rho0[i * ny + j] = (Ex[i * ny + j] - Ex[(i - 1) * ny + j] + Ey[i * ny + j] - Ey[i * ny + j - 1]) / ds +
						(phi[(i + 1) * ny + j] + phi[i * ny + j + 1] + phi[(i - 1) * ny + j] + phi[i * ny + j - 1] - 4 * phi[i * ny + j]) / (ds * ds);
			}
		}
	}

	int stepsarray[9] = {0, 0, 200, 200, 200, 200, 200, 50, 20};
	
	downscale(MG_rho0, MG_rho, log2_resolution, nx, ny);
	
	for (int fineness = 2; fineness <= log2_resolution; fineness++) {
		int nx_tmp = (1 << fineness);
		int ny_tmp = (1 << fineness);
		double gridsize = width/(1 << fineness);
		
		int poissonsteps = stepsarray[std::min(fineness, 8)];
		double alpha = (gridsize*gridsize);

		JacobiIteration(nx, ny, poissonsteps, nx_tmp, ny_tmp, alpha, fineness, computePhi,
				MG_phi1, MG_phi2, MG_rho, MG_epsx, MG_epsy, MG_eps_avg);

		if (fineness == log2_resolution)
			break;

		for (int i = 0; i < nx_tmp; i++) {
			for (int j = 0; j < ny_tmp; j++) {
				MG_phi2[i * ny + j] = MG_phi1[i * ny + j];
			}
		}
		
		for (int i = 0; i < nx_tmp; i++)
		{
			for (int j = 0; j < nx_tmp; j++) {
				MG_phi1[(2*i) * ny + 2*j] = MG_phi2[i * ny + j];
				MG_phi1[(2*i+1) * ny + 2*j] = MG_phi2[i * ny + j];
				MG_phi1[(2*i) * ny + 2*j+1] = MG_phi2[i * ny + j];
				MG_phi1[(2*i+1) * ny + 2*j+1] = MG_phi2[i * ny + j];
			}
		}
		
		/*for (int i = 1; i < 2*nx_tmp-1; i++)
		{
			for (int j = 1; j < 2*nx_tmp-1; j++) {
				MG_phi1[i][j] = this.bilinearinterp(MG_phi2, (i-0.5)/2.0, (j-0.5)/2.0);
			}
		}*/
		
		for (int i = 0; i < 2*nx_tmp; i++) {
			for (int j = 0; j < 2*ny_tmp; j++) {
				MG_phi2[i * ny + j] = MG_phi1[i * ny + j];
			}
		}
	}


	if (correctEfield) {
		for (int i = 0; i < nx - 1; i++) {
			for (int j = 0; j < ny - 1; j++) {
				Ex[i * ny + j] += (MG_phi1[(i + 1) * ny + j] - MG_phi1[i * ny + j]) / ds;
				Ey[i * ny + j] += (MG_phi1[i * ny + j + 1] - MG_phi1[i * ny + j]) / ds;
			}
		}
	}

	if (computePhi) {
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				phi[i * ny + j] += MG_phi1[i * ny + j];
			}
		}
	}
}

} // extern "C"
