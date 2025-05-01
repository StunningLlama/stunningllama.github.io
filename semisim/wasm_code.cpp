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

//void multigridSolve(bool correctEfield, bool computePhi) {
//	for (int i = 0; i < nx; i++) {
//		for (int j = 0; j < ny; j++) {
//			MG_phi1[i][j] = 0;
//			MG_rho[i][j] = 0;
//		}
//	}
//
//	if (correctEfield)
//	{
//		for (int i = 1; i < nx-1; i++) {
//			for (int j = 1; j < ny-1; j++) {
//				MG_rho0[i][j] = ((Ex[i][j]*epsx[i][j]-Ex[i-1][j]*epsx[i-1][j] + Ey[i][j]*epsy[i][j]-Ey[i][j-1]*epsy[i][j-1])/ds) - rho_free[i][j];
//			}
//		}
//	}
//	if (computePhi) {
//		for (int i = 1; i < nx-1; i++) {
//			for (int j = 1; j < ny-1; j++) {
//				MG_rho0[i][j] = (Ex[i][j]-Ex[i-1][j] + Ey[i][j]-Ey[i][j-1])/(ds) + (phi[i+1][j]+phi[i][j+1]+phi[i-1][j]+phi[i][j-1]-4*phi[i][j])/(ds*ds);
//			}
//		}
//	}
//
//	int nxint = nx-2;
//	int nyint = ny-2;
//
//	int[] stepsarray = {0, 0, 200, 200, 200, 200, 200, 50, 20, 20, 20};
//
//	int maxfineness = (int)std::floor(std::log(nx)/std::log(2));
//	for (int fineness = 2; fineness <= maxfineness; fineness++) {
//		int nxcgint = (1 << fineness) - 1;
//		int nycgint = (1 << fineness) - 1;
//		double gridsize = width/(1 << fineness);
//
//
//		downscale(MG_rho0, MG_rho, nx-1, ny-1, nxcgint+1, nycgint+1, 0, 0, ds, 0, 0, gridsize);
//
//		int poissonsteps = stepsarray[fineness];
//		double alpha = (gridsize*gridsize);
//		double beta = 4;
//
//		JacobiIteration(poissonsteps, nxcgint, nycgint, alpha, beta, gridsize, fineness, debug, computePhi);
//
//		if (fineness < maxfineness) {
//
//			for (int i = 0; i < nxcgint+2; i++) {
//				for (int j = 0; j < nycgint+2; j++) {
//					MG_phi2[i][j] = MG_phi1[i][j];
//				}
//			}
//
//			nxcgint = (1 << (fineness+1)) - 1;
//			nycgint = (1 << (fineness+1)) - 1;
//			gridsize = width/(1 << (fineness+1));
//
//
//			for (int i = 1; i < nxcgint+1; i++)
//			{
//				for (int j = 1; j < nycgint+1; j++) {
//					MG_phi1[i][j] = bilinearinterp(MG_phi2, i/2.0, j/2.0);
//				}
//			}
//
//
//			for (int i = 0; i < nxcgint+2; i++) {
//				MG_phi2[i][0] = 0;
//				MG_phi2[i][nycgint+1] = 0;
//			}
//			for (int j = 0; j < nycgint+2; j++) {
//				MG_phi2[0][j] = 0;
//				MG_phi2[nxcgint+1][j] = 0;
//			}
//		} else {
//			for (int i = 0; i < nxcgint+2; i++) {
//				for (int j = 0; j < nycgint+2; j++) {
//					MG_phi2[i][j] = MG_phi1[i][j];
//				}
//			}
//
//			for (int i = 1; i < nxint+1; i++)
//			{
//				for (int j = 1; j < nyint+1; j++) {
//					MG_phi1[i][j] = bilinearinterp(MG_phi2, i*(double)(nxcgint+1.0)/((double)nx-1.0), j*(double)(nycgint+1.0)/((double)ny-1.0));
//				}
//			}
//
//		}
//	}
//
//
//	int poissonsteps = stepsarray[maxfineness+1];
//
//	double alpha = (ds*ds);
//	double beta = 4;
//
//	for (int i = 0; i < nx; i++) {
//		for (int j = 0; j < ny; j++) {
//			MG_rho[i][j] = MG_rho0[i][j];
//			MG_epsx[maxfineness+1][i][j] = epsx[i][j];
//			MG_epsy[maxfineness+1][i][j] = epsy[i][j];
//		}
//	}
//
//	JacobiIteration(poissonsteps, nxint, nyint, alpha, beta, ds, maxfineness+1, debug, computePhi);
//
//	if (correctEfield) {
//		for (int i = 0; i < nx-1; i++)
//		{
//			for (int j = 0; j < ny-1; j++)
//			{
//				Ex[i][j] = Ex[i][j] + (MG_phi1[i+1][j]-MG_phi1[i][j])/ds;
//				Ey[i][j] = Ey[i][j] + (MG_phi1[i][j+1]-MG_phi1[i][j])/ds;
//			}
//		}
//	}
//	if (computePhi)
//	{
//		for (int i = 0; i < nx; i++)
//		{
//			for (int j = 0; j < ny; j++)
//			{
//				phi[i][j] = phi[i][j] + MG_phi1[i][j];
//			}
//		}
//	}
//}
//
//
//void JacobiIteration(int steps, int xbound, int ybound, double alpha, double beta, double gridsize, int fineness, boolean debug, boolean calcPhi) {
//
//	if (calcPhi) {
//		for (int poissonit = 0; poissonit < steps; poissonit++) {
//			for (int i = 1; i < xbound+1; i++) {
//				for (int j = 1; j < ybound+1; j++) {
//					MG_phi2[i][j] = (MG_phi1[i-1][j] + MG_phi1[i+1][j] + MG_phi1[i][j-1] + MG_phi1[i][j+1] + MG_rho[i][j]*alpha)/4.0;
//				}
//			}
//			for (int i = 1; i < xbound+1; i++) {
//				for (int j = 1; j < ybound+1; j++) {
//					MG_phi1[i][j] = ((MG_phi2[i-1][j] + MG_phi2[i+1][j] + MG_phi2[i][j-1] + MG_phi2[i][j+1]) + MG_rho[i][j]*alpha)/4.0;
//				}
//			}
//		}
//	} else {
//		for (int i = 1; i < xbound+1; i++) {
//			for (int j = 1; j < ybound+1; j++) {
//				MG_eps_avg[i][j] = (MG_epsx[fineness][i-1][j]+MG_epsx[fineness][i][j]+MG_epsy[fineness][i][j-1]+MG_epsy[fineness][i][j]);
//			}
//		}
//
//		for (int poissonit = 0; poissonit < steps; poissonit++) {
//			for (int i = 1; i < xbound+1; i++) {
//				for (int j = 1; j < ybound+1; j++) {
//					MG_phi2[i][j] = ((MG_phi1[i-1][j]*MG_epsx[fineness][i-1][j]
//																			 + MG_phi1[i+1][j]*MG_epsx[fineness][i][j]
//																													+ MG_phi1[i][j-1]*MG_epsy[fineness][i][j-1]
//																																						   + MG_phi1[i][j+1]*MG_epsy[fineness][i][j])
//							+ MG_rho[i][j]*alpha)/MG_eps_avg[i][j];
//				}
//			}
//			for (int i = 1; i < xbound+1; i++) {
//				for (int j = 1; j < ybound+1; j++) {
//					MG_phi1[i][j] = ((MG_phi2[i-1][j]*MG_epsx[fineness][i-1][j]
//																			 + MG_phi2[i+1][j]*MG_epsx[fineness][i][j]
//																													+ MG_phi2[i][j-1]*MG_epsy[fineness][i][j-1]
//																																						   + MG_phi2[i][j+1]*MG_epsy[fineness][i][j])
//							+ MG_rho[i][j]*alpha)/MG_eps_avg[i][j];
//				}
//			}
//		}
//	}
//}
//
//void downscale(double[][] source, double[][] dest, int sx, int sy, int dx, int dy, double soffsetx, double soffsety, double sspacing, double doffsetx, double doffsety, double dspacing) {
//	for (int i = 0; i <= dx; i++) {
//		for (int j = 0; j < dy; j++) {
//			dest[i][j] = 0;
//		}
//	}
//	double scalefactor = (sspacing*sspacing/(dspacing*dspacing));
//	for (int i = 0; i <= sx; i++) {
//		for (int j = 0; j <= sy; j++) {
//			double cx = (i*sspacing + soffsetx - doffsetx)/dspacing;
//			double cy = (j*sspacing + soffsety - doffsety)/dspacing;
//			int xfloor = (int)std::floor(cx);
//			int yfloor = (int)std::floor(cy);
//			double fx = cx - xfloor;
//			double fy = cy - yfloor;
//			if (xfloor < 0) {
//				xfloor = 0;
//				fx = 0.0;
//			} else if (xfloor >= dx) {
//				xfloor = dx - 1;
//				fx = 1.0;
//			}
//			if (yfloor < 0) {
//				yfloor = 0;
//				fy = 0.0;
//			} else if (yfloor >= dy) {
//				yfloor = dy - 1;
//				fy = 1.0;
//			}
//			double srcval = source[i][j]*scalefactor;
//			dest[xfloor][yfloor] += srcval*(1-fx)*(1-fy);
//			dest[xfloor+1][yfloor] += srcval*fx*(1-fy);
//			dest[xfloor][yfloor+1] += srcval*(1-fx)*fy;
//			dest[xfloor+1][yfloor+1] += srcval*fx*fy;
//		}
//	}
//}

void JacobiIteration(int nx, int ny, int steps, int xbound, int ybound, double alpha, double beta, double gridsize, int fineness, bool calcPhi,
		double* MG_phi1,
		double* MG_phi2,
		double* MG_rho,
		double* MG_epsx,
		double* MG_epsy,
		double* MG_eps_avg

) {
	if (calcPhi) {
		for (int p = 0; p < steps; p++) {
			for (int i = 1; i < xbound + 1; i++) {
				for (int j = 1; j < ybound + 1; j++) {
					int ij = i * ny + j;
					MG_phi2[ij] = (MG_phi1[(i - 1) * ny + j] + MG_phi1[(i + 1) * ny + j] +
							MG_phi1[i * ny + (j - 1)] + MG_phi1[i * ny + (j + 1)] +
							MG_rho[ij] * alpha) / 4.0;
				}
			}
			for (int i = 1; i < xbound + 1; i++) {
				for (int j = 1; j < ybound + 1; j++) {
					int ij = i * ny + j;
					MG_phi1[ij] = (MG_phi2[(i - 1) * ny + j] + MG_phi2[(i + 1) * ny + j] +
							MG_phi2[i * ny + (j - 1)] + MG_phi2[i * ny + (j + 1)] +
							MG_rho[ij] * alpha) / 4.0;
				}
			}
		}
	} else {
		for (int i = 1; i < xbound + 1; i++) {
			for (int j = 1; j < ybound + 1; j++) {
				int ij = i * ny + j;
				MG_eps_avg[ij] = MG_epsx[(fineness * nx + i - 1) * ny + j] + MG_epsx[(fineness * nx + i) * ny + j] +
						MG_epsy[(fineness * nx + i) * ny + j - 1] + MG_epsy[(fineness * nx + i) * ny + j];
			}
		}

		for (int p = 0; p < steps; p++) {
			for (int i = 1; i < xbound + 1; i++) {
				for (int j = 1; j < ybound + 1; j++) {
					int ij = i * ny + j;
					MG_phi2[ij] = (MG_phi1[(i - 1) * ny + j] * MG_epsx[(fineness * nx + i - 1) * ny + j] +
							MG_phi1[(i + 1) * ny + j] * MG_epsx[(fineness * nx + i) * ny + j] +
							MG_phi1[i * ny + (j - 1)] * MG_epsy[(fineness * nx + i) * ny + j - 1] +
							MG_phi1[i * ny + (j + 1)] * MG_epsy[(fineness * nx + i) * ny + j] +
							MG_rho[ij] * alpha) / MG_eps_avg[ij];
				}
			}

			for (int i = 1; i < xbound + 1; i++) {
				for (int j = 1; j < ybound + 1; j++) {
					int ij = i * ny + j;
					MG_phi1[ij] = (MG_phi2[(i - 1) * ny + j] * MG_epsx[(fineness * nx + i - 1) * ny + j] +
							MG_phi2[(i + 1) * ny + j] * MG_epsx[(fineness * nx + i) * ny + j] +
							MG_phi2[i * ny + (j - 1)] * MG_epsy[(fineness * nx + i) * ny + j - 1] +
							MG_phi2[i * ny + (j + 1)] * MG_epsy[(fineness * nx + i) * ny + j] +
							MG_rho[ij] * alpha) / MG_eps_avg[ij];
				}
			}
		}
	}
}

void downscale(int nx, int ny, double* source, double* dest, int sx, int sy, int dx, int dy, double soffsetx, double soffsety, double sspacing, double doffsetx, double doffsety, double dspacing) {
	for (int i = 0; i <= dx; i++) {
		for (int j = 0; j < dy; j++) {
			dest[i * ny + j] = 0;
		}
	}
	double scalefactor = (sspacing * sspacing / (dspacing * dspacing));
	for (int i = 0; i <= sx; i++) {
		for (int j = 0; j <= sy; j++) {
			double cx = (i * sspacing + soffsetx - doffsetx) / dspacing;
			double cy = (j * sspacing + soffsety - doffsety) / dspacing;
			int xfloor = (int) std::floor(cx);
			int yfloor = (int) std::floor(cy);
			double fx = cx - xfloor;
			double fy = cy - yfloor;
			if (xfloor < 0) {
				xfloor = 0;
				fx = 0.0;
			} else if (xfloor >= dx) {
				xfloor = dx - 1;
				fx = 1.0;
			}
			if (yfloor < 0) {
				yfloor = 0;
				fy = 0.0;
			} else if (yfloor >= dy) {
				yfloor = dy - 1;
				fy = 1.0;
			}
			double srcval = source[i * ny + j] * scalefactor;
			dest[xfloor * ny + yfloor] += srcval * (1 - fx) * (1 - fy);
			dest[(xfloor + 1) * ny + yfloor] += srcval * fx * (1 - fy);
			dest[xfloor * ny + (yfloor + 1)] += srcval * (1 - fx) * fy;
			dest[(xfloor + 1) * ny + (yfloor + 1)] += srcval * fx * fy;
		}
	}
}

double bilinearinterp(int nx, int ny, double* array, double x, double y) {
	int xfloor = (int)std::floor(x);
	int yfloor = (int)std::floor(y);
	double fx = x - xfloor;
	double fy = y - yfloor;
	if (std::abs(x-std::round(x)) < 1e-2 || std::abs(y-std::round(y)) < 1e-2) {
		int i = (int)std::round(x);
		int j = (int)std::round(y);
		if (i < 0) i = 0;
		if (j < 0) j = 0;
		if (i >= nx) i = nx - 1;
		if (j >= ny) j = ny - 1;
		return array[i*ny + j];
	}

	if (xfloor < 0) {
		xfloor = 0;
		fx = 0.0;
	} else if (xfloor >= nx - 1) {
		xfloor = nx - 2;
		fx = 1.0;
	}
	if (yfloor < 0) {
		yfloor = 0;
		fy = 0.0;
	} else if (yfloor >= ny - 1) {
		yfloor = ny - 2;
		fy = 1.0;
	}
	double va = array[xfloor*ny + yfloor]*(1.0-fx) + array[(xfloor+1) * ny + yfloor]*fx;
	double vb = array[xfloor * ny + yfloor+1]*(1.0-fx) + array[(xfloor+1) * ny + yfloor+1]*fx;

	return va*(1.0-fy) + vb*fy;
}


void multigridSolve(int nx, int ny, bool correctEfield, bool computePhi,
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
			MG_rho[i * ny + j] = 0;
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

	int nxint = nx - 2;
	int nyint = ny - 2;
	int stepsarray[11] = {0, 0, 200, 200, 200, 200, 200, 50, 20, 20, 20};
	int maxfineness = (int) std::floor(std::log(nx) / std::log(2));

	for (int fineness = 2; fineness <= maxfineness; fineness++) {
		int nxcgint = (1 << fineness) - 1;
		int nycgint = (1 << fineness) - 1;
		double gridsize = width / (1 << fineness);

		downscale(nx, ny, MG_rho0, MG_rho, nx - 1, ny - 1, nxcgint + 1, nycgint + 1, 0, 0, ds, 0, 0, gridsize);

		int poissonsteps = stepsarray[fineness];
		double alpha = gridsize * gridsize;
		double beta = 4;

		JacobiIteration(nx, ny, poissonsteps, nxcgint, nycgint, alpha, beta, gridsize, fineness, computePhi,
				MG_phi1, MG_phi2, MG_rho, MG_epsx, MG_epsy, MG_eps_avg);

		if (fineness < maxfineness) {
			for (int i = 0; i < nxcgint + 2; i++) {
				for (int j = 0; j < nycgint + 2; j++) {
					MG_phi2[i * ny + j] = MG_phi1[i * ny + j];
				}
			}

			nxcgint = (1 << (fineness + 1)) - 1;
			nycgint = (1 << (fineness + 1)) - 1;
			gridsize = width / (1 << (fineness + 1));

			for (int i = 1; i < nxcgint + 1; i++) {
				for (int j = 1; j < nycgint + 1; j++) {
					MG_phi1[i * ny + j] = bilinearinterp(nx, ny, MG_phi2, i / 2.0, j / 2.0);
				}
			}

			for (int i = 0; i < nxcgint + 2; i++) {
				MG_phi2[i * ny] = 0;
				MG_phi2[i * ny + nycgint + 1] = 0;
			}
			for (int j = 0; j < nycgint + 2; j++) {
				MG_phi2[j] = 0;
				MG_phi2[(nxcgint + 1) * ny + j] = 0;
			}
		} else {
			for (int i = 0; i < nxcgint + 2; i++) {
				for (int j = 0; j < nycgint + 2; j++) {
					MG_phi2[i * ny + j] = MG_phi1[i * ny + j];
				}
			}

			for (int i = 1; i < nxint + 1; i++) {
				for (int j = 1; j < nyint + 1; j++) {
					MG_phi1[i * ny + j] = bilinearinterp(nx, ny, MG_phi2, i * (double) (nxcgint + 1) / (nx - 1.0),
							j * (double) (nycgint + 1) / (ny - 1.0));
				}
			}
		}
	}

	int poissonsteps = stepsarray[maxfineness + 1];
	double alpha = ds * ds;
	double beta = 4;

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			MG_rho[i * ny + j] = MG_rho0[i * ny + j];
			MG_epsx[((maxfineness + 1) * nx + i) * ny + j] = epsx[i * ny + j];
			MG_epsy[((maxfineness + 1) * nx + i) * ny + j] = epsy[i * ny + j];
		}
	}

	JacobiIteration(nx, ny, poissonsteps, nxint, nyint, alpha, beta, ds, maxfineness + 1, computePhi,
			MG_phi1, MG_phi2, MG_rho, MG_epsx, MG_epsy, MG_eps_avg
	);

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
