// TRACKING_MP
// ===========
// Author: 		Ximenes R. Resende
// email:  		xresende@gmail.com, ximenes.resende@lnls.br
// affiliation:	LNLS - Laboratorio Nacional de Luz Sincrotron
// Date: 		2012-08-20

#include "fmapdp_mp.h"


static double  fmapdp_mp_z;
static int     fmapdp_mp_nturn;
static bool    fmapdp_mp_diffusion;

void tracking_mp_fmapdp(int cpu_id, int data_id, void *data_in, void *shm_data_out);


void fmapdp_mp(int nr_cpus, long Nbx, long Nbe, long Nbtour, double x0, double xmax,
		double emin, double emax, double z, bool diffusion)
{

	double tiny_amp = 1e-6;  // meters

	int    nr_pnts;
	double *points;
	double *result;
	const char fic[] = "fmapdp.out";

	nr_pnts = Nbx * Nbe;

	points = (double*) malloc (2 * nr_pnts * sizeof(double));

	int k = 0;
	for(int i=0; i<Nbe; ++i) {
		double dp = emin + i * (emax - emin) / (Nbe - 1);
		for(int j=0; j<Nbx; ++j) {
			double x = x0 + j * (xmax - x0) / (Nbx - 1);
			if (fabs(x) < tiny_amp) { x = tiny_amp; }
			points[2*k+0] = dp;
			points[2*k+1] = x;
			k++;
		}
	}


	result = (double*) malloc (4 * nr_pnts * sizeof(double));

	fmapdp_mp_z         = z;
	fmapdp_mp_nturn     = Nbtour;
	fmapdp_mp_diffusion = diffusion;
	mp_task_mgr(nr_cpus, nr_pnts, (void*) points, (void*) result, 4*nr_pnts*sizeof(double), tracking_mp_fmapdp);

	
	// salva resultado em arquivo
	FILE * outf;
	/* Opening file */
	if ((outf = fopen(fic, "w")) == nullptr) {
		fprintf(stdout, "fmap: error while opening file %s\n", fic);
		exit(1);
	}
	fprintf(outf,"# TRACY III Sirius-- %s --\n", fic);
	fprintf(outf,"# nu = f(x) \n");
	fprintf(outf,"#    dp[m]         x[m]           fx            fz           dfx           dfz\n");
	k = 0;
	for(int i=0; i<Nbe; ++i) {
		for(int j=0; j<Nbx; ++j) {
			double dp = points[2*k+0];
			double x  = points[2*k+1];
			double nux1 = result[4*k+0];
			double nuz1 = result[4*k+1];
			double dfx  = result[4*k+2];
			double dfz  = result[4*k+3];
			if (fmapdp_mp_diffusion) {
				fprintf(outf,"%10.6e %10.6e %10.6e %10.6e %10.6e %10.6e\n", dp, x, nux1, nuz1, dfx, dfz);
			}else{
				fprintf(outf,"%10.6e %10.6e %10.6e %10.6e\n", dp, x, nux1, nuz1);
			}
			k++;
		}
	}
	fclose(outf);



	free(points);
	free(result);

}




#define NTERM2  10
void tracking_mp_fmapdp(int cpu_id, int data_id, void *data_in, void *shm_data_out) {

	/*
	double xp = 0;
	double zp = 0;
	double ctau = 0;
	double fx[NTERM2], fz[NTERM2], fx2[NTERM2], fz2[NTERM2], dfx, dfz;
	double nux1 = 0.0, nuz1 = 0.0, nux2 = 0.0, nuz2 = 0.0;
	int    nb_freq[2] = {0, 0};
	double dp = ((double*)data_in)[data_id * 2 + 0];
	double x  = ((double*)data_in)[data_id * 2 + 1];
	double z  = fmapdp_mp_z;
	bool   status = true;
	double Tab[DIM][NTURN];
	double Tab0[DIM][NTURN];
	long nturn = fmapdp_mp_nturn;	
	
	// IF 6D Tracking diffusion turn off
	if (globval.Cavity_on == true) 	fmapdp_mp_diffusion = false;
	if (fmapdp_mp_diffusion) nturn = 2*fmapdp_mp_nturn;


	if (!globval.Cavity_on) Trac_Simple4DCOD(x,xp,z,zp,dp,ctau,nturn,Tab,&status);
	else Trac_Simple6DCOD(x,xp,z,zp,dp,ctau,nturn,Tab,&status);

	if (status) { // if trajectory is stable
		Get_NAFF(NTERM2, fmapdp_mp_nturn, Tab, fx, fz, nb_freq);
		Get_freq(fx,fz,&nux1,&nuz1);  // gets nux and nuz
		if (fmapdp_mp_diffusion) {
			Get_Tabshift(Tab,Tab0,fmapdp_mp_nturn,fmapdp_mp_nturn); // shift data for second round NAFF
			Get_NAFF(NTERM2, fmapdp_mp_nturn, Tab0, fx2, fz2, nb_freq); // gets frequency vectors
			Get_freq(fx2,fz2,&nux2,&nuz2);// gets nux and nuz
		}
	} else { //zeroing output
		nux1 = 0.0; nuz1 = 0.0;
		nux2 = 0.0; nuz2 = 0.0;
	}
	dfx = nux1 - nux2; dfz = nuz1 - nuz2;
	if (fmapdp_mp_diffusion) {
		printf("<fmapdp_mp: cpu %i, data %i> %10.6e %10.6e, %10.6e %10.6e, %10.6e %10.6e\n", cpu_id, data_id, dp, x, nux1, nuz1, dfx, dfz);
	}else{
		printf("<fmapdp_mp: cpu %i, data %i> %10.6e %10.6e, %10.6e %10.6e\n", cpu_id, data_id, dp, x, nux1, nuz1);
	}


	((double*)shm_data_out)[4*data_id + 0] = nux1;
	((double*)shm_data_out)[4*data_id + 1] = nuz1;
	((double*)shm_data_out)[4*data_id + 2] = dfx;
	((double*)shm_data_out)[4*data_id + 3] = dfz;
	*/

}
#undef NTERM2


	
