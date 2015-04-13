// TRACKING_MP
// ===========
// Author: 		Ximenes R. Resende
// email:  		xresende@gmail.com, ximenes.resende@lnls.br
// affiliation:	LNLS - Laboratorio Nacional de Luz Sincrotron
// Date: 		2012-08-20

#include "fmap_mp.h"

static double  fmap_mp_energy;
static int     fmap_mp_nturn;
static bool    fmap_mp_diffusion;

void tracking_mp_fmap(int cpu_id, int data_id, void *data_in, void *shm_data_out);


void fmap_mp(int nr_cpus, long Nbx, long Nbz, long Nbtour, double x0, double xmax,
		  double z0, double zmax, double energy, bool diffusion)
{

	double tiny_amp = 1e-6;  // meters

	int    nr_pnts;
	double *points;
	double *result;
	const char fic[] = "fmap.out";

	nr_pnts = Nbx * Nbz;

	points = (double*) malloc (2 * nr_pnts * sizeof(double));

	int k = 0;
	for(int i=0; i<Nbx; ++i) {
		double x = x0 + i * (xmax - x0) / (Nbx - 1);
		if (fabs(x) < tiny_amp) { x = tiny_amp; }
		for(int j=0; j<Nbz; ++j) {
			double z = z0 + j * (zmax - z0) / (Nbz - 1);
			if (fabs(z) < tiny_amp) { z = tiny_amp; }
			points[2*k+0] = x;
			points[2*k+1] = z;
			k++;
		}
	}


	result = (double*) malloc (4 * nr_pnts * sizeof(double));

	fmap_mp_energy    = energy;
	fmap_mp_nturn     = Nbtour;
	fmap_mp_diffusion = diffusion;
	mp_task_mgr(nr_cpus, nr_pnts, (void*) points, (void*) result, 4*nr_pnts*sizeof(double), tracking_mp_fmap);

	
	// salva resultado em arquivo
	FILE * outf;
	/* Opening file */
	if ((outf = fopen(fic, "w")) == nullptr) {
		fprintf(stdout, "fmap: error while opening file %s\n", fic);
		exit(1);
	}
	fprintf(outf,"# TRACY III Sirius-- %s --\n", fic);
	fprintf(outf,"# nu = f(x) \n");
	fprintf(outf,"#    x[m]          z[m]           fx             fz            dfx            dfz\n");
	k = 0;
	for(int i=0; i<Nbx; ++i) {
		for(int j=0; j<Nbz; ++j) {
			double x = points[2*k+0];
			double z = points[2*k+1];
			double nux1 = result[4*k+0];
			double nuz1 = result[4*k+1];
			double dfx  = result[4*k+2];
			double dfz  = result[4*k+3];
			if (fmap_mp_diffusion) {
				fprintf(outf,"%10.6e %10.6e %10.6e %10.6e %10.6e %10.6e\n", x, z, nux1, nuz1, dfx, dfz);
			}else{
				fprintf(outf,"%10.6e %10.6e %10.6e %10.6e\n", x, z, nux1, nuz1);
			}
			k++;
		}
	}
	fclose(outf);



	free(points);
	free(result);

}




#define NTERM2  10
void tracking_mp_fmap(int cpu_id, int data_id, void *data_in, void *shm_data_out) {

	/*
	double xp = 0;
	double zp = 0;
	double ctau = 0;
	double fx[NTERM2], fz[NTERM2], fx2[NTERM2], fz2[NTERM2], dfx, dfz;
	double nux1 = 0.0, nuz1 = 0.0, nux2 = 0.0, nuz2 = 0.0;
	int    nb_freq[2] = {0, 0};
	double x = ((double*)data_in)[data_id * 2 + 0];
	double z = ((double*)data_in)[data_id * 2 + 1];
	bool   status = true;
	double Tab[DIM][NTURN];
	double Tab0[DIM][NTURN];
	long nturn = fmap_mp_nturn;


	if (fmap_mp_diffusion) nturn = 2*fmap_mp_nturn;

	if (!globval.Cavity_on) Trac_Simple4DCOD(x,xp,z,zp,fmap_mp_energy,ctau,nturn,Tab,&status);
	else Trac_Simple6DCOD(x,xp,z,zp,fmap_mp_energy,ctau,nturn,Tab,&status);

	if (status) { // if trajectory is stable
		// gets frequency vectors
		Get_NAFF(NTERM2, fmap_mp_nturn, Tab, fx, fz, nb_freq);
		Get_freq(fx,fz,&nux1,&nuz1);  // gets nux and nuz
		if (fmap_mp_diffusion) {
			Get_Tabshift(Tab,Tab0,fmap_mp_nturn,fmap_mp_nturn);
			// gets frequency vectors
			Get_NAFF(NTERM2, fmap_mp_nturn, Tab0, fx2, fz2, nb_freq);
			Get_freq(fx2,fz2,&nux2,&nuz2); // gets nux and nuz
		}
	} else { //zeroing output
		nux1 = 0.0; nuz1 = 0.0;
		nux2 = 0.0; nuz2 = 0.0;
	}
	dfx = nux1 - nux2; dfz = nuz1 - nuz2;
	if (fmap_mp_diffusion) {
		printf("<fmap_mp: cpu %i, data %i> %10.6e %10.6e, %10.6e %10.6e, %10.6e %10.6e\n", cpu_id, data_id, x, z, nux1, nuz1, dfx, dfz);
	}else{
		printf("<fmap_mp: cpu %i, data %i> %10.6e %10.6e, %10.6e %10.6e\n", cpu_id, data_id, x, z, nux1, nuz1);
	}


	((double*)shm_data_out)[4*data_id + 0] = nux1;
	((double*)shm_data_out)[4*data_id + 1] = nuz1;
	((double*)shm_data_out)[4*data_id + 2] = dfx;
	((double*)shm_data_out)[4*data_id + 3] = dfz;
	*/

}
#undef NTERM2


	
