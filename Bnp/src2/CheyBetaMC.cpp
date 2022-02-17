/*--------------------------------------------------------------
	FILE: CheyBetaMC.c
	PURPOSE: Cheyette beta Monte-Carlo
	AUTHOR: Dimitri Mayevski
	DATE: 26/02/2003
  --------------------------------------------------------------*/

#include "srt_h_all.h"
#include "grf_h_all.h"
#include "math.h"
#include "RandomGen.h"
#include "CheyBetaMC.h"

static void CopyToGD(long path, GrfnDeal *gd, double ***cells, GrfnSSParam *vars)
{
	memcpy(&gd->cells[0][0], &cells[path][0][0], gd->sslength * gd->sswidth * sizeof(double));
	memcpy(&gd->ssparam, &vars[path], sizeof(GrfnSSParam));
}

static void CopyFromGD(long path, GrfnDeal *gd, double ***cells, GrfnSSParam *vars)
{
	memcpy(&cells[path][0][0], &gd->cells[0][0], gd->sslength * gd->sswidth * sizeof(double));
	memcpy(&vars[path], &gd->ssparam, sizeof(GrfnSSParam));
}

#define MAXCPN 1000

Err CheyBetaMC_Price(SProductDesc *g, int nt, long npaths, double cutcoef,
					 SCheyBeta *pmdl, int rgmode, double *pv, double *stddev)
{
	Err			err = NULL;
	int			ntimes, i, j, k, iex = -1, isig = 0;
	long		d1, d2, dT, seed = -12345678;
	double		*times = NULL;
	double		fr, xcut_u, xcut_d, phi_mid = 0.0, gamT, dffT, xj, num;
	double		dt, sqrtdt, dw, sig, beta, lam;
	double		std, cumxvar = 0.0;
	double		dffs[MAXCPN], gams[MAXCPN], dfs[MAXCPN], *dfs_ = dfs;
	SMktData	mkt_data = { &dfs_, MKT_IRSTANDARD, NULL };
	double		*x = NULL, *phi = NULL, **pvs = NULL;
	SRandomGen	rg;
	// Containers to save cells and vars from GRFN deal to keep track of path-dependence:
	double		***cells = NULL;
	GrfnSSParam	*vars = NULL;
	GrfnDeal	*gd;

	memset(&rg, 0, sizeof(SRandomGen));

	if (cutcoef > 50.0) cutcoef = 50.0;
	if (g->nccy > 1) return serror("Error (CheyBeta MC): Only one underlying allowed");

	// Fill times vector
	ntimes = g->nex + pmdl->nsig + 1;
	times = (double *) calloc(ntimes, sizeof(double));
	times[0] = 0.0;
	memcpy(times + 1, g->ex, g->nex * sizeof(double));
	memcpy(times + g->nex + 1, pmdl->sigtms, pmdl->nsig * sizeof(double));

	num_f_sort_vector(ntimes, times);
	num_f_unique_vector(&ntimes, times);
	while (times[ntimes-1] > g->ex[g->nex-1] + 1e-5) ntimes--;
	num_f_fill_vector_newalgo( &ntimes, &times,
		(int)(nt * times[ntimes-1] + 1e-5) );

	// Allocate memory
	x = (double *) calloc(npaths, sizeof(double));
	phi = (double *) calloc(npaths, sizeof(double));
	pvs = dmatrix(0, npaths-1, 0, g->ninst-1);

	if (!x || !phi || !pvs) { err = serror("Memory failure in CheyBetaMC_Price");  goto FREE_RETURN; }
	memset(x, 0, npaths * sizeof(double));
	memset(phi, 0, npaths * sizeof(double));
	memset(&pvs[0][0], 0, npaths * g->ninst * sizeof(double));

	if (g->type == PRODUCT_GRFN)
	{
		gd = ((FIRSTAllMkts *) g->spec_desc) -> gd;
		cells = f3tensor(0, npaths-1, 0, gd->sslength-1, 0, gd->sswidth-1);
		vars = (GrfnSSParam *) calloc(npaths, sizeof(GrfnSSParam));
		if (!cells || !vars) { err = serror("Memory failure in CheyBetaMC_Price");  goto FREE_RETURN; }
		for (j=0; j < npaths; j++) CopyFromGD(j, gd, cells, vars);
	}

	// Check that the final event is at the final time step
	if ( fabs(g->ex[g->nex-1] - times[ntimes-1]) > 1e-5 )
	{ err = serror("Last event and last time step do not coincide");  goto FREE_RETURN; }

	// Initialize random numbers generator
	if (rgmode == 0) err = RandomGen_Init(&rg, seed);
	else if (rgmode == 1) err = ABS_Init(&rg, seed, npaths, 1, 0);
	else if (rgmode == 2) err = ABS_Init(&rg, seed, npaths, 1, 1);
	if (err) goto FREE_RETURN;

	dT = (long)(pmdl->today + times[ntimes-1] * DAYS_IN_YEAR + 1e-5);
	lam = pmdl->lambda;

	d1 = pmdl->today;
	gamT = ( 1.0 - exp(-lam * times[ntimes-1]) ) / lam;

	// Generate all paths simultaneously
	for (i=0; i < ntimes-1; i++)
	{
		dt = times[i+1] - times[i];
		sqrtdt = sqrt(dt);
		while ( isig < pmdl->nsig - 1 && pmdl->sigtms[isig] < times[i] + 1e-5 ) isig++;
		sig = pmdl->sig[isig];
		beta = pmdl->beta[isig];

		d2 = (long)(pmdl->today + times[i+1] * DAYS_IN_YEAR + 1e-5);
		fr = swp_f_zr(d1, d2, pmdl->ycname);

		std = sqrt(cumxvar);
		xcut_u = cutcoef * fr * std;
		xcut_d = -cutcoef * fr * std;

		std = CheyBeta_Vol(0.0, phi_mid, sig, beta, gamT, fr, xcut_u, xcut_d);

		cumxvar += std * std * dt;
		phi_mid += (std * std - 2.0 * lam * phi_mid) * dt;

		// Move to the next step on all paths
		for (j=0; j < npaths; j++)
		{
			err = rg.Gauss(&rg, &dw);
			if (err) goto FREE_RETURN;
			dw *= sqrtdt;

			std = CheyBeta_Vol(x[j], phi[j], sig, beta, gamT, fr, xcut_u, xcut_d);
			x[j] += -lam * x[j] * dt + std * dw;
			phi[j] += (std * std - 2.0 * lam * phi[j]) * dt;
		}

		d1 = d2;
		gamT = ( 1.0 - exp(-lam * (times[ntimes-1] - times[i+1])) ) / lam;

		while ( iex < g->nex - 1 && g->ex[iex+1] < times[i+1] + 1e-5 ) iex++;

		// Evaluate events if any
		if ( iex >= 0 && ( (fabs(g->ex[iex] - times[i+1]) < 1e-5) || (g->am && g->am[iex]) ) )
		{
			dffT = swp_f_df(d1, dT, pmdl->ycname);

			// Precalculate dffs and gams
			for (k=0; k < g->nmat[0][iex]; k++)
			{
				d2 = (long)(pmdl->today + g->mat[0][iex][k] * DAYS_IN_YEAR + 1e-5);
				dffs[k] = swp_f_df(d1, d2, pmdl->ycname);
				gams[k] = ( 1.0 - exp(-lam * (g->mat[0][iex][k] - times[i+1])) ) / lam;
			}

			// Calculate payoff for all paths
			for (j=0; j < npaths; j++)
			{
				xj = x[j] - phi[j] * gamT;
				num = dffT * exp( -gamT * (xj + 0.5 * gamT * phi[j]) );

				for (k=0; k < g->nmat[0][iex]; k++)
					dfs[k] = dffs[k] * exp( -gams[k] * (xj + 0.5 * gams[k] * phi[j]) );

				for (k=0; k < g->ninst; k++) pvs[j][k] *= num;
				if (g->type == PRODUCT_GRFN) CopyToGD(j, gd, cells, vars);

				err = g->Payoff(g, iex, times[i+1], &mkt_data, pvs[j]);
				if (err) goto FREE_RETURN;

				if (g->type == PRODUCT_GRFN) CopyFromGD(j, gd, cells, vars);
				for (k=0; k < g->ninst; k++) pvs[j][k] /= num;
			}
		}
	}

	// Calculate pv and stddev
	memset(pv, 0, g->ninst * sizeof(double));
	memset(stddev, 0, g->ninst * sizeof(double));
	dffT = swp_f_df(pmdl->today, dT, pmdl->ycname);

	for (k=0; k < g->ninst; k++)
	{
		for (j=0; j < npaths; j++) pv[k] += pvs[j][k] * dffT / npaths;
		for (j=0; j < npaths; j++)
		{
			num = pvs[j][k] * dffT - pv[k];
			stddev[k] += num * num / (npaths-1);
		}
		stddev[k] = sqrt(stddev[k] / npaths);
	}

FREE_RETURN:
	free(times);
	free(x);
	free(phi);
	if (pvs) free_dmatrix(pvs, 0, npaths-1, 0, g->ninst-1);
	if (cells) free_f3tensor(cells, 0, npaths-1, 0, gd->sslength-1, 0, gd->sswidth-1);
	free(vars);
	if (rgmode > 0) ABS_Free(&rg);

	return err;
}
