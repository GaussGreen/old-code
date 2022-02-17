/*--------------------------------------------------------------
	FILE: srt_f_ra_calib.c
	PURPOSE: Range Accrual product SABR sigma-beta calibration
	AUTHOR: Dimitri Mayevski
	DATE: 01/08/2002
  --------------------------------------------------------------*/

#include "srt_h_all.h"
#include "math.h"
#include "num_h_zbrent.h"
#include "srt_h_range_accrual.h"
#include "srt_h_ra_calib.h"

typedef struct {
	long		d;
	double		phi_f, R_df, R_ffx, df;
	int			ndfs_d, ndfs_f;
	double		*dffs_d, *dffs_f;
	double		*gam_d, *gam_f;
	double		target;
} CalibRA_ExData;

typedef struct {
	SrtProduct		*product;
	CalibRA_ExData	*ex_data;
	int				nx;
	double			*x, *w;
	double			std;
} TryBumpData;

#define MAXNX 100
#define MAXCPN 1000

static Err try_bump(double bump, double *diff, void *static_data)
{
	Err			err = NULL;
	TryBumpData *p = (TryBumpData *) static_data;
	SrtRangeAccrualSpec *RA = (SrtRangeAccrualSpec *)(p->product->spec_desc);
	double		dfs_d[MAXCPN];
	double		*dfs[1] = {dfs_d};
	int			i, j;
	double		sum, c;

	RA->bump_data->bump = bump;

	for (i=0, sum=0.0; i < p->nx; i++)
	{
		/* Calculate all domestic dfs */
		for (j=0; j < p->ex_data->ndfs_d; j++)
		{
			c = p->ex_data->gam_d[j];
			dfs_d[j] = p->ex_data->dffs_d[j] *
				exp( - c * (p->x[i] + 0.5 * c * p->ex_data->phi_f) );
		}
		/* Calculate payoff */
		err = p->product->Payoff(p->product, p->ex_data->d, p->ex_data->d, dfs, &c);
		if (err) return err;

		sum += p->w[i] * c * INV_SQRT_TWO_PI *
			exp( -0.5 * p->x[i] * p->x[i] / p->ex_data->phi_f ) / p->std;
	}

	*diff = sum * p->ex_data->df - p->ex_data->target;
	return NULL;
}

static Err try_bump_quanto(double bump, double *diff, void *static_data)
{
	Err			err = NULL;
	TryBumpData *p = (TryBumpData *) static_data;
	SrtRangeAccrualSpec *RA = (SrtRangeAccrualSpec *)(p->product->spec_desc);
	double		dfs_d[MAXCPN], dfs_f[MAXCPN];
	double		*dfs[2] = {dfs_d, dfs_f};
	int			i, j;
	double		sum, c;

	RA->bump_data->bump = bump;

	for (i=0, sum=0.0; i < p->nx; i++)
	{
		/* Calculate integrated domestic dfs times the density */
		for (j=0; j < p->ex_data->ndfs_d; j++)
		{
			c = p->x[i] - p->ex_data->R_ffx + p->ex_data->gam_d[j] * p->ex_data->R_df;
			dfs_d[j] = p->ex_data->dffs_d[j] * INV_SQRT_TWO_PI *
				exp( - 0.5*c*c / p->ex_data->phi_f ) / p->std;
		}
		/* Calculate foreign dfs */
		for (j=0; j < p->ex_data->ndfs_f; j++)
		{
			c = p->ex_data->gam_f[j];
			dfs_f[j] = p->ex_data->dffs_f[j] *
				exp( - c * (p->x[i] + 0.5 * c * p->ex_data->phi_f) );
		}
		/* Calculate payoff already multiplied by the density */
		err = p->product->Payoff(p->product, p->ex_data->d, p->ex_data->d, dfs, &c);
		if (err) return err;

		sum += p->w[i] * c;
	}

	*diff = sum * p->ex_data->df - p->ex_data->target;
	return NULL;
}

static Err calib_bump(
		SrtProduct		*product,
		CalibRA_ExData	*ex_data,
		int				nx )
{
	Err			err = NULL;
	TryBumpData data;
	double		x[MAXNX], w[MAXNX];
	double		x1, x2, f1, f2, result;
	double		tol = 1.0e-5;
	SrtRangeAccrualSpec *RA = (SrtRangeAccrualSpec *)(product->spec_desc);
	SrtVolBumpData *bump_tmp;
	int			quanto = (ex_data->dffs_f != NULL);

	/* Calculate Gauss-Legendre abscissas */
	if (nx > MAXNX) return serror("nx should be no more than %d", MAXNX);
	data.std = sqrt(ex_data->phi_f);
	gauleg( ex_data->R_ffx - 4.0 * data.std,
			ex_data->R_ffx + 4.0 * data.std, x-1, w-1, nx );

	/* Fill in static data */
	data.product = product;
	data.ex_data = ex_data;
	data.nx = nx;
	data.x = x;
	data.w = w;

	/* Find bump structure for the current exercise date */
	if (RA->bump_data->date != ex_data->d)
	{
		RA->bump_data = RA->bumps;
		while (RA->bump_data && RA->bump_data->date != ex_data->d)
		{
			bump_tmp = RA->bump_data;
			RA->bump_data = RA->bump_data->next;
		}
		if (!RA->bump_data)		/* ex date not in the list */
		{
			bump_tmp->next = (SrtVolBumpData *) malloc(sizeof(SrtVolBumpData));
			RA->bump_data = bump_tmp->next;
			RA->bump_data->date = ex_data->d;
			RA->bump_data->bump = 0.0;
			RA->bump_data->next = NULL;
		}
	}

	x1 = 0.0;
	x2 = 0.01;

	/* Bracket the root */
	err = num_f_zbrac( (quanto ? try_bump_quanto : try_bump),
		&x1, &x2, &f1, &f2, &data );
	if (err) return err;

	/* Calibrate */
	err = num_f_zbrent( (quanto ? try_bump_quanto : try_bump),
		x1, x2, f1, f2, tol, &data, &result );
	if (err) return err;

	RA->bump_data->bump = result;
	return NULL;
}

Err srt_f_calib_range_accrual(
		char			*yc_d,
		char			*yc_f,
		SrtProduct		*product,
		int				nex,
		long			*ex_dates,
		double			*sig_d,
		double			*sig_f,
		double			*sig_fx,
		double			lam_d,
		double			lam_f,
		double			rho_df,
		double			rho_ffx,
		int				nx )
{
	Err			err = NULL;
	int			i, j, k;
	CalibRA_ExData	*ex_data = NULL;
	SrtCurvePtr	yc_ptr;
	long		today;
	long		**dates = NULL;
	int			n_unds, *n_dfs = NULL;
	double		t1, t2;
	double		a_SdSf_LdLf, a_2Sf_2Lf, a_SdSf_Lf, a_2Sf_Lf, a_SfSx_Lf, a_2Sd_2Ld;
	double		exp_Ld, exp_Lf;
	double		*dfs[2] = { NULL, NULL };
	double		*gam[2] = { NULL, NULL };
	char		*yc[2] = { yc_d, yc_f };
	double		lam[2] = { lam_d, lam_f };

	yc_ptr = lookup_curve(yc_d);
	if (!yc_ptr) return serror("Yield curve %s not found", yc_d);
	today = get_today_from_curve(yc_ptr);

	ex_data = (CalibRA_ExData *) calloc(nex, sizeof(CalibRA_ExData));
	if (!ex_data) { err = serror("Memory failure");  goto FREE_RETURN; }
	memset(ex_data, 0, nex * sizeof(CalibRA_ExData));

	/* Fill in ex_data */
	a_SdSf_LdLf = a_2Sf_2Lf = a_SdSf_Lf = a_2Sf_Lf = a_SfSx_Lf = a_2Sd_2Ld = 0.0;
	for (i=0, t1=0.0; i < nex; i++)
	{
		ex_data[i].d = ex_dates[i];
		t2 = (ex_dates[i] - today) * YEARS_IN_DAY;
		exp_Ld = exp(-lam_d * (t2-t1));

		if (yc_f)
		{
			exp_Lf = exp(-lam_f * (t2-t1));

			a_SdSf_LdLf = a_SdSf_LdLf * exp_Ld * exp_Lf +
				sig_d[i] * sig_f[i] * (1.0 - exp_Ld * exp_Lf) / (lam_d + lam_f);

			a_2Sf_2Lf = a_2Sf_2Lf * exp_Lf * exp_Lf +
				sig_f[i] * sig_f[i] * (1.0 - exp_Lf * exp_Lf) / (2.0 * lam_f);

			a_SdSf_Lf = a_SdSf_Lf * exp_Lf +
				sig_d[i] * sig_f[i] * (1.0 - exp_Lf) / lam_f;

			a_2Sf_Lf = a_2Sf_Lf * exp_Lf +
				sig_f[i] * sig_f[i] * (1.0 - exp_Lf) / lam_f;

			a_SfSx_Lf = a_SfSx_Lf * exp_Lf +
				sig_f[i] * sig_fx[i] * (1.0 - exp_Lf) / lam_f;

			ex_data[i].phi_f = a_2Sf_2Lf;
			ex_data[i].R_df = rho_df * a_SdSf_LdLf;
			ex_data[i].R_ffx = - rho_df * (a_SdSf_Lf - a_SdSf_LdLf) / lam_d +
				(a_2Sf_Lf - a_2Sf_2Lf) / lam_f - rho_ffx * a_SfSx_Lf;
		}
		else
		{
			a_2Sd_2Ld = a_2Sd_2Ld * exp_Ld * exp_Ld +
				sig_d[i] * sig_d[i] * (1.0 - exp_Ld * exp_Ld) / (2.0 * lam_d);
			ex_data[i].phi_f = a_2Sd_2Ld;
			ex_data[i].R_df = 0.0;
			ex_data[i].R_ffx = 0.0;
		}

		ex_data[i].df = swp_f_df(today, ex_dates[i], yc_d);

		err = product->RequestDfDates(product, ex_dates[i], &n_unds, &n_dfs, &dates);
		if (err) goto FREE_RETURN;
		if (yc_f && n_unds != 2)
		{ err = serror("Product is not quanto");  goto FREE_RETURN; }
		if (!yc_f && n_unds != 1)
		{ err = serror("Foreign market not specified");  goto FREE_RETURN; }

		ex_data[i].ndfs_d = n_dfs[0];
		dfs[0] = ex_data[i].dffs_d = (double *) calloc(n_dfs[0], sizeof(double));
		gam[0] = ex_data[i].gam_d = (double *) calloc(n_dfs[0], sizeof(double));
		if (yc_f)
		{
			ex_data[i].ndfs_f = n_dfs[1];
			dfs[1] = ex_data[i].dffs_f = (double *) calloc(n_dfs[1], sizeof(double));
			gam[1] = ex_data[i].gam_f = (double *) calloc(n_dfs[1], sizeof(double));
		}

		/* Get DFs needed to calculate target */
		for (j=0; j < n_unds; j++) for (k=0; k < n_dfs[j]; k++)
			dfs[j][k] = swp_f_df(today, dates[j][k], yc[j]);

		/* Calculate target */
		err = product->Payoff(product, today, ex_dates[i], dfs, &(ex_data[i].target));
		if (err) goto FREE_RETURN;

		/* Now replace dfs by forward dfs and calculate gammas */
		for (j=0; j < n_unds; j++) for (k=0; k < n_dfs[j]; k++)
		{
			dfs[j][k] = swp_f_df(ex_dates[i], dates[j][k], yc[j]);
			gam[j][k] = ( 1.0 -
				exp(-lam[j] * (dates[j][k] - ex_dates[i]) * YEARS_IN_DAY) ) / lam[j];
		}

		/* Free dates and n_dfs */
		for (j=0; j < n_unds; j++) free(dates[j]);
		free(dates);	dates = NULL;
		free(n_dfs);	n_dfs = NULL;
		t1 = t2;
	}	/* for (i=0, t1=0.0; i < nex; i++) */

	/* At this point ex_data is fully initialized */
	/* Now calibrate all vol bumps in RA */

	for (i=0; i < nex; i++)
	{
		err = calib_bump(product, &(ex_data[i]), nx);
		if (err) goto FREE_RETURN;
	}

FREE_RETURN:

	if (dates) for (j=0; j < n_unds; j++) free(dates[j]);
	free(dates);
	free(n_dfs);

	if (ex_data) for (i=0; i < nex; i++)
	{
		free(ex_data[i].dffs_d);
		free(ex_data[i].dffs_f);
		free(ex_data[i].gam_d);
		free(ex_data[i].gam_f);
	}
	free(ex_data);

	return err;
}
