/*--------------------------------------------------------------
	FILE: srt_f_range_accrual.c
	PURPOSE: Range accrual product implementation for GRFN
	AUTHOR: Dimitri Mayevski
	DATE: 12/06/2002
  --------------------------------------------------------------*/

#include "srt_h_all.h"
#include "math.h"
#include "OPFNCTNS.H"
#include "swp_h_vol.h"
#include "srt_h_range_accrual.h"

#define K_TINY 1.0e-06
#define MAXCPN 1000

static Err SrtRA_Destroy(SrtProduct *product)
{
	SrtRangeAccrualSpec *RA;
	SrtRAObservStep *stp, *tmp;
	SrtVolBumpData *bump_data, *bump_tmp;

	if (product)
	{
		RA = (SrtRangeAccrualSpec *)(product->spec_desc);
		if (RA)
		{
			if (RA->dates) free(RA->dates);
			if (RA->dates_f) free(RA->dates_f);
			if (RA->cvg_f) free(RA->cvg_f);
			if (RA->idx) free(RA->idx);
			if (RA->obs_steps)
			{
				stp = tmp = RA->obs_steps[0];
				while (stp)
				{
					stp = stp->next;
					free(tmp);
					tmp = stp;
				}
				free(RA->obs_steps);
			}
			bump_data = bump_tmp = RA->bumps;
			while (bump_data)
			{
				bump_data = bump_data->next;
				free(bump_tmp);
				bump_tmp = bump_data;
			}
			free(RA);
		}
		free(product);
	}
	return NULL;
}

static Err SrtRA_RequestDfDates(SrtProduct *product, long date,
								int *pn_unds, int **pn_dfs, long ***pdates)
{
	SrtRangeAccrualSpec *RA = (SrtRangeAccrualSpec *)(product->spec_desc);
	int i, j, ndfs_d, ndfs_f;
	long d, d1, dates_d[MAXCPN], *dates_f;
	long idx[MAXCPN];

	*pn_unds = (RA->quanto ? 2 : 1);
	*pn_dfs = (int *) calloc(*pn_unds, sizeof(int));
	*pdates = (long **) calloc(*pn_unds, sizeof(long*));

	for ( i=0; i < RA->n_periods &&
		  date > add_unit( RA->dates[i], -RA->spotlag_f,
		  SRT_BDAY, MODIFIED_SUCCEEDING );
		  i++ );

	if (i == RA->n_periods)
	{
		for (j=0; j < *pn_unds; j++)
		{
			(*pn_dfs)[j] = 0;
			(*pdates)[j] = NULL;
		}
		return NULL;
	}

	ndfs_d = RA->n_periods - i;
	dates_f = dates_d + ndfs_d;

	dates_f[j=0] = d = RA->dates[i];
	d1 = d + 15 * RA->tenor_f;
	while (1)
	{
		dates_f[++j] = bus_date_method(d1, MODIFIED_SUCCEEDING);
		d = add_unit(d, RA->tenor_f, SRT_MONTH, NO_BUSDAY_CONVENTION);
		d1 = add_unit(d1, RA->tenor_f, SRT_MONTH, NO_BUSDAY_CONVENTION);
		dates_f[++j] = bus_date_method(d, MODIFIED_SUCCEEDING);
		if (dates_f[j-2] >= RA->dates[RA->n_periods]) break;
	}
	ndfs_f = j+1;

	memcpy(dates_d, RA->dates + i+1, ndfs_d * sizeof(long));

	if (RA->quanto)
	{
		(*pn_dfs)[0] = ndfs_d;
		(*pn_dfs)[1] = ndfs_f;
		(*pdates)[0] = (long *) calloc(ndfs_d, sizeof(long));
		(*pdates)[1] = (long *) calloc(ndfs_f, sizeof(long));
		memcpy((*pdates)[0], dates_d, ndfs_d * sizeof(long));
		memcpy((*pdates)[1], dates_f, ndfs_f * sizeof(long));
	}
	else
	{
		(*pn_dfs)[0] = ndfs_d + ndfs_f;
		(*pdates)[0] = (long *) calloc(ndfs_d + ndfs_f, sizeof(long));
		indexx_ll(dates_d, idx, ndfs_d + ndfs_f);
		for (j=0; j < ndfs_d + ndfs_f; j++)	(*pdates)[0][j] = dates_d[idx[j]];
	}
	return NULL;
}

static Err SrtRA_Payoff(SrtProduct *product, long today, long date,
						double **dfs, double *payoff)
{
	Err err = NULL;
	SrtRangeAccrualSpec *RA = (SrtRangeAccrualSpec *)(product->spec_desc);
	int i, j, k, nfwds;
	long d, d1;
	double dfs_d[MAXCPN], *dfs_f;
	double fwds[MAXCPN];
	long dates_d[MAXCPN], *dates_f;
	double coef, fwd, fixing, mat, df, sum, betavol, vol, bs[4], bump_dir[4];
	SrtRAObservStep *stp;
	SrtVolBumpData *bump_tmp;

	*payoff = 0.0;

	/* If stored deterministic values are no more valid recalculate them */
	if (date != RA->last_date)
	{
		for ( i=0; i < RA->n_periods &&
			  date > add_unit( RA->dates[i], -RA->spotlag_f,
			  SRT_BDAY, MODIFIED_SUCCEEDING );
			  i++ );
		RA->last_date = date;
		RA->period_idx = i;
		if (i == RA->n_periods) return NULL;

		RA->ndfs_d = RA->n_periods - i;
		dates_f = dates_d + RA->ndfs_d;

		dates_f[j=0] = d = RA->dates[i];
		d1 = d + 15 * RA->tenor_f;
		while (1)
		{
			dates_f[++j] = bus_date_method(d1, MODIFIED_SUCCEEDING);
			d = add_unit(d, RA->tenor_f, SRT_MONTH, NO_BUSDAY_CONVENTION);
			d1 = add_unit(d1, RA->tenor_f, SRT_MONTH, NO_BUSDAY_CONVENTION);
			dates_f[++j] = bus_date_method(d, MODIFIED_SUCCEEDING);
			if (dates_f[j-2] >= RA->dates[RA->n_periods]) break;
		}
		RA->ndfs_f = j+1;

		if (!RA->quanto)
		{
			memcpy(dates_d, RA->dates + i+1, RA->ndfs_d * sizeof(long));
			RA->idx = (long *) realloc( RA->idx,
				(RA->ndfs_d + RA->ndfs_f) * sizeof(long) );
			indexx_ll(dates_d, RA->idx, RA->ndfs_d + RA->ndfs_f);
		}
		RA->dates_f = (long *) realloc(RA->dates_f, (RA->ndfs_f-2) * sizeof(long));
		RA->cvg_f = (double *) realloc(RA->cvg_f, (RA->ndfs_f-2) * sizeof(double));

		memcpy(RA->dates_f, dates_f, (RA->ndfs_f-2) * sizeof(long));
		for (j=0; j < RA->ndfs_f-2; j++)
			RA->cvg_f[j] = coverage(dates_f[j], dates_f[j+2], RA->basis_f);

	}
	if (RA->period_idx == RA->n_periods) return NULL;

	dfs_f = dfs_d + RA->ndfs_d;

	if (RA->quanto)
	{
		memcpy(dfs_d, dfs[0], RA->ndfs_d * sizeof(double));
		memcpy(dfs_f, dfs[1], RA->ndfs_f * sizeof(double));
	}
	else for (j=0; j < RA->ndfs_d + RA->ndfs_f; j++) dfs_d[RA->idx[j]] = dfs[0][j];

	/* Calculate forwards cash */
	nfwds = RA->ndfs_f - 2;
	for (j=0; j < nfwds; j++) fwds[j] = (dfs_f[j] / dfs_f[j+2] - 1.0) / RA->cvg_f[j];

	/* Find vol bump corresponding to today */
	if (RA->bump_data->date != today)
	{
		RA->bump_data = RA->bumps;
		while (RA->bump_data && RA->bump_data->date != today)
		{
			bump_tmp = RA->bump_data;
			RA->bump_data = RA->bump_data->next;
		}
		if (!RA->bump_data)		/* today not in the list */
		{
			bump_tmp->next = (SrtVolBumpData *) malloc(sizeof(SrtVolBumpData));
			RA->bump_data = bump_tmp->next;
			RA->bump_data->date = today;
			RA->bump_data->bump = 0.0;
			RA->bump_data->next = NULL;
		}
	}

	/* Calculate the sum of call spreads interpolating forwards */
	for (k=0, stp = RA->obs_steps[i = RA->period_idx]; stp->next; stp = stp->next)
	{
		if (stp->date >= RA->dates[i+1]) i++;
		for (; k < nfwds-1 && RA->dates_f[k+1] < stp->date; k++);

		coef = ((double)(stp->date - RA->dates_f[k])) /
			(RA->dates_f[k+1] - RA->dates_f[k]);
		fwd = coef * fwds[k+1] + (1.0-coef) * fwds[k];

		fixing = add_unit(stp->date, -RA->spotlag_f, SRT_BDAY, MODIFIED_SUCCEEDING);
		mat = (fixing - today) * YEARS_IN_DAY;
		df = dfs_d[i - RA->period_idx];		/* df(dates[i+1]) */

		/* Adjust forward */
		fwd += stp->spread;
		fwd *= exp(stp->correction * mat);

		/* Calculate bump_dir */
		if (fwd > stp->K[3]) bump_dir[3] = bump_dir[2] = -1.0;
		else if (fwd < stp->K[2]) bump_dir[3] = bump_dir[2] = 1.0;
		else bump_dir[3] = bump_dir[2] = 0.0;

		if (fwd < stp->K[0]) bump_dir[0] = bump_dir[1] = -1.0;
		else if (fwd > stp->K[1]) bump_dir[0] = bump_dir[1] = 1.0;
		else bump_dir[0] = bump_dir[1] = 0.0;

		/* Calculate puts */
		for (j=0; j < 4; j++)
		{
			if (stp->K[j] < K_TINY) bs[j] = 0.0;
			else if (fwd < K_TINY) bs[j] = df * (stp->K[j] - fwd);
			else
			{
				betavol = stp->betavol + RA->bump_data->bump * bump_dir[j];
				if (betavol < K_TINY) betavol = K_TINY;
				err = srt_f_optsarbvol(fwd, stp->K[j], mat, betavol, stp->alpha,
					stp->beta, stp->rho, SRT_BETAVOL, SRT_LOGNORMAL, &vol);
				if (err) return err;
				bs[j] = srt_f_optblksch( fwd, stp->K[j], vol, mat, df,
										  SRT_PUT, PREMIUM );
/*				bs[j] = srt_f_optblksch( fwd, stp->K[j], stp->vol_atK[j], mat, df,
										  SRT_PUT, PREMIUM );*/
			}
		}
		sum = (bs[3] - bs[2]) / (stp->K[3] - stp->K[2]) -
			  (bs[1] - bs[0]) / (stp->K[1] - stp->K[0]);

		*payoff += sum * stp->cpn;
	}
	return NULL;
}

SrtProduct *srt_f_init_range_accrual(
			char			*name,
			char			*yc_d,
			char			*yc_f,
			char			*volc_d,
			char			*volc_f,
			char			*refrate_d,
			char			*refrate_f,
			int				n_fxvol_dates,
			long			*fxvol_dates,
			double			*fxvol,
			double			*qtocorrel,
			int				n_periods,
			long			*dates,		/* n_periods + 1 */
			double			*cpns,
			char			*basis,
			double			*upper_barr,
			double			*lower_barr,
			char			*recpay,
			double			c_spread,
			int				obs_freq,
			double			rho_df )
{
	Err err = NULL;
	SrtProduct *product = NULL;
	SrtRangeAccrualSpec *RA;
	SrtRAObservStep *stp;
	int i, j;
	long d, d1, d2;
	SrtBasisCode basis_d, basis_f, ibasis;
	SrtCompounding freq_d, freq_f;
	double vol_fx, rho_ffx, coef, fra_f, vol_f, fra_d, cvg_d, vol_d, power;
	double quanto_corr, DRS_corr;
	SrtReceiverType irecpay;
	SrtCurvePtr ycd_ptr;

	ycd_ptr = lookup_curve(yc_d);
	if (!ycd_ptr) return NULL;

	product = (SrtProduct *) malloc(sizeof(SrtProduct));
	if (!product) return NULL;

	strcpy(product->product_name, name);
	strupper(product->product_name);
	strip_white_space(product->product_name);

	strcpy(product->product_lbl, "Range Accrual");
	product->Destroy = SrtRA_Destroy;
	product->RequestDfDates = SrtRA_RequestDfDates;
	product->Payoff = SrtRA_Payoff;

	product->spec_desc = malloc(sizeof(SrtRangeAccrualSpec));
	RA = (SrtRangeAccrualSpec *)(product->spec_desc);
	if (!RA) goto RETURN_ERR;
	memset(RA, 0, sizeof(SrtRangeAccrualSpec));

	/* Allocate memory and copy data to RA */

	RA->n_periods = n_periods;
	RA->dates = (long *) calloc(n_periods + 1, sizeof(long));
	if (!RA->dates) goto RETURN_ERR;

	RA->obs_steps =
		(SrtRAObservStep **) calloc(n_periods + 1, sizeof(SrtRAObservStep *));

 	if (!RA->obs_steps) goto RETURN_ERR;
	RA->obs_steps[0] = NULL;

	memcpy(RA->dates, dates, (n_periods + 1) * sizeof(long));
	RA->quanto = (fxvol != NULL);

	err = srt_f_get_spot_lag_from_refrate(refrate_f, &RA->spotlag_f);
	if (err) goto RETURN_ERR;

	err = swp_f_get_ref_rate_details(refrate_f, &basis_f, &freq_f);
	if (err) goto RETURN_ERR;
	RA->tenor_f = 12 / freq_f;
	RA->basis_f = basis_f;

	err = swp_f_get_ref_rate_details(refrate_d, &basis_d, &freq_d);
	if (err) goto RETURN_ERR;

	err = interp_basis(basis, &ibasis);
	if (err) goto RETURN_ERR;

	err = interp_rec_pay(recpay, &irecpay);
	if (err) goto RETURN_ERR;

	/* Initialize observation schedule */

	stp = RA->obs_steps[0] = (SrtRAObservStep *) malloc(sizeof(SrtRAObservStep));
	if (!stp) goto RETURN_ERR;
	stp->next = NULL;
	i = 0;
	stp->date = d = dates[0];

	while (i < n_periods)
	{
		stp->next = (SrtRAObservStep *) malloc(sizeof(SrtRAObservStep));
		stp = stp->next;
		if (!stp) goto RETURN_ERR;
		stp->next = NULL;

		d = add_unit(d, obs_freq, SRT_DAY, SUCCEEDING);
		if (d >= dates[i+1])
		{
			i++;
			RA->obs_steps[i] = stp;
			d = dates[i];
		}
		stp->date = d;
	}

	/* Fill in observation steps */

	j = 0;
	for (i=0, stp = RA->obs_steps[0]; stp->next; stp = stp->next)
	{
		if (stp->date >= dates[i+1]) i++;
		d = add_unit(stp->date, RA->tenor_f, SRT_MONTH, MODIFIED_SUCCEEDING);
		stp->spread = swp_f_spread(stp->date, d, refrate_f);
		stp->cpn = cpns[i] * (stp->next->date - stp->date) / (dates[i+1] - dates[i]) *
			coverage(dates[i], dates[i+1], ibasis);

		fra_f = swp_f_fra(stp->date, d, basis_f, yc_f, refrate_f);
		err = swp_f_SABRvol(volc_f, stp->date, d, fra_f, &vol_f, &power, SABR_LOGVOL);
		if (err) goto RETURN_ERR;

		if (RA->quanto)		/* Calculate quanto adjustment */
		{
			/* Interpolate FX vol and rho_ffx at d */
			for (; j < n_fxvol_dates-1 && fxvol_dates[j+1] < d; j++);
			if (fxvol_dates[j] > d || j == n_fxvol_dates-1)
			{
				vol_fx = fxvol[j];
				rho_ffx = qtocorrel[j];
			}
			else
			{
				coef = ((double)(d - fxvol_dates[j])) /
					(fxvol_dates[j+1] - fxvol_dates[j]);
				vol_fx = coef * fxvol[j+1] + (1.0-coef) * fxvol[j];
				rho_ffx = coef * qtocorrel[j+1] + (1.0-coef) * qtocorrel[j];
			}
			quanto_corr = - rho_ffx * vol_fx * vol_f;
		}
		else quanto_corr = 0.0;

		/* Calculate DRS adjustment */
		if (dates[i+1] == d) DRS_corr = 0;
		else
		{
			d1 = (dates[i+1] < d ? dates[i+1] : d);
			d2 = (dates[i+1] > d ? dates[i+1] : d);

			fra_d = swp_f_fra(d1, d2, basis_d, yc_d, refrate_d);
			cvg_d = coverage(d1, d2, basis_d);
			err = swp_f_SABRvol(volc_d, d1, d2, fra_d, &vol_d, &power, SABR_LOGVOL);
			if (err) goto RETURN_ERR;

			DRS_corr = rho_df * vol_d * vol_f * cvg_d * fra_d / (1.0 + cvg_d * fra_d);
			if (dates[i+1] > d) DRS_corr = -DRS_corr;
		}
		stp->correction = quanto_corr + DRS_corr;

		/* Calculate vols at barriers */
		if (irecpay == SRT_RECEIVER)
		{
			stp->K[0] = lower_barr[i];
			stp->K[1] = lower_barr[i] + c_spread;
			stp->K[2] = upper_barr[i] - c_spread;
			stp->K[3] = upper_barr[i];
		}
		else
		{
			stp->K[0] = lower_barr[i] - c_spread;
			stp->K[1] = lower_barr[i];
			stp->K[2] = upper_barr[i];
			stp->K[3] = upper_barr[i] + c_spread;
		}
/*		for (k=0; k < 4; k++) if (stp->K[k] >= K_TINY)
		{
			err = swp_f_SABRvol( volc_f, stp->date, d, stp->K[k],
				&stp->vol_atK[k], &power, SABR_LOGVOL );
			if (err) goto RETURN_ERR;
		}*/
		/* Get SABR parameters for the step */
		err = swp_f_SABRvol( volc_f, stp->date, d, 0.0,
				&stp->betavol, &power, SABR_BETAVOL );
		if (err) goto RETURN_ERR;
		err = swp_f_SABRvol( volc_f, stp->date, d, 0.0,
				&stp->alpha, &power, SABR_ALPHA );
		if (err) goto RETURN_ERR;
		err = swp_f_SABRvol( volc_f, stp->date, d, 0.0,
				&stp->beta, &power, SABR_BETA );
		if (err) goto RETURN_ERR;
		err = swp_f_SABRvol( volc_f, stp->date, d, 0.0,
				&stp->rho, &power, SABR_RHO );
		if (err) goto RETURN_ERR;
	}

	/* Add today with zero bump to the vol bump list */
	d = get_today_from_curve(ycd_ptr);
	RA->bumps = (SrtVolBumpData *) malloc(sizeof(SrtVolBumpData));
	if (!RA->bumps) goto RETURN_ERR;
	RA->bumps->date = d;
	RA->bumps->bump = 0.0;
	RA->bumps->next = NULL;
	RA->bump_data = RA->bumps;

	RA->last_date = -1;

	err = add_product_to_list(product);
	if (err) goto RETURN_ERR;

	return product;

RETURN_ERR:
	SrtRA_Destroy(product);
	return NULL;
}
