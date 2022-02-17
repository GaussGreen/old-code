/*--------------------------------------------------------------
	FILE: srt_f_cheybeta_diagcalib.c
	PURPOSE: Cheyette beta diag calib (new)
	AUTHOR: Dimitri Mayevski
	DATE: 10/10/2002
  --------------------------------------------------------------*/

#include "srt_h_all.h"
#include "swp_h_vol.h"
#include "opfnctns.h"
#include "num_h_zbrent.h"
#include "srt_h_cheybeta_new.h"
#include "srt_h_cheybeta_diagcalib.h"
#include "swp_h_utils.h"
#include "Product.h"
#include "CheyBetaPDE.h"
#include "CPDCalib.h"
#include "math.h"

typedef struct _CheyBeta_TrySigmaData {
	SCheyBetaPDE		*pde;
	double				*mkt_prices;
	int					i;
} CheyBeta_TrySigmaData;

static Err try_sigma(double sig, double *diff, void *static_data)
{
	Err		err = NULL;
	CheyBeta_TrySigmaData *p = (CheyBeta_TrySigmaData *) static_data;
	SCheyBetaPDE		*pde = p->pde;
	SCheyBeta			*pmdl = pde->pmdl;
	int					i = p->i;
	double				result;

	err = CheyBetaPDE_ResetFwd(pde, i ? pmdl->sigtms[i-1] : 0.0);
	if (err) return err;

	pmdl->sig[i] = sig;
	err = CheyBetaPDE_PriceSwaption(pde, i, &result);
	if (err) return err;

	*diff = result - p->mkt_prices[i];
	return NULL;
}

Err cheybeta_diagcalib_new(
	SCheyBeta		*pmdl,				/*	beta, lambda, ycname and today must be
											initialized inside. Sigmas are returned.
											To be freed by caller */
	char			*vc_name,			/*	Name of the market vol curve */
	double			vol_shift,
	int				shift_type,			/*	0:	Additive
											1:	Multiplicative */
	int				nex,
	long			*ex_dates,			/*	Exercise dates */
	long			end_date,			/*	End date for diagonal */
	char			**tenors,			/*	Instruments' tenors */
	char			*freq,
	char			*basis,
	char			*refrate,			/*	Name of the reference rate */
	double			*strikes,
	int				strike_type,		/*	0: ATM
											1: CASH
											2: SWAP
											3: STD */
	double			max_std,
	int				nt,					/*  Parameters of the		*/
	int				nx,					/*  CheyBeta forward PDE	*/
	int				nphi,				/*  pricing routine			*/
	double			cutcoef,
	Err				(*get_cash_vol)(				/*	Function to get cash vol (for LGM calib) */
						char	*vol_curve_name,	
						double	start_date, 
						double	end_date,
						double	cash_strike,
						int		zero,
						char	*ref_rate_name,
						double	*vol,
						double	*power)
)
{
	Err				err = NULL;
	SCashFlows		g;
	SProductDesc	product;
	SCheyBetaPDE	pde;
	double			*mkt_prices = NULL;
	SrtCompounding	reffreq;
	SrtBasisCode	refbasis;
	char			*sreffreq, *srefbasis;
	long			start, end;
	int				i, spot_lag, elements;
	double			fwd, lvl, strike, std, power;
	SrtCurvePtr		yc_ptr;
	SwapDP			sdp;
	double			x1, x2, f1, f2;
	CheyBeta_TrySigmaData data;

	/* ____For preliminary LGM calibration (interface with cpd_calib_diagonal_2)____ */
	char			**tenors_LGM = NULL;
	double			*strikes_LGM = NULL;
	int				*arr_of_1 = NULL;
	int				nsig_LGM;
	double			*sigtms_LGM = NULL, *sig_LGM = NULL;
	cpd_diag_calib_param param_LGM;
	/* _____________________________________________________________________________ */

	memset(&g, 0, sizeof(SCashFlows));
	memset(&pde, 0, sizeof(SCheyBetaPDE));

	yc_ptr = lookup_curve(pmdl->ycname);
	if (!yc_ptr) return serror("Yield Curve not found");
	spot_lag = get_spotlag_from_curve(yc_ptr);
	elements = SE_NOTIONALS | SE_FIXED | SE_SPREADS;

	err = swp_f_get_ref_rate_details(refrate, &refbasis, &reffreq);
	if (err) return err;

	err = translate_basis(&srefbasis, refbasis);
	if (err) return err;

	err = translate_compounding(&sreffreq, reffreq);
	if (err) return err;

	pmdl->nsig = nex;
	pmdl->sig = (double *) calloc(nex, sizeof(double));
	pmdl->sigtms = (double *) calloc(nex, sizeof(double));
	if (!pmdl->sig || !pmdl->sigtms) return serror("Memory failure");

	err = CashFlows_Init(&g, nex);
	if (err) goto FREE_RETURN;

	mkt_prices = (double *) calloc(nex, sizeof(double));
	if (!mkt_prices) { err = serror("Memory failure");  goto FREE_RETURN; }

	/* Initialize cash flows and market prices */
	for (i=0; i < nex; i++)
	{
		pmdl->sigtms[i] = g.ex[i] = (ex_dates[i] - pmdl->today) * YEARS_IN_DAY;
		start = add_unit(ex_dates[i], spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);
		if (end_date != 0) end = end_date;
		else if (tenors)
		{
			err = add_tenor(start, tenors[i], NO_BUSDAY_CONVENTION, &end);
			if (err) goto FREE_RETURN;
		}
		else
		{
			end = add_unit(start, 12 / reffreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
			freq = sreffreq;
			basis = srefbasis;
		}

		err = swp_f_initSwapDP(start, end, freq, basis, &sdp);
		if (err) goto FREE_RETURN;

		err = swp_f_ForwardRate_SwapDP(&sdp, pmdl->ycname, refrate, &fwd);
		if (err) goto FREE_RETURN;

		err = swp_f_Level_SwapDP(&sdp, pmdl->ycname, &lvl);
		if (err) goto FREE_RETURN;

		/* ATM vol */
		err = swp_f_vol(vc_name, start, end, fwd, &std, &power);
		if (err) goto FREE_RETURN;

		std += (shift_type == 1 ? std * vol_shift : vol_shift);
		if (power > 0.5) std *= fwd;
		std *= sqrt(g.ex[i]);

		/* Strike */
		if (!strikes || strike_type == 0) strike = fwd;
		else if (strike_type == 2) strike = strikes[i];
		else if (strike_type == 3) strike = fwd + strikes[i] * std;
		else { err = serror("Unknown strike type");  goto FREE_RETURN; }

		/* Apply max std */
		if (strike > fwd + max_std * std) strike = fwd + max_std * std;
		if (strike < fwd - max_std * std) strike = fwd - max_std * std;

		/* Vol at strike */
		err = swp_f_vol(vc_name, start, end, strike, &std, &power);
		if (err) goto FREE_RETURN;
		std += (shift_type == 1 ? std * vol_shift : vol_shift);

		/* Market price */
		mkt_prices[i] = (power > 0.5 ?
			srt_f_optblksch(fwd, strike, std, g.ex[i], lvl, SRT_PUT, PREMIUM) :
			srt_f_optblknrm(fwd, strike, std, g.ex[i], lvl, SRT_PUT, PREMIUM) );

		/* Cash flows */
		err = SwapElements_Init(start, end, freq, basis, refrate, "REC", strike,
			elements, pmdl->today, &(g.mat[i]), &(g.cpn[i]), &(g.nmat[i]));
		if (err) goto FREE_RETURN;

	}	/* for (i=0; i < nex; i++) */

	ProductDesc_InitSwaptions(&product, &g);

	/* _______Calibrate LGM for later use as a first guess_______________________________________ */

	param_LGM.vol_shift = vol_shift;
	param_LGM.shift_type = shift_type;
	param_LGM.strike_type = strike_type;
	param_LGM.max_std = max_std;
	param_LGM.min_time = 0.0;
	param_LGM.skip_last = 0;
	param_LGM.fx_vol_shift = 0.0;

	tenors_LGM = svector_size(0, nex-1, 20);
	strikes_LGM = (double *) calloc(nex, sizeof(double));
	arr_of_1 = (int *) calloc(nex, sizeof(int));

	if (!tenors_LGM || !strikes_LGM || !arr_of_1)
	{ err = serror("Memory failure");  goto FREE_RETURN; }

	if (tenors) for (i=0; i < nex; i++) strcpy(tenors_LGM[i], tenors[i]);
	else if (end_date != 0) for (i=0; i < nex; i++) strcpy(tenors_LGM[i], "DIAG");
	else for (i=0; i < nex; i++) sprintf(tenors_LGM[i], "%dM", 12 / reffreq);

	if (strikes) memcpy(strikes_LGM, strikes, nex * sizeof(double));
	else memset(strikes_LGM, 0, nex * sizeof(double));

	for (i=0; i < nex; i++) arr_of_1[i] = 1;

	err = cpd_calib_diagonal_2( pmdl->ycname, vc_name, refrate, get_cash_vol, freq, basis,
		nex, ex_dates, arr_of_1, tenors_LGM, end_date, strikes_LGM, pmdl->lambda, 1, 0.0, 0.0, 0.0,
		&nsig_LGM, &sigtms_LGM, &sig_LGM, &param_LGM, NULL );
	if (err) goto FREE_RETURN;

	if (nsig_LGM != nex) { err = serror("Wrong nsig after LGM calibration"); goto FREE_RETURN; }

	for (i=0; i < nex; i++)
	{
		if ( fabs(pmdl->sigtms[i] - sigtms_LGM[i]) > 1e-5 )
		{ err = serror("Wrong dates in LGM calibration"); goto FREE_RETURN; }

		pmdl->sig[i] = sig_LGM[i];
	}
	/* ________________________________________________________________________________________ */

	err = CheyBetaPDE_Init(&pde, &product, nt, nx, nphi, cutcoef, 1, pmdl);
	if (err) goto FREE_RETURN;

	data.pde = &pde;
	data.mkt_prices = mkt_prices;

	for (i=0; i < nex; i++)
	{
		data.i = i;
		x1 = pmdl->sig[i] * 0.9;
		x2 = pmdl->sig[i] * 1.1;

		err = num_f_zbrac(try_sigma, &x1, &x2, &f1, &f2, &data);
		if (err) goto FREE_RETURN;

		err = num_f_zbrent(try_sigma, x1, x2, f1, f2, 1e-6, &data, &(pmdl->sig[i]));
		if (err) goto FREE_RETURN;

		err = try_sigma(pmdl->sig[i], &f1, &data);
		if (err) goto FREE_RETURN;
	}

FREE_RETURN:

	CheyBetaPDE_Free(&pde);
	CashFlows_Free(&g);
	free(mkt_prices);

	if (tenors_LGM) free_svector_size(tenors_LGM, 0, nex-1, 20);
	free(strikes_LGM);
	free(arr_of_1);
	free(sigtms_LGM);
	free(sig_LGM);

	return err;
}
