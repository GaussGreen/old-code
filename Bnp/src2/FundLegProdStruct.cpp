#include "srt_h_all.h"
#include "opfnctns.h"
#include "srt_h_allFx3F.h"
#include "math.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"
#include "LGMSVUtil.h"
#include "FundLegProdStruct.h"

#define	BIG_BIG_NUMBER	1E40
#define MAXCPN 1000

/*	Functions for the funding leg */


Err fill_funding_leg(
/*	Coupons that started before today are disregarded */
long		today,
/*	EOD Flag */
int			eod_flag,				/*	0: I, 1: E */
int			spot_lag,
double		fund_not,
int			fund_ncpn,
long		*fund_fix,
long		*fund_start,
long		*fund_pay,
char		**fund_basis,
double		*fund_spr,
double		*fund_mrg,
FUNDING_LEG	fund_leg)
{
	int				i, j;
	SrtBasisCode	bas;
	FUNDING_CPN		cpn;
	Err				err			= NULL;

	/*	Spot Lag */
	fund_leg->spot_lag = spot_lag;

	/*	Notional */
	fund_leg->notional = fund_not;

	/*	Initialise pointers to NULL */
	fund_leg->cpn = NULL;

	/*	Skip coupons fixed before today */
	i = 0;
	while (i < fund_ncpn && fund_fix[i] < today + eod_flag)
//	while (i < fund_ncpn && add_unit( fund_start[i], -spot_lag, 
//						SRT_BDAY, MODIFIED_SUCCEEDING ) < today + eod_flag)
	{
		i++;
	}

	/*	Check that at least one coupon is left */
	if (i == fund_ncpn)
	{
		/*	err = "All funding coupons start before today in ccf_fill_fund_leg"; */
		fund_leg->num_cpn = 0;
		goto FREE_RETURN;
	}

	/*	Allocate memory */
	fund_leg->num_cpn = fund_ncpn - i;
	fund_leg->cpn = (funding_cpn*) calloc (fund_leg->num_cpn, sizeof (funding_cpn));
	if (!fund_leg->cpn)
	{
		err = "Allocation error in fill_fund_leg";
		goto FREE_RETURN;
	}

	/*	Fill coupons information */
	j = 0;
	while (i < fund_ncpn)
	{
		cpn = fund_leg->cpn + j;

		/*	Dates */
		cpn->start_date = fund_start[i];
		cpn->pay_date = fund_pay[i];
		
		/*	Times */
		cpn->start_time = (cpn->start_date - today) * YEARS_IN_DAY;
		cpn->pay_time = (cpn->pay_date - today) * YEARS_IN_DAY;
		
		/*	Coupon */
		err = interp_basis (fund_basis[i], &bas);
		if (err)
		{
			goto FREE_RETURN;
		}
		cpn->cvg = coverage (fund_start[i], fund_pay[i], bas);
		cpn->cpn = fund_not * cpn->cvg * (fund_spr[i] + fund_mrg[i]);

		i++;
		j++;
	}
	
	err = check_funding_leg (fund_leg);

FREE_RETURN:

	if (err)
	{
		free_funding_leg (fund_leg);
	}

	return err;
}

/*	Check dates consistency */
Err check_funding_leg(
FUNDING_LEG	fund_leg)
{
	int				i;
	
	/*	Check that start and pay dates are increasing */
	for (i=1; i<fund_leg->num_cpn; i++)
	{
		if (fund_leg->cpn[i].start_date < fund_leg->cpn[i-1].start_date)
		{
			return "Start dates should be increasing in funding leg";
		}

		if (fund_leg->cpn[i].pay_date < fund_leg->cpn[i-1].pay_date)
		{
			return "Pay dates should be increasing in funding leg";
		}
	}

	/*	Check that pay dates are after start dates */
	for (i=0; i<fund_leg->num_cpn; i++)
	{
		if (fund_leg->cpn[i].pay_date < fund_leg->cpn[i].start_date)
		{
			return "Pay dates should be after start dates in funding leg";
		}
	}

	/*	OK */
	return NULL;
}

/*	Free */
Err free_funding_leg(
FUNDING_LEG	fund_leg)
{
	if (fund_leg->cpn)
	{
		free (fund_leg->cpn);
		fund_leg->cpn = NULL;
	}

	return NULL;
}



Err FundLeg_RequestDfDates(funding_leg *fundleg, long date, int *n_dates, long **pdates)
{
	int i, j, ndates;
	Err err=NULL;

	for ( i=0; i < fundleg->num_cpn &&
		  date - 10 > add_unit( (fundleg->cpn[i]).start_date, -fundleg->spot_lag,
		  SRT_BDAY, MODIFIED_SUCCEEDING );
		  i++ );

	if (i == fundleg->num_cpn)
	{
		*pdates = NULL;
		*n_dates = 0;
		err = "No more funding coupon";
		smessage("No more funding coupon");
		goto FREE_RETURN;
	}

	ndates = fundleg->num_cpn - i + 1;
	*pdates = (long *) calloc(ndates, sizeof(long));
	if(!(*pdates))
	{
		smessage("Memory allocation failed in FundLeg_RequestDfDates");
		err = "Memory allocation failed in FundLeg_RequestDfDates";
		goto FREE_RETURN;
	}


	*n_dates = ndates;

	(*pdates)[0] = (fundleg->cpn[i]).start_date;
	for(j=i;j<fundleg->num_cpn;++j)
	{
		(*pdates)[j+1-i] = (fundleg->cpn[j]).pay_date;
	}

FREE_RETURN:

	if(err)
	{
		if(*pdates)
		{
			free(*pdates);
			*pdates = NULL;
		}
	}

	return err;
}


Err FundLeg_Payoff(funding_leg *fundleg, long today, long date,
						double *dfs, double *payoff)
{
	Err err = NULL;
	int i, k;

	for ( i=0; i < fundleg->num_cpn &&
		  date > add_unit( (fundleg->cpn[i]).start_date, -fundleg->spot_lag,
		  SRT_BDAY, MODIFIED_SUCCEEDING );
		  i++ );

	*payoff = dfs[0] * fundleg->notional;

	for (k=i; k<fundleg->num_cpn; ++k)
	{
		*payoff += fundleg->cpn[k].cpn * dfs[k-i+1];
	}

	*payoff -= dfs[fundleg->num_cpn-i]* fundleg->notional;

//	*payoff = dfs[fundleg->num_cpn-i]* fundleg->notional;

	return NULL;
}



