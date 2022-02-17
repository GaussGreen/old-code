#include "math.h"

#include "MBSPTCalib.h"
#include "MBSPTCaller.h"
#include "TBACalib.h"
#include "TBAProdStruct.h"


/*****  IOPO_Price  ***********************************************************
*
*       Calculate the IO price in the lattice.
*		AKA CashFlow Engine
*/
void    IOPO_Price (      double          *IO,                                    /* IO price in the lattice at the current period */
			double             *PO,                                   /* PO price in the lattice at the current period */
			double          **Zero,
			double          *Discount,                              /* One period discount factor at each node at period NbPer */
			double          *Discount1,                             /* One period discount factor without OAS */
			double          *pu,                                    /* Array of probabilities in the up state in the interest rate tree */
			double          *pd,                                    /* Array of probabilities in the down state in the interest rate tree */
			double          *p0,                                    /* Array of probabilities in the flat state in the interest rate tree */
			double          *Amort,                                 /* Array of main amortization factor in the interest rate lattice */
			double          Coupon,                                 /* Underlying coupon (annualized) */
			int             CouponNb,                               /* Coupon Number: = deal_data->Term at the end of the tree */
			int             F,                                      /* Frequency of underlying payment as an integer */
			int             CpPaid,                                 /* Is the current period a coupon payment date */
			int             NbPer,                                  /* Current time period */
			int		Lockout,                                /* Lockout date of the SFS expressed in number of nodes */
			int             ExEnd,                                  /* Maturity of the IO expressed in number of nodes */
			TREE_DATA       *tree_data,                             /* Structure of tree data */
			DEAL_DATA       *deal_data)                             /* Structure of deal data */
{
	double
		TotalAmort,							/* Total amortization factor including all effects */
		TimeBalance;                                                    /* Notional balance left at this time after natural (deterministic) repayment */
	int                                                                     /* on the amortization factor due to interest rate (Amort).             */
		Bottom,                                                         /* Bottom node index in the interest rate tree */
		Top,                                                            /* Top node index in the interest rate tree */
		i,                                                              /* Node index for the interest rate */
		j;                                                              /* Delay index */

	int	MM;								/* amortization schedule */


	if (deal_data->Balloon == 'Y') MM = deal_data->BalloonSch;
		else MM = deal_data->Term;

	Bottom  = (int) max (0., NbPer - tree_data->NodeMax);                    /* If NbPer>NodeMax the interest rate tree is cut at the top and the bottom */
	Top     = (int) min (2*NbPer, NbPer + tree_data->NodeMax);

	if (NbPer == ExEnd)                                                     /* Initialization of the IO price at the end of the tree */
	{
		for (i = Bottom; i <= Top; i ++){
			if (deal_data->invfloater == 'Y' ) /* inverse floator */	
				IO[i] = 0.0;
			else if (deal_data->invfloater == 'C' ) /* cap */	
				IO[i] = 0.0;
			else                                
				IO[i] = Coupon / F * pow(Discount[i], (double) 
					deal_data->pay_delay/30.);

			PO[i] = 100. * pow(Discount[i], (double) deal_data->pay_delay/30.);     /* taking into account payment delay */
		}
		return;
	}


	Dev (   IO,                                                             /* Discounted expected value function */
		Discount,
		pu,
		pd,
		p0,
		NbPer,
		tree_data);
	Dev(    PO,								/* Discounted expected value function */
		Discount,
		pu,
		pd,
		p0,
		NbPer,
		tree_data);

		if (deal_data->Schedule == 'N') TimeBalance = 1.0;		/* if no scheduled amortization, Balance remains at 1.*/
			else
		TimeBalance = (pow (1. + deal_data->Gwac/F/100., (double) MM - CouponNb + 1) - 1. - deal_data->Gwac/F/100.)  /* Balance is calculated so that it is 1 at onset */
			    / (pow (1. + deal_data->Gwac/F/100., (double) MM - CouponNb + 1) - 1.);                 /* of mortgage and goes to zero at maturity.      */


	if (CpPaid)                                                             /* If we are at an coupon period */
	{
		if ((NbPer >= Lockout) && (CouponNb < deal_data->Term - deal_data->Delay))      /* If we are after the lockout date and "Delay" number of coupons */
		{                                                                               /* before maturity we amortize according to the index level.      */
			for (i = Bottom; i <= Top; i ++)
			{
				TotalAmort = Amort[i];                          /* Amortization that depends on index level only */
				IO[i] *= (1. - TotalAmort);                     /* Amortization of the IO linked to the index */
				PO[i] = PO[i] * (1. - TotalAmort) + 100. * TotalAmort * Zero[deal_data->Delay][i]
					* pow(Discount[i], (double) deal_data->pay_delay / 30.);      /* If there is a delay, the par payment  is reduced, also there is a payment delay  */

				for (j = 1; j <= deal_data->Delay; j++)                 /* If there is a delay in the index amortization we should not have */
					IO[i] += Coupon / F * TotalAmort * Zero[j][i];  /* amortized the next few coupons. We correct that here.            */

			if (deal_data->invfloater == 'Y')               /* inverse floater */  
			        IO[i] +=(1. - TotalAmort) *  max(0, Coupon / F - (1./Discount1[i] - 1.) * 100. *
					(double) (tree_data->Ppy + 1.0) / F ) * 
					pow(Discount[i], (double) deal_data->pay_delay / 30. + 1.0);
			else if (deal_data->invfloater == 'C')               /* cap */  
			        IO[i] += (1. - TotalAmort) *  min(Coupon / F, (1./Discount1[i] - 1.) * 100. *
					(double) (tree_data->Ppy + 1.0) / F ) * 
					pow(Discount[i], (double) deal_data->pay_delay / 30. + 1.0);
			}  /* for i */
		}  /* if */
		

		for (i = Bottom; i <= Top; i ++)
		{
			IO[i] *= TimeBalance;                                   /* At a coupon date there always is a scheduled amortization. */

			if (deal_data->invfloater != 'Y' && deal_data->invfloater != 'C')               /* if not inverse floater or cap */  
				IO[i] += Coupon / F * pow(Discount[i], (double) deal_data->pay_delay /30.);                                    /* The coupon is paid on balance left before all type of amortization */

			PO[i] = PO[i] * TimeBalance + 100. * (1. - TimeBalance) * pow(Discount[i], (double) deal_data->pay_delay /30.);        /* At a coupon date there always is a scheduled amortization. */


		}  /* for i */

	}  /* if */

	return;

}  /* IOPO_Price */