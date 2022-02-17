
/**********************************************************
*
*	instrument.cpp
*
***********************************************************/
#include "kplatform.h"
#include "instruments.h"
#include "macros.h"
#include "gtobf.h"

#ifdef _DEBUG
// extern std::ofstream out;
// extern std::ofstream out("c:\\yy.dat");
#endif

void Instrument::update()
{
	DTimeSlice tsStock = m_tree->get_stock_price();
	DTimeSlice tsAsset = (m_tree->*m_underFunc)();
	int	dim = m_tree->get_curr_limit()[0];
	KDate expiryDate = get_expiryDate();
	KDate currDate = m_tree->get_curr_date();
//	double integratValue;								//HY4.1v, commented out in HY4.1a

	// if ts is empty allocate
	if(m_ts.dim() == 0)
	{
		m_ts.resize(m_tree);
	}

	if(m_ts_opt.dim() == 0)				//HY3.4v
	{
		m_ts_opt.resize(m_tree);		//HY3.4v
	}

	if(m_ts_vega.dim() == 0)
	{
		m_ts_vega.resize(m_tree);
	}

	if(m_ts_opt_vega.dim() == 0)				//HY3.4v
	{
		m_ts_opt_vega.resize(m_tree);			//HY3.4v
	}

	int i;

//	if(currDate == expiryDate)

	if(currDate >= expiryDate)				//Hy4.1v
	{
//		out << "\nDIM="<<dim<<"\n"<<std::endl;
		for (i = -dim; i <= dim; i++)
		{
/*			if(DefaultSwap *pt = dynamic_cast<DefaultSwap*>this)					//HY4.1v
			{
				if(m_option != NULL)												//HY4.1v
				{
					integratValue = get_integratedValue(m_ts, i, m_startLimit);		//HY4.1v
				}
			}
			else if(Bond *pt1 = dynamic_cast<Bond*>this)							//HY4.1v
			{
				if(m_option != NULL || m_cv_option != NULL)							//HY4.1v
				{
					integratValue = get_integratedValue(m_ts, i, m_startLimit);		//HY4.1v
				}
			}
			else
			{
				integratValue = m_ts[i];
			}
*/
			payoff_new(&m_ts[i],&m_ts_opt[i],&m_ts_vega[i],&m_ts_opt_vega[i],tsAsset[i],tsStock[i],m_ts[i]); //HY4.1v
//			payoff_new(&m_ts[i], &m_ts_opt[i], &m_ts_vega[i], &m_ts_opt_vega[i], tsAsset[i],tsStock[i]);	//HY3.4v

//			m_ts[i] = payoff(m_ts[i],tsAsset[i],tsStock[i],&m_ts_vega[i]);
//			m_ts_vega[i] = 0.0;
		}
//		out<<"\nDate: "<<GtoFormatDate(currDate)<<"\t dim="<<dim<<"\n"<<std::endl;
//		out<<" dimExpiry : "<<dim<<"\n"<<m_ts<<std::endl;
//		out<<" dimExpiry : "<<dim<<"\n"<<tsAsset<<std::endl;
	}
	else
	{
//		m_tree->ddv(m_ts_vega,m_ts);
//		m_tree->dev(m_ts);

		m_tree->ddv(m_ts_vega, m_ts, false);				//HY3.4v
		m_tree->dev(m_ts, false);						//HY3.4v

		m_tree->ddv(m_ts_opt_vega,m_ts_opt, true);		//HY3.4v
		m_tree->dev(m_ts_opt, true);						//HY3.4v

		jump_dev();										//jump version

//		m_tree->ddv(m_ts_cvOpt_vega,m_ts_cvOpt, true, true);		//HY3.4v
//		m_tree->dev(m_ts_cvOpt, true, true);						//HY3.4v

//		out<<"\nDate: "<<GtoFormatDate(currDate)<<"\t dim="<<dim<<"\n"<<std::endl;
//		out<<"\n dim : "<<dim<<"\n"<<m_ts<<std::endl;
//		out<<"\n dim : "<<dim<<"\n"<<tsAsset<<std::endl;
	
//		out << "\nDIM="<<dim<<"\n"<<std::endl;
		for (i = -dim; i <= dim; i++)
		{
/*			//HY4.1v
			if(DefaultSwap *pt = dynamic_cast<DefaultSwap*>(this))		//HY4.1a
			{
				if(pt->get_option_pt() != NULL)					
				{
					int startLimit = m_tree->get_startLimit();									
					int endLimit = m_tree->get_endLimit();										
					int numNodes = 1 + int((startLimit-1) * (endLimit-dim) / (endLimit-startLimit));

					integratValue = get_integratedValue(m_ts, i, numNodes);		
				}
			}
			else if(Bond *pt1 = dynamic_cast<Bond*>(this))				
			{
				if(pt1->get_option_pt() != NULL || pt1->get_option_cv_pt() != NULL)				
				{
					int startLimit = m_tree->get_startLimit();									
					int endLimit = m_tree->get_endLimit();										
					int numNodes = 1 + int((startLimit-1) * (endLimit-dim) / (endLimit-startLimit));

					integratValue = get_integratedValue(m_ts, i, numNodes);		
				}
			}
			else
			{
				integratValue = m_ts[i];
			}
			//end HY4.1v
*/

//			payoff_new(&m_ts[i], &m_ts_opt[i], &m_ts_vega[i], &m_ts_opt_vega[i],
//							tsAsset[i],tsStock[i], integratValue);					//HY4.1v

			payoff_new(&m_ts[i], &m_ts_opt[i], &m_ts_vega[i], &m_ts_opt_vega[i],
							tsAsset[i],tsStock[i], m_ts[i]);						//HY4.1a

//			m_ts[i] = payoff(m_ts[i],tsAsset[i],tsStock[i],&m_ts_vega[i]);
//			m_ts_vega[i] = vegaSet(m_ts_vega[i],tsAsset[i]);
		}
		
//		out<<"\n dimA : "<<dim<<"\n"<<m_ts<<std::endl;
//		out<<" asset : "<<"\n"<<tsAsset<<std::endl;
		//get stock for the next run
	//	DTimeSlice stock = m_tree->get_stock_price();

	}
}

void	Instrument::jump_dev()				//jump version
{
	DTimeSlice	ts_probNoJump = m_tree->get_probNoJump();
	int			i;

	int	dim  = m_tree->get_curr_limit()[0];
	for (i = -dim; i <= dim; i++)
	{
		m_ts[i] = m_ts[i] * ts_probNoJump[i] + (1.0 - ts_probNoJump[i]) * get_defaultValue();	
		m_ts_opt[i] = m_ts_opt[i] * ts_probNoJump[i];				//a bond's embedded options are worth zero on default
		m_ts_vega[i] = m_ts_vega[i] * ts_probNoJump[i];			
		m_ts_opt_vega[i] = m_ts_opt_vega[i] * ts_probNoJump[i];
	}
}

double Instrument::vegaSet(double y, double asset)
{
	KDate currDate = m_tree->get_curr_date();
	double limit1 = (get_strikes())->get_rate(currDate);

	if(asset < limit1)
	{
		return 0.0;
	}
	else
	{
		return y;
	}
}
	

double  Instrument::get_integratedValue(DTimeSlice ts, int ts_idx, int num)		//HY4.1v added function
{
	int i;
	double vollim = m_tree->get_vollim();
	DTimeSlice  stock_ts = m_tree->get_stock_price();

	if(vollim == 0.0 || IS_ALMOST_ZERO(stock_ts[ts_idx]))		// only when node is above barrier.
	{
		return ts[ts_idx];
	}
	else
	{
		double prob, prob1, prob2, temp1, temp2;

		double integratedValue = 0.0;
		double cumulativeProb = 0.0;

		int k1, k2, dim, nn;

		dim = m_tree->get_curr_limit()[0];
		k1 = ts_idx - num;
		k2 = ts_idx + num;

		if(k1 < -dim)
		{
			k1 = -dim;
		}
		if(k2 > dim)
		{
			k2 = dim;
		}

		nn=MIN(k2 - ts_idx, ts_idx - k1);
		k1 = ts_idx - nn;
		k2 = ts_idx + nn;

/*		int k2 = ts_idx + num;
		int dim = m_tree->get_curr_limit()[0];

		if(k2 > dim)
		{
			k2 = dim;
		}

		int k1 = 2*ts_idx - k2;
		if(k1 < -dim)
		{
			k1 = -dim;
		}
*/
		DTimeSlice  asset_ts = m_tree->get_asset_value();
		AssetToStockMapping1* ats = dynamic_cast<AssetToStockMapping1*>(m_tree->get_assetToStockMap());

		double level = ats->get_lim();
//		double dps = ats->get_debtPerShare();

//		double assetToBarrier = log(level*dps/asset_ts[ts_idx]);
		double assetToBarrier = log(level/asset_ts[ts_idx]);

		double assetInit_adj = asset_ts[ts_idx] * exp(vollim*vollim);		//HY4.2c

		for (i = k1; i <= k2; i++)
		{
			if(!(IS_ALMOST_ZERO(stock_ts[i])))
			{
//				double assetDist = log(asset_ts[i]/asset_ts[ts_idx]);		//HY4.1a, JP
//				temp1 = assetDist / vollim - 0.5 * vollim;
//				temp2 = (assetDist - 2.0 * assetToBarrier) / vollim - 0.5 * vollim;
//				prob1 = exp(-0.5 * temp1 * temp1);
//				prob2 = exp(-0.5 * temp2 * temp2);
//				prob = prob1 - prob2 * exp(assetToBarrier);

				//HY4.2b, GP, adjust barrier
//				temp1 = log(asset_ts[i]/asset_ts[ts_idx]);
//				temp2 = log(asset_ts[i] * asset_ts[ts_idx] / level / level) + 2.0*vollim*vollim;
//				prob1 = exp(-temp1*temp1/(2.0*vollim*vollim));
//				prob2 = exp(-temp2*temp2/(2.0*vollim*vollim));
//				prob = exp(-temp1/2.0)*(prob1 - prob2);

				//HY4.2c, GP, adjust initial asset
				temp1 = log(asset_ts[i]/assetInit_adj);
				temp2 = log(asset_ts[i] * assetInit_adj / level / level);
				prob1 = exp(-temp1*temp1/(2.0*vollim*vollim));
				prob2 = exp(-temp2*temp2/(2.0*vollim*vollim));
				prob = exp(-temp1/2.0)*(prob1 - prob2);

				integratedValue += prob * ts[i];
				cumulativeProb += prob;
			}
		}
		integratedValue = integratedValue / cumulativeProb;

		if(IS_ALMOST_ZERO(stock_ts[k1]))
		{
//			double dist = asset_ts[ts_idx]/level/dps * exp(vollim * vollim);
			double dist = asset_ts[ts_idx]/level * exp(vollim * vollim);
			temp1 = -vollim/2.0 + log(dist)/vollim;
			temp2 = -vollim/2.0 - log(dist)/vollim;

			double p_nd = GtoNormalCum(temp1) - dist * GtoNormalCum(temp2);
//			integratedValue = integratedValue * p_nd + get_defaultValue() * (1.0 - p_nd);  //get_defaultValue() gets either CDS or bond, not correct for option
			integratedValue = integratedValue * p_nd + ts[k1] * (1.0 - p_nd);		//ts[k1] over-estimated
		}
		return integratedValue;
	}
}

double  Instrument::get_spot_ts_value()		//HY4.1v
{
/*	double vollim = m_tree->get_vollim();

	if(vollim == 0.0)
	{
*/
		return m_ts[0];
/*	}
	else
	{
		int num = m_tree->get_curr_limit()[0];

		double spot_ts_value = get_integratedValue(m_ts, 0, num);
		return spot_ts_value;
	}
*/
}


double  Instrument::get_option_ts_value()		//HY4.1v
{
		return m_ts_opt[0];
}


/*
double  Instrument::get_option_ts_value()		//HY4.1a
{
	double vollim = m_tree->get_vollim();

	if(vollim == 0.0)
	{
		return m_ts_opt[0];
	}
	else
	{
		int num = m_tree->get_curr_limit()[0];

		double option_ts_value = get_integratedValue(m_ts_opt, 0, num);
		return option_ts_value;
	}

}
*/

double  Instrument::get_vega_value()		//HY4.1v
{
/*	double vollim = m_tree->get_vollim();

	if(vollim == 0.0)
	{
*/
		return m_ts_vega[0];
/*
	}
	else
	{
		int num = m_tree->get_curr_limit()[0];

		double vega = get_integratedValue(m_ts_vega, 0, num);
		return vega;
	}
*/
}


double  Instrument::get_option_vega_value()		//HY4.1v
{
		return m_ts_opt_vega[0];
}


/*
double  Instrument::get_option_vega_value()		//HY4.1a
{
	double vollim = m_tree->get_vollim();

	if(vollim == 0.0)
	{

		return m_ts_opt_vega[0];
	}
	else
	{
		int num = m_tree->get_curr_limit()[0];

		double vega = get_integratedValue(m_ts_opt_vega, 0, num);
		return vega;
	}
}
*/
	
double  Instrument::get_delta_value()
{
	double delta;
	DTimeSlice& stock = m_tree->get_stock_price();

/*	double vollim = m_tree->get_vollim();

	if(vollim == 0.0)
	{
*/
		delta= (m_ts[-1]-m_ts[0])/(stock[-1]-stock[0])+
				  (m_ts[1]-m_ts[0])/(stock[1]-stock[0]) -
				  (m_ts[-1]-m_ts[1])/(stock[-1]-stock[1]);
/*	}
	else
	{
		int num = m_tree->get_curr_limit()[0];
		double m_ts_minus = get_integratedValue(m_ts, -1, num);
		double m_ts_0 = get_integratedValue(m_ts, 0, num);
		double m_ts_plus = get_integratedValue(m_ts, 1, num);

		delta= (m_ts_minus - m_ts_0)/(stock[-1] - stock[0])+
					  (m_ts_plus - m_ts_0)/(stock[1] - stock[0]) -
					  (m_ts_minus - m_ts_plus)/(stock[-1] - stock[1]);
	}
*/
	return delta;
}

/*
double  Instrument::get_option_delta_value()		//HY4.1a
{
	double delta;
	DTimeSlice& stock = m_tree->get_stock_price();
	double vollim = m_tree->get_vollim();

	if(vollim == 0.0)
	{
		delta= (m_ts_opt[-1] - m_ts_opt[0])/(stock[-1] - stock[0])+
				  (m_ts_opt[1] - m_ts_opt[0])/(stock[1] - stock[0]) -
				  (m_ts_opt[-1] - m_ts_opt[1])/(stock[-1] - stock[1]);
	}
	else
	{
		int num = m_tree->get_curr_limit()[0];
		double m_ts_minus = get_integratedValue(m_ts_opt, -1, num);
		double m_ts_0 = get_integratedValue(m_ts_opt, 0, num);
		double m_ts_plus = get_integratedValue(m_ts_opt, 1, num);

		delta= (m_ts_minus - m_ts_0)/(stock[-1] - stock[0])+
				  (m_ts_plus - m_ts_0)/(stock[1] - stock[0]) -
				  (m_ts_minus - m_ts_plus)/(stock[-1] - stock[1]);
	}
	return delta;
}
*/


double  Instrument::get_option_delta_value()				//HY3.4v
{
	DTimeSlice& stock = m_tree->get_stock_price();

	double delta= (m_ts_opt[-1]-m_ts_opt[0])/(stock[-1]-stock[0])+
				  (m_ts_opt[1]-m_ts_opt[0])/(stock[1]-stock[0]) -
				  (m_ts_opt[-1]-m_ts_opt[1])/(stock[-1]-stock[1]);
	return delta;
}


/*
double  Instrument::get_gamma_value()
{
	DTimeSlice& stock = m_tree->get_stock_price();

	double gamma= ((m_ts[-1]-m_ts[0])/(stock[-1]-stock[0])-
				  (m_ts[0]-m_ts[1])/(stock[0]-stock[1]))/((stock[-1]-stock[1])/2);
	return gamma;
}
*/

double  Instrument::get_gamma_value()
{
	double gamma;
	DTimeSlice& stock = m_tree->get_stock_price();

/*	double vollim = m_tree->get_vollim();

	if(vollim == 0.0)
	{
*/
		gamma= ((m_ts[-1] - m_ts[0])/(stock[-1] - stock[0])-
				  (m_ts[0] - m_ts[1])/(stock[0] - stock[1]))/((stock[-1] - stock[1])/2);
/*	}
	else
	{
		int num = m_tree->get_curr_limit()[0];
		double m_ts_minus = get_integratedValue(m_ts, -1, num);
		double m_ts_0 = get_integratedValue(m_ts, 0, num);
		double m_ts_plus = get_integratedValue(m_ts, 1, num);

		gamma= ((m_ts_minus - m_ts_0)/(stock[-1] - stock[0])-
				  (m_ts_0 - m_ts_plus)/(stock[0] - stock[1]))/((stock[-1] - stock[1])/2);
	}
*/
	return gamma;
}


/*
double  Instrument::get_option_gamma_value()		//HY4.1a
{
	double gamma;
	DTimeSlice& stock = m_tree->get_stock_price();
	double vollim = m_tree->get_vollim();

	if(vollim == 0.0)
	{
		gamma= ((m_ts_opt[-1]-m_ts_opt[0])/(stock[-1]-stock[0])-
				  (m_ts_opt[0]-m_ts_opt[1])/(stock[0]-stock[1]))/((stock[-1]-stock[1])/2);
	}
	else
	{
		int num = m_tree->get_curr_limit()[0];
		double m_ts_minus = get_integratedValue(m_ts_opt, -1, num);
		double m_ts_0 = get_integratedValue(m_ts_opt, 0, num);
		double m_ts_plus = get_integratedValue(m_ts_opt, 1, num);

		gamma= ((m_ts_minus - m_ts_0)/(stock[-1] - stock[0])-
				  (m_ts_0 - m_ts_plus)/(stock[0] - stock[1]))/((stock[-1] - stock[1])/2);
	}
	return gamma;
}
*/


double  Instrument::get_option_gamma_value()
{
	DTimeSlice& stock = m_tree->get_stock_price();

	double gamma= ((m_ts_opt[-1]-m_ts_opt[0])/(stock[-1]-stock[0])-
				  (m_ts_opt[0]-m_ts_opt[1])/(stock[0]-stock[1]))/((stock[-1]-stock[1])/2);
	return gamma;
}


double DefaultSwap::payoff(double y,double asset,double stock, double *vega) const 
{
	KDate expiryDate = get_expiryDate();
	KDate currDate = m_tree->get_curr_date();
	double limit1 = (get_strikes())->get_rate(currDate);
	double limit2 = m_strikes2->get_rate(expiryDate);
	double payment;												//
	double returnValue;

	SettleInfo settlement = get_settlement();					//
	KDate settlementDate = settlement.get_settlementDate();		//

	if(currDate == expiryDate)
	{
		if(expiryDate <= settlementDate)					//
		{
			returnValue=0.0;									//
		}
		else
		{
			*vega = 0.0;
			if( asset > limit2)
			{
				returnValue = -m_feePayment.get_Payment(currDate);
			}
			else
			{
				if( asset > limit1)
				{
					returnValue = (1-m_recovery)*(limit2-asset)/(limit2-limit1)-m_feePayment.get_accruCoupon_recover(currDate);
				}
				else
				{
					returnValue = (1-m_recovery)-m_feePayment.get_accruCoupon_recover(currDate);
				}
			}
		}
	}
	else
	{
			//add option 
		double optionValue;
		
		if(m_option == NULL)
		{
			optionValue = 0.0;
		}
		else
		{
			double accruCoupon = m_feePayment.get_accruCoupon(currDate);
			optionValue = m_option->payoff(y+accruCoupon, currDate, vega);
		}
//		out << "("<<y<<", "<<optionValue<<")";
		
		y = y + optionValue;

		if (currDate <= settlementDate)						//
		{
			payment = 0.0;									//
		}
		else
		{
			payment = m_feePayment.get_Payment(currDate);	//
		}

//		double payment= m_feePayment.get_Payment(currDate);

		if(asset > limit2)
		{
//			out<<"("<<y<<", "<<payment<<")  ";
			returnValue = y-payment;				//long protection, short fee.
		}
		else
		{
			if(asset > limit1)
			{
				double feeRecover = m_feePayment.get_couponAccruPercent();
				returnValue = y-payment*(feeRecover+(1-feeRecover)/(limit2-limit1)*(asset-limit1));
			}
			else
			{
				double accrupay = m_feePayment.get_accruCoupon_recover(currDate);
	//			out<<"(accru, "<<accrupay<<")  ";
				returnValue = 1-m_recovery-accrupay;
				*vega = 0.0;
			}
		}
	}
	return returnValue;
}


void DefaultSwap::payoff_new(double *underly, double *optValue,
							 double *underVega, double *optVega,
							 double asset, double stock, double integratUnderly) 
{
	KDate expiryDate = get_expiryDate();
	KDate currDate = m_tree->get_curr_date();
	double limit1 = (get_strikes())->get_rate(currDate);
	double limit2 = m_strikes2->get_rate(expiryDate);
	double payment;												//
	double y;

	SettleInfo settlement = get_settlement();					//
	KDate settlementDate = settlement.get_settlementDate();		//
			
//	if(currDate == expiryDate)
	if(currDate >= expiryDate)				//HY4.1v
	{
		*underVega = 0.0;									
		*optVega = 0.0;									

		*optValue = 0.0;									
		if(expiryDate <= settlementDate)					//
		{
			*underly = 0.0;									
		}
		else
		{
			if( asset > limit2)
			{
				payment = -m_feePayment.get_Payment(currDate);
			}
			else
			{
				if( asset > limit1)
				{
					payment = (1-m_recovery)*(limit2-asset)/(limit2-limit1)-m_feePayment.get_accruCoupon_recover(currDate);
				}
				else
				{
					payment = (1-m_recovery)-m_feePayment.get_accruCoupon_recover(currDate);
				}
			}
			*underly = payment;									
		}
	}
	else
	{
		if (currDate <= settlementDate)						//
		{
			payment = 0.0;									//
		}
		else
		{
			payment = m_feePayment.get_Payment(currDate);	//
		}

		if(asset <= limit1)
		{
			double accrupay = m_feePayment.get_accruCoupon_recover(currDate);
//			out<<"(accru, "<<accrupay<<")  ";
			*underly = 1-m_recovery-accrupay;

			*underVega = 0.0;
			*optVega = 0.0;
			*optValue = 0.0;		//assume only putable protection
		}
		else
		{
			//add option 
			double optionValue;
			
			if(m_option == NULL)
			{
				optionValue = 0.0;
			}
			else
			{
				double accruCoupon = m_feePayment.get_accruCoupon(currDate);
//				y = *underly + *optValue + accruCoupon;
				y = integratUnderly + *optValue + accruCoupon;
				optionValue = m_option->payoff(y, currDate, optVega);

				if(*optVega == 0.0 && optionValue !=0.0)
				{
					*underVega = 0.0;
				}
			}
//			out << "("<<y<<", "<<optionValue<<")";
			*optValue += optionValue;

			//add coupon
			if(asset > limit2)
			{
//				out<<"("<<y<<", "<<payment<<")  ";
				*underly = *underly - payment;				//long protection, short fee.
			}
			else
			{
				double feeRecover = m_feePayment.get_couponAccruPercent();
				*underly = *underly - payment*(feeRecover+(1-feeRecover)/(limit2-limit1)*(asset-limit1));
			}
		}
	}
}	


double DefaultSwap::get_defaultValue()
{
	KDate currDate = m_tree->get_curr_date();

	double accrupay = m_feePayment.get_accruCoupon_recover(currDate);
	double defaultValue = 1-m_recovery-accrupay;

	return defaultValue;
}


double DSAverageLife::payoff(double y,double asset, double stock, double *vega) const 
{
	KDate expiryDate = get_expiryDate();
	KDate currDate = m_tree->get_curr_date();
	double limit1 = get_strikes()->get_rate(currDate);
	if(currDate == expiryDate)
	{
		return 0.0;
	}
	else
	{
		if(asset > limit1)
		{
			KTimePoint	&tp = m_tree->get_curr_tp();
			double     timeTillNextTp =  tp[__timeTillNextTp];
			return y+timeTillNextTp;
		}
		else
		{
			*vega = 0.0;
			return 0.0;
		}
	}
}

void DSAverageLife::payoff_new(double *underly, double *optValue,
							   double *underVega, double *optVega,
							   double asset, double stock, double integratUnderly) 
{
	*optValue = 0.0;
	*optVega = 0.0;
	*underly = payoff(*underly, asset, stock, underVega);
}	


double DSAverageLife::get_defaultValue()			//not used, needs to be defined when used
{
	return 0.0;
}


double CallOptionT::payoff(double y,double asset, double sstock, double *vega) const 
{
	KDate expiryDate = get_expiryDate();
	KDate currDate = m_tree->get_curr_date();
	double strike = get_strikes()->get_rate(currDate);
	double optionValue =0.0;

	//recover the true stock forward price
	KTimePoint	&tp = m_tree->get_curr_tp();		//HY4.2
	double debtDriftFactor = tp[__cumDebtDrift];	//HY4.2
	asset *= debtDriftFactor;						//HY4.2
	
	if(currDate == expiryDate)
	{
		optionValue = MAX(asset-strike,0.0);
		*vega = 0.0;
	}
	else
	{
		if(m_isamerican == false)
		{
			optionValue = y;
		}
		else
		{
			optionValue = MAX(asset-strike,y);
			if( (asset-strike) > y)
			{
				*vega = 0.0;
			}
		}
	}

	return optionValue;
}


void CallOptionT::payoff_new(double *underly, double *optValue,
							 double *underVega, double *optVega,
							 double asset, double stock, double integratUnderly) 
{
	*optValue = 0.0;
	*optVega = 0.0;
	*underly = payoff(*underly, asset, stock, underVega);
}	

double CallOptionT::get_defaultValue()
{
	return 0.0;
}


double PutOptionT::payoff(double y,double asset, double sstock, double *vega) const 
{
	KDate expiryDate = get_expiryDate();
	KDate currDate = m_tree->get_curr_date();
	double strike = get_strikes()->get_rate(currDate);
	double optionValue =0.0;

	//recover the true stock forward price
	KTimePoint	&tp = m_tree->get_curr_tp();		//HY4.2
	double debtDriftFactor = tp[__cumDebtDrift];	//HY4.2
	asset *= debtDriftFactor;						//HY4.2

	if(currDate == expiryDate)
	{
		optionValue = MAX(strike - asset,0.0);
	}
	else
	{
		if(m_isamerican == false)
		{
			optionValue = y;
		}
		else
		{
			optionValue = MAX(strike - asset,y);
			if( (strike - asset)>y)
			{
				*vega = 0.0;
			}
		}
	}

	return optionValue;
}


void PutOptionT::payoff_new(double *underly, double *optValue,
							double *underVega, double *optVega,
							double asset, double stock, double integratUnderly) 
{
	*optValue = 0.0;
	*optVega = 0.0;
	*underly = payoff(*underly, asset, stock, underVega);
}	
	
double PutOptionT::get_defaultValue()
{
	KDate		currDate = m_tree->get_curr_date();
	KDate		expiryDate = get_expiryDate();
	KRateCurve	rateCurve = m_tree->get_index_zero();

	double		strike, defaultValue;

	if(m_isamerican == true)
	{
		strike = get_strikes()->get_rate(currDate);
		defaultValue = strike;
	}
	else
	{
		strike = get_strikes()->get_rate(expiryDate);
		double pvExpiry = rateCurve.get_pv(expiryDate);
		double pvCurrDate = rateCurve.get_pv(currDate);

		defaultValue = strike * pvExpiry / pvCurrDate;
//		defaultValue = strike;
	}
	return defaultValue;
}


double Bond::payoff(double y, double asset, double stock, double *vega) const 
{
	KDate expiryDate = get_expiryDate();
	KDate currDate = m_tree->get_curr_date();
	double limit1 = (get_strikes())->get_rate(currDate);
	KRateCurve rateCurve = m_tree->get_index_zero();
	double payment;
	double limit2 = m_strikes2->get_rate(expiryDate);
	double returnValue;
	double optionValue = 0.0;
	double cvOptionValue = 0.0;

	SettleInfo settlement = get_settlement();					//
	KDate settlementDate = settlement.get_settlementDate();		//


	if(currDate == expiryDate)
	{
		if(expiryDate <= settlementDate)					//
		{
			returnValue=0.0;									//
		}
		else
		{
			*vega = 0.0;
			if( asset > limit2)
			{
				payment = m_bondPayment.get_Payment(currDate,&rateCurve);
			}
			else
			{
				double principalClaim = m_bondPayment.get_principal_claim(currDate);
				if( asset > limit1)
				{
					payment = (m_recovery + (1.0-m_recovery)*(asset-limit1)/(limit2-limit1) )*principalClaim
						+m_bondPayment.get_accruCoupon_recover(currDate,&rateCurve);
				}
				else
				{
					payment = m_recovery*principalClaim + 
							  m_bondPayment.get_accruCoupon_recover(currDate,&rateCurve);			
				}
//				out<<"("<<y<<", "<<principalClaim*100<<")  ";
			}
			
			// do convertion option
			if(m_cv_option != NULL)
			{
				if(stock < 1e-5)
				{
					cvOptionValue = m_cv_option->payoff(payment,currDate,vega);
				}
				else
				{
					cvOptionValue = stock*m_cv_option->payoff(payment/stock,currDate,vega);
				}
			}

//			out<<"("<<stock<<", "<<payment <<" , "<<cvOptionValue<<")  ";
			returnValue = payment+cvOptionValue;
		}
	}
	else
	{
		// add option here
		
		if(m_option == NULL)
		{
			optionValue = 0.0;
		}
		else
		{
			double accruCoupon = m_bondPayment.get_accruCoupon(currDate, true, &rateCurve);
			optionValue = m_option->payoff(y-accruCoupon, stock, currDate, vega);		
		}
		
//	out<<"("<<y<<", "<<optionValue<<")  ";
		y = y + optionValue;

		if (currDate <= settlementDate)						//
		{
			payment = 0.0;									//
		}
		else
		{
			payment = m_bondPayment.get_Payment(currDate,&rateCurve);
		}

		double principalClaim = m_bondPayment.get_principal_claim(currDate);
//		out<<"("<<y<<", "<<principalClaim*100<<")  ";
		if(asset > limit2)
		{
			returnValue = y + payment;
		}
		else
		{
			if(asset > limit1)
			{
				double couponRecover = m_bondPayment.get_couponAccruPercent();
				returnValue = y + payment*(couponRecover+(1-couponRecover)/(limit2-limit1)*(asset-limit1));
			}
			else
			{
				returnValue = m_recovery*principalClaim + m_bondPayment.get_accruCoupon_recover(currDate,&rateCurve);
				*vega = 0.0;
			}
		}

		// do convertion option
		if(m_cv_option != NULL)
		{
			if(stock < 1e-5)
			{
				cvOptionValue = m_cv_option->payoff(returnValue,currDate, vega);
			}
			else
			{
				cvOptionValue = stock*m_cv_option->payoff(returnValue/stock,currDate, vega);
			}
		}
//		out<<"("<<stock<<", "<<returnValue <<" , "<<cvOptionValue<<")  ";

		returnValue = returnValue + cvOptionValue;
	}
	return returnValue;
}



void Bond::payoff_new(double *underly,double *optValue,
					  double *underVega,double *optVega,
					  double asset,double stock, double integratUnderly)												//HY3.4v
{
	KDate expiryDate = get_expiryDate();
	KDate currDate = m_tree->get_curr_date();
	double limit1 = (get_strikes())->get_rate(currDate);
	KRateCurve rateCurve = m_tree->get_index_zero();
	double payment;
	double limit2 = m_strikes2->get_rate(expiryDate);

	double optionValue = 0.0;
	double cvOptionValue = 0.0;
	double y;

	SettleInfo settlement = get_settlement();					//
	KDate settlementDate = settlement.get_settlementDate();		//

//	if(currDate == expiryDate)
	if(currDate >= expiryDate)				//HY4.1v
	{
		*underVega = 0.0;									
		*optVega = 0.0;									

		if(expiryDate <= settlementDate)					
		{
			*underly = 0.0;									
			*optValue = 0.0;									
		}
		else
		{
			if( asset > limit2)
			{
				payment = m_bondPayment.get_Payment(currDate,&rateCurve);
			}
			else
			{
				double principalClaim = m_bondPayment.get_principal_claim(currDate);
				if( asset > limit1)
				{
					payment = (m_recovery + (1.0-m_recovery)*(asset-limit1)/(limit2-limit1) )*principalClaim
						+m_bondPayment.get_accruCoupon_recover(currDate,&rateCurve);
				}
				else
				{
					payment = m_recovery*principalClaim + 
							  m_bondPayment.get_accruCoupon_recover(currDate,&rateCurve);			
				}
//				out<<"("<<y<<", "<<principalClaim*100<<")  ";
			}
			*underly = payment;
			
			// do conversion option
			if(m_cv_option != NULL)
			{
				//recover the true stock forward price
				KTimePoint	&tp = m_tree->get_curr_tp();		//HY4.2
				double debtDriftFactor = tp[__cumDebtDrift];	//HY4.2
				stock *= debtDriftFactor;						//HY4.2

				if(stock < 1e-5)
				{
					cvOptionValue = m_cv_option->payoff(payment,currDate,optVega);	//change to integratUnderly
				}
				else
				{
					cvOptionValue = stock*m_cv_option->payoff(payment/stock,currDate,optVega);
				}
			}
			*optValue = cvOptionValue;

//			out<<"("<<stock<<", "<<payment <<" , "<<cvOptionValue<<")  ";
		}
	}
	else
	{
		if (currDate <= settlementDate)						//
		{
			payment = 0.0;									//
		}
		else
		{
			payment = m_bondPayment.get_Payment(currDate,&rateCurve);
		}

		double principalClaim = m_bondPayment.get_principal_claim(currDate);
//		out<<"("<<y<<", "<<principalClaim*100<<")  ";

		if(asset <= limit1)
		{
			*underly = m_recovery*principalClaim + m_bondPayment.get_accruCoupon_recover(currDate,&rateCurve);
			*underVega = 0.0;
			*optVega = 0.0;

			*optValue = 0.0;		//even holder's put option will be worthless
		}
		else
		{
			// add option here
			if(m_option == NULL)
			{
				optionValue = 0.0;
			}
			else
			{
				double accruCoupon = m_bondPayment.get_accruCoupon(currDate, true, &rateCurve);
//				y = *underly + *optValue - accruCoupon;
				y = integratUnderly + *optValue - accruCoupon;					//HY4.1v

				if(asset < limit2 && m_option->get_option_direction() == 1)		//HY3.4v modify strike for put btw lim1 and lim2
				{
					double alpha = (m_recovery + (1.0-m_recovery)*(asset-limit1)/(limit2-limit1));
					optionValue = m_option->payoff(y/alpha, stock, currDate, optVega);
					optionValue *= alpha;
				}
				else
				{
					optionValue = m_option->payoff(y, stock, currDate, optVega);
				}						//end, HY3.4v modify strike for put btw lim1 and lim2

				if(*optVega == 0.0 && optionValue !=0.0)
				{
					*underVega = 0.0;
				}
			}
//			out<<"("<<y<<", "<<optionValue<<")  ";
			*optValue += optionValue;

			// do conversion option
			if(m_cv_option != NULL)
			{
//				y = *underly + *optValue;
				y = integratUnderly + *optValue;		//HY4.1v

				//recover the true stock forward price
				KTimePoint	&tp = m_tree->get_curr_tp();		//HY4.2
				double debtDriftFactor = tp[__cumDebtDrift];	//HY4.2
				stock *= debtDriftFactor;						//HY4.2

				if(stock < 1e-5)
				{
					cvOptionValue = m_cv_option->payoff(y,currDate, optVega);
				}
				else
				{
					cvOptionValue = stock*m_cv_option->payoff(y/stock,currDate, optVega);
				}

				if(*optVega == 0.0 && cvOptionValue != 0.0)
				{
					*underVega = 0.0;
				}
			}
//			out<<"("<<stock<<", "<<returnValue <<" , "<<cvOptionValue<<")  ";
			*optValue += cvOptionValue;

			if(asset > limit2)		//add coupon
			{
				*underly += payment;
			}
			else
			{
				double couponRecover = m_bondPayment.get_couponAccruPercent();
				*underly += payment*(couponRecover+(1-couponRecover)/(limit2-limit1)*(asset-limit1));
			}
		}
	}
}


double Bond::get_defaultValue()
{
	KDate currDate = m_tree->get_curr_date();
	KRateCurve rateCurve = m_tree->get_index_zero();

	double principalClaim = m_bondPayment.get_principal_claim(currDate);
	double defaultValue = m_recovery*principalClaim + m_bondPayment.get_accruCoupon_recover(currDate,&rateCurve);

	return defaultValue;
}