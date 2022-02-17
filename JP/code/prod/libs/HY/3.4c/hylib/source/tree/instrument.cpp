
/**********************************************************
*
*	instrument.cpp
*
***********************************************************/
#include "kplatform.h"
#include "instruments.h"
#include "macros.h"

#ifdef _DEBUG
// extern std::ofstream out;
// extern std::ofstream out("c:\\yy.dat");
#endif

void Instrument::update()
{
	DTimeSlice tsStock = m_tree->get_stock_price();
	DTimeSlice tsAsset = (m_tree->*m_underFunc)();
	int	dim  = m_tree->get_curr_limit()[0];
	KDate expiryDate = get_expiryDate();
	KDate currDate = m_tree->get_curr_date();

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

	if(currDate == expiryDate)
	{
//		out << "\nDIM="<<dim<<"\n"<<std::endl;
		for (i = -dim; i <= dim; i++)
		{
			payoff_new(&m_ts[i], &m_ts_opt[i], &m_ts_vega[i], &m_ts_opt_vega[i], tsAsset[i],tsStock[i]);	//HY3.4v

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

//		m_tree->ddv(m_ts_cvOpt_vega,m_ts_cvOpt, true, true);		//HY3.4v
//		m_tree->dev(m_ts_cvOpt, true, true);						//HY3.4v

//		out<<"\nDate: "<<GtoFormatDate(currDate)<<"\t dim="<<dim<<"\n"<<std::endl;
//		out<<"\n dim : "<<dim<<"\n"<<m_ts<<std::endl;
//		out<<"\n dim : "<<dim<<"\n"<<tsAsset<<std::endl;

//		out << "\nDIM="<<dim<<"\n"<<std::endl;
		for (i = -dim; i <= dim; i++)
		{
			payoff_new(&m_ts[i], &m_ts_opt[i], &m_ts_vega[i], &m_ts_opt_vega[i],
						tsAsset[i],tsStock[i]);										//HY3.4v

//			m_ts[i] = payoff(m_ts[i],tsAsset[i],tsStock[i],&m_ts_vega[i]);
//			m_ts_vega[i] = vegaSet(m_ts_vega[i],tsAsset[i]);
		}
		
//		out<<"\n dimA : "<<dim<<"\n"<<m_ts<<std::endl;
//		out<<" asset : "<<"\n"<<tsAsset<<std::endl;
		//get stock for the next run
	//	DTimeSlice stock = m_tree->get_stock_price();

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
	

double  Instrument::get_delta_value()
{
	DTimeSlice& stock = m_tree->get_stock_price();

	double delta= (m_ts[-1]-m_ts[0])/(stock[-1]-stock[0])+
				  (m_ts[1]-m_ts[0])/(stock[1]-stock[0]) -
				  (m_ts[-1]-m_ts[1])/(stock[-1]-stock[1]);
	return delta;
}


double  Instrument::get_option_delta_value()				//HY3.4v
{
	DTimeSlice& stock = m_tree->get_stock_price();

	double delta= (m_ts_opt[-1]-m_ts_opt[0])/(stock[-1]-stock[0])+
				  (m_ts_opt[1]-m_ts_opt[0])/(stock[1]-stock[0]) -
				  (m_ts_opt[-1]-m_ts_opt[1])/(stock[-1]-stock[1]);
	return delta;
}


double  Instrument::get_gamma_value()
{
	DTimeSlice& stock = m_tree->get_stock_price();

	double gamma= ((m_ts[-1]-m_ts[0])/(stock[-1]-stock[0])-
				  (m_ts[0]-m_ts[1])/(stock[0]-stock[1]))/((stock[-1]-stock[1])/2);
	return gamma;
}


double  Instrument::get_option_gamma_value()				//HY3.4v
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


void DefaultSwap::payoff_new(double *underly,double *optValue,
							 double *underVega,double *optVega,
							 double asset,double stock) 
{
	KDate expiryDate = get_expiryDate();
	KDate currDate = m_tree->get_curr_date();
	double limit1 = (get_strikes())->get_rate(currDate);
	double limit2 = m_strikes2->get_rate(expiryDate);
	double payment;												//
	double y;

	SettleInfo settlement = get_settlement();					//
	KDate settlementDate = settlement.get_settlementDate();		//
			
	if(currDate == expiryDate)
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
				y = *underly + *optValue + accruCoupon;
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
							   double asset, double stock) 
{
	*optValue = 0.0;
	*optVega = 0.0;
	*underly = payoff(*underly, asset, stock, underVega);
}	


double CallOptionT::payoff(double y,double asset, double sstock, double *vega) const 
{
	KDate expiryDate = get_expiryDate();
	KDate currDate = m_tree->get_curr_date();
	double strike = get_strikes()->get_rate(currDate);
	double optionValue =0.0;

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
							 double asset, double stock) 
{
	*optValue = 0.0;
	*optVega = 0.0;
	*underly = payoff(*underly, asset, stock, underVega);
}	


double PutOptionT::payoff(double y,double asset, double sstock, double *vega) const 
{
	KDate expiryDate = get_expiryDate();
	KDate currDate = m_tree->get_curr_date();
	double strike = get_strikes()->get_rate(currDate);
	double optionValue =0.0;

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


void PutOptionT::payoff_new(double *underly,double *optValue,
							double *underVega,double *optVega,
							double asset,double stock) 
{
	*optValue = 0.0;
	*optVega = 0.0;
	*underly = payoff(*underly, asset, stock, underVega);
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
					  double asset,double stock)												//HY3.4v
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

	if(currDate == expiryDate)
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
			
			// do convertion option
			if(m_cv_option != NULL)
			{
				if(stock < 1e-5)
				{
					cvOptionValue = m_cv_option->payoff(payment,currDate,optVega);
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
				y = *underly + *optValue - accruCoupon;

				if(asset < limit2 && m_option->get_option_direction() == 1)		//HY3.4v modify strike for put btw lim1 and lim2
				{
					double alpha = (m_recovery + (1.0-m_recovery)*(asset-limit1)/(limit2-limit1));
					optionValue = m_option->payoff(y/alpha, stock, currDate, optVega);
					optionValue *= alpha;
				}
				else
				{
					optionValue = m_option->payoff(y, stock, currDate, optVega);
				}

				if(*optVega == 0.0 && optionValue !=0.0)
				{
					*underVega = 0.0;
				}							//end, HY3.4v modify strike for put btw lim1 and lim2
			}
//			out<<"("<<y<<", "<<optionValue<<")  ";
			*optValue += optionValue;

			// do convertion option
			if(m_cv_option != NULL)
			{
				y = *underly + *optValue;
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
