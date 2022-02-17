/****************************************************************************/
/*                                 */
/****************************************************************************/



#include "cashflow.h"       
#include "fltrate.h"         /* GtoForwardRateSafe */
#include "ldate.h"           /* GtoDayCountFraction */
#include "date_sup.h"         /* GtoDaySubtract */
#include "datelist.h"         /* GtoNewBusDateList */
#include "kratecurve.h" 
#include "macros.h"
      

extern std::ofstream out;


double CashFlow::get_couponRate(KRateCurve* zc, long interpType) const
{
	static char routine[] = "CashFlow::get_couponRate";
	double couponRate;

	try{
		if(get_cpnRateType() == HY_FLOAT)
		{
			
			couponRate = zc->get_rate(interpType,
						  m_accruStartDate,
						 (_TFloatRate*) &m_fRateInfo);
//	out<<"r="<<couponRate <<", "<<m_couponRate;
			couponRate += m_couponRate;
		}
		else
		{
			couponRate = m_couponRate;
		}

	}

	catch (KException& e)
	{
		throw e + KException(" failed to get floating rate.",routine);
	}

	return couponRate;
}

double CashFlow::get_couponPay(KRateCurve* zc, long interpType) const
{
	static char routine[] = "CashFlow::get_couponPay";
	double couponRate, yearFract;
	try
	{
		couponRate = get_couponRate(zc, interpType);
		yearFract = GetYearDiff(m_accruStartDate, m_accruEndDate, m_dayCountCpn);
	}


	catch(KException& e)	
	{
		throw e + KException("get_couponPay failed.",routine);
	}

	return couponRate*yearFract;
}

double CashFlow::get_cashflow(KRateCurve* zc, long interpType) const 
{
	static char routine[] = "CashFlow::get_cashflow";
	double couponPay, amortPay;

	try
	{
		couponPay = get_couponPay(zc,interpType);
	}

	catch (KException& e)
	{
		throw e + KException("get_couponPay failed.",routine);
	}

	amortPay = get_amortPay();

	return couponPay+amortPay;

}


double CashFlow::get_accruCoupon(KDate date, KRateCurve* zc, long interpType) const 
{
	static char routine[] = "CashFlow::get_accruCoupon";
	double yearFract,yearFractAccru ,couponPay, accruCoupon;

	if(date < m_accruStartDate)
	{
		accruCoupon = 0.0;
	}
	else
	{
		try
		{
			KDate tempDate = ( date < m_accruEndDate ? date : m_accruEndDate );
			yearFractAccru = GetYearDiff(m_accruStartDate,date,m_dayCountAccru);
			yearFract = GetYearDiff(m_accruStartDate,m_accruEndDate,m_dayCountAccru);
			couponPay = get_couponPay(zc,interpType);
			accruCoupon = couponPay * yearFractAccru/yearFract;
		}
	
		catch(KException& e)
		{
			throw e + KException("get_accruCoupon failed.",routine);
		}
	}

	return accruCoupon;
}

std::ostream& operator<<(std::ostream& out,const CashFlow& cf )
{
	out<<cf.get_accruStartDate()<<"\t"
	   <<cf.get_accruEndDate()<<"\t"
	   <<cf.get_amortPay()<<"\t";

	return out;
}

std::ostream & operator<<(std::ostream & out, const CashFlowList& cf )
{
	out<<"PayDate\t AccStartDate\t AccEndDate\t AmortPay\n";
	CashFlowList::const_iterator iter = cf.begin();
	for(int i=0; i<cf.size(); i++)
	{
		out<<iter->first<<"\t"<<iter->second<<"\n";
		iter++;
	}

	return out;
}






CashFlowList::CashFlowList(std::string bondSpec)
{
	static char routine[] = "CashFlowList";
	char      spaceString[80]=" _,";     /* stirng seperates tokens */
	char      *token;

	if(bondSpec.empty())
	{
		goto done;
	}

	{
		std::string bondSpecCopy = bondSpec;
		char*  bondSpecPtr = &bondSpecCopy[0];
		/* get issue date */
		token = strtok(bondSpecPtr, spaceString);
		if(token == NULL)
		{
			throw KException(" Need to have issue date.", routine);
		}
		char conponStartDate[32];
		(void)strcpy(conponStartDate,token);
		
		

		/* get maturity date */
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException(" Need to have maturity date.", routine);
		}
		char maturityDate[32];
		(void)strcpy(maturityDate,token);
		
   
		/* get coupon type */
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException("Need to have coupon type.", routine);
		}

		long rateTypeLong;
		if(toupper(token[0]) == 'F')
		{
			switch(toupper(token[1]))
			{
			  case 'I':
				rateTypeLong = HY_FIX;
				break;
			  case 'L':
				rateTypeLong = HY_FLOAT;
				break;
			  default:
				throw KException("Wrong coupon type.", routine);
			}
		}
		else
		{
			throw KException("Wrong coupon type.", routine);
		}

		/* get payment frequency */
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException("Need to have payment frequency.", routine);
  
		}
		long pmtFrequency;
		pmtFrequency = atol(token);

		/* get coupon rate */
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException("Need to have coupon rate.\n", routine);
		}

		double couponRate = atof(token);

		/* get day count convention for accrued interest */
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException("Need to have day count convention.", routine);
		}
		long dayCountAccru;
		if(GtoStringToDayCountConv(token,&dayCountAccru) == FAILURE)
		{
			throw KException("Wrong Accrual day count convention.", routine);
		}

		/* get day count convention for coupon payment*/
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException("Need to have day count convention.", routine);
		}
		long dayCountCpn;
		if(GtoStringToDayCountConv(token,&dayCountCpn) == FAILURE)
		{
			throw KException("Wrong coupon day count convention.", routine);
		}

		/* get the stub type */
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException("Need to have stub type.\n", routine);
		}
		long stubType;
		if(toupper(token[0]) == 'F')
		{
		stubType = FALSE;  /* front stub */
		}
		else
		{
		stubType = TRUE;  /* end stub */
		}

		/* get business day convention */
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException("Need to have business day convention.", routine);
		}
		long busiDayConv = (long)token[0];

		/* get principal guarantee flag */
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException("Need to have principal guarantee flag.", routine);
		}
		long principalGuar = atoi(token);

		/* get coupon accrual on default flag */
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException("Need to have coupon accrual flag.", routine);
		}

		double cpnAccruFlag = atof(token);
		m_couponAccruPercent = cpnAccruFlag;

		/* get issue date */
		bool standardFlag;
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			standardFlag =TRUE;
			goto skipIssueDate;
		}
		char issueDateC[32];
		(void)strcpy(issueDateC, token);
		standardFlag = FALSE;

		/* get date claim become par */
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException("Need to have date for claim become par.\n", routine);
		}
		char claimParDateC[32];
		(void)strcpy(claimParDateC, token);
		
		/* get issue price */
		token = strtok(NULL, spaceString);
		if(token == NULL)
		{
			throw KException("Need to have issue price.\n", routine);
		}
		m_issuePrice = atof(token)/100.0;

		/* finished reading the bond string */
skipIssueDate:
		
		/* note the GtoStringToDate function cannot be called while strtok is
			active, as it also uses strtok function */

		/* convert issue and maturity dates char to TDate */
		TDate couponStartTDate;
		if(GtoStringToDate(conponStartDate,
			   &couponStartTDate) == FAILURE)
		{
			throw KException(" not a valid conpon start date.", routine);
		}
		/* set issue date */
		if(standardFlag ==TRUE)
		{
			m_issueDate =couponStartTDate;
			m_accreteYield =0.0;
			m_issuePrice = 1.0;
		}
		else
		{
			TDate temp;
			if(GtoStringToDate(issueDateC,
			   &temp) == FAILURE)
			{
				throw KException(" not a valid issue date.", routine);
			}
			m_issueDate=temp;

			/* principal claim related accrete yield */
			TDate claimParTDate;
			if(GtoStringToDate(claimParDateC,
				   &claimParTDate) == FAILURE || m_issueDate > claimParTDate)
			{
				throw KException(" not a valid claim par date.", routine);
			}
			double years=GetYearDiff(m_issueDate,claimParTDate);
			m_claimParDate = claimParTDate;
			if(years == 0.0)
			{
				m_accreteYield =0.0;
			}
			else
			{
				m_accreteYield = log(1.00/m_issuePrice)/years;
			}
		}
		
		TDate maturityTDate;
		if(GtoStringToDate(maturityDate,
			   &maturityTDate) == FAILURE)
		{
			throw KException(" not a valid maturity date.", routine);
		}
		m_maturityDate = maturityTDate;
		
		TDateInterval payInterval;
		/* determining the payment interval */
		if(pmtFrequency == 0)
		{
			if(GtoDateSubtract(m_maturityDate,
							   couponStartTDate,
							   &payInterval) == FAILURE)
			{
				throw KException("Failed to creat payInterval.", routine);
			}
		}
		else if(pmtFrequency == 365)
		{
			if(GtoMakeDateInterval(1,
							  'D',
								&payInterval) == FAILURE)
			{
				throw KException("Failed to creat payInterval.", routine);
			}
		}
		else
		{
			if(GtoFreqAndTypeToInterval(pmtFrequency,
										'N',
										&payInterval) == FAILURE)
			{
				throw KException("Failed to creat payInterval.", routine);
			}
		}

   
   

		/* get payment date list */
		TDateList   *payDateList=NULL;
		TDateList   *accruDateList =NULL;
		payDateList = GtoNewBusDateList(couponStartTDate,
				m_maturityDate,
				&payInterval,
				stubType,
				busiDayConv,
				"NONE");
		if(payDateList == NULL)
		{
			throw KException("Failed to creat payDateList.", routine);
		}

		/* get accrual date list */
		
		accruDateList = GtoNewBusDateList(couponStartTDate,
				m_maturityDate,
				&payInterval,
				stubType,
				GTO_BAD_DAY_NONE,
				(char*)"NONE");
		if(accruDateList == NULL)
		{
			throw KException("Failed to creat accruDateList.", routine);
		}

		long nCashFlows = payDateList->fNumItems -1;  /* no payment on couponStartDate date */
		CashFlow cashflow;
		
		double principal = 0.0;
		for(int i=0; i<nCashFlows; i++)
		{
			if(i == nCashFlows-1 && principalGuar != -1)
			{
				principal = get_principal_claim(m_maturityDate);
			}
			cashflow.set_CashFlow(1.0,
								  principal,
								  couponRate,
								  payInterval,
								  dayCountCpn,
								  dayCountAccru,
								  rateTypeLong,
								  accruDateList->fArray[i],
								  accruDateList->fArray[i+1]);

			insert(value_type(KDate(payDateList->fArray[i+1]),cashflow));
    
		}

		GtoFreeDateList(payDateList);
		GtoFreeDateList(accruDateList);
	}
	
done:
	;
	
}

CashFlowList::CashFlowList(int num, KDate *payDates, KDate *accruStartDates, KDate *accruEndDates, 
		double* couponFlows, double* amortFlows,double couponAccPercent,KDate issueDate, KDate claimParDate,double issuePrice)
{
	static char routine[] = "CashFlowList";
	double currentPrincipal = 1.0;

	/* payInterval does not mattwer */
	TDateInterval payInterval;
	if(GtoFreqAndTypeToInterval(2,
								'N',
								&payInterval) == FAILURE)
	{
		throw KException("Failed to creat payInterval.", routine);
	}

	double yearFract;
	CashFlow cashflow;
	double totalAmort = 0;

	for(int i=0; i<num; i++)
	{
		yearFract = GetYearDiff(accruStartDates[i],accruEndDates[i],GTO_ACT_365);	

//		printf("year = %f",yearFract);
		cashflow.set_CashFlow(currentPrincipal,
							  amortFlows[i],
							  (yearFract<0.00000001)? 0:couponFlows[i]/yearFract,
							  payInterval,
							  GTO_ACT_365,
							  GTO_ACT_365,
							  HY_FIX,
							  accruStartDates[i],
							  accruEndDates[i]);
		
		currentPrincipal -= amortFlows[i];
		totalAmort += amortFlows[i];
		insert(value_type(KDate(payDates[i]),cashflow));

	}

	if(totalAmort != 1)
	{
//		throw KException("Total amortization is not 1", routine);
	}

	m_issueDate = issueDate;
	m_issuePrice = issuePrice;
	m_couponAccruPercent = couponAccPercent;
	m_claimParDate = claimParDate;
	
	double years=GetYearDiff(m_issueDate,claimParDate);
	if(years == 0.0)
	{
		m_accreteYield =0.0;
	}
	else
	{
		m_accreteYield = log(1.00/m_issuePrice)/years;
	}
}



void CashFlowList::reset_rate(double x)
{
	iterator iter;

	for(iter = begin(); iter != end(); iter++)
	{
		iter->second.set_couponRate(x);
	}
}


void CashFlowList::reset_rateType(long x)
{
	iterator iter;

	for(iter = begin(); iter != end(); iter++)
	{
		iter->second.set_cpnRateType(x);
	}
}


void CashFlowList::reset_amortPay(double x)
{
	iterator iter;

	for(iter = begin(); iter != end(); iter++)
	{
		iter->second.set_amortPay(x);
	}
}


void CashFlowList::shift_rate(double x)							//
{																//
	double temp;
	iterator iter;												//	
		
	for(iter = begin(); iter != end(); iter++)					//
	{															//
		temp = iter->second.get_couponRate() + x;				//
		iter->second.set_couponRate(temp);						//
	}															//
}																//	


void CashFlowList::set_resetCoupon(KDate date, double rate)
{
	static char routine[] = "CashFlowList::set_resetCoupon";

	iterator iter = lower_bound(date);
	
	if(iter == end())
	{
		throw KException("end of the flow reached.\n", routine);
	}
	else
	{
		while(true)
		{
			if((iter->second.get_cpnRateType() != HY_FIX) && (iter->second.get_accruStartDate()<date))
			{
				iter->second.set_couponRate(rate);
				iter->second.set_cpnRateType(HY_FIX);
			}
			if(iter == begin())
			{
				break;
			}
			iter--;
		}
	}
}

double CashFlowList::get_Payment(KDate currDate, KRateCurve* zc, long interpType) const
{
	static char routine[] = "CashFlowList::get_Payment";
	const_iterator currCashFlowIter;
	double payment;

	try
	{
		currCashFlowIter = find(currDate);
		if(currCashFlowIter == end())  // cashflow does not exist
		{
			payment = 0.0;
		}
		else
		{
			if(zc == NULL && currCashFlowIter->second.get_cpnRateType() == HY_FLOAT)
			{
				throw KException("Float rate needs zero curve.\n");
			}
			else
			{
				payment = currCashFlowIter->second.get_cashflow(zc, interpType);
			}

		}
	}

	catch (KException& e)
	{
		throw e + KException("get_payment failed.", routine);
	}

	return payment;
}

double CashFlowList::get_accruCoupon(KDate currDate,bool accrueFlag, KRateCurve* zc, long interpType) const
{
	static char routine[] = "CashFlowList::get_accruCoupon";
	const_iterator currCashFlowIter;
	double accruCoupon;

	try
	{
		if(accrueFlag == true)		//by default
		{
			currCashFlowIter = upper_bound(currDate);		//
		}
		else
		{
			currCashFlowIter = lower_bound(currDate);		
		}

		if(currCashFlowIter == end())
		{
			accruCoupon = 0.0;
		}
		else
		{
			accruCoupon = currCashFlowIter->second.get_accruCoupon(currDate, zc, interpType);
		}
	}

	catch (KException& e)
	{
		throw e + KException("get_accrual failed.", routine);
	}

	return accruCoupon;
}


double CashFlowList::get_accruCoupon_recover(KDate currDate, KRateCurve* zc, long interpType) const
{
	static char routine[] = "CashFlowList::get_accruCoupon_recover";

	double accruCoupon;

	try
	{
		accruCoupon = get_accruCoupon(currDate, false, zc, interpType)*m_couponAccruPercent;
	}

	catch (KException& e)
	{
		throw e + KException("get_accrual failed.", routine);
	}

	return accruCoupon;
}

double CashFlowList::get_principal_claim(KDate currDate) const
{
	static char routine[] = "CashFlowList::get_principal_claim";
	double years, prinClaim;

	years = MAX(0,GetYearDiff(m_issueDate, currDate));
	prinClaim = MIN(1.00,m_issuePrice*exp(m_accreteYield*years));

	return prinClaim;
}


KVector(KDate) CashFlowList::get_datelist() const 
{
	KVector(KDate) dates;
	const_iterator iter;

	for (iter = begin(); iter != end(); iter++)
	{
		dates.push_back( iter->first);
	}
	return dates;
}











int    GtoHYInterpTypeCToI(
    char      *interpTypeChar,         /* (I) Interp type in charactors */
    long      *interpTypeLong)         /* (O) Interp type in long */
{
    static  char routine[] = "GtoHYInterpTypeCToI";
    int  status = FAILURE;       /* Until proven SUCCESS */

    switch(toupper(interpTypeChar[0]))
    {
      case 'L':
        *interpTypeLong =GTO_LINEAR_INTERP;  /* 0 */
        break;

      case 'F':
        *interpTypeLong =GTO_FLAT_FORWARDS;  /* 124 */
        break;

      case 'S':
        *interpTypeLong =GTO_SPLINE_INTERP;  /* 1 */
        break;

      case 'C':
        *interpTypeLong =GTO_CONST_SPOT_VOL_INTERP;  /* 2 */
        break;
      default:
        KException("Invalid interpolation type. \n", 
                  routine);
        goto error;
    }
    status = SUCCESS;

  error:
    if(status == FAILURE)
    {
        KException("Failed. \n", routine);
    }
    return(status);
}
