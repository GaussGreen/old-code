/****************************************************************************/
/*      Declaration file for structure of data.                             */
/****************************************************************************/

#ifndef	_CASHFLOW_H_
#define	_CASHFLOW_H_


#include "kratecurve.h"       
#include "General/General.h"
#include "bastypes.h"        /* TCouponDates */
#include "cdate.h"           /* TDate TDateInterval */
#include "fltrate.h"  /* Float rate */
#include "kdate.h"   

#define HY_FIX   1
#define HY_FLOAT 2

class CashFlow
{
private:
    double          m_currPrincipal; /* Principal on which coupon is paid */
    double          m_amortPay;      /* Amortization payment */
    double          m_couponRate;    /* Rate for fixed or spread for float */
    TFloatRate      m_fRateInfo;    /* Floating coupon maturity */
    long            m_dayCountCpn;   /* Day count convention for cpn payment */
    long            m_dayCountAccru; /* Day count convention for accru */
    long            m_cpnRateType;   /* HY_FIX, HY_FLOAT */
	KDate           m_accruStartDate;
	KDate           m_accruEndDate;

public:
	void set_CashFlow(double        currPrincipal,
					  double        amortPay,
					  double        couponRate,
					  TDateInterval payInterval,
					  long          dayCountCpn,
					  long          dayCountAccru,
					  long          cpnRateType,
					  KDate         accruStartDate,
					  KDate         accruEndDate)
	{
		try{
			set_currPrincipal(currPrincipal);
			set_amortPay(amortPay);
			set_couponRate(couponRate);
			set_dayCountCpn(dayCountCpn);
			set_dayCountAccru(dayCountAccru);
			set_cpnRateType(cpnRateType);
			set_fRateInfo(payInterval,dayCountCpn);
			set_accruStartDate(accruStartDate);
			set_accruEndDate(accruEndDate);

		}

		catch (KException)
		{
			throw KException("Fail to set CashFlow.");
		}
	}
	void set_currPrincipal(double currPrincipal){m_currPrincipal=currPrincipal;}
	void set_amortPay(double amortPay){m_amortPay=amortPay;}
	void set_couponRate(double couponRate){m_couponRate=couponRate;}
	void set_dayCountCpn(long dayCountCpn){m_dayCountCpn=dayCountCpn;}
	void set_dayCountAccru(long dayCountAccru){m_dayCountAccru=dayCountAccru;}
	void set_cpnRateType(long cpnRateType){m_cpnRateType=cpnRateType;}
	void set_fRateInfo(TDateInterval payInterval,long dayCountCpn)
	{
		if(GtoFloatRateSet(&payInterval,
						   &payInterval,
						   dayCountCpn,
						   0,  // spotOffsetdays
						   0,  // spread
						   1,  //weight
						   &m_fRateInfo) == FAILURE)
		{
			throw KException("Fail to set up fRateInfo.");
		}
	}



	void set_accruStartDate(KDate date){m_accruStartDate = date;}
	void set_accruEndDate(KDate date){m_accruEndDate = date;}

	double get_currPrincipal() const {return m_currPrincipal;}
	double get_amortPay() const {return m_amortPay;}
	double get_couponRate() const {return m_couponRate;}
	long get_dayCountCpn() const {return m_dayCountCpn;}
	long get_dayCountAccru() const {return m_dayCountAccru;}
	long get_cpnRateType() const {return m_cpnRateType;}
	KDate get_accruStartDate() const {return m_accruStartDate;}
	KDate get_accruEndDate() const {return m_accruEndDate;}
	TFloatRate get_fRateInfo() const {return m_fRateInfo;}

	double get_couponPay(KRateCurve* zc, long interpType = GTO_FLAT_FORWARDS) const ;
	double get_couponRate(KRateCurve* zc, long interpType = GTO_FLAT_FORWARDS) const ;
	double get_cashflow(KRateCurve* zc, long interpType = GTO_FLAT_FORWARDS) const ;

	double get_accruCoupon(KDate date, KRateCurve* zc, long interpType = GTO_FLAT_FORWARDS) const ;
//	string  header();

};

std::ostream & operator<<(std::ostream & out, const CashFlow& cf );

typedef std::map<KDate,CashFlow> CashFlowMap;

class CashFlowList:public CashFlowMap, public CM::Object
{
	KDate  m_issueDate;
	KDate  m_maturityDate;
	double  m_couponAccruPercent;
	double m_issuePrice;
	KDate m_claimParDate;
	double m_accreteYield;

public:
	CashFlowList(){}
	CashFlowList(std::string bondSpec);
	CashFlowList(int num, KDate *payDates, KDate *accruStartDates, KDate *accruEndDates, 
		double* couponFlows, double* amortFlows,double couponAccPercent,KDate issueDate,KDate claimParDate, double issuePrice);
	
	void reset_rate(double x);
	void reset_rateType(long x);
	void reset_amortPay(double x);

	void shift_rate(double x);											//
	
	void set_issueDate(KDate issueDate){m_issueDate = issueDate;}
	void set_maturityDate(KDate maturityDate){m_maturityDate = maturityDate;}

	void set_resetCoupon(KDate date, double rate);

	double get_Payment(KDate date, KRateCurve* zc = NULL, long interpType=GTO_FLAT_FORWARDS) const;
	double get_accruCoupon(KDate currDate,bool accrueFlag=true, KRateCurve* zc= NULL, long interpType=GTO_FLAT_FORWARDS) const;
	double get_couponAccruPercent() const {return m_couponAccruPercent;} 
	double get_accruCoupon_recover(KDate currDate, KRateCurve* zc= NULL, long interpType=GTO_FLAT_FORWARDS) const;
	double get_principal_claim(KDate currDate) const;
	KDate get_claimParDate() const {return m_claimParDate;} 
	double get_issuePrice() const {return m_issuePrice;} 
	KVector(KDate) get_datelist() const;

};

std::ostream & operator<<(std::ostream & out, const CashFlowList& cf );



int    GtoHYInterpTypeCToI(
    char      *interpTypeChar,         /* (I) Interp type in charactors */
    long      *interpTypeLong);         /* (O) Interp type in long */

#endif

