/****************************************************************************/
/*      Declaration file for structure of data.                             */
/****************************************************************************/

#ifndef	_INSTRUMENT_H_
#define	_INSTRUMENT_H_

#include "kplatform.h"  

#include "kratecurve.h"       
#include "tree_asset.h"
#include "General/General.h"
#include "cashflow.h"   /* CashFlowList */
#include "kplatdep.h"   /* KVector */
#include "strikes.h"
#include "settleInfo.h"			//


typedef DTimeSlice& (AssetTree::*TUnderFunc)(void);

class Instrument:public CM::Object
{
	KDate   m_expiryDate;
	KRateCurve*	 m_strikes;
	SettleInfo m_settlement;			//

protected:

	AssetTree*   m_tree;
	DTimeSlice   m_ts;
	DTimeSlice   m_ts_opt;						//HY3.4v
	DTimeSlice   m_ts_vega;
	DTimeSlice   m_ts_opt_vega;					//HY3.4v

//	DTimeSlice& (AssetTree::*m_underFunc)();
	TUnderFunc  m_underFunc;
public:
	Instrument(){} 

	virtual ~Instrument(void){}

	void  set_expiryDate(KDate expiry){ m_expiryDate = expiry;} 
	KDate get_expiryDate() const { return m_expiryDate;}

	void  set_settlement(SettleInfo Settlement){ m_settlement = Settlement;}		//
	SettleInfo get_settlement() const { return m_settlement;}						//

	void	set_strikes(KRateCurve* strikes){m_strikes = strikes;}
	KRateCurve* get_strikes() const { return m_strikes;}

	virtual double payoff(double y, double asset, double stock, double *vega) const =0;		
	virtual void payoff_new(double *underly,double *optValue, double *underVega,
							double *optVega, double asset,double stock) = 0;	//HY3.4v		

	double vegaSet(double y, double asset);

	virtual KVector(KDate) get_datelist() const { KVector(KDate) dates; dates.push_back(m_expiryDate);return dates;}
	
	virtual CashFlowList get_cashflowlist() const {return CashFlowList();}
	virtual CashFlowList* get_cashflowlistRef() {return &CashFlowList();}					//

	virtual void  set_tree(AssetTree* tree) 
	{
		m_tree=tree;
		m_underFunc = &AssetTree::get_asset_value;
	}

//	void	resize(){ m_ts.resize(m_tree);}
	void update();
	double  get_spot_ts_value(){return m_ts[0];}
	double  get_option_ts_value(){return m_ts_opt[0];}					//HY3.4v
//	double  get_cvOption_ts_value(){return m_ts_cvOpt[0];}				//HY3.4v

	double  get_delta_value();
	double  get_option_delta_value();									//HY3.4v
//	double  get_cvOption_delta_value();									//HY3.4v

	double  get_gamma_value();
	double  get_option_gamma_value();									//HY3.4v
//	double  get_cvOption_gamma_value();									//HY3.4v

	double  get_vega_value(){return m_ts_vega[0];}
	double  get_option_vega_value(){return m_ts_opt_vega[0];}			//HY3.4v
//	double  get_cvOption_vega_value(){return m_ts_cvOpt_vega[0];}			//HY3.4v

	virtual double get_accrual(KDate date,KRateCurve* zc, bool accrueFlag = true){return 0.0;}

	virtual Instrument* clone() const = 0;										//

};


class DefaultSwap:public Instrument
{
//	KRateCurve*	 m_strike1;
	KRateCurve*  m_strikes2;
	double       m_recovery;
	CashFlowList m_feePayment;
	OptionContext*  m_option;

public:
	DefaultSwap(KRateCurve* strikes1, KRateCurve* strikes2, 
				double recovery, std::string feeSpec = "",
				OptionContext*  option = NULL ):
				 m_strikes2(strikes2), m_recovery(recovery) ,m_feePayment(feeSpec),m_option(option)
				{  set_strikes(strikes1); set_expiryDate(m_strikes2->endDate());}

	DefaultSwap(KRateCurve* strikes1, KRateCurve* strikes2, 
				double recovery, CashFlowList feeSpec,
				OptionContext* option = NULL):
				 m_strikes2(strikes2), m_recovery(recovery) ,m_feePayment(feeSpec),m_option(option)
				{  set_strikes(strikes1); set_expiryDate(m_strikes2->endDate());}

	~DefaultSwap(void){}
	double payoff(double y, double asset, double stock, double *vega) const;
	void payoff_new(double *underly,double *optValue, double *underVega,
							double *optVega, double asset,double stock);		//HY3.4v		

	CashFlowList get_cashflowlist() const {return m_feePayment;}
	CashFlowList* get_cashflowlistRef() {return &m_feePayment;}					//
	
	KRateCurve* get_strike2(){return m_strikes2;}

	KVector(KDate) get_datelist() const 
	{
		KVector(KDate) dates = m_feePayment.get_datelist();
		if(m_option != NULL)
		{
			KVector(KDate) datesStrike = m_option->get_datelist();
			dates.insert(dates.end(),datesStrike.begin(),datesStrike.end());
		}
		dates.push_back(get_expiryDate());

		return dates;
	}

	OptionContext* get_option_pt(){return m_option;}
	double get_recoveryRate(){return m_recovery;}
	double get_accrual(KDate date, KRateCurve* zc,bool accrueFlag = true)
	{
		return m_feePayment.get_accruCoupon(date,accrueFlag, zc);
	}

	Instrument* clone() const{return new DefaultSwap(*this);}					//

};



class CallOptionT:public Instrument
{
//	KRateCurve* m_strike;
	bool m_isamerican;

public:
	CallOptionT(KRateCurve* strikes, bool isAmerican)
	{
		set_strikes(strikes);
		set_expiryDate(get_strikes()->endDate());
		m_isamerican = isAmerican;
	}
	CallOptionT(KRateCurve* strikes, std::string isAmerican = "E" )
	{
		set_strikes(strikes); 
		set_expiryDate( (get_strikes())->endDate());
		if(toupper( isAmerican[0] ) == 'E')
		{
			m_isamerican = false;
		}
		else
		{
			m_isamerican = true;
		}
//		m_underFunc = underFunc;
	}		
	
	bool get_isamerican(){return m_isamerican;}
	void  set_tree(AssetTree* tree)
	{
		m_tree=tree;
		m_underFunc = &AssetTree::get_stock_price;
	}
	
	double payoff(double y,double under, double stock, double *vega) const;
	void payoff_new(double *underly,double *optValue, double *underVega,
							double *optVega, double asset,double stock);		//HY3.4v		

	Instrument* clone() const{return new CallOptionT(*this);}							//
};


class PutOptionT:public Instrument
{
//	KRateCurve* m_strike;
	bool m_isamerican;

public:
	PutOptionT(KRateCurve* strikes, bool isAmerican)
	{
		set_strikes(strikes);
		set_expiryDate(get_strikes()->endDate());
		m_isamerican = isAmerican;
	}

	PutOptionT(KRateCurve* strikes, std::string isAmerican = "E")
	{
		set_strikes(strikes);
		set_expiryDate(get_strikes()->endDate());
		if( toupper( isAmerican[0] ) == 'E')
		{
			m_isamerican = false;
		}
		else
		{
			m_isamerican = true;
		}
	}		
	
	bool get_isamerican(){return m_isamerican;}
	void  set_tree(AssetTree* tree)
	{
		m_tree=tree;
		m_underFunc = &AssetTree::get_stock_price;
	}
	
	double payoff(double y,double under, double stock, double *vega) const;
	void payoff_new(double *underly,double *optValue, double *underVega,
							double *optVega, double asset,double stock);		//HY3.4v		

	Instrument* clone() const{return new PutOptionT(*this);}								//
};


class DSAverageLife: public Instrument
{
	

public:
	DSAverageLife(KRateCurve* strikes){ set_strikes(strikes); set_expiryDate(get_strikes()->endDate());}

	double payoff(double y, double asset, double stock, double *vega) const;
	void payoff_new(double *underly,double *optValue, double *underVega,
							double *optVega, double asset,double stock);		//HY3.4v		

	Instrument* clone() const{return new DSAverageLife(*this);}								//
};

class Bond:public Instrument
{
	KRateCurve*  m_strikes2;
	double       m_recovery;
	CashFlowList m_bondPayment;
	OptionContext*  m_option;
	OptionContext*  m_cv_option;

public:
	Bond( KRateCurve* strikes1, KRateCurve* strikes2, 
		  double recovery, std::string bondSpec, double resetRate, 
		  OptionContext* option = NULL,
		  OptionContext* cvOption = NULL):
				 m_strikes2(strikes2), m_recovery(recovery) ,m_bondPayment(bondSpec),
				 m_option(option),m_cv_option(cvOption)
				{  set_strikes(strikes1); set_expiryDate(m_strikes2->endDate()); 
				   m_bondPayment.set_resetCoupon(m_strikes2->valueDate(),resetRate); }
	Bond(KRateCurve* strikes1, KRateCurve* strikes2, 
		 double recovery, CashFlowList bondSpec, double resetRate =0, 
		 OptionContext* option = NULL,
		 OptionContext* cvOption = NULL):
				 m_strikes2(strikes2), m_recovery(recovery) ,m_bondPayment(bondSpec),
				 m_option(option),m_cv_option(cvOption)
				{  set_strikes(strikes1); set_expiryDate(m_strikes2->endDate()); 
				   m_bondPayment.set_resetCoupon(m_strikes2->valueDate(),resetRate); }

	~Bond(void){}

	double payoff(double y, double asset, double stock, double *vega) const;
	void payoff_new(double *underly,double *optValue, double *underVega,
							double *optVega, double asset,double stock);		//HY3.4v		

	CashFlowList get_cashflowlist() const {return m_bondPayment;}
	CashFlowList* get_cashflowlistRef() {return &m_bondPayment;}							//
	
	KVector(KDate) get_datelist() const 
	{
		KVector(KDate) dates = m_bondPayment.get_datelist();
		if(m_option != NULL)
		{
			KVector(KDate) datesStrike = m_option->get_datelist();
			dates.insert(dates.end(),datesStrike.begin(),datesStrike.end());
		}
		dates.push_back(get_expiryDate());

		return dates;
	}

	KRateCurve* get_strike2(){return m_strikes2;}

	OptionContext*  get_option_pt(){return m_option;}
	OptionContext*  get_option_cv_pt(){return m_cv_option;}
	double get_recoveryRate(){return m_recovery;}

	double get_accrual(KDate date, KRateCurve* zc,bool accrueFlag = true)
	{
		return m_bondPayment.get_accruCoupon(date,accrueFlag, zc);
	}

	Instrument* clone() const{return new Bond(*this);}											//
};

#endif

