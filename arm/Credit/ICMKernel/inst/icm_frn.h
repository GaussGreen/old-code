
#if !defined(_ICM_FRN_H__)
#define _ICM_FRN_H__


#include "ICMKernel\inst\icm_leg.h"

/*********************************************************************************/
/*! \class  ICM_Frn icm_frn.h "icm_frn.h"
 *  \author 
 *	\version 1.0
 *	\date   June 2003
 *	\brief  Describes a Floating Rate Note */
/***********************************************************************************/

class ICM_Frn : public ICM_Leg  
{
private:

public:


	ICM_Frn() : ICM_Leg()
	{}

	virtual ~ICM_Frn()
	{}

	void BitwiseCopy(const ARM_Object* srcleg);

	void Copy(const ARM_Object* srcleg);

	ARM_Object* Clone(void);

/* ****************************************************************************************************************** */
/*! \fn ICM_Frn(double Spread,
			ARM_Date Int_Accrual_Date,
			ARM_Date Maturity, 
			ARM_Date First_Period_Reference_Date, 
			ARM_IRIndex* irindex,
			double InitialRate =0.,
			qPAYMENT_PREMIUM_LEG AccruedOnDefault = qACCRUED_SETTLED,
			qPAYS_ON_DEFAULT AmortOnDef = qPays_Recovery,
			qPAYS_ON_DEFAULT IntOnDef = qPays_Recovery,
			double NotionalAmount = 100., 
			int	DayCount = K30_360, 
			int	AccruedDayCount = -1000, 
			int SettlementGap   = -1000, 
			int	stubrule = K_SHORTSTART,
			ARM_Currency* discountCcy = ARM_DEFAULT_CURRENCY,
			char* resetCalName = NULL, 
			char* payCalName = NULL,
			int nxChange = K_NX_END)

 *  \brief 	 {<B> Constructeur de bond a coupon fixe PAIR-PAIR </B>} 
	<ul>	
    	<li> - Definition </li>
 		<li> - Remboursement 100%  </li> 
 		<li> - Coupon tous EGAUX = Taux / Frequency * Notional </li>
	</ul>
  \param CouponRate Taux Coupon fixe annualisé
  \param Int_Accrual_Date Date de début du calcul des intérets (Definition de Bloomberg) 
  \param Maturity Maturité du bond (UnAdjusted) 
  \param Frequency fréquence de payment des coupons
  \param First_Period_Reference_Date Date de Payment de la premiere periode. (Unadjusted) 
         <li> Par Défaut construit a partir de la maturité avec Short Start </li>
  \param NotionalAmount Notionnel du bond 
  \param Devise devise des coupons et du notionnel
  \param AccruedDayCount Base pour le calcul des courrus. 
		<li> Si non affecté il sera attribué par défaut de la devise <\li>
  \param SettlementGap   Settlement Gap en jours
		<li> Si non affecté il sera attribué par défaut de la devise <\li> 
	*/
/* ******************************************************************************************************************** */
	
	ICM_Frn(double Spread,
			const ARM_Date& Int_Accrual_Date,
			const ARM_Date& Maturity, 
			const ARM_Date* First_Period_Reference_Date,
			ARM_IRIndex* irindex,
			double InitialRate , /* =0. */ 
			double LastIndexFixing /* = CREDIT_DEFAULT_VALUE*/ ,
			qPAYMENT_PREMIUM_LEG AccruedOnDefault /* = qACCRUED_SETTLED*/ ,
			double NotionalAmount /* = 100.*/ , 
			int	DayCount ,/* = K30_360, /*< Si non affecté il sera attribué par défaut*/
			int	AccruedDayCount, /* = -1000, /*< Si non affecté il sera attribué par défaut*/
			int SettlementGap,   /* = -1000, /*<  Si non affecté il sera attribué par défaut de la devise */ 
			int	stubrule /* = K_SHORTSTART*/ ,
			const std::string& discountCcy /*= ARM_DEFAULT_CURRENCY */ ,
			const std::string& resetCalName /* = NULL*/ , 
			const std::string& payCalName /* = NULL*/ ,
			int nxChange /* = K_NX_END*/ ) : ICM_Leg()
		{
			Init();
			
			Set(Spread, 
				Int_Accrual_Date,  
				 Maturity, 
				 First_Period_Reference_Date, //Reference Date
				 irindex,
				 InitialRate , // InitialRate	
				 LastIndexFixing ,
				 AccruedOnDefault ,
				 NotionalAmount ,
				 AccruedDayCount ,
				 SettlementGap ,
				 DayCount, // dayCount
				 stubrule, //stubRule
				 discountCcy, // discountCcy
				 resetCalName,
				 payCalName,
				 nxChange);

		}


	void Set(double Spread,
				  const ARM_Date& Int_Accrual_Date,  
				  const ARM_Date& Maturity, 
				  const ARM_Date* First_Period_Reference_Date,
				  ARM_IRIndex* irindex,
				  double InitialRate ,/* =0., // InitialRate*/ 	
				  double LastIndexFixing /* = CREDIT_DEFAULT_VALUE*/ ,
				  qPAYMENT_PREMIUM_LEG AccruedOnDefault /* = qACCRUED_SETTLED*/ ,	
				  double Notional /* = 100.*/ ,
				  int AccruedDayCount /* = -1000*/ ,
				  int SettlementGap /* = -1000*/ ,
				  int DayCount, /* = K30_360, // dayCount*/ 
				  int stubrule, /* = K_SHORTSTART, //stubRule*/ 
				  const std::string&  Devise ,/* = ARM_DEFAULT_CURRENCY, // discountCcy*/ 
				  const std::string& resetCalName /* = NULL*/ ,
				  const std::string& payCalName /* = NULL*/ ,
				  int nxChange /* = K_NX_END*/ );

	virtual double Accrued(ARM_Date & settlement);

private :

	void Set(int SettlementGap, double NotionalAmount,double redemption =100.);
	void Init(void);
	virtual void CptCashFlowValues(void);
	double CouponPeriod(int index);
	
};

#endif // !defined(AFX_ICM_Frn_H__705C12EA_BBC1_4D17_AA0E_664E45678280__INCLUDED_)
