
#if !defined(_ICM_CLN_H__)
#define _ICM_CLN_H__

#include "ICMKernel\inst\icm_leg.h"

/*********************************************************************************/
/*! \class  ICM_Cln icm_cln.h "icm_cln.h"
 *  \author damien pouponneau
 *	\version 1.0
 *	\date   January 2006
 *	\brief  Definition of a <B> Credit Linked Note </B>
	The definition of the members is closed to the ISDA defintion. */
/***********************************************************************************/

class ICM_Cln : public ICM_Leg
{
private:

	double itsRedempNotional;

public:

	ICM_Cln() : ICM_Leg() {}
	virtual ~ICM_Cln() {}

	void BitwiseCopy(const ARM_Object* srcleg);
	void Copy(const ARM_Object* srcleg);
	ARM_Object* Clone(void);

	 
	ICM_Cln(const ARM_Date& startDate, 
            const ARM_Date& endDate, 
			const ARM_Date* refDate,
   			const ARM_Date* FstCpnEffDate,
			ARM_IRIndex* irIndex,
			const double& spread /* = 0.0*/ ,
			const double& notional /* = 1.e7*/ ,
			const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/ ,
			const int&				    AccruedDayCount /* = KACTUAL_365*/ ,
			const int& rcvOrPay /* = K_PAY*/ , 
			const int& dayCount/* = K30_360*/ , 
            const int& decompFreq /* = K_COMP_PROP*/ ,
            const int& stubRule /* = K_SHORTSTART*/ ,
			const int& resetgap /* = 10000*/ ,
			const std::string& resetCalName /* = NULL*/ , 
            const std::string& discountCcy /* = ARM_DEFAULT_CURRENCY*/ ,
            const std::string& payCalName /* = NULL*/ ,
            const int& nxChange /* = K_NX_NONE*/ ,
			const bool& includematurity /* = EXCLUDE_MATURITY*/ ,
			const int& adjStartDate /* = K_ADJUSTED*/ ,
			const qCredit_Leg_Type& LegType /* = qRunning_Leg*/ ,
			const double& Binary /* = CREDIT_DEFAULT_VALUE*/ ,
			const string& name /* = ISSUER_UNDEFINE*/ );

	void Set(const ARM_Date& startDate, 
            const ARM_Date& endDate, 
			const ARM_Date* refDate,
   			const ARM_Date* FstCpnEffDate,
			ARM_IRIndex* irIndex,
			const double& spread /* = 0.0*/,
			const double& notional /* = 1.e7*/,
			const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/,
			const int&				    AccruedDayCount /* = KACTUAL_365*/,
			const int& rcvOrPay /* = K_PAY*/, 
			const int& dayCount/* = K30_360*/, 
            const int& decompFreq /* = K_COMP_PROP*/,
            const int& stubRule /* = K_SHORTSTART*/,
			const int& resetgap /* = 10000*/,
			const std::string& resetCalName /* = NULL*/, 
            const std::string& discountCcy /* = ARM_DEFAULT_CURRENCY*/,
            const std::string& payCalName /* = NULL*/,
            const int& nxChange /* = K_NX_NONE*/,
			const bool& includematurity /* = EXCLUDE_MATURITY*/,
			const int& adjStartDate /* = K_ADJUSTED*/,
			const qCredit_Leg_Type& LegType /* = qRunning_Leg*/,
			const double& Binary /* = CREDIT_DEFAULT_VALUE*/,
			const string& name /* = ISSUER_UNDEFINE*/);

	inline double GetRedempNotional() { return itsRedempNotional;}

private :

	void Init(void);

};

#endif // !defined(AFX_ICM_Cln_H__705C12EA_BBC1_4D17_AA0E_664E45678280__INCLUDED_)
