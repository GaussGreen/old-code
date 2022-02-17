#include "ARMKernel\glob\firsttoinc.h"
#include "ICMKernel\inst\icm_cmcds.h"

// ----------------------------
//	Copy of members data
// ----------------------------
void ICM_Cmcds::BitwiseCopy(const ARM_Object* srcCmcds)
{
    ICM_Cmcds* cmcds = (ICM_Cmcds *) srcCmcds;

	itsCapLevel = cmcds->itsCapLevel ;
	itsFloorLevel = cmcds->itsFloorLevel ;
}

// -------------
//	Copy Method 
// -------------
void ICM_Cmcds::Copy(const ARM_Object* srcCmcds)
{
	ICM_Cds::Copy(srcCmcds);
    BitwiseCopy(srcCmcds);
}

// --------------
//	Clone Method
// --------------
ARM_Object* ICM_Cmcds::Clone(void)
{
     ICM_Cmcds* theClone = new ICM_Cmcds();

     theClone->Copy(this);
 
     return(theClone);
}

ICM_Cmcds::ICM_Cmcds(const ARM_Date& EffectiveDate,
					 const ARM_Date& MaturityDate,
					 const ARM_Date*FirstPeriodReferenceDate,
					 const ARM_Date* FstCpnEffDate,
					 const ARM_Date& ProtectionStartDate,
					 const ARM_Date& ProtectionEndDate,
					 const double& CapLevel ,
					 const double& FloorLevel ,
					 const double& FixedRate,
					 ICM_Credit_Index* Index,
					 const int& FrequencyFeeLeg,
					 const int& DayCountBasis,
					 const double& FixedPayerAmount, 
					 const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
					 const std::string& ccy, 
					 const double& FloatingPayerAmount,
					 const int& stubrule,
					 const double& CreditLag,
					 const int& FrequencyDefLeg,
					 const int& intRule,
					 const bool& includematurity,
					 const int& adjStartDate,
					 const std::string& payCalName,
					 const qSecurity_TYPE& cdstype)
{
	Init();

	Set(EffectiveDate,
		MaturityDate,
		FirstPeriodReferenceDate,
		FstCpnEffDate,
		ProtectionStartDate,
		ProtectionEndDate,
		CapLevel,
		FloorLevel,
		FixedRate,
		Index,
		FrequencyFeeLeg,
		DayCountBasis,
		FixedPayerAmount, 
		AccruedOnDefault,
		ccy, 
		FloatingPayerAmount,
		stubrule,
		CreditLag,
		FrequencyDefLeg,
		intRule,
		includematurity,
		adjStartDate,
		payCalName,
		cdstype);
}


// -------------------------------
//	Set Method for CMCDS Product
// -------------------------------

void ICM_Cmcds::Set(const ARM_Date& EffectiveDate,
					const ARM_Date& MaturityDate,
					const ARM_Date* FirstPeriodReferenceDate,
					const ARM_Date* FstCpnEffDate,
					const ARM_Date& ProtectionStartDate,
					const ARM_Date& ProtectionEndDate,
					const double& CapLevel ,
					const double& FloorLevel ,
					const double& FixedRate,
					ICM_Credit_Index* Index,
					const int& FrequencyFeeLeg,
					const int& DayCountBasis,
					const double& FixedPayerAmount, 
					const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
					const std::string&Ccy, 
					const double& FloatingPayerAmount,
					const int& stubrule,
					const double& CreditLag,
					const int& FrequencyDefLeg,
					const int& intRule,
					const bool& includematurity,
					const int& adjStartDate,
					const std::string& payCalName,
					const qSecurity_TYPE& cdstype)
{
	ICM_Cds::Set(EffectiveDate, MaturityDate, FirstPeriodReferenceDate,FstCpnEffDate,ProtectionStartDate,
				 ProtectionEndDate, FixedRate, 
				 FixedPayerAmount, FloatingPayerAmount,0,0,	// notionals 
				 Index,
				 FrequencyFeeLeg, DayCountBasis,
				  AccruedOnDefault, Ccy, 
				 stubrule, CreditLag,
				 FrequencyDefLeg, intRule, includematurity,
				 adjStartDate, payCalName, cdstype,
					 ISSUER_UNDEFINE, // const string& name /* = ISSUER_UNDEFINE*/ ,
					 CREDIT_DEFAULT_VALUE // const double& Binary /* = CREDIT_DEFAULT_VALUE*/ );
				 ) ;
	
	itsCapLevel = CapLevel ;
	itsFloorLevel = FloorLevel ;
}

