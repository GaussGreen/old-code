#ifndef _CMCDS_H
#define _CMCDS_H
 
/********************************************************************************/
/*! \file icm_cmcds.h
 * 
 *  \brief Describes a Constant Maturity CDS
 *  \author : Fakher Ben Atig
 *	\version 1.0
 *	\date   January 2005 */
/*
 *********************************************************************************/

#include "ICMKernel\inst\icm_cds.h"

class ICM_Cmcds : public ICM_Cds 
{
   private:

	   double itsCapLevel ;
	   double itsFloorLevel ;

   public:

      
// -------------------------------------------------------------------------
// Constructors For CMCDS
// -------------------------------------------------------------------------
	
	   void Init()
	   {
		   SetName(ICM_CAPFLOORCMCDS);
		   itsCapLevel = -1. ;
		   itsFloorLevel = -1. ;
	   };

	   ICM_Cmcds(void){};

	   ICM_Cmcds(const ARM_Date& EffectiveDate,
				 const ARM_Date& MaturityDate,
				 const ARM_Date* FirstPeriodReferenceDate ,
 				 const ARM_Date* FstCpnEffDate,
				 const ARM_Date& ProtectionStartDate,
				 const ARM_Date& ProtectionEndDate,
				 const double& CapLevel ,
				 const double& FloorLevel ,
				 const double& FixedRate,
				 ICM_Credit_Index* Index,
				 const int& FrequencyFeeLeg /* = K_QUARTERLY*/ ,
				 const int& DayCountBasis /* = KACTUAL_360*/ ,
				 const double& FixedPayerAmount /* = 100.*/ ,
				 const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/ ,
				 const std::string& ccy /* = ARM_DEFAULT_CURRENCY*/ , 
				 const double& FloatingPayerAmount /* = 0.0*/ ,
				 const int& stubrule ,/* , /* 1 <-> SHORTSTART */
				 const double& CreditLag /* = DEFAULT_CREDIT_LAG*/ ,
				 const int& FrequencyDefLeg /* = -1*/ ,
				 const int& intRule /* = K_ADJUSTED*/ ,
				 const bool& includematurity /* = false*/ ,
				 const int& adjStartDate /* = 1*/ ,
				 const std::string& payCalName /* = NULL*/ ,
				 const qSecurity_TYPE& cdstype /* = qCDS_INDEX*/ );

	   void Set(const ARM_Date& EffectiveDate,
				const ARM_Date& MaturityDate,
				const ARM_Date* FirstPeriodReferenceDate ,
				const ARM_Date* FstCpnEffDate,
				const ARM_Date& ProtectionStartDate,
				const ARM_Date& ProtectionEndDate,
				const double& CapLevel,
				const double& FloorLevel,
				const double& FixedRate,
				ICM_Credit_Index* Index,
				const int& FrequencyFeeLeg /* = K_QUARTERLY*/ ,
				const int& DayCountBasis /* = KACTUAL_360*/ ,
				const double& FixedPayerAmount /* = 100.*/ , 
				const qPAYMENT_PREMIUM_LEG& AccruedOnDefault /* = qACCRUED_SETTLED*/ ,
				const std::string& ccy /* = ARM_DEFAULT_CURRENCY*/ , 
				const double& FloatingPayerAmount /* = 0.0*/ ,
				const int& stubrule, /* = 1, /* 1 <-> SHORTSTART */
				const double& CreditLag /* = DEFAULT_CREDIT_LAG*/ ,
				const int& FrequencyDefLeg /* = -1*/ ,
				const int& intRule /* = K_ADJUSTED*/ ,
				const bool& includematurity /* = false*/ ,
				const int& adjStartDate /* = 1*/ ,
				const std::string& payCalName /* = NULL*/ ,
				const qSecurity_TYPE& cdstype /* = qCDS_INDEX*/ );

	   
	   virtual ~ICM_Cmcds(void) {};

	   void BitwiseCopy(const ARM_Object* srcCmcds);

	   void Copy(const ARM_Object* srcCmcds) ;

	   ARM_Object* Clone(void) ;

	   void Print(void)
	   {
		   printf("\n\n ===> ICM_Cmcds\n");
		   //...
	   }

       double GetCapLevel() { return itsCapLevel ; }
	   double GetFloorLevel() { return itsFloorLevel ; }

	   void SetCapLevel(const double& value) { itsCapLevel = value ; }
	   void SetFloorLevel(const double& value) { itsFloorLevel = value ; }
	   
};



#endif
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
