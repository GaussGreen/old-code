#pragma warning (disable:4786)
#include "ICMKernel\inst\ICM_Cln.h"

void ICM_Cln::BitwiseCopy(const ARM_Object* srcleg)
{    ICM_Cln* leg = (ICM_Cln *) srcleg;  }

void ICM_Cln::Copy(const ARM_Object* srcleg)
{    ICM_Leg::Copy(srcleg);
     BitwiseCopy(srcleg);  }

ARM_Object* ICM_Cln::Clone(void)
{    ICM_Cln* theClone = new ICM_Cln();
     theClone->Copy(this);
     return(theClone);  }

void ICM_Cln::Init(void)
{	SetName(ICM_CLN); 
	itsRedempNotional = 0.; }

ICM_Cln::ICM_Cln(const ARM_Date& startDate, 
					const ARM_Date& endDate, 
					const ARM_Date* refDate,
   					const ARM_Date* FstCpnEffDate,
					ARM_IRIndex* irIndex,
					const double& spread,
					const double& notional,
					const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
					const int&	AccruedDayCount,
					const int& rcvOrPay, 
					const int& dayCount, 
					const int& decompFreq,
					const int& stubRule,
					const int& resetgap,
					const std::string&  resetCalName, 
					const std::string&  discountCcy,
					const std::string&  payCalName,
					const int& nxChange,
					const bool& includematurity,
					const int& adjStartDate,
					const qCredit_Leg_Type& LegType,
					const double& Binary,
					const string& name)
{
	Init();

	Set(startDate,endDate,refDate,FstCpnEffDate,irIndex,spread,notional,AccruedOnDefault,
					AccruedDayCount,rcvOrPay,dayCount,decompFreq,stubRule,
					resetgap,resetCalName,discountCcy,payCalName,nxChange,
					includematurity,adjStartDate,LegType,Binary,name);

}	


void ICM_Cln:: Set(const ARM_Date& startDate, 
					const ARM_Date& endDate, 
					const ARM_Date* refDate,
   					const ARM_Date* FstCpnEffDate,
					ARM_IRIndex* irIndex,
					const double& spread,
					const double& notional,
					const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
					const int&	AccruedDayCount,
					const int& rcvOrPay, 
					const int& dayCount, 
					const int& decompFreq,
					const int& stubRule,
					const int& resetgap,
					const std::string&  resetCalName, 
					const std::string&  discountCcy,
					const std::string&  payCalName,
					const int& nxChange,
					const bool& includematurity,
					const int& adjStartDate,
					const qCredit_Leg_Type& LegType,
					const double& Binary,
					const string& name)
{
	ARM_IRIndex* IR_INDEX=NULL;
	ARM_IRIndex FlatIndex;

	if (irIndex==NULL)
	{IR_INDEX =&FlatIndex;}
	else 
	{IR_INDEX = irIndex;};

	ICM_Leg::Set(startDate,endDate,refDate,FstCpnEffDate,IR_INDEX,spread ,AccruedOnDefault,
					AccruedDayCount,0.,0.,rcvOrPay,dayCount, 
					decompFreq,stubRule,resetgap,resetCalName,discountCcy,
					payCalName,nxChange,includematurity,adjStartDate,LegType,Binary,
					name); 

	itsRedempNotional = notional;

	ARM_ReferenceValue refvalfix(notional , 1 /* price */, 0 /* K_CONSTANT */);
	SetAmount(&refvalfix,100.);

//	ARM_ReferenceValue refvalspread(spread*100., 1 /* price */, 0 /* K_CONSTANT */);
//	SetVariableSpread(&refvalspread);	
}





