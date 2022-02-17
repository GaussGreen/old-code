
#include "firsttoinc.h"
#include "ICMKernel\util\icm_schedule_info.h"
#include "ICMKernel\inst\icm_leg.h"
#include "ARMKernel\crv\zerocurv.h"
#include "ICMKernel\util\icm_utils.h"


ICM_Schedule_Info::ICM_Schedule_Info(const ARM_Date&	aEffectiveDate,
							const ARM_Date&	aMaturityDate,			
							 int	aPayFrequency,
							 int aResetFreq ,
							 int	aDayCount,
							 int	aStubrule,
							 int	aIntRule,
							const std::string& aPayCalName,
							 int aPayTiming,
							 int aResetTiming,
							 int aFwdRule,
							 bool	aIncludeMaturity,
							 int aAdj,
							 int	aStartAdj,
							 int aAccDayCount,
							const ARM_Date*	apReferenceDate,
							const ARM_Date*	apFirstCpnEffDate,
							 qCDS_ADJ CDSAdj)
{
	Set(aEffectiveDate,
		aMaturityDate,			
		aPayFrequency,
		aResetFreq ,
		aDayCount,
		aStubrule,
		aIntRule,
		aPayCalName,
		aPayTiming,
		aResetTiming,
		aFwdRule,
		aIncludeMaturity,
		aAdj,
		aStartAdj,
		aAccDayCount,
		apReferenceDate,
		apFirstCpnEffDate,
		CDSAdj);
}

// FIXMEFRED: mig.vc8 (28/05/2007 15:05:13): cast
void ICM_Schedule_Info::Set(const ARM_Date&	aEffectiveDate,
						const ARM_Date&	aMaturityDate,			
						 int	aPayFrequency,
						 int aResetFreq ,
						 int	aDayCount,
						 int	aStubrule,
						 int	aIntRule,
						const std::string& aPayCalName,
						 int aPayTiming,
						 int aResetTiming,
						 int aFwdRule,
						 bool	aIncludeMaturity,
						 int aAdj,
						 int	aStartAdj,
						 int aAccDayCount,
						const ARM_Date*	apReferenceDate,
						const ARM_Date*	apFirstCpnEffDate,
						 qCDS_ADJ CDSAdj)
{
	Init();
	SetName(ICM_SCHEDULE_INFO);
	itsEffectiveDate = aEffectiveDate;
	itsMaturityDate = aMaturityDate;			
	itsPayFrequency = aPayFrequency;
	itsResetFreq = aResetFreq;
	itsDayCount = aDayCount;
	itsStubrule = aStubrule;
	itsIntRule = aIntRule;
	itsPayCalName = aPayCalName;
	itsPayTiming = aPayTiming;
	itsResetTiming = aResetTiming;
	itsFwdRule = aFwdRule;
	itsIncludeMaturity = aIncludeMaturity;
	itsAdj = aAdj;
	itsStartAdj = aStartAdj;
	itsAccDayCount = aAccDayCount;
	itsCDSAdj = CDSAdj;
	// clone for *
	if(itspReferenceDate) delete itspReferenceDate;
	itspReferenceDate = NULL;
	if (apReferenceDate)
		itspReferenceDate = (ARM_Date*)(apReferenceDate)->Clone();
	
	if (itspFirstCpnEffDate) delete itspFirstCpnEffDate;
	itspFirstCpnEffDate = NULL;
	if (apFirstCpnEffDate)
		itspFirstCpnEffDate = (ARM_Date*)(apFirstCpnEffDate)->Clone();
}
void ICM_Schedule_Info::Init(void) 
{
	// from ARM_SwapLeg ,ARM Index and ICM_LEG			
	itsPayFrequency = K_QUARTERLY;
	itsResetFreq = K_QUARTERLY;
	itsDayCount = K30_360;
	itsStubrule = K_SHORTSTART;
	itsIntRule = K_ADJUSTED; // from ICM_LEG
	itsPayCalName = ""; //tmpPayCal   = ccy->GetPayCalName(K_INDEX_TYPE);
	itsPayTiming = K_ARREARS;
	itsResetTiming = K_ADVANCE;
	itsFwdRule = K_MOD_FOLLOWING;
	itsIncludeMaturity = false;
	itsAdj = K_ADJUSTED;
	itsStartAdj = K_ADJUSTED;
	itsAccDayCount = KACTUAL_365;
	itspReferenceDate = NULL;
	itspFirstCpnEffDate = NULL;
	itsCDSAdj = qCredit_Adjust20; // STDCDS
}


void ICM_Schedule_Info::BitwiseCopy(const ARM_Object* srcS)
{
	
	ICM_Schedule_Info* src = (ICM_Schedule_Info *) srcS;
	itsEffectiveDate = src->itsEffectiveDate;
	itsMaturityDate = src->itsMaturityDate;			
	itsPayFrequency = src->itsPayFrequency;
	itsResetFreq = src->itsResetFreq;
	itsDayCount = src->itsDayCount;
	itsStubrule = src->itsStubrule;
	itsIntRule = src->itsIntRule;
	itsPayCalName = src->itsPayCalName;
	itsPayTiming = src->itsPayTiming;
	itsResetTiming = src->itsResetTiming;
	itsFwdRule = src->itsFwdRule;
	itsIncludeMaturity = src->itsIncludeMaturity;
	itsAdj = src->itsAdj;
	itsStartAdj = src->itsStartAdj;
	itsAccDayCount = src->itsAccDayCount;
	itsCDSAdj = src->itsCDSAdj;
	if(itspReferenceDate) delete itspReferenceDate;
	itspReferenceDate = NULL;
	if (src->itspReferenceDate)
		itspReferenceDate = (ARM_Date*)(src->itspReferenceDate)->Clone();
	
	if (itspFirstCpnEffDate) delete itspFirstCpnEffDate;
	itspFirstCpnEffDate = NULL;
	if (src->itspFirstCpnEffDate)
		itspFirstCpnEffDate = (ARM_Date*)(src->itspFirstCpnEffDate)->Clone();

}
void ICM_Schedule_Info::Copy(const ARM_Object* srcInfo)
{
	 ARM_Object::Copy(srcInfo);
     BitwiseCopy(srcInfo);
}



ARM_Object* ICM_Schedule_Info::Clone(void)
{
	ICM_Schedule_Info* theClone = new ICM_Schedule_Info();

     theClone->Copy(this);
 
     return(theClone);
}

void ICM_Schedule_Info::View(char* id , FILE* ficOut)
{
	FILE* fOut;
    char fOutName[200];
	bool status = true;
   
 
    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
 
       (void) unlink(fOutName);

       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

	char ds[20];
    char de[20];
    char dr[20];
    char df[20];

	fprintf(fOut, "\t\t\t ----------------- SCHEDULE INFORMATIONS ----------------- \n");
	((ARM_Date) itsEffectiveDate).JulianToStrDateDay(ds);
	fprintf(fOut, "StartDate : %14s\n", ds);
    ((ARM_Date) itsMaturityDate).JulianToStrDateDay(de);
	fprintf(fOut, "EndDate : %14s\n", de);
	if (itspReferenceDate){
		((ARM_Date*) itspReferenceDate)->JulianToStrDateDay(dr);
		fprintf(fOut, "RollDate : %14s\n", dr);
	} else {
		fprintf(fOut, "RollDate : not given\n");
	}
	if(itspFirstCpnEffDate){
		((ARM_Date*) itspFirstCpnEffDate)->JulianToStrDateDay(df);
		fprintf(fOut, "FirstCouponEffectiveDate : %14s\n", df);	
	}else {
		fprintf(fOut, "FirstCouponEffectiveDate : not given\n");
	}

	fprintf(fOut, "\n Reset Frequency\t: %s",  ARM_ParamView::GetMappingName(S_FREQUENCY, itsResetFreq));
    fprintf(fOut, "\n Pay Frequency\t: %s",  ARM_ParamView::GetMappingName(S_FREQUENCY, itsPayFrequency));
	fprintf(fOut, "\n DayCount \t      : %s ", ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsDayCount));
	fprintf(fOut, "\n Acc DayCount \t      : %s ", ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsAccDayCount));
	fprintf(fOut, "\n Int Rule \t\t: %s", ARM_ParamView::GetMappingName(S_INTEREST_RULES, itsIntRule)); 
	
	fprintf(fOut, "\n Pay Timing\t: %s", ARM_ParamView::GetMappingName(S_TIMING_MOD, itsPayTiming));
	fprintf(fOut, "\n Reset Timing\t: %s", ARM_ParamView::GetMappingName(S_TIMING_MOD, itsResetTiming));  
	fprintf(fOut, "\n Forward Rule\t: %s",  ARM_ParamView::GetMappingName(S_FORWARD_RULES, itsFwdRule));
	fprintf(fOut,"\n Stub Rule\t: %s \n", ARM_ParamView::GetMappingName(S_STUB_RULES, itsStubrule));
	fprintf(fOut,"\n Calendar name\t: %s \n",itsPayCalName.c_str());
	string cdsAdj = "";
	ICM_EnumsCnv::toString(itsCDSAdj, cdsAdj);
	fprintf(fOut," CDS ADJ \t: %s \n",cdsAdj.c_str());

	if (itsIncludeMaturity) fprintf(fOut,"IncludeMaturity:OK\n"); else fprintf(fOut,"IncludeMaturity:NO\n");
	if ( itsStartAdj == 1) fprintf(fOut,"Adj StartSate:OK\n"); else fprintf(fOut,"Adj StartSate:NO\n");
	
    if ( ficOut == NULL )
    {
       fclose(fOut);
    }

}


ICM_Schedule_Info::~ICM_Schedule_Info() {
	if (itspReferenceDate) delete itspReferenceDate;
	itspReferenceDate = NULL;

	if (itspFirstCpnEffDate) delete itspReferenceDate;
	itspFirstCpnEffDate = NULL;

}

