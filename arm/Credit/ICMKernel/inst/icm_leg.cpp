#include "firsttoinc.h"
#include "ICMKernel\inst\icm_leg.h"
#include "ICMKernel\pricer\icm_pricer.h"
#include "ICMKernel/glob/icm_mktdatamng.h"
#include "ICMKernel/mod/modelmulticurves.h"
#include "ICMKernel/inst/icm_credit_index.h"


ICM_Leg::~ICM_Leg()
	{
		if (itsCreditInfos)
			delete itsCreditInfos;
		itsCreditInfos = NULL;

		if (itsNotExchange)
			delete itsNotExchange;
		itsNotExchange = NULL;

		if (itsCreditIndex)
			delete itsCreditIndex;
		itsCreditIndex = NULL;

		if (itsFwdFixedDate)
			delete itsFwdFixedDate;
		itsFwdFixedDate = NULL;

		if (itsFwdCalcTypes)
			delete itsFwdCalcTypes;
		itsFwdCalcTypes = NULL;

		if (itsRefRiskyType)
			delete itsRefRiskyType;
		itsRefRiskyType = NULL;

		if (itsRefRiskyDate)
			delete itsRefRiskyDate;
		itsRefRiskyDate = NULL;

		if (itsSwapLeg)
			delete itsSwapLeg;
		itsSwapLeg = NULL;

		if (itsRefPartRate)
			delete itsRefPartRate;
		itsRefPartRate = NULL;

		if (itsFstCpnEffDate) 
			delete itsFstCpnEffDate;
		itsFstCpnEffDate = NULL; 

	}

void ICM_Leg::SetCreditIndex(ICM_Credit_Index* index) 
{
	if (itsCreditIndex)
		delete itsCreditIndex;
	itsCreditIndex = index ;

	ARM_SwapLeg::SetIRIndex(index);
	// ResetCreditInfos(); 
}

void ICM_Leg::SetCreditLegType(const qCredit_Leg_Type& type) 
{
	itsCreditLegType=type;

	switch (itsCreditLegType)
	{
		case qNone_Leg :
		{itsCreditLegStyle = qStyle_None_Leg;
		break;}
		case qZeroCoupon_Leg :
		case qInfine_Leg :
		case qFunded_Leg :
		case qRunning_Leg :
		case qCorridor_Leg :
		{itsCreditLegStyle = qStyle_Premium_Leg;
		break;}
		case qInflation_Leg :
		case qSwapLeg :
		case qForward_IRrates_Leg :
		case qCMCDS_Leg :
		case qConstant_Maturity_Leg :
		{itsCreditLegStyle = qStyle_Premium_Leg;
		SetLegType(K_FLOATING_LEG);	
		break;}
		case qStandart_Recovery_Leg :
		{itsCreditLegStyle = qStyle_Recovery_Leg;
		break;}
		default :
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : Unknown LegType ");
	}
}


void ICM_Leg::BitwiseCopy(const ARM_Object* srcleg)
{
    ICM_Leg* leg = (ICM_Leg *) srcleg;

	itsAccruedOnDefault = leg->itsAccruedOnDefault;
	itsCreditLegStyle = leg->itsCreditLegStyle;
	itsCreditLegType = leg->itsCreditLegType;
	itsRefSpread = leg->itsRefSpread ;

	if (itsCreditIndex) delete itsCreditIndex;

	if (leg->itsCreditIndex)
		itsCreditIndex = (ICM_Credit_Index*) leg->itsCreditIndex->Clone();

	if (itsCreditInfos)
		delete itsCreditInfos;
	itsCreditInfos = NULL;

	if (leg->itsCreditInfos)
	{
	itsCreditInfos = (ICM_Security*) (leg->itsCreditInfos)->Clone(); 
	}

	itsIncludeMaturity = leg->itsIncludeMaturity;

	if (itsNotExchange)
		delete itsNotExchange;
	itsNotExchange = NULL;

	if (leg->itsNotExchange)
	{
	itsNotExchange = (ARM_ReferenceValue*) (leg->itsNotExchange)->Clone(); 
	}

	SetFwdFixedDate(leg->itsFwdFixedDate); 

	itsRatesPricer = leg->itsRatesPricer;
	itsCreditLag = leg->itsCreditLag;
	itsBinary = leg->itsBinary;
	itsSingleName = leg->itsSingleName;

	if (itsRefRiskyDate)
		delete itsRefRiskyDate;
	itsRefRiskyDate = NULL;

	if (itsRefRiskyType)
		delete itsRefRiskyType;
	itsRefRiskyType = NULL;

	if (leg->itsRefRiskyDate)
	{itsRefRiskyDate = (ARM_ReferenceValue*) (leg->itsRefRiskyDate)->Clone();}

	if (leg->itsRefRiskyType)
	{itsRefRiskyType = (ARM_ReferenceValue*) (leg->itsRefRiskyType)->Clone();}

	if (itsFwdCalcTypes)
		delete itsFwdCalcTypes;
	itsFwdCalcTypes = NULL;

	if (leg->itsFwdCalcTypes)
	{itsFwdCalcTypes = (ARM_ReferenceValue*) (leg->itsFwdCalcTypes)->Clone();}

	if (leg->itsSwapLeg)
	{itsSwapLeg = (ARM_SwapLeg*) (leg->itsSwapLeg)->Clone();}
	itsIntRule = leg->itsIntRule ;			

	SetFstCpnEffDate(leg->itsFstCpnEffDate);
	itsPastFixings = leg->itsPastFixings ; 
}

void ICM_Leg::Copy(const ARM_Object* srcleg)
{
     ARM_SwapLeg::Copy(srcleg);
 
     BitwiseCopy(srcleg);
}

void ICM_Leg::Rebuild(void)
{
    // compute expected forward rates and flow values after setting the model
    ARM_SwapLeg::Rebuild();
}

ARM_Object* ICM_Leg::Clone()
{
     ICM_Leg* theClone = new ICM_Leg();

     theClone->Copy(this);
 
     return(theClone);
}

ARM_Object* ICM_Leg::Clone() const
{
     ICM_Leg* theClone = new ICM_Leg();
     theClone->Copy(this);
     return(theClone);
}

void ICM_Leg::Init(void) // initialization with default values
{
	SetName(ICM_LEG);

	SetModel(NULL);

	itsAccruedOnDefault=qACCRUED_SETTLED;
	itsCreditInfos = NULL;
	itsIncludeMaturity = false;
	itsNotExchange = NULL;
	itsCreditLegStyle = qStyle_None_Leg;
	itsCreditLegType = qNone_Leg;
	itsCreditIndex = NULL;
	itsRefSpread = CREDIT_DEFAULT_VALUE;
	itsFwdFixedDate = 0; // = new char[5];
	
	itsFwdCalcTypes = NULL;
	itsRatesPricer = NULL;
	itsCreditLag = 0;
	itsBinary = CREDIT_DEFAULT_VALUE;
	itsSingleName = "undefine";
	itsRefRiskyDate = NULL;
	itsRefRiskyType = NULL;
	itsSwapLeg = NULL;
	itsRefPartRate = NULL;
	itsIntRule = K_ADJUSTED ;
	itsFstCpnEffDate = NULL;
}

// ---------------------------------------------------------------------
// Set of fix Leg with ARM_Date refDate
// ---------------------------------------------------------------------

 void ICM_Leg::Set(const ARM_Date& startDate, 
				 const ARM_Date& endDate, 
				 const ARM_Date* refDate,
         		 const ARM_Date* FstCpnEffDate,
				 const double& fixedRate,
				 const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
				 const int&    AccruedDayCount,
				 const double&	LastIndexFixing,
				 const int& rcvOrPay , 
				 const int& freq, 
				 const int& dayCount, 
				 const int& decompFreq,
				 const int& payTiming,
				 const int& intRule,
				 const int& stubRule,
				 const std::string& discountCcy,
				 const std::string& payCalName,
				 const int& nxChange,
				 const bool& includematurity,
				 const int& adjStartDate,
				 const qCredit_Leg_Type& LegType,
				 const double& Binary,
				 const string& name)
{

	SetName(ICM_LEG);

	char buff [50] ; 
	if (refDate) unconst(*refDate).JulianToStrDate(buff); 
	
	

	itsIntRule=intRule ; 
	ARM_Currency ccy(discountCcy.c_str()); 
	ARM_SwapLeg::Set(unconst(startDate), unconst(endDate), (fixedRate*100.) , rcvOrPay, freq, dayCount, 
		decompFreq, payTiming, intRule==K_MATUNADJUSTED?K_ADJUSTED:intRule, stubRule, &ccy,
		payCalName.empty()?0:(char*)payCalName.c_str(), nxChange, refDate?buff:0 ,adjStartDate);
	
	SetFstCpnEffDate(FstCpnEffDate);
	if (itsFstCpnEffDate) 
		CustomizeFirstPeriod(unconst(*itsFstCpnEffDate)); 
	
	itsAccruedOnDefault = AccruedOnDefault; 
	itsIncludeMaturity = includematurity;	
	SetCreditLegType(LegType);
	SetBinary(Binary);
	SetSingleName(name);
	// should: SetFwdCalcTypes(qCPN_FIXED); 
}	

// ---------------------------------------------------------------------
// Set of floating Leg with char* refDate
// ---------------------------------------------------------------------

void ICM_Leg::Set(const ARM_Date& startDate, 
				   const ARM_Date& endDate, 
				   const ARM_Date* refDate, 
         		   const ARM_Date* FstCpnEffDate,
				   ARM_IRIndex* irIndex,
				   const double& spread,
				   const qPAYMENT_PREMIUM_LEG& AccruedOnDefault,
				   const int&				    AccruedDayCount,
				   const double&				InitialRate,
				   const double&				LastIndexFixing,
				   const int& rcvOrPay, 
				   const int& dayCount, 
				   const int& decompFreq,
				   const int& stubRule,
				   const int& resetgap,
				   const std::string& resetCalName, 
				   const std::string& discountCcy,
				   const std::string& payCalName,
				   const int& nxChange,
				   const bool& includematurity,
				   const int& adjStartDate,
				   const qCredit_Leg_Type& LegType,
				   const double& Binary,
				   const string& name )
{

	if (startDate.GetJulian()>=endDate.GetJulian())
      throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
 	  "EndDate >= StartDate");


	char buff [50] ; 
	if (refDate) unconst(*refDate).JulianToStrDate(buff); 

	ARM_Currency ccy(discountCcy.c_str()); 
	ARM_SwapLeg::Set(unconst(startDate), unconst(endDate), irIndex, rcvOrPay, (spread*100.), stubRule, 
					 decompFreq, &ccy, dayCount, resetgap, 
					 resetCalName.empty()?0:(char*)payCalName.c_str(), 
					 payCalName.empty()?0:(char*)payCalName.c_str(), 
					 1, nxChange, refDate?buff:0,adjStartDate);

	SetFstCpnEffDate(FstCpnEffDate);
	if (itsFstCpnEffDate) 
		CustomizeFirstPeriod(unconst(*itsFstCpnEffDate)); 

	itsAccruedOnDefault = AccruedOnDefault; 
	itsIncludeMaturity = includematurity;	
	SetCreditLegType(LegType);
	SetBinary(Binary);
	SetSingleName(name);
	// should: SetFwdCalcTypes(qCPNFWD_FIXED); 
}

// ---------------------------------------------------------------------
// Set of floating Leg 
// ---------------------------------------------------------------------

void ICM_Leg::Set(const ARM_Date& StartDate, 
            const ARM_Date& EndDate, 
			const ARM_Date* refDate,
         	const ARM_Date* FstCpnEffDate,
			const double& FixedRate,
			const double& FixedNotional,
			ARM_ReferenceValue* Notional,
			ARM_ReferenceValue* Rates,
			ARM_ReferenceValue* Exchange,
			const int& frequency, 
            const int& dayCount, 
            const int& payTiming,
            const int& intRule,
            const int& stubRule,
            const std::string&  discountCcy,
            const std::string& payCalName,
			const qCredit_Leg_Type& LegType,
			const bool& includematurity,
			const int& adjStartDate,
			ARM_IRIndex* irIndex,
			const double& Binary,
			const string& name,
			const int& NXchange,
			const qPAYMENT_PREMIUM_LEG& AccruedOnDefault)
{
	SetName(ICM_LEG);

	char buff [50] ; 
	if (refDate) unconst(*refDate).JulianToStrDate(buff); 
	itsIntRule=intRule; 

	SetFstCpnEffDate(FstCpnEffDate);

	ARM_Currency ccy(discountCcy.c_str()); 
	ARM_SwapLeg::Set(unconst(StartDate), unconst(EndDate), (FixedRate*100.), K_RCV, frequency, dayCount, 
					K_COMP_PROP, payTiming, intRule==K_MATUNADJUSTED?K_ADJUSTED:intRule, stubRule, 
					&ccy,
					payCalName.empty()?0:(char*)payCalName.c_str(), 
					NXchange, refDate?buff:0,adjStartDate);


	if (irIndex) SetIRIndex(irIndex);
	CptCashFlowDates();
	
	if (itsFstCpnEffDate) 
		CustomizeFirstPeriod(unconst(*itsFstCpnEffDate)); 

	if (Rates) SetVariableSpread((ARM_ReferenceValue*)Rates->Clone());
			
	if (Notional) 
		SetAmount((ARM_ReferenceValue*)Notional->Clone());
	else
	{
		ARM_ReferenceValue refvalfix(FixedNotional , 1 /* price */, 0 /* K_CONSTANT */);
		SetAmount(&refvalfix,100.0);
	}	

	if (Exchange) SetNotExchange((ARM_ReferenceValue*)Exchange->Clone());
	SetCreditLegType(LegType);

	//if(irIndex != NULL)
	//	SetCreditIndex((ICM_Credit_Index*)irIndex->Clone()) ;

	if (irIndex!=NULL) 
		this->SetFwdCalcTypes(&ARM_ReferenceValue(qCPN_FIXED_AND_FWD)) ;

	itsIncludeMaturity = includematurity;	
	itsAccruedOnDefault = AccruedOnDefault;
	SetRefSpread(10000*FixedRate) ;
	SetBinary(Binary);
	SetSingleName(name);
	// should: SetFwdCalcTypes(qCPNFWD_FIXED); 
}

// ---------------------------------------------------------------------
// Set Libor & Euribor risky Leg 
// ---------------------------------------------------------------------

void ICM_Leg::Set(const ARM_Date& StartDate, 
            const ARM_Date& EndDate, 
			const ARM_Date* refDate,
            const ARM_Date* FstCpnEffDate,
			const ARM_INDEX_TYPE& liborType,
			const double& spread, 
			const int& resetFreq, 
			const int& payFreq, 
            const int& resetTiming, 
			const int& payTiming,
            const std::string&  discountCcy,
            const int& intRule,
            const int& resetgap,
            const std::string& resetCalName,
            const std::string& payCalName,
            const int& stubRule,
            const int& adjStartDate,
			const int& dayCount,
			const bool& includematurity,
			const double& Binary,
			const string& name)
{
	SetName(ICM_LEG);

	char buff [50] ; 
	if (refDate) unconst(*refDate).JulianToStrDate(buff); 

	SetFstCpnEffDate(FstCpnEffDate);
	ARM_Currency ccy(discountCcy.c_str()); 

	ARM_SwapLeg::Set(unconst(StartDate),unconst(EndDate),liborType, K_RCV,
                    spread,resetFreq,payFreq,resetTiming,payTiming,
                    &ccy,resetgap,
					resetCalName.empty()?0:(char*)resetCalName.c_str(),
                    payCalName.empty()?0:(char*)payCalName.c_str(),
					K_COMP_PROP,K_NX_NONE,stubRule,
                    refDate?buff:0,adjStartDate);
	if (itsFstCpnEffDate) 
		CustomizeFirstPeriod(unconst(*itsFstCpnEffDate)); 

	SetCreditLegType(qForward_IRrates_Leg);
	itsIncludeMaturity = includematurity;	
	SetRefSpread(10000*spread) ;
	SetBinary(Binary);
	SetSingleName(name);
	// should: SetFwdCalcTypes(qCPNFWD_FIXED); 
}


// ---------------------------------------------------------------------
// Summit Leg 
// ---------------------------------------------------------------------

void ICM_Leg::Set(ARM_Vector* flowStartDates, 
		    ARM_Vector* flowEndDates,
            ARM_Vector* paymentDates, 
			ARM_Vector* resetDates, 
            ARM_Vector* intDays, 
			ARM_Vector* fwdRates,
            ARM_ReferenceValue* Notional, 
			ARM_IRIndex* irIndex,
            const int& rcvOrPay, 
			const double& spread,
            const double& fixRate, 
            const std::string&  discountCcy,
            const int& NxId,
			const int& decompPricingFlag,
			const qCredit_Leg_Type& LegType,
			const double& Binary,
			const string& name)
{
	SetName(ICM_LEG);

	ARM_Currency ccy(discountCcy.c_str()); 
	ARM_SwapLeg::Set(flowStartDates,flowEndDates,paymentDates,resetDates, 
                 intDays,fwdRates,Notional,irIndex,rcvOrPay,spread, 
                 fixRate,&ccy,NxId,decompPricingFlag);

	SetCreditLegType(LegType);
	SetBinary(Binary);
	SetSingleName(name);
}

// -------------------------------------------------------------------------------------------
//	Retrieves the current period index such that index = max{i in N/ startdate[i]<= AsofDate}
// -------------------------------------------------------------------------------------------

int ICM_Leg::PeriodIndex(const ARM_Date& AsofDate)   
{
	int numflows = GetNumFlows();
	int index = 0;
	ARM_Vector* startDates = GetFlowStartDates();
	ARM_Vector* endDates = GetFlowEndDates();

	if (AsofDate.GetJulian() > endDates->Elt(numflows-1))
		return (-1);

	while (startDates->Elt(index) <= AsofDate.GetJulian()) 
	{index++;
	 if (index >= numflows) break;}

	index--;
	return (index);
}

/* ****************************************************************************************************************** */
/*!	\fn void View()
/* ****************************************************************************************************************** */
void ICM_Leg::View(char* id, FILE* ficOut)
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

	fprintf(fOut, "\t\t\t ----------------- Leg Informations ----------------- \n");
	if (itsIncludeMaturity) fprintf(fOut,"IncludeMaturity:OK\n"); else fprintf(fOut,"IncludeMaturity:NO\n");
	if ( GetAdjStartDateFlag()) fprintf(fOut,"Adj StartSate:OK\n"); else fprintf(fOut,"Adj StartSate:NO\n");
	switch (itsCreditLegStyle)
	{
		case qStyle_Premium_Leg:
			fprintf(fOut,"Leg Style:Premium Leg\n");
			break;
		case qStyle_Recovery_Leg:
			fprintf(fOut,"Leg Style:Recovery Leg\n");
			break;
		case qStyle_None_Leg:
		default :
			fprintf(fOut,"Leg Style:Gohst Leg\n");
	}

	switch (itsCreditLegType)
	{
		case qCorridor_Leg:
			fprintf(fOut,"Leg Type:Corridor Leg\n");
			break;
		case qRunning_Leg:
			fprintf(fOut,"Leg Type:Running Leg\n");
			break;
		case qFunded_Leg:
			fprintf(fOut,"Leg Type:Funded Leg\n");
			break;
		case qInfine_Leg:
			fprintf(fOut,"Leg Type:Infine Leg\n");
			break;
		case qZeroCoupon_Leg:
			fprintf(fOut,"Leg Type:Zero Coupon Leg\n");
			break;
		case qConstant_Maturity_Leg:
			fprintf(fOut,"Leg Type:CM Leg\n");
			break;
		case qCMCDS_Leg:
			fprintf(fOut,"Leg Type:CM-Cds Leg\n");
			break;
		case qForward_IRrates_Leg:
			fprintf(fOut,"Leg Type:IR Fwd Leg\n");
			break;
		case qSwapLeg:
			fprintf(fOut,"Leg Type:Swap Leg\n");
			break;
		case qStandart_Recovery_Leg:
			fprintf(fOut,"Leg Type:Recovery Leg\n");
			break;
		case qInflation_Leg:
			fprintf(fOut,"Leg Type:Inflation Leg\n");
			break;
		case qNone_Leg: 
			fprintf(fOut,"LegType: None_Leg\n"); 
			break ;
		default :
			fprintf(fOut,"Leg Type:Unknown Leg\n");
			status = false;
	}
	
	fprintf(fOut,"Credit Lag:%f\n",(double)itsCreditLag);
	if (itsBinary == CREDIT_DEFAULT_VALUE) fprintf(fOut,"Binary:NO\n"); else fprintf(fOut,"Binary:%f\n",itsBinary);
	fprintf(fOut,"Single Issuer Name :%s\n",itsSingleName.c_str());
	fprintf(fOut,"RefSpread :%f\n",itsRefSpread );
	fprintf(fOut,"PastFixings size: %d\n",itsPastFixings.size())  ;
	if (itsPastFixings.size()!=0) 
	{
		std::map<double,double>::const_iterator it=itsPastFixings.begin() ;
		while (it!=itsPastFixings.end()) 
		{
			char d[20]; 
			ARM_Date(it->first).JulianToStrDate(d); 
			fprintf(fOut,"%s : %10.4lf\n",d,it->second) ;
			++it; 
		}
	}

	if (GetIRIndex()) {GetIRIndex()->View(id, fOut);}
	
	if (status){
	ICM_Security* security = GetCreditInfos();
	security->View(id, fOut);}
	
    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}


void ICM_Leg::GenScheduleDates(int datesType, int viewInitialExch, int modfoll, int creditgap, char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[200];

    ARM_Vector* julDates;
    int vectSize;
	int prev = 0;

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

    switch(datesType)
    {
        case qACC_START_DATE :
        {
            julDates = GetFlowStartDates();
        };
        break;
        case qACC_END_DATE :
        {
            julDates = GetFlowEndDates();
        };
        break;
        case qPAY_DATE :
        {
            julDates = GetPaymentDates();
        };
        break;
        case qOBS_START_DATE :
        {
            julDates = GetFlowStartDates();
        };
        break;
        case qOBS_END_DATE :
        {
            julDates = GetFlowEndDates();
			prev = (-1)*creditgap;
        };
        break;

        default :
        {
        };
        break;
    }

    if ( julDates == NULL )
    {
       vectSize = 0;
    }
    else
    {
       vectSize = julDates->GetSize();
    }

    // Now Print out Vect values in date format

    char d1[20];

// TMP : MA+JPR 17 Oct 2003
    if ((viewInitialExch == K_YES) && ( datesType != K_RESET_DATES ) && (GetAmortization()) && ( julDates != NULL ) 
        && 
        (( GetNxFlag() == K_NX_START ) 
         || 
         ( GetNxFlag() == K_NX_BOTH )
        ) 
       )
    {
       ARM_Vector* dates = GetAmortization()->GetDiscreteDates();
		   
       fprintf(fOut,"%d\n", vectSize+1);

       ((ARM_Date) (dates->Elt(0))).JulianToStrDate(d1);

       fprintf(fOut,"%s\n", d1);
    }
    else
    {
       fprintf(fOut,"%d\n", vectSize);
    }

    for (int i = 0;  i< vectSize; i++)
    {
	   ARM_Date TmpDate = (ARM_Date) julDates->Elt(i);

	   if (prev) 
	   {
		   for (int il=0;il<prev; il++)
			   TmpDate.PreviousBusinessDay(1,GetCurrencyUnit()->GetCcyName());
	   }	

	   if (datesType == qPAY_DATE) 
	   {				
		TmpDate.AdjustToBusDate(GetCurrencyUnit()->GetCcyName(),modfoll);
	   }	
	   
       TmpDate.JulianToStrDate(d1);
       fprintf(fOut,"%s\n", d1);
    }

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}

//	---------------------------------------------------------------------------------
void 
ICM_Leg::CptCashFlowDatesCredit(void)
{
	ARM_SwapLeg::CptCashFlowDates();

	ResetCreditInfos(); 
	//if (itsCreditInfos)
	//	delete itsCreditInfos;
	//itsCreditInfos = NULL;
	
	CptCreditInfos();
}

// -------------------------------------------------------------------------
// --						 CptCreditInfos								 ---
// -------------------------------------------------------------------------

void ICM_Leg::CptCreditInfos(void) const 
{
	// ICMQUANTIFYER("ICM_Leg::CptCreditInfos"); 

	if (itsCreditInfos)
		delete itsCreditInfos;
	itsCreditInfos = NULL;

	int NbFlows = GetNumFlows();

// 	if (itsFstCpnEffDate) 
// 		CustomizeFirstPeriod(unconst(*itsFstCpnEffDate)); 
	
	ARM_Vector CouponRates(NbFlows,0.);
	//to avoid bad Fixed Rate in Floating Case
	//values for Fixed Rate : GetFixedRate() :
	//#define K_FLOAT_RATE    -1000000.0
	//#define K_FIXED_RATE    0.0
	//#define K_MARKET_RAT    -1.0
	//#define K_ATMF_CAPLETS  -2.0
	if (GetFixedRate()>0.) {CouponRates+=GetFixedRate();}

	ARM_Vector Notionals(NbFlows,0.);
	ARM_Vector NotionalXchange(NbFlows,0.);
	ARM_Vector ParticipationRates(NbFlows,1.);

	if ((!IsFixedLeg()) && GetVariableSpread())
	{for (int i=0;i<NbFlows;i++)
			CouponRates.Elt(i) = GetVariableSpread()->CptReferenceValue((ARM_Date)GetFlowStartDates()->Elt(i));	}

	if (GetRefPartRate())
	{for (int i=0;i<NbFlows;i++)
			ParticipationRates.Elt(i) = GetRefPartRate()->CptReferenceValue((ARM_Date)GetFlowStartDates()->Elt(i));	}

	if (GetAmount())
	{
		// ARM TAUX convention: notionals "Reference value" is indexed on Pay Dates. 
		for (int i=0;i<NbFlows;i++)	Notionals.Elt(i) = GetAmount()->CptReferenceValue((ARM_Date)GetPaymentDates()->Elt(i));
	}

	if (itsNotExchange)	{
	for (int i=0;i<NbFlows;i++) NotionalXchange.Elt(i) = itsNotExchange->CptReferenceValue((ARM_Date)GetFlowStartDates()->Elt(i));}

	// InterestDays are computed here.. 
	ARM_Vector InterestDays(*GetInterestDays()) ;
	ARM_Vector EndDates(*GetFlowEndDates());  
	ARM_Vector StartDates(*GetFlowStartDates());  
	//	For K_MATURITYUNADJ , the last period adjustment is not handled by ARM
	//	so we do it here
	if (itsIntRule==K_MATUNADJUSTED && ARM_SwapLeg::GetEndDateNA() != ARM_SwapLeg::GetEndDate() )
	{
		int payDates = GetPaymentDates()->size() ;
		if ( (NbFlows>=2) && 
			 (GetPaymentDates()->Elt(NbFlows-1)==GetPaymentDates()->Elt(NbFlows-2) ) )
		{
			EndDates.Elt(NbFlows-2)= ARM_SwapLeg::GetEndDateNA().GetJulian(); 
			InterestDays.Elt(NbFlows-2) =   DaysBetweenDates(GetDayCount(),
				GetFlowStartDates()->Elt(NbFlows-2),
				EndDates.Elt(NbFlows-2)
				); 
			EndDates.Elt(NbFlows-1)= ARM_SwapLeg::GetEndDateNA().GetJulian(); 
			StartDates.Elt(NbFlows-1)= ARM_SwapLeg::GetEndDateNA().GetJulian(); 
		}
		else {
		// ICMLOG("WARNING: Handling K_MATUNADJUSTED leg.. "); 
			EndDates.Elt(NbFlows-1)= ARM_SwapLeg::GetEndDateNA().GetJulian(); 
		InterestDays.Elt(NbFlows-1) =   DaysBetweenDates(GetDayCount(),
			GetFlowStartDates()->Elt(NbFlows-1),
			EndDates.Elt(NbFlows-1)
			); 
		}
	}

	if (itsIncludeMaturity) InterestDays.Elt(NbFlows-1)++ ; 

	itsCreditInfos = new ICM_Security(ARM_SwapLeg::GetStartDateNA(),
										ARM_SwapLeg::GetEndDateNA(),
										&StartDates, // GetFlowStartDates(),
										&EndDates, // GetFlowEndDates(),
										&InterestDays,
										GetPaymentDates(),
										&CouponRates,
										&Notionals,
										&NotionalXchange,
										itsIncludeMaturity,
										GetNotionalAmount(),
										// 0,
										GetDayCount(),
										// GetDayCount(),
										GetCurrencyUnit()->GetCcyName());

	itsCreditInfos->setPastFixings(itsPastFixings); 

	//
	//	Some specifities by leg type 
	//
	itsCreditInfos->SetPaymentFreq(ARM_SwapLeg::GetPaymentFreq());
	switch (itsCreditLegType)
	{
		case qCorridor_Leg:
		{
			itsCreditInfos->SetFwdStartDates(*GetFwdRateStartDates());	// cloning.
			itsCreditInfos->SetFwdEndDates(*GetFwdRateEndDates());
			itsCreditInfos->SetResetDates(*GetResetDates());
		}
		break ;	
		case qConstant_Maturity_Leg: 
		{
			// JLA. Was done at pricing time... 
			if ( GetFwdCalcTypes()) itsCreditInfos->SetFwdType(GetFwdCalcTypes());
				else itsCreditInfos->SetFwdType(qCPN_FWD);

			ICM_Credit_Index* index = GetCreditIndex(); 
			if (!index) 
				ICMTHROW(ERR_INVALID_ARGUMENT,"Can't create CMCDS reset schedule: no index"); 

			ARM_Vector ResetDates(NbFlows,0.);
			ARM_Vector FwdStartDates(NbFlows,0.);
			ARM_Vector FwdEndDates(NbFlows,0.);
			for(int i=0;i<NbFlows;i++)
			{
				double refDate ; 
				if (index->GetResetTiming()==K_ADVANCE) refDate = StartDates[i]; // GetFlowStartDates()->Elt(i) ;
				else refDate = EndDates[i]  ; //  GetFlowEndDates()->Elt(i) ;
				int resetDay = ARM_Date(refDate).GetDayOfWeek() ;
				//	reset date	= "second next friday after accrual start date" , adjusted K_FOLLOWING.
				//	or			= "second prev friday before ..... ",  adjusted K_FOLLOWING
				//	fwd start	= reset +2D business according to the current market (taken from ARM_SwapLeg::CptCashFlowDates )
				//	fwd end		= might have special 20th maturity date
				double resetDateNA ; 
				if (index->GetCM_resetOccur()>=0) 
					resetDateNA = refDate + MODULO((index->GetCM_resetWeekDay()-resetDay-1),7)+1+(index->GetCM_resetOccur()-1)*7 ; 
				else 
					resetDateNA = refDate + MODULO((index->GetCM_resetWeekDay()-resetDay),7)+index->GetCM_resetOccur()*7 ; 

				double resetDate = ARM_Date(resetDateNA).GapBusinessDay(0,GetCurrencyUnit()->GetCcyName()).GetJulian(); 
				//	
				int fwdGap ; 
				if( strcmp( index->GetCurrencyUnit()->GetCcyName(), "AUD") == 0 )fwdGap = 0;
				else if( strcmp( index->GetCurrencyUnit()->GetCcyName(), "CAD") == 0 )fwdGap= 2;
				else fwdGap= index->GetCurrencyUnit()->GetSpotDays();
				double fwdStartDate = ARM_Date(resetDate).NextBusinessDay(fwdGap,GetCurrencyUnit()->GetCcyName()).GetJulian(); 			
				char TEMP[4];
				if(index->GetIsYT()){
					int nb = floor(index->GetYearTerm());
					sprintf(TEMP,"%iY",nb);
					FwdEndDates.Elt(i,0)=AddPeriod(ARM_Date(resetDate),TEMP,GetCurrencyUnit()->GetCcyName(),false,index->GetAdjForTenor()).GetJulian() ;
				} else {
					 ICMTHROW(ERR_INVALID_ARGUMENT,"No YearTerm in Index."); 
				}
				ResetDates.Elt(i)= resetDate; 
				FwdStartDates.Elt(i)= fwdStartDate ; 
			}
			itsCreditInfos->SetFwdStartDates(FwdStartDates);	// cloning.
			itsCreditInfos->SetFwdEndDates(FwdEndDates);
			itsCreditInfos->SetResetDates(ResetDates);
		}
		break ;	// qConstant_Maturity_Leg 
	default:
		{
	
			//JLA. Was done at pricing time ... 
			if ( GetFwdCalcTypes()) itsCreditInfos->SetFwdType(GetFwdCalcTypes());
				else itsCreditInfos->SetFwdType(qCPN_FWD); 

			ARM_Vector FwdEndDates(NbFlows,0.);
			for (int i=0;i<NbFlows;i++){FwdEndDates.Elt(i) = GetFwdRateEndDates()->Elt(i);}
			if (GetCreditIndex())
			if ((GetFwdRateEndDates())&&
				(GetCreditIndex()->GetYearTerm())&&
				(GetCreditIndex()->GetName()==ICM_CREDIT_INDEX))
			{
				if (GetCreditIndex()->GetAdjForTenor()==qCredit_Adjust20)
				{
				char TEMP[4];
				int nb = floor(GetCreditIndex()->GetYearTerm());
				sprintf(TEMP,"%iY",nb);

				for (int i=0;i<NbFlows;i++){
					FwdEndDates.Elt(i) = AddPeriod((ARM_Date)(GetFwdRateStartDates()->Elt(i)),
						TEMP,GetCurrencyUnit()->GetCcyName(),false,qCredit_Adjust20).GetJulian();}
				}
				else
				{
				for (int i=0;i<NbFlows;i++){
					FwdEndDates.Elt(i) = GetFwdRateStartDates()->Elt(i) + 365.*GetCreditIndex()->GetYearTerm();}
				}
			}
			if (GetFwdRateStartDates()) itsCreditInfos->SetFwdStartDates(*GetFwdRateStartDates());
			if (GetResetDates()) itsCreditInfos->SetResetDates(*GetResetDates());
			itsCreditInfos->SetFwdEndDates(FwdEndDates);
		} ;
	}	

	//in case of participation rate
	if (GetRefPartRate()){itsCreditInfos->SetPartRates(ParticipationRates);}

	//Calculation des risky start & end dates
	if (itsRefRiskyType) 
	{
		//JLA: here we take the appropriates (eventually adjusted) dates
		ARM_Vector StartRiskDates(itsCreditInfos->GetAccStartDates()); 
		//JLA: fix : ARM_Vector EndRiskDates(itsCreditInfos->GetAccEndDate());
		ARM_Vector EndRiskDates(itsCreditInfos->GetAccEndDates());

		bool frozenflg = true;
		ARM_Date frozenmaturity;
		for (int i=0;i<NbFlows;i++)
		{
			qRISK_TYPE type = (qRISK_TYPE)(int) itsRefRiskyType->CptReferenceValue((ARM_Date) StartDates[i] /** GetFlowStartDates()->Elt(i) **/ );
			switch (type)
			{
			case qNoRisk:
				StartRiskDates.Elt(i) = 0.;
				EndRiskDates.Elt(i) = 0.;
				break;
			case qRiskToDate:
				StartRiskDates.Elt(i) = itsRefRiskyDate->CptReferenceValue((ARM_Date) StartDates[i] /** GetFlowStartDates()->Elt(i) **/ );
				EndRiskDates.Elt(i) = StartRiskDates.Elt(i);
				if (frozenflg) {frozenmaturity=(ARM_Date)StartRiskDates.Elt(i); frozenflg=false;}
				break;
			case qStdRisk:
				StartRiskDates.Elt(i) = StartDates[i]; // GetFlowStartDates()->Elt(i);
				EndRiskDates.Elt(i) =	EndDates[i] ; // GetFlowEndDates()->Elt(i);
			}
		}
		itsCreditInfos->SetStartRiskDates(StartRiskDates);
		itsCreditInfos->SetEndRiskDates(EndRiskDates);
		//Frozen Maturity case
		if (!frozenflg) {itsCreditInfos->SetFrozenMaturity(frozenmaturity);}
	}
}

// --------------------------------------------
// Reset Forward Coupons 
// --------------------------------------------
void ICM_Leg::ResetCoupons(ARM_Model* model)
{
	ARM_Date AsOf = model->GetStartDate();
	ICM_Security* security = GetCreditInfos();
	security->ComputeYF(AsOf);
	security->SetFwdSpreads(0.);
}

// --------------------------------------------
// Compute Forward Coupons 
// --------------------------------------------
void CptCouponsForFeeLeg(ICM_Leg& leg, ARM_Object* object ,const ARM_Date& AsOf)
{
	ARM_ZeroCurve* ircurve = NULL;
	ARM_ZeroCurve* cpn_ircurve = NULL;
	ARM::ARM_InfCurv* infcrv = NULL;
	ARM_VolCurve* pVolcurve = NULL;

	double YFStartDateFwd = 0.;
	double YFEndDateFwd = 0.;

	if (object->GetName() == ICM_DEFAULT_CURVE_MODEL) 
	{	ICM_DefaultCurveModel* model = dynamic_cast<ICM_DefaultCurveModel*>(object);
		ircurve = model->GetZeroCurve(); 
		cpn_ircurve = model->GetZeroCurve(); 
		pVolcurve = model->GetVolCurve();
		//AsOf = ircurve->GetAsOfDate();
	}
	else if (object->GetName() == ICM_MODELMULTICURVES) 
	{	ICM_ModelMultiCurves* model = dynamic_cast<ICM_ModelMultiCurves*>(object);
		ircurve = model->GetZeroCurve(); 
		cpn_ircurve = model->GetCpnIRCurve(); 
		infcrv = model->GetCpnInfCurve();
		pVolcurve = model->GetVolCurve();
		//AsOf = ircurve->GetAsOfDate();
	}
	else{ ICM_MktDataMng* Mng = dynamic_cast<ICM_MktDataMng*>(object);
		string zc_name = GetZeroCurveName(leg.GetCreditInfos()->GetCcy(),AsOf);
		ircurve = dynamic_cast<ARM_ZeroCurve*>(Mng->find(zc_name)); 
		cpn_ircurve = dynamic_cast<ARM_ZeroCurve*>(Mng->find(zc_name)); 
	}
	
	//ARM_Date AsOf = ircurve->GetAsOfDate();
	ICM_Security* security = leg.GetCreditInfos();
	security->ComputeYF(AsOf);
	int size = security->GetAccStartDates().GetSize();
	int il = 0,findex=-1;
	double StartDateFwd=0.,EndDateFwd=0.,spread=0./** ,PartRate=1. */ , ResetDate =0., PayDate = 0., fwdpv01 = 1.0, unadjspread = 0.;
	ARM_Date FwdFixedDate;
	qCredit_Leg_Type CreditLegType = leg.GetCreditLegType();

	//use of external pricer for rates computing
	if (leg.GetRatesPricer()) 
	{
		ARM_Vector rates;
		leg.GetRatesPricer()->GenerateRates(rates);
		ARM_IRIndex* index = leg.GetIRIndex();

		if (size != rates.GetSize()) throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,"ERROR:Incompatible rates vector size");
		if (index==NULL) 
		{ 
			for (int il=0; il<size; il++) 
			{
				//JLA : this is no index given: 
				// fixed spread, fwd spread should be 0
				security->GetCouponRates().Elt(il)= rates.Elt(il) ;
			}
		}
		else if (!index->IsFixedIndex())
		{ 
			// JLA: this is pay only fwd spread, computed from the index 
			//	might pay fixed + spread, accordingly to GetFwdCalcTypes. 
			//	otherwise forced to qCPN_FWD...
		  if (leg.GetFwdCalcTypes()) security->SetFwdType(leg.GetFwdCalcTypes());
		  else security->SetFwdType(qCPN_FWD);

		  for (int il=0; il<size; il++)	{security->GetFwdSpreads().Elt(il)= rates.Elt(il);}
		}
		else
		{ 
			// this is a fixed index: use it (like fixed )
			for (int il=0; il<size; il++)	{security->GetCouponRates().Elt(il)= rates.Elt(il);}
		}
		security->updateFwdSpreadsFromPastFixings() ;
		return;
	}

	switch (CreditLegType)
	{
	case qNone_Leg :
	case qRunning_Leg :
	case qFunded_Leg :
	case qInfine_Leg :
	case qZeroCoupon_Leg :
	break;
	case qSwapLeg :
		{
		//ICM_ModelMultiCurves* mmc = dynamic_cast<ICM_ModelMultiCurves*>(model);
		//ARM_ZeroCurve* CpnIRcrv = mmc->GetCpnIRCurve(); 
		ARM_ZeroCurve* CpnIRcrv = cpn_ircurve;
		ARM_IRIndex* index = leg.GetIRIndex();
		if (index==NULL) break;
		if (CpnIRcrv)
		{	ARM_YCModel ycmodel(CpnIRcrv);
			leg.SetModelVariable(NULL); // bugfix : two successive calls can lead to the pointer value for ycmodel
			leg.SetModel(&ycmodel);

			if (size != leg.GetDecompRates()->GetSize()) throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,"ERROR:Incompatible rates vector size");

			if (!index->IsFixedIndex())
			{
				// JLA: this is now done when the CreditInfo is created. 
				// security->SetFwdType(qCPN_FIXED_AND_FWD);
			for (int il=0; il<size; il++) 
			{
				security->GetFwdSpreads().Elt(il)= leg.GetDecompRates()->Elt(il);
			}
			}
 			//else
			//{for (int il=0; il<size; il++) {security->GetCouponRates()->Elt(il)= leg.GetDecompRates()->Elt(il);}}
		}

		break;
		}
	case qForward_IRrates_Leg :
		{
		if (leg.GetFwdCalcTypes()) security->SetFwdType(leg.GetFwdCalcTypes());
			else security->SetFwdType(qCPN_FWD);

		ARM_ZeroCurve* CpnIRcrv = cpn_ircurve; 

		const ARM_Vector& YFFwdStartDates = security->GetYFFwdStartDates();
		const ARM_Vector& YFFwdAccEndDates = security->GetYFFwdEndDates();
		double FwdRate = 0.;

		for (int il=0; il<size; il++)
		{
			StartDateFwd = YFFwdStartDates.Elt(il);
			EndDateFwd = YFFwdAccEndDates.Elt(il);
			
			if (StartDateFwd>=EndDateFwd) continue;

			spread = CpnIRcrv->ForwardYield(StartDateFwd,EndDateFwd);
			security->GetFwdSpreads().Elt(il)=spread;
		}

		break;
		}
	case qConstant_Maturity_Leg :
		{
		
		//	JLA This is now done when the CreditInfo is created.... 
		//	
		// if (leg.GetFwdCalcTypes()) security->SetFwdType(leg.GetFwdCalcTypes());
		// 	else security->SetFwdType(qCPN_FWD);

		const ARM_Vector& YFFwdStartDates = security->GetYFFwdStartDates();
		const ARM_Vector& YFFwdEndDates = security->GetYFFwdEndDates();

		const ARM_Vector& FwdStartDates = security->GetFwdStartDates();
		const ARM_Vector& FwdEndDates = security->GetFwdEndDates();
		//leg.GetCreditIndex()->CptDefaultCurve(dynamic_cast<ICM_ModelMultiCurves*>(model));
		std::auto_ptr<const ICM_DefaultCurve> DefCurveIndex(leg.GetCreditIndex()->CptDefaultCurve(dynamic_cast<ICM_ModelMultiCurves*>(object)));
		// JLA Participation Rate is not part of forward spread
		// PartRate = leg.GetRefSpread();
		if (!pVolcurve) {
			// find in MarketDataMng
			ICM_MktDataMng* Mng = dynamic_cast<ICM_MktDataMng*>(object);
			if(!Mng) ICMTHROW(ERR_INVALID_ARGUMENT,"No Vol Curve for index"); 
			string IndexName = leg.GetCreditIndex()->GetIndexName();
			pVolcurve = Mng->GetVolCurve(IndexName, leg.GetCreditIndex()->GetCurrencyUnit()->GetCcyName(), AsOf);
			if (!pVolcurve) ICMTHROW(ERR_INVALID_ARGUMENT,"No Vol Curve for index"); 
		}
		if (leg.GetFwdFixedDate())
		{
			// FwdFixedDate = (ARM_Date)leg.GetFwdFixedDate();
			findex = leg.PeriodIndex(*leg.GetFwdFixedDate());
		}

		for (int il=0; il<size; il++)
		{
			if ((il>findex)&&(findex>=0)) {security->GetFwdSpreads().Elt(il) = security->GetFwdSpreads().Elt(findex);continue;} 

			YFStartDateFwd = YFFwdStartDates.Elt(il);
			YFEndDateFwd = YFFwdEndDates.Elt(il);

			if ( YFEndDateFwd <0 ) {
				if( leg.GetCreditIndex()->GetIsYT()) {
					YFEndDateFwd =  YFStartDateFwd + leg.GetCreditIndex()->GetYearTerm();			
				} else {
					YFEndDateFwd = FwdEndDates.Elt(il)/365;
				}
			}
			spread = DefCurveIndex->FwdSpread(YFStartDateFwd,YFEndDateFwd);
			spread += DefCurveIndex->AjustConvexity(YFStartDateFwd,YFEndDateFwd,spread,pVolcurve);
			security->GetFwdSpreads().Elt(il)=/** PartRate***/ spread/100.;
				
		}
		break;
	}
	case qCMCDS_Leg :
		{
			if (leg.GetFwdCalcTypes()) security->SetFwdType(leg.GetFwdCalcTypes());
			else security->SetFwdType(qCPN_FWD);

			const ARM_Vector& YFFwdStartDates = security->GetYFFwdStartDates();
			const ARM_Vector& YFFwdEndDates = security->GetYFFwdEndDates();
			const ARM_Vector& FwdStartDates = security->GetFwdStartDates();
			const ARM_Vector& FwdEndDates = security->GetFwdEndDates();

			const ARM_Vector& ResetDates = security->GetResetDates();
			const ARM_Vector& PayDates = security->GetPayDates();
			
			//leg.GetCreditIndex()->CptDefaultCurve(dynamic_cast<ICM_ModelMultiCurves*>(model));
			 std::auto_ptr<const ICM_DefaultCurve> DefCurveIndex(leg.GetCreditIndex()->CptDefaultCurve(dynamic_cast<ICM_ModelMultiCurves*>(object)));
			// PartRate = leg.GetRefSpread();

			if (leg.GetFwdFixedDate())
			{
				findex = leg.PeriodIndex(*leg.GetFwdFixedDate());
			}

			for (int il=0; il<size; il++)
			{
				if ((il>findex)&&(findex>=0))
				{
					security->GetAdjFwdSpreads().Elt(il) = security->GetAdjFwdSpreads().Elt(findex);
					security->GetFwdSpreads().Elt(il) = security->GetFwdSpreads().Elt(findex);
					security->GetUnAdjFwdSpreads().Elt(il) = security->GetUnAdjFwdSpreads().Elt(findex);
					security->GetFwdPV01s().Elt(il) = security->GetFwdPV01s().Elt(findex);
					continue;
				}

				YFStartDateFwd = YFFwdStartDates.Elt(il);
				YFEndDateFwd = YFFwdEndDates.Elt(il);

				if (YFEndDateFwd<0 ) {
					if( leg.GetCreditIndex()->GetIsYT()) {
						YFEndDateFwd =  YFStartDateFwd + leg.GetCreditIndex()->GetYearTerm();			
					} else {
						YFEndDateFwd = FwdEndDates.Elt(il)/365;
					}
				}
				spread = DefCurveIndex->FwdSpread(YFStartDateFwd,YFEndDateFwd);
				security->GetFwdSpreads().Elt(il)=/** PartRate***/ spread/100.;
				spread += DefCurveIndex->AjustConvexity(YFStartDateFwd,YFEndDateFwd,spread,pVolcurve);
				
				ResetDate = ResetDates.Elt(il);
				PayDate =	PayDates.Elt(il);

				if(ResetDate==PayDate)
					fwdpv01 = 1.0 ;
				else
				{
					fwdpv01 = DefCurveIndex->RiskyPV01(ResetDate,PayDate) ;
					spread = unadjspread + DefCurveIndex->PayLagAdjustment(StartDateFwd,EndDateFwd,unadjspread, pVolcurve);					
				}
				
				security->GetFwdPV01s().Elt(il) = fwdpv01 ;

				security->GetAdjFwdSpreads().Elt(il)=/** PartRate**/ spread/100.;
				security->GetFwdSpreads().Elt(il)=/** PartRate**/ spread/100.;
			}
			break;
		}
	case qInflation_Leg :
		{

		if (leg.GetFwdCalcTypes()) security->SetFwdType(leg.GetFwdCalcTypes());
			else security->SetFwdType(qCPN_FWD);

		const ARM_Vector& YFFwdStartDates = security->GetYFFwdStartDates();
		const ARM_Vector& YFFwdEndDates = security->GetYFFwdEndDates();
		ARM_Vector DiscPrices(size);
		double zcInterpol = 0.;
		int il=0,start=0;

		//recherche du 1er indice correspondant à la date suivant la date AsoF
		for (il=0; il<size; il++)
		{security->GetFwdSpreads().Elt(il)= 0.;
		 if (YFFwdEndDates.Elt(il)>0.) {start=il;break;}}

		for (il=start; il<size; il++)
		{
			zcInterpol = infcrv->ZCRateInterpolate((ARM_Date)(AsOf.GetJulian()+365.*YFFwdEndDates.Elt(il)))/100.;
			DiscPrices.Elt(il)=1./pow(1.+ zcInterpol,YFFwdEndDates.Elt(il));
		}

		for (il=start+1; il<size; il++)
		{security->GetFwdSpreads().Elt(il)= 100.*(DiscPrices.Elt(il-1)/DiscPrices.Elt(il) - 1.);}
		security->GetFwdSpreads().Elt(start)= 100.*(1./DiscPrices.Elt(start) - 1.);
		
		break;
		}
	case qStandart_Recovery_Leg :
		break;
	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : Unknown LegType ");
	}
	security->updateFwdSpreadsFromPastFixings(/*PartRate*/) ;
}


//	-----------------------------------------------------------------------
void 
ICM_Leg::setPastFixing(const ARM_Date&resetDate,double value)
{
	itsPastFixings[resetDate.GetJulian()]=value; 
	ResetCreditInfos(); 
}
void ICM_Leg::Set(const qCredit_Leg_Type& LegType ,
			 ICM_Credit_Index*	CreditIndex )
	{
		SetCreditLegType(LegType);
		
		if (CreditIndex)
		{
		if (itsCreditIndex)
			delete itsCreditIndex;
		itsCreditIndex = (ICM_Credit_Index*)CreditIndex->Clone();
		}
		ResetCreditInfos();
	}
//	-----------------------------------------------------------------------
// virtual 
void ICM_Leg::SetVariableSpread(ARM_ReferenceValue* spreads)
{
	SetLegType(K_FLOATING_LEG);
	ARM_SwapLeg::SetVariableSpread(spreads);

	ICM_Security* security = GetCreditInfos(); 
	int NbFlows = GetNumFlows();
	ARM_Vector CouponRates(NbFlows,0.);
	for (int i=0;i<NbFlows;i++)
	{
		CouponRates.Elt(i) = spreads->CptReferenceValue((ARM_Date)GetFlowStartDates()->Elt(i));
	}
	security->SetCouponRates(CouponRates);
}	
