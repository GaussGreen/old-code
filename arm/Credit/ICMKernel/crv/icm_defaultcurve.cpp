
#include "ARMKernel/glob/firsttoinc.h" 
#include "ICMKernel\crv\icm_defaultcurve.h"
#include "ICMKernel\pricer\icm_pricer_cds.h"
#include "ICMKernel\util\icm_brentsolver.h"
#include "ARMKernel\util\interpol.h"
#include "ARMKernel\crv\volcube.h"
#include "ICMKernel\mod\icm_defcurvemodel.h"
#include "ICMKernel\util\icm_rootfinder1D.h"

#include <set>
// ----------------------------------------------------------------------------------------------------------------
// Converting Summit Adjustment maturity code to ARM adjustment maturity
// ----------------------------------------------------------------------------------------------------------------

qCDS_ADJ FromSummitAdjCDSToARM(const char* AdjCode) 
{
    char temp[20];
    strcpy(temp,AdjCode);
    char* tmp = ARM_UPPERCASE(temp);

	std::string list ;
	bool ok; 
	qCDS_ADJ ret; 
	ICM_EnumsCnv::cnv(tmp,ret,ok,list); 
	if (!ok) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Summit "<<AdjCode<<" not in " << list ) ;
	return ret; 

}

// ----------------------------------------------------------------------------------------------------------------
// Ajoute une période de type '5Y' à une date avec un ajustement Business (AdjorNot) et un type de calendrier (adj)
// ----------------------------------------------------------------------------------------------------------------
ARM_Date AddPeriod(const ARM_Date& date, char*Matu , const std::string& ccy , bool AdjorNot, qCDS_ADJ adj)
{
	ARM_Date datretour;
	int nbperiod=0,i=0;
	char TypeMatu;
	bool flg = true;

	if (!Matu) {return datretour;}

	sscanf(Matu, "%d%c", &nbperiod, &TypeMatu);
    TypeMatu = toupper(TypeMatu);

	datretour = date;

	//Ajustement calendrier CDS
	if (adj!=qCredit_Default) 
	{ 
		flg=false;
		datretour = DateRoll(datretour,adj);  
	}

    if ( TypeMatu == 'D' ) // Ex : "1D"
    {
		datretour.AddDays(nbperiod);
    }
    else if ( TypeMatu == 'W' )  //  Ex : "1W"   
    {   
		datretour.AddDays(7*nbperiod);
    }
    else if ( TypeMatu == 'M' ) //  Ex : "9M"
    {  
		datretour.AddMonths(nbperiod);
    }
    else if ( TypeMatu == 'Y')  // Ex : "2Y"
    {   
		datretour.AddYears(nbperiod);
	}
	else
		ICMTHROW(ERR_INVALID_ARGUMENT,"AddPeriod::Invalid Plot "<<TypeMatu);
	
	if ((AdjorNot) && (flg))
		datretour.AdjustToBusDate((char*)ccy.c_str(),1 /* next business day */);

	return datretour;
}


//------------------------------------------------
// Initialization of datas members
//------------------------------------------------
void ICM_DefaultCurve::Init(void)
{
	SetName(ICM_DEFAULTCURVE);

	itsZeroCurve = NULL;
 	
	itsDates = NULL;
	itsYearTerms = NULL;
	itsNbDays = NULL;
	itsRates = NULL;
	itsLambdas = NULL;
	itsLambdas = NULL;
	itsSurvivalProba = NULL;

	itsRecovery = 0.;
	itsCurrency = "";  
	its_pricer = NULL;
	its_Current_indice = 0;

	itsLagStartDate = ARM_Currency("EUR").GetCreditStartDateLag();		//Default Summit Case
	itsLagProtectStartDate = itsLagStartDate; //Default Summit Case

	itsAdjBusiness = true;
	itsCdsAdj= qCredit_Default;

	its_STDCDS_Frequency = K_QUARTERLY;
	its_STDCDS_Basis = KACTUAL_360;
	its_STDCDS_Notional = DEFAULT_NOMINAL;
	its_STDCDS_Currency= "";  
	its_STDCDS_Stub = K_SHORTSTART;
	its_STDCDS_CreditLag = 30;
	its_STDCDS_AccOnDef = qACCRUED_SETTLED;

	its_STDCDS_AdjustedStartDateOnBusinessDay = K_ADJUSTED;
	its_STDCDS_IncludeMaturity = true;
	its_STDCDS_Adjusted = K_ADJUSTED;

	itsIsSummitCurve = true;
	itsIsNameInDefault = false;

	itsVolCurve = NULL;

	itsDirectionalLambdasShift	=	NULL;


	std::vector<ICM_Cds*>::iterator it = itsCalibrationInstruments.begin() ;
	while (it != itsCalibrationInstruments.end()) { delete *it ;++it; }
	itsCalibrationInstruments.clear(); 
	itsCalibrationAlgo= qDEFCURVE_DICHO ;
	itsCalibrationData="STD"; 
	itsInterpolYF.Resize(0);	
	itsInterpolLambda.Resize(0); 
	itsInterpolDates.Resize(0); 
	itsInterpolRates.Resize(0); 
	itsInterpolSP.Resize(0); 
	itsUpFronts = NULL;
	itsAmortRefValues.clear();
}
 
void ICM_DefaultCurve::SetAmortRefValues(vector<ARM_ReferenceValue*>&	AmortRefValues)
{
	int i=0;

	if (itsAmortRefValues.size()>0)
	{
		for (i=0;i<itsAmortRefValues.size();i++)
		{
		if (itsAmortRefValues[i])
			delete itsAmortRefValues[i];
		itsAmortRefValues[i]=NULL;
		}
	}
	itsAmortRefValues.resize(AmortRefValues.size());
	for (i=0;i<itsAmortRefValues.size();i++)
	{
		if (AmortRefValues[i])
			itsAmortRefValues[i]=(ARM_ReferenceValue*) AmortRefValues[i]->Clone();
	}

}

void ICM_DefaultCurve::SetUpFronts(ARM_Vector* UpFront)
{
	if (itsUpFronts)
		delete itsUpFronts;
	itsUpFronts = (ARM_Vector*) UpFront->Clone();
}

ICM_DefaultCurve::ICM_DefaultCurve(const ICM_DefaultCurve&ref) : ARM_Object(ref) 
{

  Init();	

  itsAsOfDate = ref.itsAsOfDate;
  itsZeroCurve=NULL; 
  if (ref.itsZeroCurve)
     itsZeroCurve = (ARM_ZeroCurve*) ref.itsZeroCurve->Clone();

  itsTerms=ref.itsTerms; 
  itsRecovery = ref.itsRecovery;
  itsCurrency=ref.itsCurrency; 
  itsCdsAdj = ref.itsCdsAdj; 
  itsAdjBusiness = ref.itsAdjBusiness;
  itsLagStartDate = ref.itsLagStartDate;
  itsLagProtectStartDate = ref.itsLagProtectStartDate;
  SetLabel(ref.itsLabel);
  itsDates = NULL;
  if (ref.itsDates)
     itsDates = (ARM_Vector*) ref.itsDates->Clone();
  itsNbDays = NULL;
  if (ref.itsNbDays)
     itsNbDays = (ARM_Vector*) ref.itsNbDays->Clone();
  itsYearTerms = NULL;
  if (ref.itsYearTerms)
     itsYearTerms = (ARM_Vector*) ref.itsYearTerms->Clone();
  itsRates = NULL;
  if (ref.itsRates)
     itsRates = (ARM_Vector*) ref.itsRates->Clone();
  itsLambdas = NULL;
  if (ref.itsLambdas)
     itsLambdas = (ARM_Vector*) ref.itsLambdas->Clone();
  itsSurvivalProba = NULL;
  if (ref.itsSurvivalProba)
     itsSurvivalProba = (ARM_Vector*) ref.itsSurvivalProba->Clone();
  its_Current_indice = ref.its_Current_indice;
  its_pricer = ref.its_pricer;

  // Members Standard Cds
  its_STDCDS_StartDate = ref.its_STDCDS_StartDate;
  its_STDCDS_Frequency = ref.its_STDCDS_Frequency;
  its_STDCDS_Basis = ref.its_STDCDS_Basis;
  its_STDCDS_Notional = ref.its_STDCDS_Notional;
  its_STDCDS_Currency =ref.its_STDCDS_Currency ;
  its_STDCDS_Stub = ref.its_STDCDS_Stub;
  its_STDCDS_CreditLag = ref.its_STDCDS_CreditLag;
  its_STDCDS_Adjusted = ref.its_STDCDS_Adjusted;
  its_STDCDS_AccOnDef = ref.its_STDCDS_AccOnDef;
  itsIsSummitCurve = ref.itsIsSummitCurve;

  itsVolCurve=NULL; 
  if (ref.itsVolCurve)
	itsVolCurve = (ARM_VolCurve*) ref.itsVolCurve->Clone();

  its_STDCDS_ProtectionStartDate = ref.its_STDCDS_ProtectionStartDate;	
  its_STDCDS_IncludeMaturity = ref.its_STDCDS_IncludeMaturity;
  its_STDCDS_AdjustedStartDateOnBusinessDay = ref.its_STDCDS_AdjustedStartDateOnBusinessDay;
  itsIsNameInDefault = ref.itsIsNameInDefault;

  itsDirectionalLambdasShift = NULL;

  if (ref.itsDirectionalLambdasShift)
     itsDirectionalLambdasShift = (ARM_Vector*) ref.itsDirectionalLambdasShift->Clone();


	itsCalibrationInstruments.resize(ref.itsCalibrationInstruments.size()); 
	std::vector<ICM_Cds*>::const_iterator it = ref.itsCalibrationInstruments.begin(); 
	int i=0; 
	while (it != ref.itsCalibrationInstruments.end()) 
	{
		itsCalibrationInstruments[i] = dynamic_cast<ICM_Cds*>((*it)->Clone()); 
		++it; 
		++i; 
	}
  
	itsInterpolYF=ref.itsInterpolYF;	
	itsInterpolLambda=ref.itsInterpolLambda; 
	itsInterpolDates=ref.itsInterpolDates; 
	itsInterpolRates=ref.itsInterpolRates; 
	itsInterpolSP=ref.itsInterpolSP; 
	itsCalibrationData=ref.itsCalibrationData; // temp string
	itsCalibrationAlgo =ref.itsCalibrationAlgo; 

	if (ref.itsUpFronts)
     itsUpFronts = (ARM_Vector*) ref.itsUpFronts->Clone();

	itsAmortRefValues.resize(ref.itsAmortRefValues.size());
	for (i=0;i<ref.itsAmortRefValues.size();i++)
	{itsAmortRefValues[i] = (ARM_ReferenceValue*) ref.itsAmortRefValues[i]->Clone();}
	itsParameters=ref.itsParameters; 
}


void ICM_DefaultCurve::ResetLambda(const int& indice,const double& lambda)
{
	ICMTHROW(ERR_UNIMP_METHOD_CALL,"Unimplemented ResetLambda for"<<typeid(*this).name()); 
 		
}

// ---------------------------------------------------------------------
// Evaluate NPV of default pricer & default Index for an intensity value
// ---------------------------------------------------------------------

double ICM_DefaultCurve::Evaluate(const double& x)
{
	// Variation du lambda
	ResetLambda(its_Current_indice, x );

	// NPV before variation
	its_pricer->ResetPricer();
	double npv = its_pricer->Price(qCMPPRICE);

	//calibration avec upfront
	if (itsUpFronts)
	{npv += itsUpFronts->Elt(its_Current_indice)*its_STDCDS_Notional;}
			
	return (npv);
}
 
//--------------------------------------------
//	View Method
//--------------------------------------------
void ICM_DefaultCurve::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[40];
 
    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

   if (itsZeroCurve) itsZeroCurve->View(id, fOut);

   fprintf(fOut, "\n\n------------------------------------------------------------------\n");
   fprintf(fOut, "---------------- Basis Default Curve -----------------------------\n");
   fprintf(fOut, "------------------------------------------------------------------\n\n");
   fprintf(fOut, "DefCurve Label :%s\n\n",itsLabel.c_str());	

   fprintf(fOut, "Date AsOf : %ld.%ld.%ld\n\n", itsAsOfDate.GetDay(), itsAsOfDate.GetMonth(), itsAsOfDate.GetYear());

   fprintf(fOut, "STDCDS StartDate : %ld.%ld.%ld\n", its_STDCDS_StartDate.GetDay(), its_STDCDS_StartDate.GetMonth(), its_STDCDS_StartDate.GetYear());
   fprintf(fOut, "STDCDS Protect StartDate : %ld.%ld.%ld\n", its_STDCDS_ProtectionStartDate.GetDay(), its_STDCDS_ProtectionStartDate.GetMonth(), its_STDCDS_ProtectionStartDate.GetYear());
   fprintf(fOut, "STDCDS Frequency : %ld.\n", its_STDCDS_Frequency);
   fprintf(fOut, "STDCDS Basis : %ld.\n", its_STDCDS_Basis);
   fprintf(fOut, "STDCDS Currency : %s\n", its_STDCDS_Currency.c_str() );
   fprintf(fOut, "STDCDS Stub : %ld.\n", its_STDCDS_Stub);
   fprintf(fOut, "STDCDS CreditLag : %ld.\n", its_STDCDS_CreditLag);
   fprintf(fOut, "STDCDS Adjusted : %ld.\n", its_STDCDS_Adjusted);

   string sCdsAdj;
   ICM_EnumsCnv::toString(itsCdsAdj, sCdsAdj);
   fprintf(fOut, "STDCDS Calendard Type : %s\n", sCdsAdj.c_str());
   if (its_STDCDS_IncludeMaturity)
	fprintf(fOut, "STDCDS Include Maturity : Include\n");
   else
	fprintf(fOut, "STDCDS Include Maturity : Not Include\n");

   if (its_STDCDS_AdjustedStartDateOnBusinessDay)
	fprintf(fOut, "STDCDS Adjusted Start Date onto business day: Yes\n");
   else
	fprintf(fOut, "STDCDS Adjusted Start Date onto business day: No\n");

   fprintf(fOut, "Recovery : %f\n", itsRecovery);
   if (GetIsNameInDefault()) fprintf(fOut, "IsNameIsDefault: TRUE\n");
   else fprintf(fOut, "IsNameIsDefault: FALSE\n");

   fprintf(fOut,"Calibration %s [Data:%s]\n",ICM_EnumsCnv::toString(itsCalibrationAlgo).c_str(),itsCalibrationData.c_str()) ; 
   fprintf(fOut, "\n\nDate\t\tNbDays\tYearTerms\tRate (bp)\t\tIntensity\tSurvivalProba\n");



   for (int i = 0; i < itsRates->GetSize() ; i++)
   {
	fprintf(fOut, "\n");	
	fprintf(fOut, "%ld.%ld.%ld\t", ((ARM_Date)GetDate(i)).GetDay(), ((ARM_Date)GetDate(i)).GetMonth(), ((ARM_Date)GetDate(i)).GetYear());
	fprintf(fOut, "%ld\t", (int)GetNbDay(i));
	fprintf(fOut, "%f\t", GetYearTerm(i));
	fprintf(fOut, "%f\t", GetRate(i)*10000.);
	fprintf(fOut, "%f\t", GetLambda(i));
	fprintf(fOut, "%f\t", GetSurvivalProba(i));
   }		 
	fprintf(fOut, "\n"); 

   for ( i = 0; i < itsInterpolDates.GetSize() ; i++)
   {
	fprintf(fOut, "\n");	
	fprintf(fOut, "%ld.%ld.%ld\t", ((ARM_Date)itsInterpolDates.Elt(i)).GetDay(), ((ARM_Date)itsInterpolDates.Elt(i)).GetMonth(), ((ARM_Date)itsInterpolDates.Elt(i)).GetYear());
	fprintf(fOut, "-- \t");
	fprintf(fOut, "%f\t", itsInterpolYF.Elt(i));
	fprintf(fOut, "%f\t", itsInterpolRates.Elt(i)*10000.);
	fprintf(fOut, "%f\t", itsInterpolLambda.Elt(i));
	fprintf(fOut, "%f\t", itsInterpolSP.Elt(i));
   }		 
	if ( ficOut == NULL )
    {
       fclose(fOut);
    }
 
}


void ICM_DefaultCurve::Set (const ARM_Date& asOf,
							const ARM_Vector& Dates,
							const ARM_Vector& Intensity,
							double  Recovery,
							ARM_ZeroCurve* zc,  	
							int intRule, 
							int adjStartDate,
							qCDS_ADJ adj ,  
							const std::string& ccy,
							const std::string& label,
							const ARM_VolCurve* VolCurve,
							qDEFCURVE_CALIB_ALGO calibAlgo,
							const std::string &calibrationData,
							int LagStartDate,
							const ICM_Parameters& params)
{
	itsCalibrationData=calibrationData; 
	itsCalibrationAlgo=calibAlgo; 
	std::vector<ICM_Cds*>::iterator it = itsCalibrationInstruments.begin() ;
	while (it != itsCalibrationInstruments.end()) { delete *it ;++it; }
	itsCalibrationInstruments.clear(); 

	itsAsOfDate = asOf;
	if (itsZeroCurve) delete itsZeroCurve;
	itsZeroCurve= (ARM_ZeroCurve*) zc->Clone();

	// Default Parameters for Cds
	int lag = LagStartDate; //ARM_Currency(ccy.c_str()).GetCreditStartDateLag();
	itsLagProtectStartDate = lag;

	its_STDCDS_ProtectionStartDate = itsAsOfDate;
	its_STDCDS_ProtectionStartDate.AddDays(itsLagProtectStartDate)/*,ccy /*ccy->GetCcyName() */;

	itsLagStartDate = LagStartDate;
	its_STDCDS_StartDate = itsAsOfDate;
	its_STDCDS_StartDate.AddDays(itsLagStartDate) /*,ccy  )*/;

	its_STDCDS_Frequency = K_QUARTERLY;
	its_STDCDS_Basis = KACTUAL_360;
	its_STDCDS_Notional = DEFAULT_NOMINAL;
	its_STDCDS_Currency =ccy ;
	its_STDCDS_Stub = K_SHORTSTART;
	its_STDCDS_CreditLag = 30;
	its_STDCDS_Adjusted = intRule; // K_ADJUSTED;
	its_STDCDS_AdjustedStartDateOnBusinessDay = adjStartDate; 
	itsCdsAdj = adj ;

	int size = Dates.GetSize();

	if (itsDates) delete itsDates;
	itsDates = new ARM_Vector(Dates.GetSize()+1,0.);
	memcpy(itsDates->GetElt()+1,unconst(Dates).GetElt(),sizeof(double)*Dates.GetSize());

	if (itsLambdas) delete itsLambdas;


	itsLambdas = new ARM_Vector(Dates.GetSize()+1,0.);
	memcpy(itsLambdas->GetElt(),unconst(Intensity).GetElt(),sizeof(double)*Intensity.GetSize());
	itsLambdas->Elt(size)=Intensity.Elt(size-1);

	if (itsLambdas->Elt(0)>=1000.) itsIsNameInDefault = true;
	

	if (itsNbDays)
		delete itsNbDays;
	itsNbDays = new ARM_Vector(size+1,0.);

	if (itsYearTerms)
		delete itsYearTerms;
	itsYearTerms = new ARM_Vector(size+1,0.);

	if (itsRates)
		delete itsRates;
	itsRates = new ARM_Vector(size+1,CREDIT_DEFAULT_VALUE);

	if (itsSurvivalProba)
		delete itsSurvivalProba;
	itsSurvivalProba = new ARM_Vector(size+1,0.);

	//Création du point AsOf de la courbe
	itsDates->Elt(0) =itsAsOfDate.GetJulian();
	itsNbDays->Elt(0)=0.;
	itsYearTerms->Elt(0)=0.;
	itsSurvivalProba->Elt(0)=1.;

	itsRecovery = Recovery;

	itsCurrency =ccy; 

	SetLabel(label);
	itsCdsAdj = adj;

	//Construction des autres points 
	for (int i = 1; i<size+1; i++)
	{
		ARM_Date Date	= GetDate(i);
		int Nbday		= (int) (Date - itsAsOfDate); 
		double YearTerm= (Date.GetJulian() - itsAsOfDate.GetJulian())/K_YEAR_LEN;

		itsNbDays->Elt(i)=(double)Nbday;
		itsYearTerms->Elt(i)=YearTerm;
	}


	if (VolCurve)
		itsVolCurve = (ARM_VolCurve*) unconst(*VolCurve).Clone();

	// In this case (where no calibration is done) .. 	
	itsInterpolYF = *itsYearTerms;	
	itsInterpolLambda =*itsLambdas; 
	itsInterpolDates =*itsDates; 
	itsInterpolRates =*itsRates; 
	itsInterpolSP =*itsSurvivalProba; 

	CptTermsSurvivalProba();

	*itsSurvivalProba=itsInterpolSP; 
	*itsLambdas=itsInterpolLambda; 
	itsParameters=params; 
}


void ICM_DefaultCurve::Set (const ARM_Date& asOf,
							const ARM_Vector* Dates,
							const ARM_Vector* DefCurve,
							double& Recovery,
							ARM_ZeroCurve* zc,  							
							int intRule,
							int adjStartDate,
							qCDS_ADJ adj ,  
							bool IsDefCurveSource,	
							const std::string& ccy,
							const std::string& label,
							// ICM_DefaultCurve* defcurve,
							const ARM_VolCurve* VolCurve,
							qDEFCURVE_CALIB_ALGO calibAlgo,
							const std::string& calibrationData,
							int LagStartDate)
{
	itsCalibrationData=calibrationData; 
	itsCalibrationAlgo=calibAlgo; 
	std::vector<ICM_Cds*>::iterator it = itsCalibrationInstruments.begin() ;
	while (it != itsCalibrationInstruments.end()) { delete *it ;++it; }
	itsCalibrationInstruments.clear(); 

	itsAsOfDate = asOf;
	if (itsZeroCurve) delete itsZeroCurve;
	itsZeroCurve= (ARM_ZeroCurve*) zc->Clone();

	// Default Parameters for Cds
	int lag =  LagStartDate ; //ARM_Currency(ccy.c_str()).GetCreditStartDateLag();
	itsLagProtectStartDate = lag;

	its_STDCDS_ProtectionStartDate = itsAsOfDate;
	its_STDCDS_ProtectionStartDate.AddDays(itsLagProtectStartDate)/*,ccy)*/;

	itsLagStartDate = LagStartDate;
	its_STDCDS_StartDate = itsAsOfDate;
	its_STDCDS_StartDate.AddDays(itsLagStartDate)/*,ccy )*/;

	its_STDCDS_Frequency = K_QUARTERLY;
	its_STDCDS_Basis = KACTUAL_360;
	its_STDCDS_Notional = DEFAULT_NOMINAL;
	its_STDCDS_Currency =ccy; 
	its_STDCDS_Stub = K_SHORTSTART;
	its_STDCDS_CreditLag = 30;
	its_STDCDS_Adjusted = intRule; // K_ADJUSTED;
	its_STDCDS_AdjustedStartDateOnBusinessDay=adjStartDate ;
	itsCdsAdj = qCredit_Adjust20;

	int size = Dates->GetSize();

	if (itsDates) delete itsDates;
	itsDates = new ARM_Vector(Dates->GetSize()+1,0.);
	memcpy(itsDates->GetElt()+1,Dates->GetElt(),sizeof(double)*Dates->GetSize());

	if (itsLambdas) delete itsLambdas;
	itsLambdas = new ARM_Vector(size+1,CREDIT_DEFAULT_VALUE);


	if (itsNbDays)
		delete itsNbDays;
	itsNbDays = new ARM_Vector(size+1,0.);

	if (itsYearTerms)
		delete itsYearTerms;
	itsYearTerms = new ARM_Vector(size+1,0.);

	if (itsRates)
		delete itsRates;
	itsRates = new ARM_Vector(size+1,CREDIT_DEFAULT_VALUE);

	if (itsSurvivalProba)
		delete itsSurvivalProba;
	itsSurvivalProba = new ARM_Vector(size+1,0.);

	//Création du point AsOf de la courbe
	itsDates->Elt(0) =itsAsOfDate.GetJulian();
	itsNbDays->Elt(0)=0.;
	itsYearTerms->Elt(0)=0.;
	itsSurvivalProba->Elt(0)=1.;

	itsRecovery = Recovery;

	itsCurrency=ccy; 

	SetLabel(label);
	itsCdsAdj = adj;

	//Construction des autres points 
	for (int i = 1; i<size+1; i++)
	{
		ARM_Date Date	= GetDate(i);
		int Nbday		= (int) (Date - itsAsOfDate); 
		double YearTerm= (Date.GetJulian() - itsAsOfDate.GetJulian())/K_YEAR_LEN;
		
		itsLambdas->Elt(i)=CREDIT_DEFAULT_VALUE;
		itsNbDays->Elt(i)=(double)Nbday;
		itsYearTerms->Elt(i)=YearTerm;
	}

	if (IsDefCurveSource)
	{//Construction des autres points pour la DefCurve 
		memcpy(itsSurvivalProba->GetElt()+1,DefCurve->GetElt(),sizeof(double)*DefCurve->GetSize());
	}

	if (VolCurve)
		itsVolCurve = (ARM_VolCurve*) unconst(*VolCurve).Clone();

	// In this case (where no calibration is done) .. 	
	itsInterpolYF = *itsYearTerms;	
	itsInterpolLambda =*itsLambdas; 
	itsInterpolDates =*itsDates; 
	itsInterpolRates =*itsRates; 
	itsInterpolSP =*itsSurvivalProba; 

	CptTermsSurvivalProba();

	*itsSurvivalProba=itsInterpolSP; 
	*itsLambdas=itsInterpolLambda; 
}


void ICM_DefaultCurve::Set (const ARM_Date& asOf,
							double* YearTerms,
							ARM_Vector* Intensity,
							double& Recovery,
							ARM_ZeroCurve* zc,  
							int intRule,
							int adjStartDate,
							qCDS_ADJ adj ,  
							const std::string& Ccy,
							const std::string&  label,
							const ARM_VolCurve* VolCurve,
							qDEFCURVE_CALIB_ALGO calibAlgo,
							const std::string&calibrationData,
							int LagStartDate)
{
	itsCalibrationData=calibrationData; 
	itsCalibrationAlgo=calibAlgo; 	
	std::vector<ICM_Cds*>::iterator it = itsCalibrationInstruments.begin() ;
	while (it != itsCalibrationInstruments.end()) { delete *it ;++it; }
	itsCalibrationInstruments.clear(); 

	double YearTerm = 0.;
	itsAsOfDate = asOf;
	if (itsZeroCurve) delete itsZeroCurve;
	itsZeroCurve= (ARM_ZeroCurve*) zc->Clone();

	// Default Parameters for Cds
	int lag = LagStartDate;//ARM_Currency(Ccy.c_str()).GetCreditStartDateLag();
	itsLagProtectStartDate = lag;
	itsLagStartDate = lag;

	its_STDCDS_StartDate = itsAsOfDate;
	its_STDCDS_StartDate.AddDays(itsLagStartDate)/*,Ccy)*/;

	its_STDCDS_ProtectionStartDate = itsAsOfDate;
	its_STDCDS_ProtectionStartDate.AddDays(itsLagProtectStartDate)/*,Ccy )*/;

	its_STDCDS_Adjusted=intRule; 
	its_STDCDS_AdjustedStartDateOnBusinessDay=adjStartDate; 

	its_STDCDS_Currency =Ccy; 

	int size = Intensity->GetSize();

	if (itsDates) delete itsDates;
		itsDates = new ARM_Vector(size+1,asOf.GetJulian());

	if (itsLambdas) delete itsLambdas;
		itsLambdas = new ARM_Vector(size+1,0.);
	memcpy(itsLambdas->GetElt(),Intensity->GetElt(),sizeof(double)*Intensity->GetSize());
	itsLambdas->Elt(size)=Intensity->Elt(size-1);

	if (itsLambdas->Elt(0)>=1000.) itsIsNameInDefault = true;

	if (itsNbDays)
		delete itsNbDays;
	itsNbDays = new ARM_Vector(size+1,0.);

	if (itsYearTerms)
		delete itsYearTerms;
	itsYearTerms = new ARM_Vector(size+1,0.);

	if (itsRates)
		delete itsRates;
	itsRates = new ARM_Vector(size+1,CREDIT_DEFAULT_VALUE);

	if (itsSurvivalProba)
		delete itsSurvivalProba;
	itsSurvivalProba = new ARM_Vector(size+1,0.);

	//Création du point AsOf de la courbe
	itsDates->Elt(0)=itsAsOfDate.GetJulian();
	itsNbDays->Elt(0)=0.;
	itsYearTerms->Elt(0)=0.;
	itsSurvivalProba->Elt(0)=1.;

	itsRecovery = Recovery;

	itsCurrency=Ccy; 

	SetLabel(label);
	itsCdsAdj = adj;

	//Construction des autres points 
	for (int i =1; i<size+1; i++)
	{
		YearTerm= YearTerms[i-1];
		itsDates->Elt(i) =YearTerm*365.+itsAsOfDate.GetJulian();

		ARM_Date Date	= GetDate(i);
		int Nbday		= (int) (Date.GetJulian() - itsAsOfDate.GetJulian()); 

		itsNbDays->Elt(i)=(double)Nbday;
		
		itsYearTerms->Elt(i)=YearTerm;
	}


	if (VolCurve)
		itsVolCurve = (ARM_VolCurve*) unconst(*VolCurve).Clone();

	// In this case (where no calibration is done) .. 	
	itsInterpolYF = *itsYearTerms;	
	itsInterpolLambda =*itsLambdas; 
	itsInterpolDates =*itsDates; 
	itsInterpolRates =*itsRates; 
	itsInterpolSP =*itsSurvivalProba; 

	CptTermsSurvivalProba();

	*itsSurvivalProba=itsInterpolSP; 
	*itsLambdas=itsInterpolLambda; 
}


// -------------------------------------------------------------------
// Constructor with terms & Spread levels
// -------------------------------------------------------------------
/** 
ICM_DefaultCurve::ICM_DefaultCurve(const ARM_Date& asOf,
						 const std::vector<std::string>& terms,
						 ARM_Vector* rates,
						 double& Recovery,
						 ARM_ZeroCurve* zc,
						 int intRule,
						 int adjStartDate,
						 qCDS_ADJ adj,  	
						 const std::string& Ccy,
						 const std::string&  label,
						 bool issummitcurve,
						 const ARM_VolCurve* VolCurve,
						 long PayFrq,
						qDEFCURVE_CALIB_ALGO calibAlgo,						 
						 const std::string& calibrationMethod,
						 int LagStartDate,
						 const ICM_Parameters&params,
						 ARM_Vector* UpFronts,
						vector<ARM_ReferenceValue*>& AmortRefValues)
{

	Init();

	Set(asOf,terms,rates,Recovery,zc,intRule,adjStartDate,adj,Ccy,label,issummitcurve,VolCurve,PayFrq, 
		calibAlgo,calibrationMethod,LagStartDate,params,UpFronts,AmortRefValues);
		
}
**/ 
void ICM_DefaultCurve::Set (const ARM_Date& asOf,
					   const std::vector<std::string>& terms,
                       ARM_Vector* rates,
					   double Recovery,
					   ARM_ZeroCurve* zc,
					   int intRule,
					   int adjStartDate,
					   qCDS_ADJ adj,  	
                       const std::string& Ccy,
					   const std::string&  label,
					   bool issummitcurve,
					   const ARM_VolCurve* VolCurve,
					   long PayFrq,
						qDEFCURVE_CALIB_ALGO calibAlgo,					   
					   const std::string &calibrationData,
					   int LagStartDate,
					   const ICM_Parameters&params,
					   ARM_Vector* UpFronts,
					   vector<ARM_ReferenceValue*>& AmortRefValues)
{
	int i=0;

	itsCalibrationData=calibrationData; 
	itsCalibrationAlgo=calibAlgo; 
	bool charTermKn = true;
	std::vector<ICM_Cds*>::iterator it = itsCalibrationInstruments.begin() ;
	while (it != itsCalibrationInstruments.end()) { delete *it ;++it; }
	itsCalibrationInstruments.clear(); 

	itsAsOfDate = asOf;
	if (itsZeroCurve) delete itsZeroCurve;
	itsZeroCurve= (ARM_ZeroCurve*) zc->Clone();

	// Default Parameters for Cds
	int lag = LagStartDate; //ARM_Currency(Ccy.c_str()).GetCreditStartDateLag();
	itsLagProtectStartDate = lag;
	itsLagStartDate = lag;
	its_STDCDS_ProtectionStartDate = itsAsOfDate;
	its_STDCDS_ProtectionStartDate.AddDays(itsLagProtectStartDate)/*,Ccy)*/;

	its_STDCDS_StartDate = itsAsOfDate;
	its_STDCDS_StartDate.AddDays(itsLagStartDate)/*,Ccy)*/;

	its_STDCDS_Adjusted = intRule ; 
	its_STDCDS_AdjustedStartDateOnBusinessDay = adjStartDate; 

	its_STDCDS_Frequency = PayFrq;
	
	itsTerms=terms; 

	if (itsRates) delete itsRates;
	itsRates = new ARM_Vector(rates->GetSize()+1,0.);
	memcpy(itsRates->GetElt()+1,rates->GetElt(),sizeof(double)*rates->GetSize());
	itsRates->Elt(0)=rates->Elt(0);

	if (itsRates->Elt(0)>=10.) itsIsNameInDefault = true; //1e5 Bp eq to default

	// ----------------------------------------
	// gestion de l'upfront
	// ----------------------------------------
	if (UpFronts)
	{
		if (itsUpFronts) delete itsUpFronts;
		itsUpFronts = new ARM_Vector(UpFronts->GetSize()+1,0.);
		memcpy(itsUpFronts->GetElt()+1,UpFronts->GetElt(),sizeof(double)*UpFronts->GetSize());
		itsUpFronts->Elt(0)=UpFronts->Elt(0);
	}

	// ----------------------------------------
	// gestion de l'armortissement
	// ----------------------------------------
	if (AmortRefValues.size()>0)
	{
		if (itsAmortRefValues.size()>0)
		{
		for (i=0;i<itsAmortRefValues.size();i++)
		{
		if (itsAmortRefValues[i])
			delete itsAmortRefValues[i];
		itsAmortRefValues[i]=NULL;
		}
		}
		itsAmortRefValues.resize(AmortRefValues.size()+1);
		for (i=0;i<AmortRefValues.size();i++)
		{
		if (AmortRefValues[i])
			itsAmortRefValues[i+1]=(ARM_ReferenceValue*) AmortRefValues[i]->Clone();
		}
		itsAmortRefValues[0]=(ARM_ReferenceValue*) AmortRefValues[0]->Clone();
	}

	itsRecovery = Recovery;

	itsCurrency=Ccy;  
	
	its_STDCDS_Currency =itsCurrency ;
	SetLabel(label);

	itsCdsAdj = adj;
	itsIsSummitCurve = issummitcurve;

	if (VolCurve)
		itsVolCurve = (ARM_VolCurve*) unconst(*VolCurve).Clone();

	itsParameters=params; 
	CptSetCreditMarketData(charTermKn);
	Calibrate();
	CptTermsSurvivalProba();
}

// ------------------------------------------------------------
// Computes Dates & YeraTerms
// ------------------------------------------------------------
void ICM_DefaultCurve::CptSetCreditMarketData(bool charTermKn)
{
	int i = 0;
	int size = itsRates->GetSize();
	bool dateVectorKn = false;
	
	if (!itsDates) {
		itsDates = new ARM_Vector(size,0.);
		dateVectorKn = false;
	}
	if (!itsYearTerms) { 
		itsYearTerms = new ARM_Vector(size,0.);
		dateVectorKn = true; // because if no YearTerms, it's DateVector or Chr** terms
	}

	if (itsLambdas)
		delete itsLambdas;
	itsLambdas = new ARM_Vector(size,0.);

	if (itsSurvivalProba)
		delete itsSurvivalProba;
	itsSurvivalProba = new ARM_Vector(size,0.);
	if (itsNbDays)
		delete itsNbDays;
	itsNbDays = new ARM_Vector(size,0.);


	//Création du point AsOf de la courbe
	itsDates->Elt(0)=itsAsOfDate.GetJulian();
	itsNbDays->Elt(0)=0.;
	itsYearTerms->Elt(0)=0.;
	itsLambdas->Elt(0)=0.;
	itsSurvivalProba->Elt(0)=1.;

	bool doAdjust(itsAdjBusiness); 
	if (its_STDCDS_Adjusted==K_MATUNADJUSTED || its_STDCDS_Adjusted==K_UNADJUSTED) doAdjust=false;
	else doAdjust=true; 

	//Construction des autres points du set
	ARM_Date Date;
	int Nbday=0;
	double YearTerm=0.;


	if (!charTermKn)
		AdjustDateYearsTerms(dateVectorKn);
	else { // case of adj from itsTerms	
		if (itsDates) 
			delete itsDates;
		itsDates = new ARM_Vector(size,0.);
		itsDates->Elt(0)=itsAsOfDate.GetJulian();
		
		if (itsYearTerms) 
			delete itsYearTerms;
		itsYearTerms = new ARM_Vector(size,0.);
		itsYearTerms->Elt(0)=0.;

		for (i = 1; i<size; i++)
		{
			if (itsCdsAdj == qCredit_Default) Date = AddPeriod(its_STDCDS_StartDate, itsTerms[i-1].c_str(),itsCurrency,doAdjust,itsCdsAdj);
			else Date = AddPeriod(itsAsOfDate, itsTerms[i-1].c_str(),itsCurrency,doAdjust,itsCdsAdj);
			
			Nbday		= (int) (Date - itsAsOfDate); 
			YearTerm= (Date.GetJulian() - itsAsOfDate.GetJulian())/K_YEAR_LEN;

			itsDates->Elt(i)=Date.GetJulian();
			itsNbDays->Elt(i)=(double)Nbday;
			itsYearTerms->Elt(i)=YearTerm;
		}
	}


}



//virtual 
void ICM_DefaultCurve::Calibrate()
{
	switch (itsCalibrationAlgo)
	{
	case qDEFCURVE_DICHO: 
	case qDEFCURVE_NEWTON:
		 OldCalibrate(); 
		break; 
	case qDEFCURVE_BRENT:
		Calibrate_Stress_Test_Guess_Brent();  
		break ; 
	default:
		ICMTHROW(ERR_INVALID_ARGUMENT,"Can't calibrate "<<ICM_EnumsCnv::toString(itsCalibrationAlgo)); 
	}
}

// virtual 
void ICM_DefaultCurve::OldCalibrate()
{
	if ((itsRecovery > 1.0) || (itsRecovery < 0.0))
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_DefaultCurve::Calibrate: Recov not in [0,1] :"<<itsRecovery); 
	
	//	JLA 
	//
	//		Calibration code is working over itsInterpolXXX vectors
	if (itsCalibrationData=="STD") 
	{
		// old case : make sure we have the correct values
		itsInterpolYF = *itsYearTerms;	
		itsInterpolLambda=*itsLambdas; 
		itsInterpolDates=*itsDates; 
		itsInterpolRates=*itsRates; 
		itsInterpolSP=*itsSurvivalProba ;
	}
	else 
	{
		// here we add all intermediates points with term interval. 
		string maturity( "yYMmwWdD" );
		string::size_type pos = itsCalibrationData.find_first_of( maturity );
		if( pos == string::npos )
			ICMTHROW(ERR_INVALID_ARGUMENT,"Can't understand "<<itsCalibrationData); 
		string nbString    = itsCalibrationData.substr( 0, pos );
		int nb            = atoi( nbString.c_str() );
		string matuString    = itsCalibrationData.substr( pos, string::npos );
		if( matuString.size() != 1 )
			ICMTHROW(ERR_INVALID_ARGUMENT,"Can't understand "<<itsCalibrationData); 
		ARM_Date to; 
		int iPeriod=1 ;
		std::set<double> aggregation ;
		// to = itsAsOfDate ;
		to = itsDates->Elt(1); 
		bool doAdjust(itsAdjBusiness); 
		if (its_STDCDS_Adjusted==K_MATUNADJUSTED || its_STDCDS_Adjusted==K_UNADJUSTED) doAdjust=false;
		else doAdjust=true; 
		ARM_Date maxto = AddPeriod(its_STDCDS_StartDate,"10Y",itsCurrency,doAdjust,itsCdsAdj);
		while (to.GetJulian() < itsDates->Elt( itsDates->GetSize()-1) && to < maxto ) 
		{
			if (to.GetJulian()>itsDates->Elt(1)) aggregation.insert(to.GetJulian()) ; 
			std::stringstream sstr; sstr<<iPeriod*nb<<matuString; 
			if (itsCdsAdj == qCredit_Default)
				to = AddPeriod(its_STDCDS_StartDate,sstr.str(),itsCurrency,doAdjust,itsCdsAdj);
			else 		
				to = AddPeriod(itsAsOfDate,sstr.str(),itsCurrency,doAdjust,itsCdsAdj); 	
			iPeriod++;   
		} 
		for(unsigned int i =0; i<itsDates->GetSize(); i++) 
			aggregation.insert(itsDates->Elt(i)); 
		//
		unsigned int newsize = aggregation.size() ;
		itsInterpolDates.Resize(newsize); 
		itsInterpolYF.Resize(newsize) ; 
		itsInterpolRates.Resize(newsize); 
		std::set<double>::const_iterator it ; 
		for(i=0, it = aggregation.begin() ;i<newsize; i++,it++) 
		{
			itsInterpolDates.Elt(i)=*it; 
			itsInterpolYF.Elt(i)=(itsInterpolDates.Elt(i) - itsAsOfDate.GetJulian())/K_YEAR_LEN ;
			itsInterpolRates.Elt(i)= linInterpol(itsDates,itsInterpolDates.Elt(i),itsRates); 
		}
		// output: 
		itsInterpolLambda.Resize(newsize); 
		itsInterpolSP.Resize(newsize); 
		itsInterpolSP.Elt(0)=1; 
	}
	
	//
	//	here is the place where the itsInterpolSearch is constructed
	double Nominal = DEFAULT_NOMINAL;
	double Nominal0 = its_STDCDS_Notional;
	// double rate = itsRates->Elt(1);
	double rate = itsInterpolRates.Elt(1);

	// Init Lambda
	double Lambda;
	Lambda	=	rate;
	if (1.0-itsRecovery)
		Lambda	/=	(1.0-itsRecovery);

	double result = 1.;

	// itsLambdas->Elt(0)=Lambda;
	SetInterpolLambda(0,Lambda); 
	ICM_DefaultCurveModel* model = NULL;

	model = new ICM_DefaultCurveModel(this,GetZeroCurve(),NULL,false);

	double _inf = 0.;
	double _sup = 0.;

	int	nb_iter=0;

	// En sortie de cette boucle on a calibre tous les lambdas sauf le dernier
	// for (int it=1; it<itsDates->GetSize(); it++)
	for (int it=1; it<itsInterpolDates.GetSize(); it++)
	{
		its_Current_indice = it;

		if (itsIsNameInDefault) {ResetLambda(its_Current_indice,1000.);continue;}

		ARM_Date EndDate = itsInterpolDates.Elt(it); // GetDate(it); 
 		rate = itsInterpolRates.Elt(it);				// GetRate(it);

		if (itsAmortRefValues.size()>0)
		{its_STDCDS_Notional = itsAmortRefValues[0]->CptReferenceValue(0.);}

 		ICM_Cds cds(its_STDCDS_StartDate,
							   EndDate,
							   0,
							   0,
							   its_STDCDS_ProtectionStartDate,
							   EndDate,
							   rate,
							   ( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[it] : ARM_ReferenceValue(its_STDCDS_Notional)) ,
							   ( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[it] : ARM_ReferenceValue(its_STDCDS_Notional)) ,
							   its_STDCDS_Frequency,
							   its_STDCDS_Basis,
							   its_STDCDS_AccOnDef,
							   its_STDCDS_Currency ,
							   its_STDCDS_Stub,
							   its_STDCDS_CreditLag,
							   -1,	
							   its_STDCDS_Adjusted,
							   its_STDCDS_IncludeMaturity,
							   its_STDCDS_AdjustedStartDateOnBusinessDay,
							   std::string(),
							   qRunning_Leg,
							   qStandart_Recovery_Leg,
								ISSUER_UNDEFINE,
								CREDIT_DEFAULT_VALUE 
							   );	


		ICM_Pricer_Cds pricer; pricer.Set(&cds,model,itsParameters,itsAsOfDate);

		pricer.SetFaster(true);
		its_pricer = &pricer;

		_inf = MAX(Lambda-1.E-2,0.0);
		_sup = Lambda+1.E-2;

		result = RootFinder1D(ff1::mem_call(&ICM_DefaultCurve::Evaluate,(*this))).ZeroBracketDecreasing(_inf,_sup);	
		if (!result) { 
			Lambda = 0.0;
			/*	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_DefaultCurve::Calibrate: can't bracket "
			<< itsLabel
			<< "from inf="<<_inf<<" sup="<<_sup); */
		
		} else {
		// since we can't check "dichotomy failure" vs "other exception" we catch any exception and consider this a failure
			Lambda=-1; 
			try {
				switch(itsCalibrationAlgo)
				{
				case qDEFCURVE_DICHO:
					Lambda = RootFinder1D(ff1::mem_call(&ICM_DefaultCurve::Evaluate,(*this))).Dichotomy(_inf,_sup,nb_iter,100,1.E-5,1.E-2);
					break; 
				case qDEFCURVE_NEWTON:
					Lambda = RootFinder1D(ff1::mem_call(&ICM_DefaultCurve::Evaluate,(*this))).NewtonRaphsonWithBisection(_inf,_sup,100,0.01,1.E-5,1.E-2);
					
					break ;
				default:
					ICMTHROW(ERR_INVALID_ARGUMENT,"Can't calibrate with "<<itsCalibrationAlgo); 
				}
				
				if ((Lambda<0.0 ) || (Lambda ==_BAD_RET_VALUE_)) 
					ICMLOG("ICM_DefaultCurve::Calibrate failed: "<< itsLabel<<" "<<EndDate<<", with rate : "<< rate << " -> Inf : " <<_inf << " => Lambda="<<Lambda); 
			}
			catch (std::exception&e)
			{
				ICMMSG(WARN,"ICM_DefaultCurve::Calibrate failed: "<<itsLabel<<" "<<EndDate<<",Lambda="<<Lambda
					<<"["<<e.what()<<"]"); 
			}
			catch (Exception&e)
			{
				ICMMSG(WARN,"ICM_DefaultCurve::Calibrate failed: "<<itsLabel<<" "<<EndDate<<",Lambda="<<Lambda
					<<"["<<e.GetMessage()<<"]"); 
			}
			catch(...) {
				ICMMSG(WARN,"ICM_DefaultCurve::Calibrate failed: "<<itsLabel<<" "<<EndDate<<",Lambda="<<Lambda
					<<"[unknown exception]"); 
			}
		}
		ResetLambda(its_Current_indice,Lambda);
	}

	if (model) 
		delete model;
	model = NULL;
	its_STDCDS_Notional=Nominal0;

	//	At the end of the calibration step :
	//
	if (itsCalibrationData=="STD") 
	{	
		*itsLambdas = itsInterpolLambda	; 
		*itsSurvivalProba = itsInterpolSP; 
	}
	else 
	{
		unsigned int oldsize = itsLambdas->size(); 
		for(unsigned int i=0;i<oldsize;i++) 
		{
			itsLambdas->Elt(i) = GetPWCIntensity( itsYearTerms->Elt(i) ) ;  
			itsSurvivalProba->Elt(i) = SurvivalProba( itsYearTerms->Elt(i) ) ; 
		}
	}
}




// -------------------------------------------------------------
// Computes Risky PV01 Given 2 dates
// -------------------------------------------------------------

double ICM_DefaultCurve::RiskyPV01(const ARM_Date& t1, const ARM_Date& t2) const 
{
	return RiskyDuration(t2) - RiskyDuration(t1); 

}
// -------------------------------------------------------------
// Computes Risky PV01 Given 2 dates
// -------------------------------------------------------------

double ICM_DefaultCurve::RiskyPV01(const ARM_Date& t1, const std::string& tenor) const 
{
	ARM_Date t2 = AddPeriod(itsAsOfDate,tenor,itsCurrency,itsAdjBusiness,itsCdsAdj);
	return RiskyDuration(t2) - RiskyDuration(t1); 
}

// -------------------------------------------------------------
// Computes Risky PV01 Given 2 dates
// -------------------------------------------------------------

double ICM_DefaultCurve::DefLegPV(const ARM_Date& t1, const ARM_Date& t2) const 
{
	double PV01 = 0;;
	double Nominal = 1.;
	double DefLegPV1 = 0.;
	double DefLegPV2 = 0.;

	double ImpliedSpread = ImpliedSpreadInterpol(t1 );
	ImpliedSpread /= 10000.;

	ICM_DefaultCurveModel* model = new ICM_DefaultCurveModel(this,GetZeroCurve());

	if (its_STDCDS_StartDate.GetJulian()<t1.GetJulian())
	{
	ICM_Cds cds(its_STDCDS_StartDate,t1,0,0 //ARM_DATE
				,its_STDCDS_ProtectionStartDate,t1,
				ImpliedSpread,
				( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(Nominal)),
				//ARM_ReferenceValue(Nominal),
				( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(Nominal)),
				//ARM_ReferenceValue(Nominal),
				its_STDCDS_Frequency,its_STDCDS_Basis,// Nominal,
				its_STDCDS_AccOnDef,
				its_STDCDS_Currency,// Nominal,
				its_STDCDS_Stub,its_STDCDS_CreditLag,-1,its_STDCDS_Adjusted,
				its_STDCDS_IncludeMaturity,its_STDCDS_AdjustedStartDateOnBusinessDay,
				std::string(),
				qRunning_Leg,
				qStandart_Recovery_Leg,
				ISSUER_UNDEFINE,CREDIT_DEFAULT_VALUE);	


	ICM_Pricer_Cds pricer; pricer.Set(&cds,model,itsParameters,itsAsOfDate);
	pricer.SetFaster(true);
	
	// DefLegPV1 = pricer.DefLegPV();
	DefLegPV1 = pricer.Price(qCMPDEFLEGPV);
	}

	ImpliedSpread = ImpliedSpreadInterpol(t2);
	ImpliedSpread /= 10000.;

	if (its_STDCDS_StartDate.GetJulian()<t2.GetJulian())
	{
	ICM_Cds cds(its_STDCDS_StartDate,t2,
		0,0, // ARM_Date
		its_STDCDS_ProtectionStartDate,t2,
			    ImpliedSpread,// "NULL",
				( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(Nominal)),
				//ARM_ReferenceValue(Nominal),
				( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(Nominal)),
				//ARM_ReferenceValue(Nominal),
				its_STDCDS_Frequency,its_STDCDS_Basis,its_STDCDS_AccOnDef,
				its_STDCDS_Currency,its_STDCDS_Stub,its_STDCDS_CreditLag,
				-1,its_STDCDS_Adjusted,its_STDCDS_IncludeMaturity,its_STDCDS_AdjustedStartDateOnBusinessDay,
				std::string(),
				qRunning_Leg,
				qStandart_Recovery_Leg,
				ISSUER_UNDEFINE,CREDIT_DEFAULT_VALUE);

	ICM_Pricer_Cds pricer; pricer.Set(&cds,model,itsParameters,itsAsOfDate);
	pricer.SetFaster(true);
	
	DefLegPV2 = pricer.Price(qCMPDEFLEGPV);
	}

	if (model) delete model;

	PV01 = DefLegPV2 - DefLegPV1;
	
	return (PV01);
}

// -------------------------------------------------------------
// Computes Forward Spread between 2 dates
// -------------------------------------------------------------
double ICM_DefaultCurve::FwdSpread(const ARM_Date& t1, const ARM_Date& t2) const 
{
	double fwdsprd = 0;;
	double Nominal = 1.;
	double rate = 0.;
	int indice = 0;
	bool equal = false;
	double Feepv1 = 0.;
	double Duration1 = 0.;
	double Feepv2 = 0.;
	double Duration2 = 0.;

	if ((t1 <= its_STDCDS_StartDate) || (t1 <= its_STDCDS_ProtectionStartDate))
		return ImpliedSpreadInterpol(t2);

	double ImpliedSpread = ImpliedSpreadInterpol(t1);
	ImpliedSpread /= 10000.;

	auto_ptr<ICM_DefaultCurveModel> model(new ICM_DefaultCurveModel(this,GetZeroCurve()));

	if (its_STDCDS_StartDate.GetJulian()<t1.GetJulian())
	{
	
	ICM_Cds cds(its_STDCDS_StartDate,t1,
		0,0, // ARM_Date 
		its_STDCDS_ProtectionStartDate,t1,
				ImpliedSpread,// "NULL",
				( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(Nominal)),
				//ARM_ReferenceValue(Nominal),
				( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(Nominal)),
				//ARM_ReferenceValue(Nominal),
				// 0,0,
				its_STDCDS_Frequency,its_STDCDS_Basis,its_STDCDS_AccOnDef,
				its_STDCDS_Currency,
				its_STDCDS_Stub,its_STDCDS_CreditLag,-1,its_STDCDS_Adjusted,
				its_STDCDS_IncludeMaturity,its_STDCDS_AdjustedStartDateOnBusinessDay,
							   std::string(),// payCalName = NULL ,
							   qRunning_Leg,
				qStandart_Recovery_Leg,
								ISSUER_UNDEFINE,//name = ISSUER_UNDEFINE ,
								CREDIT_DEFAULT_VALUE// Binary = CREDIT_DEFAULT_VALUE 
				);	

	// ICM_Pricer_Cds pricer(&cds,model.get());
	ICM_Pricer_Cds pricer; pricer.Set(&cds,model.get(),itsParameters,itsAsOfDate);
	pricer.SetFaster(true);

	Feepv1 = pricer.Price(qCMPFEELEGPV);
	Duration1 = pricer.Price(qCMPDURATION);
	}

	ImpliedSpread = ImpliedSpreadInterpol(t2);
	ImpliedSpread /= 10000.;

	if (its_STDCDS_StartDate.GetJulian()<t2.GetJulian())
	{
	ICM_Cds cds(its_STDCDS_StartDate,t2,0,0,// ARM_Date 
				its_STDCDS_ProtectionStartDate,t2,
				ImpliedSpread,// "NULL",
				( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(Nominal)),
				//ARM_ReferenceValue(Nominal),
				( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(Nominal)),
				//ARM_ReferenceValue(Nominal),// 0,0,
				its_STDCDS_Frequency,its_STDCDS_Basis,its_STDCDS_AccOnDef,
				its_STDCDS_Currency,its_STDCDS_Stub,its_STDCDS_CreditLag,
				-1,its_STDCDS_Adjusted,its_STDCDS_IncludeMaturity,its_STDCDS_AdjustedStartDateOnBusinessDay,
							   std::string(),// payCalName = NULL ,
							   qRunning_Leg,
				qStandart_Recovery_Leg,
								ISSUER_UNDEFINE,//name = ISSUER_UNDEFINE ,
								CREDIT_DEFAULT_VALUE// Binary = CREDIT_DEFAULT_VALUE 
				);

	// ICM_Pricer_Cds pricer(&cds,model.get());
	ICM_Pricer_Cds pricer; pricer.Set(&cds,model.get(),itsParameters,itsAsOfDate);
	pricer.SetFaster(true);
	
	Feepv2 = pricer.Price(qCMPFEELEGPV);
	Duration2 = pricer.Price(qCMPDURATION);
	}


	fwdsprd = (Feepv2 - Feepv1)/(Duration2 - Duration1);
	fwdsprd *= 10000.;

	return (fwdsprd);
}

// -------------------------------------------------------------
// Computes Forward Spread between 2 yearterms
// -------------------------------------------------------------
double ICM_DefaultCurve::FwdSpread(const double& yt1, const double& yt2) const 
{
	double t1 = itsAsOfDate.GetJulian() + yt1*365.;
	double t2 = itsAsOfDate.GetJulian() + yt2*365.;

	ARM_Date	Date_t1(t1);
	ARM_Date	Date_t2(t2);

	return FwdSpread(Date_t1, Date_t2);
}

// -------------------------------------------------------------------------
// Compute Convexity ajustement between 2 yearterms
// -------------------------------------------------------------------------
double ICM_DefaultCurve::AjustConvexity(const double& yt1,const double& yt2, const double& FwdSpread, ARM_VolCurve* vol) const 
{ 
	double AjustConvexity = 0.;
	double SpreadFwd = FwdSpread/10000.;		// Pour ramener le Spread à sa vraie valeur
	double VolCurve = 0.;
	double term = 0.0;
	// On regarde d'abord si on a passé une vol en parametre,
	// si ce n'est pas le cas on prend celle de la courbe si elle existe
	ARM_VolCurve* lVolCurve = NULL;

	if (vol) lVolCurve = vol;
	else if(itsVolCurve) lVolCurve = itsVolCurve;
	else return AjustConvexity;

		
	ARM_VolCurve* VolATM = NULL;
	switch (lVolCurve->GetName())
	{
		case ARM_VOL_CUBE :   
		{
			ARM_VolCube* vCube = dynamic_cast<ARM_VolCube*> (lVolCurve);
			VolATM = vCube->GetATMVol();
			if( ! VolATM->hasIndex() )
					ICMTHROW(ERR_INVALID_ARGUMENT," no index in VolCurve ATM, cannot compute adjust convexity");
			term = VolATM->GetIndex().GetYearTerm();
		}
		break;
		case ARM_VOL_CURVE :
		case ARM_VOL_LIN_INTERPOL : 
		{
			VolATM = dynamic_cast<ARM_VolCurve*>(lVolCurve);
			if( ! VolATM->hasIndex() )
					ICMTHROW(ERR_INVALID_ARGUMENT," no index in VolCurve ATM, cannot compute adjust convexity");
			term = VolATM->GetIndex().GetYearTerm();
		}
		break;
		default :
		//case ARM_VOL_FLAT : 
			term = yt2-yt1;
	}

	VolCurve = (1/100.)*lVolCurve->ComputeVolatility(yt1,0.0, term);
	AjustConvexity = (yt2-yt1)/(2*(1-itsRecovery))*SpreadFwd*SpreadFwd*(exp(VolCurve*VolCurve*yt1) - 1);
			
	AjustConvexity *= 10000.;	

	return AjustConvexity;
}

// -------------------------------------------------------------------------
// Compute  Payment Lag Adjustment between 2 yearterms
// -------------------------------------------------------------------------

double ICM_DefaultCurve::PayLagAdjustment(const double& yt1, const double& yt2, const double& FwdSpread, ARM_VolCurve* volcurve) const 
{
	double PayLagAdjust = 0.;
	double SpreadFwd = FwdSpread/10000.;		// Pour ramener le Spread à sa vraie valeur
	double vol = 0.;

	ARM_ZeroCurve* InitShortCurve = (ARM_ZeroCurve*) GetZeroCurve();
	double ir_fwdshortrate = InitShortCurve->ForwardYield(yt1, yt2,0)/100.;
	double delta = 1.0/GetPayFrequency() ;
	double x = delta*(ir_fwdshortrate+SpreadFwd/(1.0-GetRecovery()));

	int N = (int)((yt2-yt1)/delta)+1;
	double cste = delta /(1.0-GetRecovery());

	double aux1 = 0., aux2=0., aux3 = 0;
	
	aux1 = cste * (1/(exp(x)-1) - N/(exp(N*x)-1) ) ;
// FIXMEFRED: mig.vc8 (28/05/2007 10:20:33):cast
	aux2 = pow(cste,2) * (1/(1-exp(x)) - 2*N/((exp(x)-1)*(exp(N*x)-1)) + pow(static_cast<double>(N),2)*(1+exp(N*x))/pow((exp(N*x)-1),2));

	aux3 = pow(cste,3)* (3*N*pow(exp(N*x)-1,2)+pow(exp(N*x)-1,3)+3*pow(static_cast<double>(N),2)*(exp(2*N*x)-1)-pow(static_cast<double>(N),3)*(exp(x)-1)*(1+4*exp(N*x)+exp(2*N*x)))/((exp(x)-1)*pow((exp(N*x)-1),3));
	

	// On regarde d'abord si on a passé une vol en parametre,
	// si ce n'est pas le cas on prend celle de la courbe si elle existe
	if (volcurve)
	{
		vol = (1/100.)*volcurve->ComputeVolatility(yt1, FwdSpread);

		PayLagAdjust = (aux1+SpreadFwd*aux2/2.) * SpreadFwd * SpreadFwd * (exp(vol*vol*yt1)-1.0);
	
		PayLagAdjust *= 10000.;	// Pour ramener le Spread en bp.

		return PayLagAdjust;
	}
	else
	{
		if (itsVolCurve)
		{
			//Prise en compte des courbes de vol 
			vol = (1/100.)*itsVolCurve->ComputeVolatility(yt1,FwdSpread);

			PayLagAdjust = (aux1+SpreadFwd*aux2/2.) * SpreadFwd * SpreadFwd * (exp(vol*vol*yt1)-1.0);
			
			PayLagAdjust *= 10000.;	

			return PayLagAdjust;
		}
		else
			return PayLagAdjust;
	}

}

// -------------------------------------------------------------------------
// Compute Convexity ajustement between 2 yearterms
// -------------------------------------------------------------------------
double ICM_DefaultCurve::AjustConvexity(const double& yt1,const double& yt2, const double& FwdSpread, double VolValue)
{ 
	double AjustConvexity = 0.;
	double SpreadFwd = FwdSpread/10000.;		// Pour ramener le Spread à sa vraie valeur
	double VolCurve = VolValue;
	
	AjustConvexity = (yt2-yt1)/(2*(1-itsRecovery))*SpreadFwd*SpreadFwd*(exp(VolCurve*VolCurve*yt1) - 1);
	
	AjustConvexity *= 10000.;	// Pour ramener le Spread en bp.

	return AjustConvexity;
}

// -------------------------------------------------------------------------
// Compute  Payment Lag Adjustment between 2 yearterms
// -------------------------------------------------------------------------

double ICM_DefaultCurve::PayLagAdjustment(const double& yt1, const double& yt2, const double& FwdSpread, double VolValue)
{
	if (CHECK_EQUAL(yt1, yt2))
		return 0.0;

	double PayLagAdjust = 0.;

	double SpreadFwd = FwdSpread/10000.;		// Pour ramener le Spread à sa vraie valeur
	double vol = 0.;

	ARM_ZeroCurve* InitShortCurve = (ARM_ZeroCurve*) GetZeroCurve();
	double ir_fwdshortrate = InitShortCurve->ForwardYield(yt1, yt2,0)/100.;
	double delta = 1.0/GetPayFrequency() ;
	double x = delta*(ir_fwdshortrate+SpreadFwd/(1.0-GetRecovery()));

	int N = (int)((yt2-yt1)/delta)+1;
	double cste = delta /(1.0-GetRecovery());

	double aux1 = 0., aux2=0., aux3 = 0;
	
	aux1 = cste * (1/(exp(x)-1) - N/(exp(N*x)-1) ) ;
// FIXMEFRED: mig.vc8 (28/05/2007 10:21:36):cast
	aux2 = pow(cste,2) * (1/(1-exp(x)) - 2*N/((exp(x)-1)*(exp(N*x)-1)) + pow(static_cast<double>(N),2)*(1+exp(N*x))/pow((exp(N*x)-1),2));

	aux3 = pow(cste,3)* (3*N*pow(exp(N*x)-1,2)+pow(exp(N*x)-1,3)+3*pow(static_cast<double>(N),2)*(exp(2*N*x)-1)-pow(static_cast<double>(N),3)*(exp(x)-1)*(1+4*exp(N*x)+exp(2*N*x)))/((exp(x)-1)*pow((exp(N*x)-1),3));
	
	PayLagAdjust = (aux1+SpreadFwd*aux2/2.) * SpreadFwd * SpreadFwd * (exp(VolValue*VolValue*yt1)-1.0);

	PayLagAdjust *= 10000.;	// Pour ramener le Spread en bp.

	return PayLagAdjust;
}

// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Bump a package of Default Curves according a given bump profile
// -------------------------------------------------------------------------
ICM_DefaultCurve** BumpPortfolio_ (const ICM_DefaultCurve** portfolio, 
									int nbcurves,
									qSENSITIVITY_TYPE typesensi,
									char* plot, 
									char* label,
									int& outNoCurve,
									double epsvalue)
{
    double sensitivity =0.;
	int NoCurve = -1;
	int i=0,h=0;

	ICM_DefaultCurve** ModifVectorDefCurves = new ICM_DefaultCurve*[nbcurves];	
	for (i=0;i<nbcurves;i++) ModifVectorDefCurves[i] = NULL;

	double result = 0.0;
	bool Parallelshift = false; 

	vector<string> Term;
	
	char TermZC[ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS];
	memset(TermZC,'\0',sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);

	// ARM_Vector* epsilon = NULL;

	double EPSL = 0.01;		//Bump sur les taux fixé à 1bp	
	double EPSL_DEF = 0.01; //Bump sur le spread fixé à 1 bp.
	double EPSL_REC = 0.1;  //Bump fixe sur le recovery de 0.1
	double EPSL_SPREL = 0.1 ;//Bump relatif des spreads de 0.1

	if (epsvalue != CREDIT_DEFAULT_VALUE)
		EPSL = EPSL_DEF = EPSL_REC = EPSL_SPREL = epsvalue;

	const ICM_DefaultCurve* InitDefCurve = NULL;
	ICM_DefaultCurve* ModifDefCurve = NULL;

	ARM_ZeroCurve* InitShortCurve = NULL;
	ARM_ZeroCurve* ModifShortCurve = NULL;

	ARM_ZeroCurve* InitShortCurve2 = NULL;
	ARM_ZeroCurve* ModifShortCurve2 = NULL;	

	ARM_ReferenceValue* Recovery = NULL;
	ARM_Vector* DiscreteValues = NULL;
	ARM_Vector* DiscreteDates = NULL;

	ARM_Date Date;

	ARM_Vector epsilon (1,epsvalue);

	if (!strcmp(plot,"NONE"))  //Detection d'un parallel Shift
	{
		Term.push_back("ALL");
		Parallelshift = true;
	}
	else
	{
		Term.push_back(plot);
		strcpy(TermZC[0],plot);
		switch (typesensi)
		{
			case ICMRECOVERY_TYPE :
			{
				epsilon.Elt(0)=EPSL_REC;
			}
			break;
			case ICMIRCURVE_TYPE :
			{
				epsilon.Elt(0)=EPSL;
			}
			break;
			case ICMSPRELSHIFT_TYPE :
			{
				epsilon.Elt(0)=EPSL_SPREL;
			}
			break;
			default:
			case ICMSPREAD_TYPE :
			{
				epsilon.Elt(0)=EPSL_DEF;
			}
			break;

		}
	}

		

	if (strcmp(label,"NONE") != NULL)  //Detection d'un label particulier
	{
		for (h=0; h<nbcurves; h++)
		{
			if (strcmp(label,portfolio[h]->GetLabel().c_str()) == NULL)
			{
			NoCurve = h;
			break;
			}
		}

		if (NoCurve == -1)
			ICMTHROW(ERR_INVALID_ARGUMENT,"BumpPortfolio_: nout found  "<<label); 

	for (h=0; h<nbcurves; h++)
		ModifVectorDefCurves[h] = (ICM_DefaultCurve*) portfolio[h]->Clone();

	}

	outNoCurve = NoCurve;

    
    {

		switch (typesensi)
		{
			case ICMRECOVERY_TYPE :
			{

				if (NoCurve>=0)
				{
					if (ModifVectorDefCurves[NoCurve])
						delete ModifVectorDefCurves[NoCurve];

					// if (Parallelshift)
					//	ModifVectorDefCurves[NoCurve] = (ICM_DefaultCurve*) portfolio[NoCurve]->GenerateShiftCurve(EPSL_REC,ICMRECOVERY_TYPE);
					// else
					ModifVectorDefCurves[NoCurve] = (ICM_DefaultCurve*) portfolio[NoCurve]->GenerateShiftCurve(Term,epsilon,ICMRECOVERY_TYPE);
					ModifVectorDefCurves[NoCurve]->SetLabel(portfolio[NoCurve]->GetLabel());

				}
				else
				{
					for (i = 0; i <nbcurves; i++)
					{
					if (ModifVectorDefCurves[i])
						delete ModifVectorDefCurves[i];

					// if (Parallelshift)
					// 	ModifVectorDefCurves[i] = (ICM_DefaultCurve*) portfolio[i]->GenerateShiftCurve(EPSL_REC,ICMRECOVERY_TYPE);
					// else
					ModifVectorDefCurves[i] = (ICM_DefaultCurve*) portfolio[i]->GenerateShiftCurve(Term,epsilon,ICMRECOVERY_TYPE);

					ModifVectorDefCurves[i]->SetLabel(portfolio[i]->GetLabel());
					}
				}

			}
			break;
			case ICMIRCURVE_TYPE :
			{

				if (NoCurve>=0)
				{
					if (ModifVectorDefCurves[NoCurve])
						delete ModifVectorDefCurves[NoCurve];

					ModifVectorDefCurves[NoCurve] = (ICM_DefaultCurve*) portfolio[NoCurve]->Clone();
					InitShortCurve = (ARM_ZeroCurve*) portfolio[NoCurve]->GetZeroCurve()->Clone();

					if (Parallelshift)
					{	
						ModifShortCurve = (ARM_ZeroCurve*)(InitShortCurve->Clone());
						ModifShortCurve->ParallelShift(EPSL);
					}
					else
						ModifShortCurve = (ARM_ZeroCurve*) InitShortCurve->GenerateShiftCurve(TermZC,&epsilon);


					// Impact sur la DefProbCurve  ********************************
					ModifVectorDefCurves[NoCurve]->SetZeroCurve((ARM_ZeroCurve*)ModifShortCurve->Clone());

					ModifVectorDefCurves[NoCurve]->Calibrate();
					
					ModifVectorDefCurves[NoCurve]->SetLabel(portfolio[NoCurve]->GetLabel());

					if (InitShortCurve)
						delete InitShortCurve;
					InitShortCurve = NULL;

				}
				else
				{
					for (i = 0; i <nbcurves; i++)
					{

					if (ModifVectorDefCurves[i])
						delete ModifVectorDefCurves[i];

					ModifVectorDefCurves[i] = (ICM_DefaultCurve*) portfolio[i]->Clone();
					InitShortCurve = (ARM_ZeroCurve*) ModifVectorDefCurves[i]->GetZeroCurve()->Clone();

					if (Parallelshift)
					{
						ModifShortCurve = (ARM_ZeroCurve*) InitShortCurve->Clone();
						ModifShortCurve->ParallelShift(EPSL);
					}
					else
						ModifShortCurve = (ARM_ZeroCurve*) InitShortCurve->GenerateShiftCurve(TermZC,&epsilon);


					// Impact sur la DefProbCurve  ********************************
					ModifVectorDefCurves[i]->SetZeroCurve((ARM_ZeroCurve*) ModifShortCurve->Clone());

					ModifVectorDefCurves[NoCurve]->Calibrate();
					
					ModifVectorDefCurves[i]->SetLabel(portfolio[i]->GetLabel());

					if (InitShortCurve)
						delete InitShortCurve;
					InitShortCurve = NULL;

					}
				}

			}
			break;
			case ICMSPRELSHIFT_TYPE :
			{

				if (NoCurve>=0)
				{
					if (ModifVectorDefCurves[NoCurve])
						delete ModifVectorDefCurves[NoCurve];

					InitDefCurve = (ICM_DefaultCurve*)portfolio[NoCurve]->Clone();

					// if (Parallelshift)
					// 	ModifVectorDefCurves[NoCurve] = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(EPSL_DEF,ICMSPRELSHIFT_TYPE);
					// else
					ModifVectorDefCurves[NoCurve] = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(Term,epsilon,ICMSPRELSHIFT_TYPE);

					ModifVectorDefCurves[NoCurve]->SetLabel(portfolio[NoCurve]->GetLabel());

					if (InitDefCurve)
						delete InitDefCurve;
					InitDefCurve = NULL;

				}
				else
				{
					for (i = 0; i <nbcurves; i++)
					{
					if (ModifVectorDefCurves[i])
						delete ModifVectorDefCurves[i];

					InitDefCurve = (ICM_DefaultCurve*)portfolio[i]->Clone();

					// if (Parallelshift)
					// 	ModifVectorDefCurves[i] = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(EPSL_DEF,ICMSPRELSHIFT_TYPE);
					// else
					ModifVectorDefCurves[i] = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(Term,epsilon,ICMSPRELSHIFT_TYPE);

					ModifVectorDefCurves[i]->SetLabel(portfolio[i]->GetLabel());	

					if (InitDefCurve)
						delete InitDefCurve;
					InitDefCurve = NULL;
					}
				}
				
			}
			break;
			case ICMSPREAD_TYPE :
			{

				if (NoCurve>=0)
				{
					if (ModifVectorDefCurves[NoCurve])
						delete ModifVectorDefCurves[NoCurve];


					InitDefCurve = (ICM_DefaultCurve*)portfolio[NoCurve]->Clone();

					// if (Parallelshift)
					// 	ModifVectorDefCurves[NoCurve] = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(EPSL_DEF,ICMSPREAD_TYPE);
					// else
					ModifVectorDefCurves[NoCurve] = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(Term,epsilon,ICMSPREAD_TYPE);

					ModifVectorDefCurves[NoCurve]->SetLabel(portfolio[NoCurve]->GetLabel());

					if (InitDefCurve)
						delete InitDefCurve;
					InitDefCurve = NULL;

				}
				else
				{
					for (i = 0; i <nbcurves; i++)
					{
					if (ModifVectorDefCurves[i])
						delete ModifVectorDefCurves[i];

					InitDefCurve = portfolio[i];

					// if (Parallelshift)
					// ModifVectorDefCurves[i] = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(EPSL_DEF,ICMSPREAD_TYPE);
					// else
					ModifVectorDefCurves[i] = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(Term,epsilon,ICMSPREAD_TYPE);

					ModifVectorDefCurves[i]->SetLabel(portfolio[i]->GetLabel());	

					InitDefCurve = NULL;

					}
				}
				
			}
			break;
			default :
			return (NULL);
		}

	// if (epsilon)
	// 	delete epsilon;
	// epsilon = NULL;

	if (InitDefCurve)
		delete InitDefCurve;
	InitDefCurve = NULL;

	if (ModifDefCurve)
		delete ModifDefCurve;
	ModifDefCurve = NULL;

	if (InitShortCurve)
		delete InitShortCurve;
	InitShortCurve = NULL;

	if (ModifShortCurve)
		delete ModifShortCurve;
	ModifShortCurve = NULL;

	if (InitShortCurve2)
		delete InitShortCurve2;
	InitShortCurve2 = NULL;

	if (ModifShortCurve2)
		delete ModifShortCurve2;
	ModifShortCurve2 = NULL;

	}
 
	return (ModifVectorDefCurves);

}


//	----------------------------------------------------------------------------------------------------
double ICM_DefaultCurve::ImpliedSpreadInterpol(const std::string& plot) const 
{
	return ImpliedSpreadInterpol(AddPeriod(itsAsOfDate,plot,itsCurrency,itsAdjBusiness,itsCdsAdj));
}
//	----------------------------------------------------------------------------------------------------
double ICM_DefaultCurve::ImpliedSpreadInterpol(const ARM_Date& date) const
{
	ICM_Cds cds(its_STDCDS_StartDate,
			    date,
				0,0, // ARM_Date
			    its_STDCDS_ProtectionStartDate,
			    date,	
			    0.01,
				( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(its_STDCDS_Notional)),
				//ARM_ReferenceValue(its_STDCDS_Notional),
				( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(its_STDCDS_Notional)),
				//ARM_ReferenceValue(its_STDCDS_Notional),// 0,0,
			    its_STDCDS_Frequency,
			    its_STDCDS_Basis,
			    its_STDCDS_AccOnDef,
			    its_STDCDS_Currency,
			    its_STDCDS_Stub,
			    its_STDCDS_CreditLag,
			    -1,
			    its_STDCDS_Adjusted,
			    its_STDCDS_IncludeMaturity,
			    its_STDCDS_AdjustedStartDateOnBusinessDay,
							   std::string(),// payCalName = NULL ,
							   qRunning_Leg,
				qStandart_Recovery_Leg,
								ISSUER_UNDEFINE,//name = ISSUER_UNDEFINE ,
								CREDIT_DEFAULT_VALUE// Binary = CREDIT_DEFAULT_VALUE 
				);

	ICM_DefaultCurveModel DefProbModel(&unconst(*this),GetZeroCurve());
	ICM_Pricer_Cds DefPricer; DefPricer.Set(&cds,&DefProbModel,itsParameters,itsAsOfDate);
	DefPricer.SetFaster(true);
	double result=DefPricer.ComputeSpread(0.);
	return (result);
}

 
// ----------------------------------------------------------------------------------
// Returns RiskyDuration for a given date
// ----------------------------------------------------------------------------------
double ICM_DefaultCurve::RiskyDuration(const ARM_Date& date) const 
{
	double Nominal = 1.;
	double Duration = 0.;

	// double Rate = ImpliedSpreadInterpol(date);
	// Rate /= 10000.;
	double Rate = 1./10000. ; // 1 bp 

	ICM_DefaultCurveModel model(&unconst(*this),GetZeroCurve());

	// if (its_STDCDS_StartDate.GetJulian()>= date.GetJulian())
	// 	ICMTHROW(ERR_INVALID_ARGUMENT,"AsOfDate + Adj > Maturity");  
       
	if (date>=itsAsOfDate && date<=its_STDCDS_StartDate) return 0; 
	if (date<itsAsOfDate ) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_DefaultCurve::RiskyDuration:"<<date<<"<"<<itsAsOfDate);  

	ICM_Cds cds(its_STDCDS_StartDate,date,0,0,its_STDCDS_ProtectionStartDate,
				  date,Rate,// "NULL",
				  ( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(Nominal)),
				  ( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(Nominal)),
				  //ARM_ReferenceValue(Nominal),ARM_ReferenceValue(Nominal),// 0,0,
				  its_STDCDS_Frequency,its_STDCDS_Basis,its_STDCDS_AccOnDef,
				  its_STDCDS_Currency,its_STDCDS_Stub,
				  its_STDCDS_CreditLag,-1,its_STDCDS_Adjusted,its_STDCDS_IncludeMaturity,
				  its_STDCDS_AdjustedStartDateOnBusinessDay,
							   std::string(),// payCalName = NULL ,
							   qRunning_Leg,
				qStandart_Recovery_Leg,
								ISSUER_UNDEFINE,//name = ISSUER_UNDEFINE ,
								CREDIT_DEFAULT_VALUE// Binary = CREDIT_DEFAULT_VALUE 
				  );	

	// ICM_Pricer_Cds pricer(&cds,&model);
	ICM_Pricer_Cds pricer; pricer.Set(&cds,&model,itsParameters,itsAsOfDate);
	pricer.SetFaster(true);

	Duration = pricer.Price(qCMPDURATION);

	return (Duration);
}


// ----------------------------------------------------------------------------------
// Returns RiskyDuration for a given date
// ----------------------------------------------------------------------------------
double ICM_DefaultCurve::RiskyDuration(const std::string & Tenor) const 
{
	ARM_Date date = AddPeriod(itsAsOfDate,Tenor,itsCurrency,itsAdjBusiness,itsCdsAdj); 
	return RiskyDuration(date); 
 
}

ARM_Security* 
ICM_DefaultCurve::stdCDS(const ARM_Date& t1, const ARM_Date& t2,double spread,double notional,bool useimpliedstartdate) const
{
	ICM_DefaultCurveModel model(this,GetZeroCurve());
	ICM_Cds* cds = NULL;
	
	ARM_Date Start = t1;

	if (useimpliedstartdate) { Start = its_STDCDS_StartDate; }

	cds = new ICM_Cds(Start,t2,0,0,Start,t2,
					  spread,// "NULL",
					  ( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(notional)),
					  ( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(notional)),
					  //ARM_ReferenceValue(notional),ARM_ReferenceValue(notional),//0,0,
					  its_STDCDS_Frequency,
					  its_STDCDS_Basis,its_STDCDS_AccOnDef,
					  its_STDCDS_Currency,its_STDCDS_Stub,
					  its_STDCDS_CreditLag,-1,its_STDCDS_Adjusted,its_STDCDS_IncludeMaturity,
					  its_STDCDS_AdjustedStartDateOnBusinessDay,
							   std::string(),// payCalName = NULL ,
							   qRunning_Leg,
								qStandart_Recovery_Leg,
								ISSUER_UNDEFINE,//name = ISSUER_UNDEFINE ,
								CREDIT_DEFAULT_VALUE// Binary = CREDIT_DEFAULT_VALUE 
					  );

	return ((ARM_Security*)cds) ;

}

ICM_DefaultCurve* ICM_DefaultCurve::createDefCurveFromBase(const ICM_DefaultCurve* DefCurveIndex ,const  ARM_Vector& vBase) const 
{
	std::auto_ptr<ICM_DefaultCurve> newDefCurve(dyn_clone(DefCurveIndex));	
	int sizeIndex  = newDefCurve->GetRates()->size();
	if(sizeIndex != (vBase.size()+1)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Vector of base is different from the nb of spread in the defCurveIndex "<<newDefCurve->GetRates()->size());
	ARM_Vector vDates = *DefCurveIndex->GetDates();
		
	ARM_Vector newSpread(sizeIndex);
	for(int i=1; i<sizeIndex; i++)	{
			newSpread[i] = this->ImpliedSpreadInterpol(vDates[i]);
			if (vBase[i-1] == 0) 
				ICMTHROW(ERR_INVALID_ARGUMENT,"element %i of Vector of base is Null "<<i-1);
			newSpread[i] = newSpread[i]/( vBase[i-1] *10000.0);
	}
	newSpread[0] = newSpread[1];
	newDefCurve->setRates(newSpread);
	newDefCurve->SetLabel(this->GetLabel());
	newDefCurve->Calibrate();
	newDefCurve->CptTermsSurvivalProba();
	
	return newDefCurve.release();
}

// -------------------------------------------------------------------
// Constructor with Dates/YF & Spread levels
// -------------------------------------------------------------------
/**
ICM_DefaultCurve::ICM_DefaultCurve(const ARM_Date& asOf,
						 const ARM_Vector& YForDates,
						 const ARM_Vector& rates,
						 const double& Recovery,
						 ARM_ZeroCurve* zc,
						 int intRule,
						 int adjStartDate,
						 qCDS_ADJ adj ,
						 const std::string& Ccy,
						 const std::string&  label,
						 bool issummitcurve,
 						 const ARM_VolCurve* VolCurve,
						 const bool& isdatesininput,
						 long PayFrq,
 						qDEFCURVE_CALIB_ALGO calibAlgo,
						 const std::string &calibData,
						 int LagStartDate)
{

	Init();
	Set(asOf,YForDates,rates,Recovery,zc,intRule,adjStartDate,adj, Ccy,label,issummitcurve,
		VolCurve,isdatesininput,PayFrq, calibAlgo,calibData, LagStartDate);
		
}
**/ 
void ICM_DefaultCurve::Set (const ARM_Date& asOf,
					    const ARM_Vector& YForDates,
						const ARM_Vector& rates,
					   const double& Recovery,
					   ARM_ZeroCurve* zc,
					   int intRule,
					   int adjStartDate,
					   qCDS_ADJ adj ,
                       const std::string& Ccy,
					   const std::string&  label,
					   bool issummitcurve,
					   //2 ICM_DefaultCurve* defcurve,
					   const ARM_VolCurve* VolCurve,
					   const bool& isdatesininput,
					   long PayFrq,
					   // bool	isBrentCalib,
					   qDEFCURVE_CALIB_ALGO calibAlgo,
					   const std::string& calibData,
					   int LagStartDate)
{
	itsCalibrationData=calibData;
	itsCalibrationAlgo=calibAlgo;
	bool charTermKn = false;
	std::vector<ICM_Cds*>::iterator it = itsCalibrationInstruments.begin() ;
	while (it != itsCalibrationInstruments.end()) { delete *it ;++it; }
	itsCalibrationInstruments.clear(); 

	int i =0;
	itsAsOfDate = asOf;
	if (itsZeroCurve) delete itsZeroCurve;
	itsZeroCurve= (ARM_ZeroCurve*) zc->Clone();

	// Default Parameters for Cds
	int lag = LagStartDate; //ARM_Currency(Ccy.c_str()).GetCreditStartDateLag();
	itsLagProtectStartDate = lag;
	itsLagStartDate = lag;
	its_STDCDS_ProtectionStartDate = itsAsOfDate;
	//its_STDCDS_ProtectionStartDate.AddDays(itsLagProtectStartDate);
	its_STDCDS_ProtectionStartDate.AddDays(itsLagProtectStartDate)/*,Ccy)*/;

	its_STDCDS_StartDate = itsAsOfDate;
	//its_STDCDS_StartDate.AddDays(itsLagStartDate);
	its_STDCDS_StartDate.AddDays(itsLagStartDate)/*,Ccy)*/;

	its_STDCDS_Frequency = PayFrq;
//	itsConstructorType = qDEFCURVE_CSTR_WITH_DATES_OR_YF;

	if (itsRates) delete itsRates;
	itsRates = new ARM_Vector(rates.size()+1,0.);
	for (i=0;i<rates.size();i++) {itsRates->Elt(i+1)=rates[i];}
	itsRates->Elt(0)=rates[0];

	if (itsRates->Elt(0)>=DEFAULT_LEVEL_IN_BP) itsIsNameInDefault = true; //1e5 Bp eq to default
	int sizeRates = rates.size();
	if (isdatesininput) 
	{
		if (itsDates) delete itsDates; 
		itsDates = new ARM_Vector((ARM_Vector*)&YForDates, sizeRates+1, 0, sizeRates-1, 1);
		itsDates->Elt(0)=itsAsOfDate.GetJulian();
	}
	else 
	{
		if (itsYearTerms) delete itsYearTerms; 
		itsYearTerms = new ARM_Vector((ARM_Vector*)&YForDates,  sizeRates+1, 0, sizeRates-1, 1);
		itsYearTerms->Elt(0)=0.;
	}

	itsRecovery = Recovery;
	itsCurrency=Ccy; 
	its_STDCDS_Currency =itsCurrency; 

	SetLabel(label);
	itsCdsAdj = adj;

	itsAdjBusiness = false;
	its_STDCDS_AdjustedStartDateOnBusinessDay = adjStartDate; // K_ADJUSTED;
	its_STDCDS_Adjusted	=	intRule; // K_ADJUSTED;
	itsIsSummitCurve = issummitcurve;

	if (VolCurve)
		itsVolCurve = (ARM_VolCurve*) unconst(*VolCurve).Clone();

// 	if (isBrentCalib)
// 		SetBrentCalibration();
// 	else
// 		SetBrentDichotomy();
	
	CptSetCreditMarketData(charTermKn);
	Calibrate();
	CptTermsSurvivalProba();
}

// --------------------------------------------------------------------------------------
// compute everything for implied Cds
// --------------------------------------------------------------------------------------
double ICM_DefaultCurve::NPV_Implicit_Cds(ARM_Date& Maturity,
										const string& Tenor,	
										const qCMPMETH& mode) const 
{

	double result =0.;

	// ARM_Currency* Ccy = itsCurrency;
	// ARM_Currency dfccy = *Ccy;

	ARM_Date StartDate = its_STDCDS_StartDate; 
	ARM_Date EndDate;

	if (Tenor == "NULL") 
		EndDate = (ARM_Date) Maturity;
	else
		EndDate = AddPeriod(itsAsOfDate,Tenor,itsCurrency,itsAdjBusiness,itsCdsAdj);
		

	if (EndDate<=its_STDCDS_StartDate) return 0.;
	
	ICM_Cds cds(its_STDCDS_StartDate,
			   EndDate,
			   0,
			   0,
			   its_STDCDS_ProtectionStartDate,
			   EndDate,	
			   0.01,
			   ( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(its_STDCDS_Notional)),
			   ( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(its_STDCDS_Notional)),
			   //ARM_ReferenceValue(its_STDCDS_Notional),ARM_ReferenceValue(its_STDCDS_Notional),// 0,0,
			   // "NULL",
			   its_STDCDS_Frequency,
			   its_STDCDS_Basis,
			   its_STDCDS_AccOnDef,
			   its_STDCDS_Currency,
			   its_STDCDS_Stub,
			   its_STDCDS_CreditLag,
			   -1,
			   its_STDCDS_Adjusted,
			   its_STDCDS_IncludeMaturity,
			   its_STDCDS_AdjustedStartDateOnBusinessDay,
							   std::string(),// payCalName = NULL ,
							   qRunning_Leg,
				qStandart_Recovery_Leg,
								ISSUER_UNDEFINE,//name = ISSUER_UNDEFINE ,
								CREDIT_DEFAULT_VALUE// Binary = CREDIT_DEFAULT_VALUE 
			   );

	ICM_DefaultCurveModel DefProbModel(this,GetZeroCurve());
// 	ICM_Pricer_Cds DefPricer(&cds,&DefProbModel);
	ICM_Pricer_Cds DefPricer; DefPricer.Set(&cds,&DefProbModel,itsParameters,itsAsOfDate);
	DefPricer.SetFaster(true);

	result=DefPricer.ComputePrice(mode);
	
	return (result);
}


extern "C" unsigned short CDSCalibrationFunction(void* Param, double X, double& Result)
{
	ICM_DefaultCurve	*TheDefaultCurve;

	TheDefaultCurve	=	(ICM_DefaultCurve*)((*((vector<void*>*)Param))[0]);

	Result = TheDefaultCurve->Evaluate(X);

	return	RetOk;	// 1
}


// -------------------------------------------------------------------------
// TRY WITH A BRENT SOLVER
// -------------------------------------------------------------------------
// DefProb Curves Calibration for Intensities with modified Initial Guess
// -------------------------------------------------------------------------
void ICM_DefaultCurve::Calibrate_Stress_Test_Guess_Brent()
{
	if ((itsRecovery > 1.0) || (itsRecovery < 0.0))
		ICMTHROW(ERR_INVALID_ARGUMENT,"Recov not in [0,1] :"<<itsRecovery); 

	if (itsCalibrationData!="STD") 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Cant use BRENT with CalibrationMethod "<<itsCalibrationData); 

	double Nominal = DEFAULT_NOMINAL;
	double rate = itsRates->Elt(1);


	// Init Lambda
	double Lambda;
	Lambda	=	rate;
	if (1.0-itsRecovery)
		Lambda	/=	(1.0-itsRecovery);


	double	Sum;

	double result = 1.;

	// itsLambdas->Elt(0)=Lambda;
	SetInterpolLambda(0,Lambda); 
	ICM_DefaultCurveModel* model = NULL;

	itsInterpolYF = *itsYearTerms;	
	itsInterpolLambda=*itsLambdas; 
	itsInterpolDates=*itsDates; 
	itsInterpolRates=*itsRates; 
	itsInterpolSP=*itsSurvivalProba ;

	//
	//	here is the place where the itsInterpolSearch is constructed
	/**
	ARM_Vector* vector = new ARM_Vector(itsInterpolYF.GetSize(),0.);
	for (unsigned  i=0; i<vector->GetSize(); i++)
		vector->Elt(i)= i;
	itsInterpolSearch=ARM_ReferenceValue((ARM_Vector*)itsInterpolYF.Clone(),vector,K_STEPUP_LEFT,0);
	itsInterpolSearch.SetCalcMethod(K_STEPUP_LEFT) ;
	**/ 

	model = new ICM_DefaultCurveModel(this,GetZeroCurve(),NULL,false);

	double _inf = 0.;
	double _sup = 0.;

	int	nb_iter;

	bool	CalibSuccessFlag;
	bool	ErrorRecoveringFlag;

	Sum	=	0.0;

	int	nbdates;
	// nbdates	=	itsDates->GetSize();
	nbdates	=	itsInterpolDates.GetSize();

	double	LambdaInit;
	double	Prev_rate;
	Prev_rate	=	0.0;

	ARM_Date	ValDate;
	ARM_Date	PrevDate;
	ARM_Date	EndDate;

	ValDate		=	GetAsOfDate();
	PrevDate	=	ValDate;

	// SOLVER DECLARATION
	BrentSolver Slv;
	std::vector<void*>	TheVector;

	int	NbMaxIter;
	Slv.GetNbMaxIter(NbMaxIter);

	double	Coeff_Up;

	Coeff_Up	=	2.0;

	TheVector.push_back(this);

	// En sortie de cette boucle on a calibre tous les lambdas sauf le dernier
	for (int it=1; it<itsDates->GetSize(); it++)
	{
		CalibSuccessFlag	=	true;
		ErrorRecoveringFlag	=	false;
		
		nb_iter	=	0;
		its_Current_indice = it;
	
		// Default
		if (itsIsNameInDefault) {ResetLambda(its_Current_indice,1000.);continue;}

		// get Initial Guess
		// EndDate = GetDate(it);
		EndDate = itsInterpolDates.Elt(it) ; // GetDate(it);
		// rate = GetRate(it);
		rate = itsInterpolRates.Elt(it);

		if (it == 1)
		{
			if (1.0-itsRecovery)
				LambdaInit	=	rate / (1.0 - itsRecovery);
			else
				LambdaInit	=	1e-2;	// ???

			Lambda	=	LambdaInit;
		}
		else
		{
			Lambda	=	Sum * (rate / Prev_rate * (EndDate - ValDate) / (PrevDate - ValDate) - 1.0);
			Lambda	/=	(EndDate - PrevDate);

		}
		
		double Guess_Lambda	=	Lambda;

 		double LamdaCalibrationBrentMin = Guess_Lambda * 0.5;
		double LamdaCalibrationBrentMax = Guess_Lambda * Coeff_Up;	// may be pushed higher for steep curves, otherwise 
		double LamdaCalibrationBrentTol = 1e-5;

		ICM_Cds cds(its_STDCDS_StartDate,
							   EndDate,
							   0,
							   0,
							   its_STDCDS_ProtectionStartDate,
							   EndDate,
							   rate,
							   ( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(its_STDCDS_Notional)),
							   ( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(its_STDCDS_Notional)),
							   //ARM_ReferenceValue(its_STDCDS_Notional),
							   //ARM_ReferenceValue(its_STDCDS_Notional),// 0,0,
							   its_STDCDS_Frequency,
							   its_STDCDS_Basis,
							   its_STDCDS_AccOnDef,
							   its_STDCDS_Currency,
							   its_STDCDS_Stub,
							   its_STDCDS_CreditLag,
							   -1,
							   its_STDCDS_Adjusted,
							   its_STDCDS_IncludeMaturity,
							   its_STDCDS_AdjustedStartDateOnBusinessDay,
							   std::string(),// payCalName = NULL ,
							   qRunning_Leg,
				qStandart_Recovery_Leg,
								ISSUER_UNDEFINE,//name = ISSUER_UNDEFINE ,
								CREDIT_DEFAULT_VALUE// Binary = CREDIT_DEFAULT_VALUE 
							   );	


		// ICM_Pricer_Cds pricer(&cds,(ARM_Model*)model);
		ICM_Pricer_Cds pricer; pricer.Set(&cds,model,itsParameters,itsAsOfDate);

		pricer.SetFaster(true);
		its_pricer = &pricer;

		try
		{
			// SOLVER
			Slv.Solve(LamdaCalibrationBrentMin, LamdaCalibrationBrentMax, LamdaCalibrationBrentTol, nb_iter, CDSCalibrationFunction, &TheVector, Lambda);
			if ((nb_iter > NbMaxIter) || (Lambda <= 0.0))
			{
				// maybe one more try with a bigger range!
				LamdaCalibrationBrentMin	=	0.0;
				LamdaCalibrationBrentMax	=	100.0;
				ErrorRecoveringFlag	=	true;
				Slv.Solve(LamdaCalibrationBrentMin, LamdaCalibrationBrentMax, LamdaCalibrationBrentTol, nb_iter, CDSCalibrationFunction, &TheVector, Lambda);
				if ((nb_iter > NbMaxIter) || (Lambda <= 0.0))
					CalibSuccessFlag	=	false;
				// do not forget!
				nb_iter += NbMaxIter;
			}
		}

		catch (...)
		{
			CalibSuccessFlag	=	false;
		}

		// FINALLY
		if (!CalibSuccessFlag)
		{
			for (int lt=its_Current_indice; lt<nbdates; lt++)
				ResetLambda(lt, 1e25);

			if (model) 
				delete model;
			model = NULL;

			ICMLOG("ICM_DefaultCurve::Calibrate: failed for "
			<< itsLabel 
			<< " with init: " << Guess_Lambda << " - "
			<<" Lambda="<<Lambda<<" setting to 1e25"); 
			return;
		}
		else
		{
			Sum	+=	Lambda * (EndDate - PrevDate);

			PrevDate	=	EndDate;
			Prev_rate	=	rate;
		}

		ResetLambda(its_Current_indice,Lambda);

	}

	if (model) 
		delete model;
	model = NULL;

	*itsSurvivalProba=itsInterpolSP; 
	*itsLambdas = itsInterpolLambda; 
}

// virtual 
ICM_DefaultCurve::~ICM_DefaultCurve(void)
{
	if (itsZeroCurve)	delete itsZeroCurve;
	itsZeroCurve = NULL;
	if (itsDates)	delete itsDates;
	itsDates = NULL;
	if (itsNbDays)	delete itsNbDays;
	itsNbDays = NULL;
	if (itsYearTerms)		delete itsYearTerms;
	itsYearTerms = NULL;
	if (itsRates)		delete itsRates;
	itsRates = NULL;
	if (itsLambdas)		delete itsLambdas;
	itsLambdas = NULL;
	if (itsSurvivalProba)		delete itsSurvivalProba;
	itsSurvivalProba = NULL;
	// if (itsSearch)		delete itsSearch;
	// itsSearch = NULL;
	if (itsVolCurve)		delete itsVolCurve;
	itsVolCurve = NULL;
	if (itsDirectionalLambdasShift)		delete itsDirectionalLambdasShift;
	itsDirectionalLambdasShift	=	NULL;
	std::vector<ICM_Cds*>::iterator it = itsCalibrationInstruments.begin() ;
	while (it != itsCalibrationInstruments.end()) { delete *it ;++it; }
	itsCalibrationInstruments.clear(); 
	if (itsUpFronts) delete itsUpFronts;

	for (int i=0;i<itsAmortRefValues.size();i++)
	{if (itsAmortRefValues[i])
		delete itsAmortRefValues[i];
	 itsAmortRefValues[i]=NULL;
	}
	itsAmortRefValues.clear();

}

// -------------------------------------------------------------
// Computes Forward Spread between 2 dates
// -------------------------------------------------------------
double ICM_DefaultCurve::FwdSpread_AsIndex(const ARM_Date& t1, const ARM_Date& t2, double& FlatRbp_fwd, double& Fwd_RPV01) const 
{
	double fwdsprd = 0;;
	double Nominal = 1.;
	double rate = 0.;
	int indice = 0;
	bool equal = false;
	double Feepv1 = 0.;
	double Duration1 = 0.;
	double Feepv2 = 0.;
	double Duration2 = 0.;
	double	VirtualCDS_FwdSpread = 0.;
	double ImpliedSpread = 0.;

	ICM_DefaultCurveModel model(&unconst(*this),GetZeroCurve());
	// cas trés proche de la maturity.
	if ((t1 <= its_STDCDS_StartDate) || (t1 <= its_STDCDS_ProtectionStartDate))
	{
		fwdsprd  = ImpliedSpreadInterpol(t2);
		//curve flat en t2
		std::auto_ptr<ICM_DefaultCurve> pDefCurveFlat(createFlatCurve(t2));
		FlatRbp_fwd = pDefCurveFlat->RiskyDuration(t2);
		return fwdsprd;
	}
	// ----------------------------------------------------------------------
	// STEP 1: COMPUTE VIRTUAL CDS FWD SPREAD
	// ----------------------------------------------------------------------

	if (its_STDCDS_StartDate.GetJulian()<t1.GetJulian())
	{
		ImpliedSpread = ImpliedSpreadInterpol(t1);
		ImpliedSpread /= 10000.;
		ICM_Cds cds(its_STDCDS_StartDate,t1,
					0,0, // ARM_Date 
					its_STDCDS_ProtectionStartDate,t1,
					ImpliedSpread,// "NULL",
					( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(Nominal)),
					( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(Nominal)),
					//ARM_ReferenceValue(Nominal),
					//ARM_ReferenceValue(Nominal),// 0,0,
					its_STDCDS_Frequency,its_STDCDS_Basis,its_STDCDS_AccOnDef,
					its_STDCDS_Currency,
					its_STDCDS_Stub,its_STDCDS_CreditLag,-1,its_STDCDS_Adjusted,
					its_STDCDS_IncludeMaturity,its_STDCDS_AdjustedStartDateOnBusinessDay,
								   std::string(),// payCalName = NULL ,
								   qRunning_Leg,
					qStandart_Recovery_Leg,
									ISSUER_UNDEFINE,//name = ISSUER_UNDEFINE ,
									CREDIT_DEFAULT_VALUE// Binary = CREDIT_DEFAULT_VALUE 
					);	
		// creation DefCurveFlat
		std::auto_ptr<ICM_DefaultCurve> DefCurve(this->createFlatCurve(t1));
		model.SetDefaultCurve(DefCurve.get());
		ICM_Pricer_Cds pricer; pricer.Set(&cds,&model,itsParameters,itsAsOfDate);
		pricer.SetFaster(true);
		Feepv1 = pricer.Price(qCMPFEELEGPV);
		Duration1 = pricer.Price(qCMPDURATION);
	}

	ImpliedSpread = ImpliedSpreadInterpol(t2);
	ImpliedSpread /= 10000.;

	if (its_STDCDS_StartDate.GetJulian()<t2.GetJulian())
	{
		ICM_Cds cds(its_STDCDS_StartDate,t2,0,0,// ARM_Date 
					its_STDCDS_ProtectionStartDate,t2,
					ImpliedSpread,// "NULL",
					( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(Nominal)),
					( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(Nominal)),
					//ARM_ReferenceValue(Nominal),ARM_ReferenceValue(Nominal), // 0,0,
					its_STDCDS_Frequency,its_STDCDS_Basis,its_STDCDS_AccOnDef,
					its_STDCDS_Currency,its_STDCDS_Stub,its_STDCDS_CreditLag,
					-1,its_STDCDS_Adjusted,its_STDCDS_IncludeMaturity,its_STDCDS_AdjustedStartDateOnBusinessDay,
							   std::string(),// payCalName = NULL ,
							   qRunning_Leg,
					qStandart_Recovery_Leg,
								ISSUER_UNDEFINE,//name = ISSUER_UNDEFINE ,
								CREDIT_DEFAULT_VALUE// Binary = CREDIT_DEFAULT_VALUE 
					);

		// creation DefCurveFlat
		std::auto_ptr<ICM_DefaultCurve> DefCurve(this->createFlatCurve(t2));
		model.SetDefaultCurve(DefCurve.get());
		ICM_Pricer_Cds pricer; pricer.Set(&cds,&model,itsParameters,itsAsOfDate);
		pricer.SetFaster(true);	
		Feepv2 = pricer.Price(qCMPFEELEGPV);
		Duration2 = pricer.Price(qCMPDURATION);
	}

	Fwd_RPV01	=	Duration2 - Duration1;
	fwdsprd = (Feepv2 - Feepv1) / Fwd_RPV01;
	fwdsprd *= 10000.;
	VirtualCDS_FwdSpread	=	fwdsprd;
	ARM_Date	Mty;
	ARM_Date	ExpiryDate;

	Mty	=	t1;	ExpiryDate	=	t2;
	
	//Construction de la Def Curve du CDS Virtuel à partir des spreads des noms de l'indice
	std::auto_ptr<ICM_DefaultCurve> DefCurve((ICM_DefaultCurve*) this->Clone()) ;
	
	// Modify DefCurve
	ARM_Vector FlatRates(DefCurve->GetRates()->size(),VirtualCDS_FwdSpread/10000.); 
	// FlatRates.Fill(VirtualCDS_FwdSpread/10000.) ;
	DefCurve->setRates(FlatRates); // because Set Rate : == pointers
	

	if (!DefCurve->GetTerms().empty()) DefCurve->CptSetCreditMarketData(true);
	else DefCurve->CptSetCreditMarketData(false);
	DefCurve->Calibrate();
	DefCurve->CptTermsSurvivalProba();
		
	ICM_DefaultCurveModel DefModel_(DefCurve.get(), GetZeroCurve());

	std::auto_ptr<ICM_Cds> newcds((ICM_Cds*)DefCurve->stdCDS(Mty, ExpiryDate));

	// ICM_Pricer_Cds PricerCds(newcds,&DefModel_);
	ICM_Pricer_Cds PricerCds ;PricerCds.Set(newcds.get(),&DefModel_,itsParameters,itsAsOfDate);
	PricerCds.SetFaster(true);
	
	FlatRbp_fwd = PricerCds.Price(qCMPDURATION);
	// ----------------------------------------------------------------------
	// STEP 5: FINALLY
	// ----------------------------------------------------------------------
	// Fwd Spread de l'indice	
	fwdsprd = VirtualCDS_FwdSpread + 10000.* Feepv1/FlatRbp_fwd ;
	// ----------------------------------------------------------------------
	return (fwdsprd);
}

// -------------------------------------------------------------
// Computes Forward Spread between 2 yearterms
// -------------------------------------------------------------
double ICM_DefaultCurve::FwdSpread_AsIndex(const double& yt1, const double& yt2, double& FlatRbp_fwd, double& Fwd_RPV01) const 
{
	double t1 = itsAsOfDate.GetJulian() + yt1*365.;
	double t2 = itsAsOfDate.GetJulian() + yt2*365.;

	ARM_Date	Date_t1(t1);
	ARM_Date	Date_t2(t2);

	return FwdSpread_AsIndex(Date_t1, Date_t2, FlatRbp_fwd, Fwd_RPV01);
}

// ----------------------------------------------------------------------------------
// Returns RiskyPV01 as a sensitivity
// ----------------------------------------------------------------------------------
double ICM_DefaultCurve::RiskyPV01AsSensitivity(const std::string & Tenor) const 
{
	double Nominal = 1.;
	
	bool status = true;
	int size = itsRates->GetSize();
	int i = 0;

	for (i=0;i<size;i++) // if (strcmp(itsTerms[i].c_str(),Tenor)==NULL) {status = false;break;}
		if (itsTerms[i]==Tenor) {status = false;break;}

	// if (status) return CREDIT_DEFAULT_VALUE;
	if (status) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_DefaultCurve::RiskyPV01AsSensitivity: "
			<< Tenor<< "not found"); 
		

	i++;

	ARM_Date Maturity = (ARM_Date) GetDate(i);

	ICM_DefaultCurveModel model(this,GetZeroCurve());
	ICM_Cds cds(its_STDCDS_StartDate,Maturity,0,0,its_STDCDS_ProtectionStartDate,
				  Maturity,1.e-4,// "NULL",
				  ( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(Nominal)),
				  ( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(Nominal)),
				  //ARM_ReferenceValue(Nominal),ARM_ReferenceValue(Nominal),// 0,0,
				  its_STDCDS_Frequency,
				  its_STDCDS_Basis,its_STDCDS_AccOnDef,
				  its_STDCDS_Currency,its_STDCDS_Stub,
				  its_STDCDS_CreditLag,-1,its_STDCDS_Adjusted,its_STDCDS_IncludeMaturity,
				  its_STDCDS_AdjustedStartDateOnBusinessDay,
							   std::string(),// payCalName = NULL ,
							   qRunning_Leg,
				qStandart_Recovery_Leg,
								ISSUER_UNDEFINE,//name = ISSUER_UNDEFINE ,
								CREDIT_DEFAULT_VALUE// Binary = CREDIT_DEFAULT_VALUE 
				  );	

	// ICM_Pricer_Cds pricer(&cds,&model);
	ICM_Pricer_Cds pricer; pricer.Set(&cds,&model,itsParameters,itsAsOfDate);
	pricer.SetFaster(true);

	// double PV01 = pricer.ComputeSensitivity(ICMSPREAD_TYPE,Tenor)/1.e-4;
	double PV01 = pricer.Hedge(ICMSPREAD_TYPE,Tenor,"",CREDIT_DEFAULT_VALUE)/1.e-4;

	return (PV01);
}

void ICM_DefaultCurve::AdjustDateYearsTerms(bool dateVectorKn){
	int size = itsRates->GetSize();

	if (dateVectorKn){ // case of cpt YFTerms and NbDays
		for (int i=1; i< size; i++){
			int Nbday		= (int) (itsDates->Elt(i) - itsAsOfDate.GetJulian()); 
			double YearTerm = (itsDates->Elt(i) - itsAsOfDate.GetJulian())/K_YEAR_LEN;
			itsNbDays->Elt(i)=(double)Nbday;
			itsYearTerms->Elt(i)=YearTerm;
		}
	} else { // case of  cpt dateTerms and Nbdays
		for (int i=1; i< size; i++){
			ARM_Date  Date = itsAsOfDate.GetJulian() + K_YEAR_LEN*itsYearTerms->Elt(i);
			int Nbday		= (int) (Date - itsAsOfDate); 
			itsNbDays->Elt(i)=(double)Nbday;
			itsDates->Elt(i)= Date.GetJulian();
		}
	}
}

double ICM_DefaultCurve::ImpliedSpreadInterpolOLD(const std::string& plot, 
											   const double& slope, 
											   const double& ExactDate, // only when the date is a date of the Defaultcurve else use InterDate.
											   const double& InterDate)
{
	double result =0.;
	double Spread = 0.;
	int indice = 0;
	bool equal = false;
	double yt = 0.;

	if (ExactDate != 0)
	{
		yt = (((ARM_Date)ExactDate).GetJulian() - itsAsOfDate.GetJulian())/K_YEAR_LEN;
		/** Inf_Equal(yt,indice,equal); **/ 
		indice=position(yt); 
		//if( equal) {
		return 10000.*GetRate(indice);
		//}
	}


	ARM_Date StartDate = its_STDCDS_StartDate; 
	ARM_Date EndDate;

	if (InterDate != 0) 
		EndDate = (ARM_Date) InterDate;
	else
		EndDate = AddPeriod(itsAsOfDate,plot,itsCurrency,itsAdjBusiness,itsCdsAdj);
		

	yt = (EndDate.GetJulian() - itsAsOfDate.GetJulian())/K_YEAR_LEN;

	

	if ((itsDates->GetSize() == 2) && (slope))
	{
		double yt = (EndDate.GetJulian() - itsAsOfDate.GetJulian())/K_YEAR_LEN;
		double alpha = yt - itsYearTerms->Elt(1);
		result = 10000.*itsRates->Elt(0)*exp(alpha*log(1.+slope));
		return (result);
	}

	if (EndDate<=its_STDCDS_StartDate) return 0.;
	
	ICM_Cds cds(its_STDCDS_StartDate,
			    EndDate,
				0,0, // ARM_Date
			    its_STDCDS_ProtectionStartDate,
			    EndDate,	
			    0.01,
				( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(its_STDCDS_Notional)),
				( (itsAmortRefValues.size()>0) ? *itsAmortRefValues[0] : ARM_ReferenceValue(its_STDCDS_Notional)),
				//ARM_ReferenceValue(its_STDCDS_Notional),ARM_ReferenceValue(its_STDCDS_Notional),// 0,0,
			    its_STDCDS_Frequency,
			    its_STDCDS_Basis,
			    its_STDCDS_AccOnDef,
			    its_STDCDS_Currency,
			    its_STDCDS_Stub,
			    its_STDCDS_CreditLag,
			    -1,
			    its_STDCDS_Adjusted,
			    its_STDCDS_IncludeMaturity,
			    its_STDCDS_AdjustedStartDateOnBusinessDay,
							   std::string(),// payCalName = NULL ,
							   qRunning_Leg,
				qStandart_Recovery_Leg,
								ISSUER_UNDEFINE,//name = ISSUER_UNDEFINE ,
								CREDIT_DEFAULT_VALUE// Binary = CREDIT_DEFAULT_VALUE 
				);

	ICM_DefaultCurveModel DefProbModel(this,GetZeroCurve());
	// ICM_Pricer_Cds DefPricer(&cds,&DefProbModel);
	ICM_Pricer_Cds DefPricer; DefPricer.Set(&cds,&DefProbModel,itsParameters,itsAsOfDate);
	DefPricer.SetFaster(true);

	result=DefPricer.ComputeSpread(0.);

	return (result);
}
void ICM_DefaultCurve::SetVolCurve(const ARM_VolCurve* value)
		{
			if (itsVolCurve)
				delete itsVolCurve;
			itsVolCurve = NULL;

			itsVolCurve = (ARM_VolCurve*) unconst(*value).Clone();
		}
void ICM_DefaultCurve::SetZeroCurve(ARM_ZeroCurve* zc) 
		{
			if (itsZeroCurve) 
				delete itsZeroCurve;
			itsZeroCurve = zc; 
		};

//	------------------------------------------------------------------------------------------------
ICM_DefaultCurve* 
ICM_DefaultCurve::createFlatCurve(const std::string& tenor) const
{
	ARM_Date date; 
	bool doAdjust(itsAdjBusiness); 
	if (its_STDCDS_Adjusted==K_MATUNADJUSTED || its_STDCDS_Adjusted==K_UNADJUSTED) doAdjust=false;
	else doAdjust=true; 
	if (itsCdsAdj == qCredit_Default) date = AddPeriod(its_STDCDS_StartDate, tenor,itsCurrency,doAdjust,itsCdsAdj);
	else date = AddPeriod(itsAsOfDate, tenor,itsCurrency,doAdjust,itsCdsAdj);
	return createFlatCurve(date); 
}
//	------------------------------------------------------------------------------------------------
ICM_DefaultCurve* 
ICM_DefaultCurve::createFlatCurve(const ARM_Date& date) const
{
	double spread = this->ImpliedSpreadInterpol(date); 
	ICM_DefaultCurve* newDefCurve = dyn_clone(this);
	
	// rates :
	if (newDefCurve->itsRates) delete newDefCurve->itsRates;
	newDefCurve->itsRates = new ARM_Vector(2,0.);
	newDefCurve->itsRates->Elt(0) = newDefCurve->itsRates->Elt(1) = spread/10000;
	// dates
	if (newDefCurve->itsDates) delete newDefCurve->itsDates; 
	newDefCurve->itsDates = new ARM_Vector(2,0.);
	newDefCurve->itsDates->Elt(0)=itsAsOfDate.GetJulian();
	newDefCurve->itsDates->Elt(1)=date.GetJulian();
	if (newDefCurve->itsYearTerms) 
		delete newDefCurve->itsYearTerms;
	newDefCurve->itsYearTerms = NULL;
	newDefCurve->CptSetCreditMarketData(false);
	newDefCurve->Calibrate();
#ifdef _DEBUG
	FILE *pFile =NULL;
	pFile = fopen("C:\\temp\\flatCurve.txt","w");
	newDefCurve->View("",pFile);
	if(pFile) fclose(pFile);
#endif
	return newDefCurve; 
}
//	------------------------------------------------------------------------------------------------
void ICM_DefaultCurve::SetRates(ARM_Vector* item) 
{
	if (!item) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_DefaultCurve::SetRates: null argument"); 
	if (itsRates)
		delete itsRates;
	itsRates = item;
}
