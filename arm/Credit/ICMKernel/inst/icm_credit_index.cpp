#include "firsttoinc.h"
#include "ICMKernel\inst\icm_credit_index.h"
#include "ICMKernel\mod\modelmulticurves.h"
#include "ICMKernel\inst\icm_cds.h"
#include "ICMKernel\pricer\icm_pricer_cds.h"
#include "ICMKernel\crv\ICM_Constant_Piecewise.h"
#include <algorithm>



void ICM_Credit_Index::SetMaturities(const ARM_Vector& Maturities){
	if (itsRunning.size() != Maturities.size()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Maturities Size is not of the same size of spread "); 
	itsMaturities = Maturities;
	if 	(itsMaturities.size() >0) {
		ARM_Date matu;
		try {
			matu = ARM_Date(itsMaturities.Elt(0));
			itsIsYT = true;
		} catch(...){
			itsIsYT = false;
			// nothing. Maturity in ARM_IRIndex is not set.
		}
	}
}
void ICM_Credit_Index::SetRunning(const ARM_Vector& Spread) {
	if (itsMaturities.size() != Spread.size()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Spread Size is not of the same size of Maturities "); 
	itsRunning = Spread;
}
//	-------------------------------------------------------------------------------------
void 
ICM_Credit_Index::Init()
{
	SetName(ICM_CREDIT_INDEX);
	itsLabels.resize(0);
	itsCreditIndexName	= "";
	itsMethod			= qAVERAGE;
	itsDefaultCurve		= NULL;
	itsForcedCurve		= NULL;
//	itsVol				= NULL;
	itsAdjForTenor		= qCredit_Default;
	itsDefaultCurveflg	= false;
	itsCM_resetWeekDay	=5; 
	itsCM_resetOccur	=2;
	itsIsHomogeneous	= false;
	itsRunning.clear();
	itsMaturities.clear();
	itsIsYT = false;
}

//	-------------------------------------------------------------------------------------
// virtual 
ICM_Credit_Index::~ICM_Credit_Index()
{
	if (itsDefaultCurve)
		delete itsDefaultCurve;
	itsDefaultCurve = NULL;

	if (itsForcedCurve)
		delete itsForcedCurve;
	itsForcedCurve = NULL;

/*	if (itsVol)
		delete itsVol;
	itsVol = NULL;
*/
}	

// virtual 
void ICM_Credit_Index::Copy(const ARM_Object* src)
{
  ARM_IRIndex::Copy(src);

  BitwiseCopy(src);
}

// virtual 
ARM_Object* ICM_Credit_Index::Clone(void)
{
  ICM_Credit_Index* theClone = new ICM_Credit_Index();

  theClone->Copy(this);

  return(theClone);
}
//	-------------------------------------------------------------------------------------
void 
ICM_Credit_Index::BitwiseCopy(const ARM_Object* src)
{
	int i = 0;
	ICM_Credit_Index* ind = (ICM_Credit_Index*) src;
	itsMethod = ind->itsMethod;
	itsLabels=ind->itsLabels; 
	itsIsHomogeneous = ind->itsIsHomogeneous;
	
	if (ind->itsDefaultCurve)
		// itsDefaultCurve = (ICM_DefaultCurve*) ind->itsDefaultCurve->Clone();
		itsDefaultCurve = dyn_clone( ind->itsDefaultCurve ); 

	if (ind->itsForcedCurve)
		itsForcedCurve = dyn_clone(ind->itsForcedCurve) ;

/*	if (ind->itsVol)
		itsVol = (ARM_VolCurve*) ind->itsVol->Clone();*/
	itsIsYT = ind->itsIsYT;
	itsMaturities = ind->itsMaturities;
	itsRunning = ind->itsRunning;
	itsAdjForTenor = ind->itsAdjForTenor;
	itsCM_resetWeekDay=ind->itsCM_resetWeekDay; 
	itsCM_resetOccur=ind->itsCM_resetOccur; 
	itsCreditIndexName = ind->itsCreditIndexName;
}
// --------------------------------------------------------------------
// Méthode View
// --------------------------------------------------------------------
void ICM_Credit_Index::View(char* id, FILE* ficOut)
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

   fprintf(fOut, "\n\n--------------------------------------------------------------------------------\n");
   fprintf(fOut, "--------------------------------- Credit Index ---------------------------------\n");
   fprintf(fOut, "--------------------------------------------------------------------------------\n\n");

   fprintf(fOut,"Index Name :\t %s\n\n",itsCreditIndexName.c_str());
   fprintf(fOut, "Issuer Labels\n");
   fprintf(fOut, "-------------\n");
	for (int i=0;i<itsLabels.size();i++)
		fprintf(fOut,"%s\n",itsLabels[i].c_str());
   

   fprintf(fOut,"\n");
   fprintf(fOut,"Coupons : ");
   if(itsIsYT == false){
	   fprintf(fOut,"In Date : \n");
	   for ( int i=0; i<itsMaturities.size(); i++){

			ARM_Date lDate((ARM_Date)itsMaturities.Elt(i));
		    string mat(lDate.GetStrDate());
			fprintf(fOut,"date %s ->\t %lf \n", mat.c_str(), itsRunning.Elt(i));
	   }
   } else {
		fprintf(fOut,"In YT: \n");
		for ( int i=0; i<itsMaturities.size(); i++){
			fprintf(fOut,"YT %4.2lf ->\t %lf \n", itsMaturities.Elt(i), itsRunning.Elt(i));
	   }
   }

   ARM_IRIndex::View(id, fOut);
   fprintf(fOut,"CM_resetWeekDay:%d\n",itsCM_resetWeekDay); 
   fprintf(fOut,"CM_resetOccur:%d\n",itsCM_resetOccur); 

//   if (itsVol) itsVol->View(id,fOut);
   if (itsDefaultCurve) itsDefaultCurve->View(id, fOut);
   
   if ( ficOut == NULL )
   {
      fclose(fOut);
   }

}


ICM_Credit_Index::ICM_Credit_Index(const int& dayCount, 
								   const int& resetFreq, 
								   const int& payFreq, 
								   const ARM_Vector& maturity,
								   const string &ccy /*ARM_Currency* ccyName*/,
								   const std::string& IndexName ,
								   const std::vector<std::string> &labels,
								   const qINDEX_CMPT_METHOD& Method,
								   const ARM_Vector& Spread,								   
								   const ICM_DefaultCurve* ForcedCurve,
								   //const ARM_VolCurve* volflat,
								   const int& fwdRule,
								   const int& intRule,	
								   const int& resetTiming,	
								   const int& resetGap,	
								   const int& payTiming,	
								   const int& payGap,
								   const qCDS_ADJ& adj,
								   int cm_resetWeekDay,
								   int cm_resetOccur)
{
	Init();


	Set(dayCount,resetFreq,payFreq,maturity,ccy,IndexName, labels,
		Method,Spread,// maturity,
		ForcedCurve,//volflat,
		fwdRule,intRule,resetTiming,resetGap,payTiming,payGap,adj,cm_resetWeekDay,cm_resetOccur);
}


void ICM_Credit_Index::Set(const int& dayCount, 
						   const int& resetFreq, 
						   const int& payFreq, 
						   const ARM_Vector& maturity,
						   const std::string &ccy /*ARM_Currency* ccyName*/,
						   const std::string& IndexName ,
						   const std::vector<std::string>& labels,
						   const qINDEX_CMPT_METHOD& Method,
						   const ARM_Vector& Spread,
						   const ICM_DefaultCurve* ForcedCurve,
						   //const ARM_VolCurve* volflat,
						   const int& fwdRule,
						   const int& intRule,	
						   const int& resetTiming,	
						   const int& resetGap,	
						   const int& payTiming,	
						   const int& payGap,
						   const qCDS_ADJ& adj,
						   int cm_resetWeekDay,
						   int cm_resetOccur)
{
	if ( (cm_resetWeekDay<0) || (cm_resetWeekDay>6) )
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Credit_Index::Set: cm_resetWeekDay should be [0..6] current="<<cm_resetWeekDay) ; 
	itsCM_resetWeekDay=cm_resetWeekDay;
	itsCM_resetOccur=cm_resetOccur;
	double Yt =0.0;

	int i = 0;
	auto_ptr<ARM_Currency> pCcy(new ARM_Currency(ccy.c_str()));
	

	itsMethod = Method ;
	itsRunning  = Spread ;
	itsMaturities = maturity;

	SetTerm(K_ANNUAL);
	
	ARM_IRIndex::Set(dayCount,
					 resetFreq,
					 payFreq,
					 maturity[0],
					 0 /* compMeth */,
					 fwdRule,
					 intRule,
					 resetTiming,
					 resetGap,
					 payTiming,
					 payGap,
					 pCcy.get());

	if (maturity.size()>0 ) { 
		// verrifiaction date 
		ARM_Date matu;
		try {
			matu = ARM_Date(maturity[0]);
			SetMaturity(ARM_Date(maturity[0]));
			itsIsYT = false;
			Yt = -1;
		} catch(...){
			itsIsYT = true;
			Yt = maturity[0];
			// matu
		}
	}
	itsCreditIndexName = IndexName;
	itsLabels.resize(0);	
	if(labels.size()!=0) itsLabels = labels; 
	else itsLabels.push_back(itsCreditIndexName);

	if (ForcedCurve)
		itsForcedCurve = dyn_clone(ForcedCurve) ;

/*	if (volflat)
		itsVol = dyn_clone( volflat );
*/	

	SetAdjForTenor(adj);
	CptFullHomog();
}					 

// --------------------------------------------------------------------
// Calcul du spread Fwd 
// --------------------------------------------------------------------
/*double ICM_Credit_Index::FwdSpread(const double& yt1,const double& yt2) const
{
	if (itsDefaultCurve == NULL) return CREDIT_DEFAULT_VALUE;
	
	double yt2_m =yt2;
	double FwdSpread = 0.;
	
	if (yt2_m<0) yt2_m =  yt1 + GetYearTerm();

	FwdSpread = itsDefaultCurve->FwdSpread(yt1,yt2_m);

	return (FwdSpread);
}

// --------------------------------------------------------------------
// Calcul de l ajustement de convexité
// --------------------------------------------------------------------
double ICM_Credit_Index::AjustConvexity(const double& yt1,const double& yt2,const double& SpreadFwd) const 
{
	if (itsDefaultCurve == NULL) return CREDIT_DEFAULT_VALUE;
	
	double yt2_m =yt2;
	double SpreadFwd_m = SpreadFwd;
	double AjustConvex = 0.;
	ARM_VolCurve* vol = NULL;
	if (itsVol)
		vol = itsVol;
	
	if (yt2_m<0) yt2_m =  yt1 + GetYearTerm();
	if (SpreadFwd_m <0) SpreadFwd_m = FwdSpread(yt1,yt2_m);

	AjustConvex = itsDefaultCurve->AjustConvexity(yt1,yt2_m,SpreadFwd_m,vol);

	return (AjustConvex);
}

// --------------------------------------------------------------------
// Payment Lag for CMCDS
// --------------------------------------------------------------------

double ICM_Credit_Index::PaymentLagAjust (const double& yt1, const double& yt2 , const double& SpreadFwd) const 
{
	if (itsDefaultCurve == NULL) return CREDIT_DEFAULT_VALUE;
	
	double yt2_m =yt2;
	double SpreadFwd_m = SpreadFwd;
	double PayLagAjust = 0.;
	ARM_VolCurve* vol = NULL;
	if (itsVol)
		vol = itsVol;
	
	if (yt2_m<0) yt2_m =  yt1 + GetYearTerm();
	if (SpreadFwd_m <0) SpreadFwd_m = FwdSpread(yt1,yt2_m);

	PayLagAjust = itsDefaultCurve->PayLagAdjustment(yt1,yt2_m,SpreadFwd_m,vol);

	return (PayLagAjust);

}

*/
// --------------------------------------------------------------------
// Calcul de la Default Curve Method 1
// --------------------------------------------------------------------
const ICM_DefaultCurve* ICM_Credit_Index::CptDefaultCurve(ICM_ModelMultiCurves* model)
{
	//if (itsDefaultCurve) return itsDefaultCurve;

	if (itsForcedCurve)
		return (ICM_DefaultCurve*)GetForcedCurve()->Clone();
		///adoptDefCurve((ICM_DefaultCurve*)GetForcedCurve()->Clone());
	ICM_DefaultCurve* pDefCurve= NULL;
	//else
	{
		double IndexSpread = GetSpread() ;
		ARM_Vector VectorYF;
		int size = 0,  i=0,j=0;
		model->GetUnionYFForCollateral(itsLabels,VectorYF);
		size = VectorYF.size();
		//ARM_Vector* spreads = new ARM_Vector(size,0.);
		ARM_Vector Rates(size);

		double spread = 0.;
		double recovery = 0.;

		qCDS_ADJ adj = model->GetDefaultCurve(itsLabels[0])->GetCdsAdj();  	
		std::string ccy = model->GetDefaultCurve(itsLabels[0])->GetCurrency();  	
		bool issummitcurve = model->GetDefaultCurve(itsLabels[0])->GetIsSummitCurve();
		std::string calibrationData =model->GetDefaultCurve(itsLabels[0])->GetCalibrationData();
		qDEFCURVE_CALIB_ALGO calibrationAlgo = model->GetDefaultCurve(itsLabels[0])->GetCalibrationAlgo();
		int lag = model->GetDefaultCurve(itsLabels[0])->GetLagStartDate();

		double Aux = 0.;
			
		switch (itsMethod)
		{
			case qAVERAGE :
			{
				for (i=0;i<size;i++)
				{
					spread = 0.;
					recovery = 0.;
					for (j=0; j<itsLabels.size(); j++)
					{
						double interpDate = model->GetDefaultCurve(itsLabels[j])->GetAsOfDate().GetJulian() + 365.*VectorYF[i] ;
						spread += model->GetDefaultCurve(itsLabels[j])->ImpliedSpreadInterpol(ARM_Date(interpDate) ); 
						recovery += model->GetDefaultCurve(itsLabels[j])->GetRecovery();
					}

					spread/=(double)itsLabels.size()*10000.;
					recovery/=(double)itsLabels.size();
					// Rates[i] = spreads->Elt(i)=spread;			
					Rates[i] =spread;			
				}
			}
			break;
			default:
			;
		}

		pDefCurve =  new ICM_Constant_Piecewise(model->GetStartDate(),
														VectorYF,
														Rates,
														recovery,
														model->GetZeroCurve(),
														K_ADJUSTED,	// intRule
														K_ADJUSTED,	// adjStartDate
														adj, // adjqCredit_Special_None_YF, // ,,
														ccy,
														"INDEX",
														issummitcurve,
														NULL,
														false,
														K_QUARTERLY,
														calibrationAlgo,
														calibrationData , 
														lag);
		// if (spreads)
		// 	delete spreads;
		// spreads = NULL;
	}
	return pDefCurve;
}


// -----------------------------------------------------------------------------------
// Calcul de la Default Curve Method 2 utilisé dans le pricing des index options
// -----------------------------------------------------------------------------------

const ICM_DefaultCurve* ICM_Credit_Index::CptDefaultCurveIndex(ICM_ModelMultiCurves* model, ICM_Cds* cds)
{
	if (GetDefaultCurveflg()) return itsDefaultCurve;

	double IndexSpread = GetSpread() ;
	int size = 0,i=0,j=0;
	ARM_Vector VectorYF;
	model->GetUnionYFForCollateral(itsLabels,VectorYF);
	size = VectorYF.size(); 
	ARM_Vector spreads(size); 
	double spread = 0.;
	double recovery = 0.;

	qCDS_ADJ adj = model->GetDefaultCurve(itsLabels[0])->GetCdsAdj();  	
	string ccy = model->GetDefaultCurve(itsLabels[0])->GetCurrency();  	
	bool issummitcurve = model->GetDefaultCurve(itsLabels[0])->GetIsSummitCurve();
	std::string calibrationData = model->GetDefaultCurve(itsLabels[0])->GetCalibrationData();
	qDEFCURVE_CALIB_ALGO calibrationAlgo = model->GetDefaultCurve(itsLabels[0])->GetCalibrationAlgo();
	int lag = model->GetDefaultCurve(itsLabels[0])->GetLagStartDate();

	double Aux = 0.;

	// string name = string("C:\\temp\\resDefaultCurveIndex_");
	// name += model->GetStartDate().GetStrDate() ;
	// name += string(".txt");	
	switch (itsMethod)
	{
		case qAVERAGE :
		{
			for (i=0;i<size;i++)
			{
				spread = 0.;
				recovery = 0.;
				for (j=0; j<itsLabels.size(); j++)
				{
					double interpDate = model->GetDefaultCurve(itsLabels[j])->GetAsOfDate().GetJulian() + 365.*VectorYF[i] ;
					spread += model->GetDefaultCurve(itsLabels[j])->ImpliedSpreadInterpol(interpDate ); 
					recovery += model->GetDefaultCurve(itsLabels[j])->GetRecovery();
				}

				spread/=(double)itsLabels.size()*10000.;
				recovery/=(double)itsLabels.size();
				spreads[i]=spread;			
			}
			
			itsDefaultCurve = new ICM_Constant_Piecewise(model->GetStartDate(),
														 VectorYF, 
														 spreads,
														 recovery,
														 model->GetZeroCurve(),
														K_ADJUSTED,	// intRule
														K_ADJUSTED,	// adjStartDate
														 adj, 
														 ccy,
														 "INDEX",
														 issummitcurve,
														 NULL,false,
														 K_QUARTERLY,calibrationAlgo,
														 calibrationData,
														 lag
														 );	
				
		}
		break;
		case qHOMOTHETIE :
		{
			for (i=0;i<size;i++)
			{
			
				spread = 0.;
				recovery = 0.;
				for (j=0; j<itsLabels.size(); j++)
				{
					double interpDate = model->GetDefaultCurve(itsLabels[j])->GetAsOfDate().GetJulian() + 365.*VectorYF[i] ;
					spread += model->GetDefaultCurve(itsLabels[j])->ImpliedSpreadInterpol(interpDate ); 
					recovery += model->GetDefaultCurve(itsLabels[j])->GetRecovery();
				}


				spread/=(double)itsLabels.size()*10000.;
				recovery/=(double)itsLabels.size();
				spreads[i]=spread;
			}
		

			itsDefaultCurve = new ICM_Constant_Piecewise(model->GetStartDate(),
														  VectorYF, 
														 spreads,
														 recovery,
														 model->GetZeroCurve(),
														K_ADJUSTED,	// intRule
														K_ADJUSTED,	// adjStartDate
														 adj, 
														 ccy,
														 "INDEX",
														 issummitcurve,
														 NULL,false,//2 NULL,
														 K_QUARTERLY,
														 calibrationAlgo,
														 calibrationData, lag);

	
			//Obtention du spread du CDS Virtuel

			ARM_Date EndDate ;

			if (GetExpiryDate().GetJulian() >2450000.)
				EndDate =GetExpiryDate();
			else
			{
				char TEMP[4];
				int nb = floor(GetYearTerm());
				sprintf(TEMP,"%iY",nb);
				EndDate = AddPeriod(model->GetStartDate(),TEMP,cds->GetFeeLeg()->GetCurrencyUnit()->GetCcyName(),false,adj /*qCredit_Adjust20*/);
			}

			double NewSpread = itsDefaultCurve->ImpliedSpreadInterpol(EndDate) ;
			
			Aux = IndexSpread/NewSpread ;
			
			if(itsDefaultCurve)
				delete itsDefaultCurve ;
			itsDefaultCurve = NULL ;

			for (i=0;i<size;i++)
			{
				spreads[i]=spreads[i]*Aux;
			}

			itsDefaultCurve = new ICM_Constant_Piecewise(model->GetStartDate(),
														 VectorYF,  
														 spreads,
														 recovery,
														 model->GetZeroCurve(),
														K_ADJUSTED,	// intRule
														K_ADJUSTED,	// adjStartDate
														 adj, 
														 ccy,
														 "INDEX",
														 issummitcurve,
														 NULL,false,//2 NULL,
														 K_QUARTERLY,calibrationAlgo,
														 calibrationData, lag);
			
		}
		break ;
		default:
		;
	}
	return itsDefaultCurve;
}

void ICM_Credit_Index::CptFullHomog(void) 
{
	itsIsHomogeneous=false;
	for (int i=1;i<itsLabels.size();i++)
	{ 
		if (itsLabels[i-1]!=itsLabels[i] ) return; 
	}
	itsIsHomogeneous=true;
}
void ICM_Credit_Index::adoptDefCurve(ICM_DefaultCurve* value)
{
	if (itsDefaultCurve)
		delete itsDefaultCurve;
	itsDefaultCurve = value;
	itsDefaultCurveflg = true;
}
void ICM_Credit_Index::SetForcedCurve(const ICM_DefaultCurve* value)
{
	if (!value) ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Credit_Index::SetForcedCurve: null argument"); 
	if (itsForcedCurve)
		delete itsForcedCurve;
	itsForcedCurve=NULL; 
	itsForcedCurve = dyn_clone(value) ;
}
/*void ICM_Credit_Index::SetVolatility(const ARM_VolCurve* value)
{
	if (!value) ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Credit_Index::SetVolatility: null argument"); 
	if (itsVol)delete itsVol;
	itsVol=NULL; 
	itsVol = dyn_clone(value); 
}*/
void ICM_Credit_Index::ResetDefCurve()
{
	if (GetName()!=ICM_CREDIT_INDEX) return;

	if (itsDefaultCurve)
		delete itsDefaultCurve;
	itsDefaultCurve = NULL;
	itsDefaultCurveflg = false;
}

