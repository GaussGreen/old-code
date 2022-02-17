#include "firsttoinc.h"
#include "ICMKernel/pricer/icm_pricer_index.h"
#include "ICMKernel/inst/icm_credit_index.h"
#include "ICMKernel/mod/modelmulticurves.h"


ICM_Pricer_CDSIndex::ICM_Pricer_CDSIndex(void) 
{ 
	SetName(ICM_PRICER_CDS_INDEX);
	itsDefCurveIndex=NULL; 
	Init();
}
ICM_Pricer_CDSIndex::~ICM_Pricer_CDSIndex() 
{
	if (itsDefCurveIndex) delete itsDefCurveIndex; 
	itsDefCurveIndex=0; 
}


//*******************************************************************************************
// Init Method
//*******************************************************************************************

void ICM_Pricer_CDSIndex::Init()
{
	ICM_Pricer_Cds::Init();
	itsDefCurveIndex = NULL;
	itsFlatCalculation =  true;	
}

//*******************************************************************************************
// Reset Method
//*******************************************************************************************

void ICM_Pricer_CDSIndex::Reset(void)
{
	ICM_Pricer_Cds::Reset();
	if (itsDefCurveIndex) delete itsDefCurveIndex; itsDefCurveIndex=0; 
	itsFlatCalculation = true;
}
//*******************************************************************************************
// Set Method
//*******************************************************************************************

void ICM_Pricer_CDSIndex::Set(ARM_Security *option, ARM_Model *mod,const ICM_Parameters& params,const ARM_Date&asof)
{
	if ((mod->GetName()!=ICM_MODELMULTICURVES))			
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : pricer not compatible with this security or this model ");
	ICM_Pricer_Cds::Set(option, mod,params,asof);
	if (! params.empty())
	{
		int res = params.getLong("FLAT_COMPUTE");
		if (res == 0)
			itsFlatCalculation= false;
		else
			itsFlatCalculation = true;
	}

}

//*******************************************************************************************
// Compute UpFront Payment for CDS Index Contract
//*******************************************************************************************

void ICM_Pricer_CDSIndex::ComputeUpfrontPay(qCMPMETH measure )
{
	//ICMQUANTIFYER("ICM_Pricer_CDSIndex::ComputeUpfrontPay"); 
	//if (itsFlatCalculation) 
	CptPricingMeFlat(measure);
	/*else 
	{
	ICM_Cds* cds = (ICM_Cds*) GetSecurity() ;
	ICM_Leg* Feeleg = cds->GetFeeLeg();
	const ICM_Credit_Index* index = Feeleg->GetCreditIndex() ; // (ICM_Credit_Index*) Feeleg->GetCreditIndex()->Clone() ;
	// 14514 double Notional = cds->GetInitialNotional() ;
	double Notional = cds->GetFeeLeg()->GetCreditInfosRef().GetNotionals().Elt(0) ;
	double RefSpread = Feeleg->GetRefSpread() ;
	int NbCurves = index->GetNbCurves() ;
	qINDEX_CMPT_METHOD  method = index->GetMethod() ;

	double FixedLeg = 0. ;
	double FloatingLeg = 0. ;
	double UpFrontPay = 0. ;
	double Spread = 0. ;

	switch(method)
	{
	case qAVERAGE:
		{		
			ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();

			ARM_ZeroCurve* InitShortCurve = NULL;
			InitShortCurve = (ARM_ZeroCurve*) model->GetZeroCurve()->Clone();
			ICM_DefaultCurveModel* DefModel = NULL;
			
			const ICM_DefaultCurve* DefCurve = NULL ;
			ARM_Date StartDate = model->GetStartDate().AddDays(1);
			ARM_Date EndDate = cds->GetEndDateNA();

			ICM_Pricer_Cds* PricerCds = NULL ;
					
			double Rbp = 0. ;
			
			for (int i = 0; i<NbCurves ; i ++ )
			{
				DefCurve = model->GetDefaultCurve(index->GetLabels()[i]) ;			
				Rbp = DefCurve->RiskyDuration(EndDate.GetJulian());
				Spread = DefCurve->ImpliedSpreadInterpol(EndDate );

				FloatingLeg += (Spread/10000) * Rbp ;
				FixedLeg += Rbp ;
				
			}

			Spread = 10000*FloatingLeg/FixedLeg ;
			FloatingLeg *= Notional;
			FixedLeg *= RefSpread * Notional ;
			UpFrontPay = FixedLeg- FloatingLeg ;
			SetSpread(Spread) ;

			if (InitShortCurve)
				delete InitShortCurve;
			InitShortCurve = NULL;

		} break; 
	case qHOMOTHETIE: 
		{
			double FlatRbp = ComputeDuration() ; // FlatRbp_Spot
			Spread = GetSpread() ;
			
			FloatingLeg = (Spread/10000.) * FlatRbp * NbCurves * Notional;
			FixedLeg = RefSpread * FlatRbp * NbCurves * Notional ;
			UpFrontPay = FixedLeg - FloatingLeg;
		} break ;
	default:
		ICMTHROW(ERR_INVALID_ARGUMENT,"Unhandled Index Method "<<method); 
	}
	SetFeeLegPrice(FixedLeg);
	SetDefLegPrice(FloatingLeg) ;
	SetPrice(UpFrontPay) ;
	}*/
}

//*******************************************************************************************
// Compute Spread
//*******************************************************************************************

double ICM_Pricer_CDSIndex::ComputeSpread(const double& MtM)
{
	if (GetSpreadFlg() && MtM == 0)return  GetSpread() ;// return BE
//	------------------------------------- idem CptPricingMeFlat ----------
	ICM_Pricer_Cds PricerCds;
	ICM_Cds* cds = (ICM_Cds*) GetSecurity() ;
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ARM_ZeroCurve* InitShortCurve = model->GetZeroCurve();
	
	// const ICM_DefaultCurve** ppDefCurve = new  const ICM_DefaultCurve*[1];	
	std::vector<const ICM_DefaultCurve*> ppDefCurve (1); 
	// const ICM_DefaultCurve* ppDefCurve [1];
	std::auto_ptr<ICM_DefaultCurve> DefCurve((ICM_DefaultCurve*) GetDefCurveIndex().Clone());
	ARM_Vector vRecov(1);
	vRecov[0] = DefCurve->GetRecovery();		
	ppDefCurve[0] = DefCurve.get();
	ICM_ModelMultiCurves DefModel( ppDefCurve, InitShortCurve, vRecov);
	
	PricerCds.Set(cds,&DefModel,ICM_Parameters(),GetAsOfDate());
	PricerCds.SetFaster(true);
//-------------------------------- idem CptPricingMeFlat -------------------
	double value = PricerCds.ComputeSpread(MtM);

	if(MtM == 0)
		setValue(qCMPSPREAD, value);	

	return value;
}

//*******************************************************************************************
// Compute Fwd Spread
//*******************************************************************************************

double ICM_Pricer_CDSIndex::Compute_Fwd_Spread(const ARM_Date &Mty, const ARM_Date & ExpiryDate, double& FlatRbp_Fwd)
{
	double FwdSpread = 0. ;
	double	Rbp_fwd = 0.0;
	
	FwdSpread	=((ICM_ModelMultiCurves*) GetModel())->GetDefaultCurves(0)->FwdSpread_AsIndex(Mty, ExpiryDate, FlatRbp_Fwd, Rbp_fwd);
	return FwdSpread ;
}

//*******************************************************************************************
// Compute le spread relatif au CDS virtuel
//*******************************************************************************************

double ICM_Pricer_CDSIndex::GetVirtualCDSSpread(const ARM_Date & Mty)
{
	return GetDefCurveIndex().ImpliedSpreadInterpol(Mty);	
};


void ICM_Pricer_CDSIndex::View(char* id, FILE* ficOut)
{
	FILE* fOut;
	char  fOutName[200];

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
	
	fprintf(fOut, "\n");
	ICM_Pricer::View(id, fOut);
	fprintf(fOut, "\t\t\t ----------------- Index Pricer ----------------- \n\n");
	fprintf(fOut, "\t\t\t ----------------- DEFAULT CURVE INDEX ----------------- \n");
	fprintf(fOut,"DEFAULT CURVE INDEX \t: \n") ;
	unconst(GetDefCurveIndex()).View(id,fOut)  ;

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}
//*******************************************************************************************
// Compute Duration (FlatRbp_Spot)
//*******************************************************************************************
double ICM_Pricer_CDSIndex::ComputeDuration(void)
{
	//ICMQUANTIFYER("ICM_Pricer_CDSIndex::ComputeDuration"); 
	//if (itsFlatCalculation) {
		CptPricingMeFlat(qCMPDURATION);
		return getValue(qCMPDURATION);
	/*} else {
		ICM_Pricer_Cds PricerCds;
		ICM_Cds* cds = (ICM_Cds*) GetSecurity() ;
		ICM_Cds * newcds = NULL;
		ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
		ARM_ZeroCurve* InitShortCurve = model->GetZeroCurve();
		ICM_DefaultCurveModel* pDefModel = NULL;	
		ARM_Date StartDate = model->GetStartDate().AddDays(1);
		ARM_Date EndDate = cds->GetEndDateNA();
		const ICM_DefaultCurve** ppDefCurve = new  const ICM_DefaultCurve*[1];
		
		std::auto_ptr<ICM_DefaultCurve> DefCurve((ICM_DefaultCurve*) GetDefCurveIndex().Clone());
		//StartDate.NextBusinessDay(DefCurve->GetLagStartDate(), DefCurve->GetCurrency());

		double NewIndexSpread = DefCurve->ImpliedSpreadInterpol(cds->GetEndDateNA()) ;
		int size2 = DefCurve->GetRates()->GetSize() ;
		//ARM_Vector = FlatRates(size2,NewIndexSpread/10000.) ;
		ARM_Vector* FlatRates = new ARM_Vector(size2);
		for (int i = 0; i < size2; i++)
			FlatRates->Elt(i)= NewIndexSpread/10000.;
			
		DefCurve->SetRates(FlatRates) ;
		DefCurve->Calibrate() ;

		pDefModel = new ICM_DefaultCurveModel (DefCurve.get(), InitShortCurve);
		newcds = (ICM_Cds*)((ICM_Cds*)DefCurve->stdCDS(StartDate, EndDate)); // Clone dans stdCDS
			
		// ICM_Pricer_Cds PricerCds(newcds,&DefModel);
		PricerCds.Set(newcds,pDefModel,ICM_Parameters(),GetAsOfDate());
		PricerCds.SetFaster(true);

		double FlatRbp = PricerCds.Price(qCMPDURATION);
		SetDuration(FlatRbp) ;
		double Spread = PricerCds.ComputeSpread() ;
		SetSpread(Spread);
		
		if (newcds) delete newcds;
		newcds = NULL;
		if (pDefModel) delete pDefModel;
			pDefModel = NULL;
		if (ppDefCurve) delete [] ppDefCurve;
		ppDefCurve = NULL;
		return FlatRbp;
	}*/
}

// ------------------------------------------- Pricing Measures ----------------------------------------------
double ICM_Pricer_CDSIndex::FeeLegPV () 
{ 
	if (!GetFeeLegPriceFlg()) 
		ComputeUpfrontPay(qCMPPRICE) ;
	
	return GetFeeLegPrice();
}
double ICM_Pricer_CDSIndex::DefLegPV () {
	if (!GetDefLegPriceFlg()) 
		ComputeUpfrontPay(qCMPPRICE) ;
	return GetDefLegPrice();
}


//*******************************************************************************************
// ComputeDefCurveIndex
//*******************************************************************************************
const ICM_DefaultCurve& 
ICM_Pricer_CDSIndex::GetDefCurveIndex()
{
	//ICMQUANTIFYER("ICM_Pricer_CDSIndex::GetDefCurveIndex"); 
	if (itsDefCurveIndex) 
		return  *itsDefCurveIndex; 
	
	const ICM_DefaultCurve* crv = NULL;
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ICM_Cds* cds = (ICM_Cds*) GetSecurity() ;
	//ICM_Credit_Index* index = (ICM_Credit_Index*) cds->GetFeeLeg()->GetCreditIndex() ;
	crv = model->GetDefaultCurve(0);
	if (itsFlatCalculation) {		
		itsDefCurveIndex = crv->createFlatCurve(cds->GetEndDateNA());
	} else {	
		itsDefCurveIndex = dyn_clone(crv);
		/*ICM_Credit_Index* index = (ICM_Credit_Index*) cds->GetFeeLeg()->GetCreditIndex() ;
		index->CptDefaultCurveIndex((ICM_ModelMultiCurves*)model, cds) ;
		// DefCurve appatient à l'index
		itsDefCurveIndex = dyn_clone(index->GetDefCurve()) ;
		if (!itsDefCurveIndex) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_CDSIndex::GetDefCurveIndex: can't compute"); */
	}
	return *itsDefCurveIndex ;

}
// cpt all pricing measure, with a flat curve
void ICM_Pricer_CDSIndex::CptPricingMeFlat(qCMPMETH measure)
{
	ICM_Pricer_Cds PricerCds;
	ICM_Cds* cds = (ICM_Cds*) GetSecurity() ;
	ICM_ModelMultiCurves* model = (ICM_ModelMultiCurves*) GetModel();
	ARM_ZeroCurve* InitShortCurve = model->GetZeroCurve();
	
	// const ICM_DefaultCurve** ppDefCurve = new  const ICM_DefaultCurve*[1];	
	// const ICM_DefaultCurve* ppDefCurve [1];
	std::vector<const ICM_DefaultCurve*> ppDefCurve (1); 
	std::auto_ptr<ICM_DefaultCurve> DefCurve((ICM_DefaultCurve*) GetDefCurveIndex().Clone());
	ARM_Vector vRecov(1);
	vRecov[0] = DefCurve->GetRecovery();		
	ppDefCurve[0] = DefCurve.get();
	ICM_ModelMultiCurves DefModel( ppDefCurve, InitShortCurve, vRecov);
	
	PricerCds.Set(cds,&DefModel,ICM_Parameters(),GetAsOfDate());
	PricerCds.SetFaster(true);
	double value = PricerCds.Price(measure);
	setValue(measure, value);
	if ( measure == qCMPPRICE) {
		double FeeLegPv = PricerCds.Price(qCMPFEELEGPV);
		setValue(qCMPFEELEGPV, FeeLegPv);
		double DefLegPv = PricerCds.Price(qCMPDEFLEGPV);
		setValue(qCMPDEFLEGPV, DefLegPv);
	}

}
