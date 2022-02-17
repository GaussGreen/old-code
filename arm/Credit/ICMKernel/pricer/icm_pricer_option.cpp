#include "firsttoinc.h"

#include "ICMKernel/pricer/icm_pricer_option.h"
#include "ICMKernel\mod\icm_defcurvemodel.h"

void ICM_Pricer_Option::Init()
{
	SetName(ICM_PRICER_OPTION);
	unsetFlgs();
	its_diffusion_type = qDIFFUSION_SPREADS_LN;
	its_diffusion=NULL;
}

void ICM_Pricer_Option::Reset(void)
{
	ICM_Pricer::Reset();	
	unsetFlgs();
}

ICM_Pricer_Option::ICM_Pricer_Option( ARM_Security *option, ARM_Model *mod,
									 const ICM_Parameters&params,const ARM_Date&asof)
{
	Init();
	Set(option,mod,params,asof);
}
//*******************************************************************************************
// Set method
//*******************************************************************************************
void ICM_Pricer_Option::Set(ARM_Security *option, ARM_Model *mod,
							const ICM_Parameters&params,const ARM_Date&asof)
{
	ICM_Pricer::Set(option, mod,params,&asof);
	ICM_Option* Option = (ICM_Option*) option;
}

// *************************************************************
// Get the Greeks in B&S Framework
// *************************************************************

double ICM_Pricer_Option::ComputeBSGreeks (const int& type)
{
	qSENSITIVITY_TYPE Type = (qSENSITIVITY_TYPE) type;
	double Greek = 0. ;	
	switch (Type)
	{
		case ICM_GREEK_DELTA_TYPE :
		{
			Greek = GetBSDelta();
		}	break;
		case ICM_GREEK_GAMMA_TYPE :
		{
			Greek = GetBSGamma();
		}	break;
		case ICM_GREEK_VEGA_TYPE :
		{
			Greek = GetBSVega();
		}	break;
		case ICM_GREEK_THETA_TYPE :
		{
			Greek = GetBSTheta();
		}	break;
		case ICM_GREEK_RHO_TYPE :
		{
			Greek = GetBSRho();
		}	break;
		default :
			Greek = GetBSDelta();
	}
	return (Greek) ;
}

//*******************************************************************************************
// Calcule la volatilité implicite d'une option
//*******************************************************************************************

double ICM_Pricer_Option::ComputeImpliedVol(const double& Price)
{
	double aux = 0.;
	double ImplVol = 0. ;
    int Nbiter = 0;
	double sup = 500.;
	double inf = 0.;
	double mid = (sup);
	double fc0 =0.;

	ICM_Option* option = dynamic_cast<ICM_Option*>(GetSecurity());
	if(!option) ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Option::ComputeImpliedVol : Security of Pricer is not an option"); 
	ARM_Date startDate = GetModel()->GetStartDate();
	ICM_DefaultCurveModel* model = dynamic_cast<ICM_DefaultCurveModel*>(GetModel());
	if(!model) ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Option::ComputeImpliedVol : Model of Pricer is not a DefCurveModel"); 
	ARM_VolCurve*  pInitialVol = dynamic_cast<ARM_VolCurve*>((model->GetVolCurve())) ;
	if(!pInitialVol) ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Option::ComputeImpliedVol : Vol of Pricer is not a ARM_VolCurve"); 
	ARM_VolFlat aVolFlat(startDate, 0.);

	// case of Clone Of Vol in Model 
	if (model->GetFlgClone()) {
		pInitialVol = dyn_clone(pInitialVol);		
		model->SetVolCurve(&aVolFlat);
		unsetFlg(qCMPPRICE); 
		fc0 = ComputePriceBS();
		if(Price < fc0)
		{
			model->SetVolCurve(pInitialVol) ;
			unsetFlg(qCMPPRICE); 
			aux = ComputePriceBS();
			if(model->GetFlgClone()) delete pInitialVol;
			return 0.0;
		}
		aVolFlat.SetVolatility( 500.);
		model->SetVolCurve(&aVolFlat);
		unsetFlg(qCMPPRICE); 
		aux = ComputePriceBS();
		if (Price < aux )
		{
			while ((fabs(aux - Price) > 0.05) && (Nbiter < 100))
			{
				mid = (sup + inf )/ 2.0;
				aVolFlat.SetVolatility(mid);	
				model->SetVolCurve(&aVolFlat) ;
				unsetFlg(qCMPPRICE); 	
				aux = ComputePriceBS();
				if (aux > Price)
					sup = mid;
				else
					inf = mid;
				Nbiter++;
			}
			ImplVol = mid;
		}
		model->SetVolCurve(pInitialVol) ;
		unsetFlg(qCMPPRICE); 
		aux = ComputePriceBS();
		delete pInitialVol;
		return (ImplVol);
	} else {
	// Case Of No Clone In Model
		model->SetVolCurve(&aVolFlat);
		unsetFlg(qCMPPRICE); 
		fc0 = ComputePriceBS();		
		if (Price < fc0)
		{
			model->SetVolCurve(pInitialVol) ;
			unsetFlg(qCMPPRICE); 
			aux = ComputePriceBS();
			return 0.0;
		}
	
		aVolFlat.SetVolatility(500.);
		unsetFlg(qCMPPRICE); 
		aux = ComputePriceBS();
	
		if (Price < aux )
		{
			while ((fabs(aux - Price) > 0.05) && (Nbiter < 100))
			{
				mid = (sup + inf )/ 2.0;			
				aVolFlat.SetVolatility(mid);	
				unsetFlg(qCMPPRICE); 	
				aux = ComputePriceBS();
				if (aux > Price)
					sup = mid;
				else
					inf = mid;
				Nbiter++;
			}
			ImplVol = mid;
		}

		model->SetVolCurve(pInitialVol) ;
		unsetFlg(qCMPPRICE); 
		aux = ComputePriceBS();
		return (ImplVol);
	}
}

//*******************************************************************************************
// Compute le Spread Fwd du ss-jacent ainsi que la mté fwd (courbe flat en cas de CDS Index)
//*******************************************************************************************
/*
void ICM_Pricer_Option::Compute_Fwd_Spread(void)
{
	if(GetFwdSpreadFlg()) return GetFwdSpread() ;

	ICM_Option* option = (ICM_Option*) GetSecurity() ;
	ARM_Date maturity = option->GetExpiry();
	ARM_Date cdsexpirydate = (ICM_Cds*) option->GetUnderlying()->GetSecurity()->GetEndDateNA();

	ICM_Pricer_CDSIndex* Underlying = (ICM_Pricer_CDSIndex*) option->GetUnderlying() ;
	
	ICM_Cds* cds = (ICM_Cds*) Underlying->GetSecurity()->Clone() ;
	qSecurity_TYPE CdsType = cds->GetSecurityType() ;

	double FwdSpread = 0. ;
	double PV01 = 0. ;
	
	if (CdsType == qCDS_INDEX)
	{
		Underlying = (ICM_Pricer_CDSIndex*) Underlying ;
		FwdSpread = Underlying->Compute_Fwd_Spread(Mty, CDS_ExpiryDate, PV01);
	}

	else
	{
		ICM_DefaultCurveModel* DefModel = (ICM_DefaultCurveModel*) GetModel();
		ICM_DefaultCurve* DefCurve = (ICM_DefaultCurve*) DefModel->GetDefaultCurve();

		FwdSpread = DefCurve->FwdSpread(Mty, CDS_ExpiryDate) ;
		PV01 = DefCurve->RiskyPV01(Mty, CDS_ExpiryDate);
	}

	SetFwdSpread(FwdSpread);
	SetFwdPV01(PV01) ;


	if(cds)
		delete cds ;
	cds = NULL;

	return FwdSpread ;
}
*/
void ICM_Pricer_Option::View(char* id, FILE* ficOut)
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
	fprintf(fOut, "\t\t\t ----------------- Option Pricer ----------------- \n\n");
	fprintf(fOut,"BSDelta \t: %.2lf (%s)\n",getFlgVector("ICM_GREEK_DELTA_TYPE")?getValue("ICM_GREEK_DELTA_TYPE"):0,getFlgVector("ICM_GREEK_DELTA_TYPE")?"OK":"NOK") ;
	fprintf(fOut,"BSGamma \t: %.2lf (%s)\n",getFlgVector("ICM_GREEK_GAMMA_TYPE")?getValue("ICM_GREEK_GAMMA_TYPE"):0,getFlgVector("ICM_GREEK_GAMMA_TYPE")?"OK":"NOK") ;
	fprintf(fOut,"BSVega \t: %.2lf (%s)\n",getFlgVector("ICM_GREEK_VEGA_TYPE")?getValue("ICM_GREEK_VEGA_TYPE"):0,getFlgVector("ICM_GREEK_VEGA_TYPE")?"OK":"NOK") ;
	fprintf(fOut,"BSTheta \t: %.2lf (%s)\n",getFlgVector("ICM_GREEK_THETA_TYPE")?getValue("ICM_GREEK_THETA_TYPE"):0,getFlgVector("ICM_GREEK_THETA_TYPE")?"OK":"NOK") ;
	fprintf(fOut,"BSRho \t: %.2lf (%s)\n",getFlgVector("ICM_GREEK_RHO_TYPE")?getValue("ICM_GREEK_RHO_TYPE"):0,getFlgVector("ICM_GREEK_RHO_TYPE")?"OK":"NOK") ;
	fprintf(fOut,"FwdSpread v: %.2lf (%s)\n",getFlg(qCMPFWDSPREAD)?getValue(qCMPFWDSPREAD):0,getFlg(qCMPFWDSPREAD)?"OK":"NOK") ;
	fprintf(fOut,"FwdPV01 \t: %.2lf \n",getValue(qCMPFWDDURATION)); 

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}

