#include "firsttoinc.h"

#include "ICMKernel/pricer/icm_pricer_indexoption.h"
#include "ICMKernel/mod/modelmulticurves.h"
#include "ICMKernel/pricer/icm_pricer_adviser.h"
#include "ICMKernel\util\icm_integrator.h"
#include "ICMKernel\mod\icm_vasicek_diffusion.h"


void ICM_Pricer_IndexOption::ICM_Pricer_IndexOptionvoid(ARM_Security *option, 
								ARM_Model *mod,const ICM_Parameters&params,
								const ARM_Date&asof)
{
	Init();
	Set(option, mod,params,asof);
}
//*******************************************************************************************
// Set method
//*******************************************************************************************
void ICM_Pricer_IndexOption::Set(ARM_Security *option, ARM_Model *mod,const ICM_Parameters&params,
								 const ARM_Date&asof)
{
	SetName(ICM_PRICER_INDEXOPTION);
	if((!option) || (!mod))
		ICMTHROW (ERR_INVALID_ARGUMENT,"ICM_Pricer_IndexOption:: Security or Model is null"); 
	ICM_Pricer_Option::Set(option, mod,params,asof);
}

// -------------------------------------------------------------------------------------------
// Compute Price and Greeks in B&S Framework 
//	Price + Delta, Gamma, Vega
// -------------------------------------------------------------------------------------------
double ICM_Pricer_IndexOption::ComputePriceBS()
{
	qDIFFUSION_TYPE difftype;

	ARM_Vector* DIFFUSION_TYPE = GetParameters().GetColVect("DIFFUSION_TYPE");
	if (DIFFUSION_TYPE)
	{
		difftype = (qDIFFUSION_TYPE) (int) DIFFUSION_TYPE->Elt(0);
		SetDiffusionType(difftype);
	}

	switch (difftype)
	{
		case qDIFFUSION_SPREADS_LN:
		case qDIFFUSION_NONE:
		default :
			{return ComputePriceBS_LN_SPREADS();}
		case qDIFFUSION_OU:
			{return ComputePriceBS_OU();}
	}


}

extern "C" void 
func_Payoff_OU(void* Param, double x, double& res)
{
	// get parameters from stack
	ICM_Pricer_IndexOption*	TheModel;
	TheModel	=	(ICM_Pricer_IndexOption*)((*(AddressVector*)Param)[0]);

	res	=	TheModel->Payoff_OU(x,Param);

}

double ICM_Pricer_IndexOption::Payoff_OU(double x,void* params)
{

	// HERMITE
	//if (itsIntegrationMethod	==	qGAUSS_HERMITE)
	x	*=	sqrt(2.);

	ICM_Pricer_IndexOption*	thepricer = NULL;
	thepricer	=	(ICM_Pricer_IndexOption*)((*(AddressVector*)params)[0]);

	//double Xt=

	return 0.;
}

// -------------------------------------------------------------------------------------------
//	Orstein-Uhlenbeck diffusion for intensity
// -------------------------------------------------------------------------------------------
double ICM_Pricer_IndexOption::ComputePriceBS_OU()
{
	ICM_VasicekDiffusion VasicekDiff(GetParameters());
	SetDiffusion(&VasicekDiff);

	double result =0.;
	int nointegrationcoef=-1;
	AddressVector	TheParameterVector;
	TheParameterVector.Append(this);
	TheParameterVector.Append(&nointegrationcoef);

	ICM_Integrator	TheIntegrator;
	TheIntegrator.SetIntegrationStep(40);
	TheIntegrator.SetIntegrationType(qGAUSS_HERMITE);
	TheIntegrator.Integrate(-6.0, 6.0, func_Payoff_OU, &TheParameterVector,result);

	return (result);
}


// *************************************************************
// Compute Price and Greeks in B&S Framework 
//	Price + Delta, Gamma, Vega
// *************************************************************

double ICM_Pricer_IndexOption::ComputePriceBS_LN_SPREADS()
{
	if(GetPriceFlg()) return GetPrice() ;
	ICM_Option* option = (ICM_Option*) GetSecurity() ;
		
	double Strike = option->GetStrike();
	double Price  = 0. ;

	ARM_Date maturity = option->GetExpiry();
	ARM_Date cdsexpirydate = option->GetUnderMaturityDate();
	ICM_ModelMultiCurves* model = dynamic_cast<ICM_ModelMultiCurves*>(GetModel());
	if(!model) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_IndexOption::ComputePriceBS : DefModel of Pricer is not a ICM_ModelMultiCurves"); 	
	double FwdSpread = 0. ;
	double PV01 = 0. ;
	ARM_Date Asof = GetModel()->GetStartDate();
	double dur = 0.;
	Compute_Fwd_Value() ;
	FwdSpread = getValue(qCMPFWDSPREAD); // en bp
	PV01 = getValue(qCMPFWDDURATION);
	
	double Mty =  CountYears(KACTUAL_360, Asof, maturity);
//	double Tenor = CountYears(KACTUAL_360, maturity, cdsexpirydate) ;

	ARM_VolCurve*  Volatility = model->GetVolCurve() ;
	if( !Volatility)
		ICMTHROW (ERR_INVALID_ARGUMENT,"ICM_Pricer_CdsOption::ComputePriceBS : BAD Volatility for pricing"); 
	double Vol = (1/100.)*Volatility->ComputeVolatility(Mty, Strike);
	
	double time_sqrt = sqrt(Mty);
	if(time_sqrt == 0){
		//ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_CdsOption::ComputeBS : Gamma is infini because of timeToMaturity = 0"); 
		time_sqrt = 0.000001;
	}
    double d1 = (log(FwdSpread/Strike)+ Vol*Vol*Mty /2.)/(Vol*time_sqrt); 
    double d2 = d1-(Vol*time_sqrt);

	double Nd1 = 0., Nd2 = 0., N_d2 = 0., N_d1 = 0. ;

	Nd1 = cdfNormal(d1) ;
	Nd2 = cdfNormal(d2) ;
	N_d1 = cdfNormal(-d1) ;
	N_d2 = cdfNormal(-d2) ;

	double Delta = 0, Gamma = 0., Vega = 0., Theta = 0., Rho = 0.;

	if (option->IsCall())
	{
		Price = (FwdSpread * Nd1 - Strike * Nd2) *  PV01 /10000 *option->GetNotional(); // /10 000 for bp-> Value
		Delta = Nd1 ;
	}
	else
	{
		Price = (Strike * N_d2 - FwdSpread * N_d1) *  PV01 /10000 *option->GetNotional(); // /10 000 for bp-> Value
		Delta = Nd1 - 1 ;
	}
	Delta = Delta;
	SetPrice(Price) ;
	setValue(qCMPPREMIUM, Price);
	Vega  = FwdSpread/10000 * time_sqrt * (1/sqrt(2*PI))*exp(-(d1*d1)/2.0) * PV01/100. ;
	setValue (qCMP_OPT_AN_VEGA,  Vega) ;
	setValue (qCMP_OPT_AN_DELTA, Delta) ;
		
	Gamma = (1/sqrt(2*PI))*exp(-d1*d1/2.0)/(FwdSpread/10000*Vol*time_sqrt);
	setValue (qCMP_OPT_AN_GAMMA, Gamma) ;
	return Price ;
}

//*******************************************************************************************
// Compute Sensitivity
//*******************************************************************************************

double ICM_Pricer_IndexOption::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
												  const std::string&  plot, 
												  const std::string & label, 
												  double  epsvalue,  double epsilonGamma //useless
												  )
{
    double sensitivity =0.;

	ICM_ModelMultiCurves* DefModel = (ICM_ModelMultiCurves*) GetModel();
	int nbcurves = DefModel->GetNbDefCurves();
	int nocurve = -1;
	double theta = 0.0;
	double rho = 0.0;
	double nbdays = 1.0 ;
	string SensiString = "";
	ICM_Pricer_Advisor recupPricer;
	ARM_Date ExecutionDate = GetModel()->GetStartDate();
	std::auto_ptr<ICM_ModelMultiCurves> MktModelShift; // autoPTR !!
	string formatedDate = ExecutionDate.toString('.'); // for format in ICM_Pricer
	std::auto_ptr<ICM_DefaultCurve> defaultCurveShifted;
	std::auto_ptr<ARM_ZeroLInterpol> newARM_ZeroLInterpol;

	int i = 0;
	double result = 0.,initialprice = 0.,modifiedprice = 0.;

	if (typesensi == ICMRECOVERY_BAR_TYPE)
	{
		ICMTHROW(ERR_UNIMP_METHOD_CALL,"Unimplemented ICMRECOVERY_BAR_TYPE method");
	}
		// pricing
		ComputePrice(qCMPPREMIUM);
		initialprice = getValue(qCMPPREMIUM);
		switch (typesensi)
		{
			//case ICMRECOVERY_TYPE :
			case ICMIRCURVE_TYPE : //rho
			case ICMIRCURVE_WITHOUT_DEFCURVE_TYPE :
			//case ICMSPREAD_TYPE :
			//case ICMSPRELSHIFT_TYPE : 
			case ICM_GREEK_RHO_TYPE :
 			{		
				MktModelShift = auto_ptr<ICM_ModelMultiCurves>(dynamic_cast<ICM_ModelMultiCurves*>
																(DefModel->GenerateShiftModel(typesensi,
																								plot, 
																								label,
																								nocurve,
																								epsvalue)));
			}
			break;
			case ICM_THETA_Type :
			{
				MktModelShift = auto_ptr<ICM_ModelMultiCurves>(dynamic_cast<ICM_ModelMultiCurves*>(DefModel->Clone()));
				ARM_ZeroCurve* InitShortCurve=dynamic_cast<ARM_ZeroCurve*>( MktModelShift->GetZeroCurve());						
				ARM_VolCurve* Initialvolcurve=dynamic_cast<ARM_VolCurve*>(MktModelShift->GetVolCurve());
				std::auto_ptr<ICM_DefaultCurve> defaultCurve(dyn_clone(MktModelShift->GetDefaultCurve(0))) ;

				ARM_Date NextBusDay = ExecutionDate;
					NextBusDay.NextBusinessDay(1) ;
				ARM_Date CalculusDay = NextBusDay;					
				double DayCount = GetSecurity()->GetDayCount() ;
				nbdays = DaysBetweenDates(DayCount,ExecutionDate, NextBusDay) ;
					
				// restripping de la ZC curve
				ARM_ZeroLInterpol* pZCInterpol = dynamic_cast<ARM_ZeroLInterpol  *>( MktModelShift->GetZeroCurve()) ;
				if ( !pZCInterpol ){
						InitShortCurve->SetAsOfDate(CalculusDay) ;
						ICMLOG("ICM_Pricer_IndexOption::ComputeSensitivity : compute THETA : the ZC curve is not restripped a the new AsOf "); 
				} else {
					if (!pZCInterpol->GetMktData()) {
							InitShortCurve->SetAsOfDate(CalculusDay) ;
							ICMLOG("ICM_Pricer_IndexOption::ComputeSensitivity : compute THETA : the ZC curve is not construct for a restripping a new AsOf "); 
					}else {
						try{
							newARM_ZeroLInterpol = std::auto_ptr<ARM_ZeroLInterpol>( new ARM_ZeroLInterpol(CalculusDay, pZCInterpol->GetMktData()->itsMktTerms,
																		pZCInterpol->GetMktData()->itsMktValue,
																		pZCInterpol->GetMktData()->itsMMVsFut,
																		pZCInterpol->GetMktData()->itsSwapVsFut,
																		pZCInterpol->GetMktData()->itsraw,
																		pZCInterpol->GetMktData()->itsCont_Lin,
																		0 /*pZCInterpol->getgetlastBucketInt()*/,
																		pZCInterpol->GetCurrencyUnit(),
																		K_DEF_FREQ /*swapFrqId*/,
																		KNOBASE /*fixDayCount */));
							defaultCurve->SetZeroCurve(dyn_clone(newARM_ZeroLInterpol.get()));
							MktModelShift->SetZeroCurve(newARM_ZeroLInterpol.get());
						}catch(...){
							// can be because of an expired futur point.
							std::auto_ptr<ARM_ZeroCurve> cZC(dyn_clone(InitShortCurve));
							cZC->SetAsOfDate(CalculusDay) ;
							defaultCurve->SetZeroCurve(dyn_clone(cZC.get()));
							MktModelShift->SetZeroCurve(cZC.get());
						}
							
					}
                       
					
				}
					
				Initialvolcurve->SetAsOfDate(CalculusDay) ;		
				MktModelShift->SetStartDate(CalculusDay) ;
				//// fait la calibration de la defCurve qui servira aux calcul du FwdSprad dans
				// les cas =! qCredit_Adjust20, qCredit_CDSDTRX, qCredit_CDSDIND , qCredit_CDSINDZ
				ARM_Vector epsilon(1,1.);
				vector<string> Term;
				ICM_DefaultCurve* defaultCurveShifted = defaultCurve->GenerateShiftCurve(Term, epsilon ,
																   ICM_GREEK_THETA_TYPE);	// New
				MktModelShift->SetDefaultCurve(0,defaultCurveShifted); // adoption will detruct defaultCurveShifted
				ARM_Date NewExecutionDate = MktModelShift->GetStartDate();
				formatedDate = NewExecutionDate.toString('.'); 
			}
			break;
			default :
					ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_IndexOption::ComputeSensitivity : bad sensi type"); 
		}
#ifdef _DEBUG
		FILE* pFile = NULL;
		pFile = fopen("C:\\temp\\ModelShifted.txt","w");
		MktModelShift->View("",pFile);
		if(pFile) fclose(pFile);
#endif
		auto_ptr<ICM_Pricer_IndexOption> newPricerShift(dynamic_cast<ICM_Pricer_IndexOption*>(recupPricer.GeneratePricer(GetSecurity(),
																			(ARM_Model*)MktModelShift.get(), // CRAD !! mais obligatoire
																			GetName(),
																			CREDIT_DEFAULT_VALUE, // because nbpaths not used 
																			&GetParameters(),
																			(char*)formatedDate.c_str()) ));

	modifiedprice = newPricerShift->ComputePrice(qCMPPREMIUM);
	if (!result) result=(modifiedprice - initialprice)/nbdays;
	this->setValue(SensiString,result);
	return (result);
}


//*******************************************************************************************
// Compute le Spread Fwd du ss-jacent ainsi que la mté fwd (courbe flat en cas de CDS Index)
//*******************************************************************************************

void ICM_Pricer_IndexOption::Compute_Fwd_Value(void) //const ARM_Date &Mty, const ARM_Date &CDS_ExpiryDate,double& PV01)
{
	double FwdSpread = 0. ,  PV01 = 0., fwRbp = 0.;
	if(getFlg(qCMPFWDSPREAD)) return;
	ICM_Option* option = (ICM_Option*) GetSecurity() ;
	ARM_Date expiry = option->GetExpiry();
	ARM_Date cdsexpirydate = option->GetUnderMaturityDate();
	ICM_ModelMultiCurves* lModel = dynamic_cast<ICM_ModelMultiCurves*>(GetModel());
	FwdSpread = lModel->GetDefaultCurve(0)->FwdSpread_AsIndex(expiry, cdsexpirydate, PV01, fwRbp);
	setValue(qCMPFWDSPREAD, FwdSpread);
	setValue(qCMPFWDDURATION, PV01) ;
}


// virtual 
void ICM_Pricer_IndexOption::Reset(void)
{
	ICM_Pricer_CdsOption::Reset();
}
