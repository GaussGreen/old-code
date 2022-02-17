#include "firsttoinc.h"

#include "ICMKernel/pricer/icm_pricer_cdsoption.h"
#include "ICMKernel\mod\icm_defcurvemodel.h"
#include "ICMKernel/pricer/icm_pricer_adviser.h"

ICM_Pricer_CdsOption::ICM_Pricer_CdsOption(ARM_Security *option,  ARM_Model *mod,
										   const ICM_Parameters&params,const ARM_Date&asof)
{
	Init();
	Set(option, mod, params,asof);
}
//*******************************************************************************************
// Set method
//*******************************************************************************************
void ICM_Pricer_CdsOption::Set( ARM_Security *option,  ARM_Model *mod,const ICM_Parameters&params,const ARM_Date&asof)
{
	if((!option) || (!mod))
		ICMTHROW (ERR_INVALID_ARGUMENT,"ICM_Pricer_CdsOption:: Security or Model is null"); 
	ICM_Pricer_Option::Set(option, mod,params,asof);
}

// *************************************************************
// Compute Price and Greeks in B&S Framework
// *************************************************************

double ICM_Pricer_CdsOption::ComputePriceBS()
{
	if(GetPriceFlg()) return GetPrice() ;

	ICM_Option* option = (ICM_Option*) GetSecurity() ;		
	double Strike = option->GetStrike() ;
	double Price  = 0. ;
	ARM_Date maturity = option->GetExpiry();
	ARM_Date cdsexpirydate = option->GetUnderMaturityDate() ;	
	ICM_DefaultCurveModel* model = (ICM_DefaultCurveModel*) GetModel();
	
	int  Ko = option->GetKoType();

	double FwdSpread = 0. ;
	double PV01 = 0. ;
	ARM_Date Asof = model->GetStartDate();
	double dur =0.;
	Compute_Fwd_Values();
	FwdSpread = getValue(qCMPFWDSPREAD); // en bp
	PV01 = getValue(qCMPFWDDURATION);

	double Mty =  CountYears(KACTUAL_360, Asof, maturity); 
	double Tenor = CountYears(KACTUAL_360, maturity, cdsexpirydate) ;
	ARM_VolCurve*  Volatility = dynamic_cast<ARM_VolCurve*> (model->GetVolCurve()) ;
	if( !Volatility)
		ICMTHROW (ERR_INVALID_ARGUMENT,"ICM_Pricer_CdsOption::ComputePriceBS : BAD Volatility for pricing"); 
	double Vol = (1/100.)*Volatility->ComputeVolatility(Mty, Strike);

	double time_sqrt = sqrt(Mty);
	if(time_sqrt == 0){
		//ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_CdsOption::ComputeBS : Gamma is infini because of timeToMaturity = 0"); 
		time_sqrt = 0.000001;
	}
    double d1 = (log(FwdSpread/Strike)+ Vol*Vol*Mty /2.)/(Vol*time_sqrt); 
    double d2 = d1- (Vol*time_sqrt);

	double Nd1 = 0., Nd2 = 0., N_d2 = 0., N_d1 = 0. ;


	Nd1 = cdfNormal(d1) ;
	Nd2 = cdfNormal(d2) ;
	N_d1 = cdfNormal(-d1) ;
	N_d2 = cdfNormal(-d2) ;

	double Delta = 0, Gamma = 0., Vega = 0., Theta = 0., Rho = 0.;
	if (option->IsCall())
	{
		Price = ((FwdSpread * Nd1 - Strike * Nd2) *  PV01)/10000 *option->GetNotional(); // /10 000 for bp-> Value
		double DL = 0. ;
		const ICM_DefaultCurve* DefCurve = (ICM_DefaultCurve*) model->GetDefaultCurve();
		ARM_ZeroCurve* InitShortCurve = (ARM_ZeroCurve*) DefCurve->GetZeroCurve();
			
		ICM_DefaultCurveModel DefModel(DefCurve, InitShortCurve);

		//ICM_Cds* Shortcds = (ICM_Cds*) cds->Clone();
		ARM_Date StartDate = model->GetStartDate(); //.AddDays(1);
		if(StartDate < maturity) {
			if(Ko == 0 ) // No KO and with no Acceleration
			{	
				std::auto_ptr<ICM_Cds> pShortcds((ICM_Cds*)(DefCurve->stdCDS(StartDate, maturity, 0.1, option->GetNotional())));
				ICM_Leg* DefLeg = pShortcds->GetDefLeg() ;
				DefLeg->SetPaymentFreq(K_ZEROCOUPON);
				DefLeg->GetIRIndex()->SetPayFrequency(K_ZEROCOUPON);
				DefLeg->GetIRIndex()->SetResetFrequency(K_ZEROCOUPON);
				DefLeg->CptCashFlowDatesCredit();

				ICM_Pricer_Cds  PricerCds ; 
				PricerCds.Set(pShortcds.get(),&DefModel,ICM_Parameters(),GetAsOfDate());
				PricerCds.SetFaster(true);

				ARM_Date ExecutionDate = DefModel.GetStartDate();
				
				// DL = PricerCds.DefLegPV() ;
				DL = PricerCds.Price(qCMPDEFLEGPV); 
			}
			if (Ko == -1) // No Ko with Acceleration
			{
				std::auto_ptr<ICM_Cds> pShortcds((ICM_Cds*)(DefCurve->stdCDS(StartDate, maturity, 0.1, option->GetNotional())));
				pShortcds->SetCreditLag(0.);
				
				// ICM_Pricer_Cds* PricerCds = new ICM_Pricer_Cds(Shortcds,DefModel);
				ICM_Pricer_Cds  PricerCds ; 
				PricerCds.Set(pShortcds.get(),&DefModel,ICM_Parameters(),GetAsOfDate());
				PricerCds.SetFaster(true);
				ARM_Date ExecutionDate = DefModel.GetStartDate();
				
				// DL = PricerCds.DefLegPV() ;
				DL = PricerCds.Price(qCMPDEFLEGPV) ;
			}
		}
		Price += DL ;
		Delta = Nd1 ;
	}
	else
	{
		Price = (Strike * N_d2 - FwdSpread * N_d1)  *  PV01 /10000 *option->GetNotional() ; // /10 000 for bp-> Value;
		Delta = (Nd1 - 1);
	}

	Delta = Delta;
	SetPrice(Price) ;
	setValue (qCMPPREMIUM,Price);
	Vega  = FwdSpread/10000 * time_sqrt * (1/sqrt(2*PI))*exp(-(d1*d1)/2.0) * PV01/100. * option->GetNotional();
	setValue (qCMP_OPT_AN_DELTA, Delta) ;
	setValue (qCMP_OPT_AN_VEGA,  Vega) ;
	if(time_sqrt == 0) {
		//ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_CdsOption::ComputeBS : Gamma is infini because of timeToMaturity = 0"); 
		time_sqrt = 0.000001;
	}
	Gamma = (1/sqrt(2*PI))*exp(-d1*d1/2.0)/(FwdSpread/10000*Vol*time_sqrt) / option->GetNotional();
	setValue (qCMP_OPT_AN_GAMMA, Gamma) ;
	return Price ;
}

//New Version
double ICM_Pricer_CdsOption::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
												const std::string&  plot, 
												const std::string&  label, 
												double  epsvalue, double epsilonGamma // useless
												)
{
	ICM_DefaultCurveModel* DefModel = (ICM_DefaultCurveModel*) GetModel();
	string SensiString = "";
	ARM_Date ExecutionDate = DefModel->GetStartDate();

	double result = 0.0;
	double initialprice = 0.0;
	double modifiedprice = 0.0;
	double theta = 0.0;
	double rho = 0.0;
	int NbDays = 1;
	ICM_Pricer_Advisor recupPricer;
	std::auto_ptr<ICM_DefaultCurveModel> MktModelShift ;
	string formatedDate = ExecutionDate.toString('.'); // for format in ICM_Pricer

	if (typesensi == ICM_ISSUER_DEFAULT) return CptDefaultPL(label, epsvalue);			
	
	ComputePrice(qCMPPREMIUM);
	initialprice = getValue(qCMPPREMIUM);
		switch (typesensi)
		{
		//case ICMRECOVERY_TYPE :
		case ICMIRCURVE_TYPE :
		case ICMIRCURVE_WITHOUT_DEFCURVE_TYPE :
		//case ICMSPREAD_TYPE :
		//case ICMSPRELSHIFT_TYPE :
		case ICM_GREEK_RHO_TYPE :
		{
			 MktModelShift = auto_ptr<ICM_DefaultCurveModel>(dynamic_cast<ICM_DefaultCurveModel*>(DefModel->GenerateShiftModel(typesensi,plot,epsvalue)));
		}
		break;
		case ICM_THETA_Type :
		{			
			MktModelShift = auto_ptr<ICM_DefaultCurveModel> (dynamic_cast<ICM_DefaultCurveModel*>(DefModel->Clone()));
			ARM_ZeroCurve* ShortCurve = dynamic_cast<ARM_ZeroCurve*>(MktModelShift->GetZeroCurve());	
			if (!ShortCurve){ ICMTHROW(ERR_SWAPLEG_PB, "MktModelShift->GetZeroCurve() is not of type ARM_ZeroCurve");}
			
			ARM_VolCurve* volcurve=dynamic_cast<ARM_VolCurve*>(MktModelShift->GetVolCurve());
			
			std::auto_ptr<ICM_DefaultCurve> DefCurve ( dyn_clone(MktModelShift->GetDefaultCurve() ) ) ;
			ARM_Date NextBusDay = ExecutionDate;
			NextBusDay.NextBusinessDay(1) ;		
			ARM_ZeroLInterpol * newARM_ZeroLInterpol=NULL;

			ARM_Vector epsilon(1,1.);
			//	ARM_CRV_TERMS Term;
			vector<string> Term;
			ShortCurve->SetAsOfDate(NextBusDay) ;
			volcurve->SetAsOfDate(NextBusDay) ;
			MktModelShift->SetStartDate(NextBusDay) ;

			// restripping de la ZC curve
			ARM_ZeroLInterpol* pZCInterpol = dynamic_cast<ARM_ZeroLInterpol  *>( MktModelShift->GetZeroCurve()) ;
			if ( !pZCInterpol ){
					ShortCurve->SetAsOfDate(NextBusDay) ;
					ICMLOG("ICM_Pricer_IndexOption::ComputeSensitivity : compute THETA : the ZC curve is not restripped a the new AsOf "); 
			} else {
					if (!pZCInterpol->GetMktData()) {
						ShortCurve->SetAsOfDate(NextBusDay) ;
						ICMLOG("ICM_Pricer_IndexOption::ComputeSensitivity : compute THETA : the ZC curve is not construct for a restripping a new AsOf "); 
					}else {
						try {
						newARM_ZeroLInterpol = new ARM_ZeroLInterpol(NextBusDay, pZCInterpol->GetMktData()->itsMktTerms,
																		pZCInterpol->GetMktData()->itsMktValue,
																		pZCInterpol->GetMktData()->itsMMVsFut,
																		pZCInterpol->GetMktData()->itsSwapVsFut,
																		pZCInterpol->GetMktData()->itsraw,
																		pZCInterpol->GetMktData()->itsCont_Lin,
																		0 , //pZCInterpol->getgetlastBucketInt()
																		pZCInterpol->GetCurrencyUnit(),
																		K_DEF_FREQ, //swapFrqId,
																		KNOBASE //fixDayCount 
																		);
							DefCurve->SetZeroCurve(newARM_ZeroLInterpol);
							MktModelShift->SetZeroCurve(newARM_ZeroLInterpol);
						}catch(...){
							if(newARM_ZeroLInterpol) 
								delete newARM_ZeroLInterpol;
							// can be because of an expired futur point. No restriping
							ARM_ZeroCurve* newZC = dyn_clone(ShortCurve);
							newZC->SetAsOfDate(NextBusDay) ;
							DefCurve->SetZeroCurve(newZC);
							MktModelShift->SetZeroCurve(newZC);
						}
						
					}		
			}
	
			std::auto_ptr<ICM_DefaultCurve> DefCurveShifted(DefCurve->GenerateShiftCurve(Term, epsilon ,
																   ICM_GREEK_THETA_TYPE));
			MktModelShift->SetDefaultCurve(DefCurveShifted.get()); 
			
			ARM_Date NewExecutionDate = MktModelShift->GetStartDate();
			formatedDate = NewExecutionDate.toString('.'); 
			int dayCount = GetSecurity()->GetDayCount();
			NbDays = DaysBetweenDates(dayCount,ExecutionDate, NextBusDay) ;
		}
		break;
		default :
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_CdsOption::ComputeSensitivity : bad sensi type"); 
		}
#ifdef _DEBUG
		FILE* pFile = NULL;
		pFile = fopen("C:\\temp\\ModelShiftedCDSOption.txt","w");
		MktModelShift->View("",pFile);
		if(pFile) fclose(pFile);
#endif
		auto_ptr<ICM_Pricer_CdsOption> newPricerShift(dynamic_cast<ICM_Pricer_CdsOption*>(recupPricer.GeneratePricer(GetSecurity(),
																			(ARM_Model*)MktModelShift.get(), // CRAD !! mais obligatoire
																			GetName(),
																			CREDIT_DEFAULT_VALUE, // because nbpaths not used 
																			&GetParameters(),
																			(char*)formatedDate.c_str()) ));

		modifiedprice = newPricerShift->ComputePrice(qCMPPREMIUM);
		if (!result)
			result = (modifiedprice - initialprice)/NbDays;
		this->setValue(SensiString,result);
	return (result);
}


double ICM_Pricer_CdsOption::CptDefaultPL(const std::string& label, double epsvalue)
{
	double result = 0. ;
	double initialprice = 0.0;
	double modifiedprice = 0.0;
	
	{
		if (GetInitialPriceFlg())
			initialprice = GetInitialPrice();
		else
			initialprice = ComputePriceBS();	
		Reset();
		if (!result)
			result = modifiedprice - initialprice;
	}
	return (result);
}


//*******************************************************************************************
// Compute le Spread Fwd du ss-jacent ainsi que la mté fwd (courbe flat en cas de CDS Index)
//*******************************************************************************************

void ICM_Pricer_CdsOption::Compute_Fwd_Values(void ) //const ARM_Date &Mty, const ARM_Date &CDS_ExpiryDate, double& dur)
{
	double FwdSpread =0., PV01 = 0.;
	if(getFlg(qCMPFWDSPREAD)) return;

	ARM_Date expiry = ((ICM_Option*)GetSecurity())->GetExpiry();
	ARM_Date cdsMaturity = ((ICM_Option*)GetSecurity())->GetUnderMaturityDate();
	
	ICM_DefaultCurveModel* DefModel = dynamic_cast<ICM_DefaultCurveModel*>(GetModel());
	if(!DefModel) ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_CdsOption::Compute_Fwd_Values : DefModel of Pricer is not a DefCurveModel"); 
	
	const ICM_DefaultCurve* DefCurve = dynamic_cast<const ICM_DefaultCurve*>(DefModel->GetDefaultCurve());
	if(!DefCurve) ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_CdsOption::Compute_Fwd_Values : DefCurve of Model is not a DefCurve"); 
	
	FwdSpread = DefCurve->FwdSpread(expiry, cdsMaturity) ;
	PV01 = DefCurve->RiskyPV01(expiry, cdsMaturity);
	setValue(qCMPFWDSPREAD,FwdSpread);
	setValue(qCMPFWDDURATION,PV01);
}

