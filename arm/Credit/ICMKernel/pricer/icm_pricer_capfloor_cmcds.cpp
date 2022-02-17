#include "firsttoinc.h"
#include "ICMKernel\pricer\icm_pricer_capfloor_cmcds.h"
#include "ICMKernel\mod\icm_defcurvemodel.h"
#include "ICMKernel\inst\icm_credit_index.h"


/************************************************************************************************************/
/*
/*	Set Method
/*
/***********************************************************************************************************/


void ICM_Pricer_CapFloorCMCDS::Set(ARM_Security *sec, ARM_Model *mod,const ICM_Parameters&params,const ARM_Date&asof)
{
	ICM_Pricer::Set(sec, mod,params,&asof);
}


/************************************************************************************************************/
/*
/*	Computes the Fwd spreads with Convexity Adjustments
/*
/***********************************************************************************************************/

const ARM_Vector& ICM_Pricer_CapFloorCMCDS::GetAdjFwdSpreads(void) 
{
	if(AdjFwdSpreadsFlg)
		return itsAdjFwdSpreads ;
	
	ICM_Leg* Feeleg = (ICM_Leg*) ((ICM_Cmcds*)GetSecurity())->GetFeeLeg();

	if(!GetFeeLegPriceFlg())
		CptCouponsForFeeLeg(*Feeleg, GetModel(),GetAsOfDate());
	
	ICM_Security* security = Feeleg->GetCreditInfos();
	SetFwdSpreads(security->GetUnAdjFwdSpreads()) ;
	SetAdjFwdSpreads(security->GetAdjFwdSpreads()) ;
	SetFwdPV01s(security->GetFwdPV01s());

	return security->GetAdjFwdSpreads() ;
}

/************************************************************************************************************/
/*
/*	Computes the Fwd spreads without any Adjustments
/*
/***********************************************************************************************************/

const ARM_Vector& ICM_Pricer_CapFloorCMCDS::GetFwdSpreads(void)
{
	if(FwdSpreadsFlg)
		return itsFwdSpreads ;

	ICM_Leg* Feeleg = (ICM_Leg*) ((ICM_Cmcds*)GetSecurity())->GetFeeLeg();

	if(!GetFeeLegPriceFlg())
		CptCouponsForFeeLeg(*Feeleg, GetModel(),GetAsOfDate());
	
	ICM_Security* security = Feeleg->GetCreditInfos();
	SetFwdSpreads(security->GetUnAdjFwdSpreads()) ;
	SetAdjFwdSpreads(security->GetAdjFwdSpreads()) ;
	SetFwdPV01s(security->GetFwdPV01s());

	return security->GetUnAdjFwdSpreads() ;
}

/************************************************************************************************************/
/*
/*	Computes the Fwd PV01s (T1 = Reset Date , T2 =  the following Payment Date.
/* If T1 = T2 we use the convention PV01 = 1 => Be carefull when using it!!
/*
/***********************************************************************************************************/

const ARM_Vector& ICM_Pricer_CapFloorCMCDS::GetFwdPV01s(void)
{
	if(FwdPV01sFlg)
		return itsFwdPV01s ;

	ICM_Leg* Feeleg = (ICM_Leg*) ((ICM_Cmcds*)GetSecurity())->GetFeeLeg();

	if(!GetFeeLegPriceFlg())
		CptCouponsForFeeLeg(*Feeleg, GetModel(),GetAsOfDate());
	
	ICM_Security* security = Feeleg->GetCreditInfos();
	SetFwdSpreads(security->GetUnAdjFwdSpreads()) ;
	SetAdjFwdSpreads(security->GetAdjFwdSpreads()) ;
	SetFwdPV01s(security->GetFwdPV01s());

	return security->GetFwdPV01s() ;
}


/************************************************************************************************************/
/*	Computes the Cap and the Floor Values 
/*	Cap = Sum(i=0,..,NumFlows, Caplets)
/*	Floor = Sum(i=0,..,NumFlows, Floorlets)
/************************************************************************************************************/

void ICM_Pricer_CapFloorCMCDS::CapFloorBSPrices() 
{
	ICM_Cmcds* Sec = (ICM_Cmcds*) GetSecurity();
	double CapLevel, FloorLevel ;
	CapLevel = Sec->GetCapLevel() ;
	FloorLevel = Sec->GetFloorLevel() ;
	
	
	if (CapLevel == -1. && FloorLevel == -1.)
	{
		SetCapValue(0.);
		SetFloorValue(0.);
	}
	else
	{
		ICM_Leg* Feeleg =(ICM_Leg*) Sec->GetFeeLeg();
		ICM_Security* security = (ICM_Security*) Feeleg->GetCreditInfos(); 
		ICM_DefaultCurveModel* DefModel = (ICM_DefaultCurveModel*) GetModel();
		ICM_DefaultCurve* DefCurve = (ICM_DefaultCurve*) DefModel->GetDefaultCurve();
		ARM_VolCurve*  Volatility = DefModel->GetVolCurve() ;
		ARM_VolCurve*  CapVolatility = (ARM_VolCurve* )Volatility->Clone();
		ARM_VolCurve*  FloorVolatility = (ARM_VolCurve* ) Volatility->Clone();
		
		//double Vol = Volatility->ComputeVolatility(1, Sec->GetCapLevel())/100; // ???
		double VolCap,VolFloor;
		
		
		const ARM_Vector& PayDates = security->GetPayDates();
		int	NumFlows = PayDates.size() ; // 14514 security->GetNumFlows();
		ARM_Date Asof = DefModel->GetStartDate();
		security->ComputeYF(DefModel->GetStartDate());
		
		double CapValue = 0. , FloorValue = 0. ;
		double Mty , time_sqrt, FwdSpread, AccruedPeriod, PV01;
		
		double Nd1 = 0., Nd2 = 0., N_d2 = 0., N_d1 = 0. , d1, d2 ;

		const ARM_Vector& AdjFwdSpreads = GetAdjFwdSpreads();
		const ARM_Vector& FwdPV01s = GetFwdPV01s();

		const ARM_Vector& YFFwdStartDates = security->GetYFFwdStartDates();
		const ARM_Vector& YFFwdEndDates = security->GetYFFwdEndDates();

		for (int i = 0; i < NumFlows; i++)
		{	
			AccruedPeriod = (YFFwdStartDates.Elt(i)-YFFwdEndDates.Elt(i)) ;
			Mty = YFFwdStartDates.Elt(i);
			time_sqrt = sqrt(Mty);
			
			FwdSpread = 100*AdjFwdSpreads.Elt(i);
			PV01 = FwdPV01s.Elt(i) ;

			VolCap = CapVolatility->ComputeVolatility(Mty, CapLevel)/100.;
			VolFloor = FloorVolatility->ComputeVolatility(Mty, FloorLevel)/100.;

			// Caplet Value Computation

			if (CapLevel != -1.)
			{
				d1 = (log((FwdSpread)/(CapLevel))+ VolCap*VolCap*Mty /2.)/(VolCap*time_sqrt) ;
				d2 = d1-(VolCap*time_sqrt) ;
				Nd1 = cdfNormal(d1) ;
				Nd2 = cdfNormal(d2) ;
				
				CapValue += (FwdSpread * Nd1 - CapLevel * Nd2) * PV01 ;
			}
			// Floorlet Value Computation
			if (FloorLevel != -1.)
			{
				d1 = (log((FwdSpread)/(FloorLevel))+ VolFloor*VolFloor*Mty /2.)/(VolFloor*time_sqrt) ; 
				d2 = d1-(VolFloor*time_sqrt) ;
				N_d1 = cdfNormal(-d1) ;
				N_d2 = cdfNormal(-d2) ;

				FloorValue += (FloorLevel * N_d2 - FwdSpread * N_d1) * PV01 ;
			}
		}
		
		delete CapVolatility;
		delete FloorVolatility;

		SetCapValue(CapValue);
		SetFloorValue(FloorValue);
	}


}
/************************************************************************************************************/
//	
//Computes the Participation Rate
//	
/************************************************************************************************************/

double ICM_Pricer_CapFloorCMCDS::ComputeSpread(const double& MtM)
{
	if (GetSpreadFlg()) return GetSpread();

	ICM_Cmcds* Sec = (ICM_Cmcds*) GetSecurity();
	ICM_Leg* Feeleg =(ICM_Leg*) Sec->GetFeeLeg();
	double PartRate = Feeleg->GetRefSpread();
	ICM_DefaultCurveModel* DefModel = (ICM_DefaultCurveModel*) GetModel();
	ARM_Date Asof = DefModel->GetStartDate();

	double FeeLegPrice = FeeLegPV();
	double DefLegPrice = DefLegPV();
	double Result  = DefLegPrice/(FeeLegPrice / PartRate) ;

	SetSpread(Result);
	return Result ;
}

/************************************************************************************************************/
//	
//Computes the FeeLeg Price
//	
/************************************************************************************************************/

double ICM_Pricer_CapFloorCMCDS::FeeLegPV()
{
	if(GetFeeLegPriceFlg()) return GetFeeLegPrice() ;

	ICM_Cds* Sec = (ICM_Cds*) GetSecurity();
	ICM_DefaultCurveModel* DefModel = (ICM_DefaultCurveModel*) GetModel();
// 	ICM_Pricer_Cds* pricer = new ICM_Pricer_Cds(Sec, DefModel) ;
	ICM_Pricer_Cds  pricer ; pricer.Set(Sec, DefModel,ICM_Parameters(),GetAsOfDate()) ;
	double CapValue = GetCapValue() ;
	double FloorValue = GetFloorValue() ;
	// ARM_Date Asof = DefModel->GetStartDate();
	double Notional = Sec->GetFeeLeg()->GetNotionalAmount() ;
	double PartRate = ((ICM_Cmcds*)GetSecurity())->GetFeeLeg()->GetRefSpread() ;
	
	// double  PL_CDS = pricer.FeeLegPV() ;
	double  PL_CDS = pricer.Price(qCMPFEELEGPV) ;

	double Result = PL_CDS - PartRate*(CapValue - FloorValue)* Notional/10000;

	SetFeeLegPrice(Result) ;

	// if(pricer)
	// 	delete pricer ;
	// pricer = NULL ;

	return  (Result);
}

/************************************************************************************************************/
//	
//Computes the DefLeg Price
//	
/************************************************************************************************************/

double ICM_Pricer_CapFloorCMCDS::DefLegPV ()
{	
	if(GetDefLegPriceFlg()) return GetDefLegPrice() ;

	ICM_Cds* Sec = (ICM_Cds*) GetSecurity();
	ICM_DefaultCurveModel* DefModel = (ICM_DefaultCurveModel*) GetModel();
	// ICM_Pricer_Cds* pricer = new ICM_Pricer_Cds(Sec, DefModel) ;
	ICM_Pricer_Cds  pricer ; pricer .Set(Sec, DefModel,ICM_Parameters(),GetAsOfDate()) ;
	// ARM_Date Asof = DefModel->GetStartDate();

	// double Result = pricer.DefLegPV() ;
	double Result = pricer.Price(qCMPDEFLEGPV) ;

	SetDefLegPrice(Result) ;

	// if(pricer)
	// 	delete pricer ;
	// pricer = NULL ;
	
	return  (Result) ;
}


/************************************************************************************************************/
//
// Compute Sensitivity
//	
/************************************************************************************************************/

double ICM_Pricer_CapFloorCMCDS::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsvalue, double epsilonGamma // useless
		)
{
	ICM_DefaultCurveModel* DefModel = (ICM_DefaultCurveModel*) GetModel();
	ARM_Date ExecutionDate = DefModel->GetStartDate();
	ICM_Leg* Feeleg = (ICM_Leg*) ((ICM_Cmcds*)GetSecurity())->GetFeeLeg();
	ICM_Credit_Index* CreditIndex = (ICM_Credit_Index*) Feeleg->GetCreditIndex()->Clone();

	double result = 0.0;
	double initialprice = 0.0;
	double modifiedprice = 0.0;

	 
	{
		if (GetInitialPriceFlg())
			initialprice = GetInitialPrice();
		else
			initialprice = ComputePrice(qCMPPRICE);

		ResetPricer();
		
		switch (typesensi)
		{
		case ICMRECOVERY_TYPE :
		case ICMIRCURVE_TYPE :
		case ICMIRCURVE_WITHOUT_DEFCURVE_TYPE :
		case ICMSPREAD_TYPE :
		case ICMSPRELSHIFT_TYPE :
		{
				ICM_DefaultCurveModel* ModelDef = DefModel->GenerateShiftModel(typesensi,plot,epsvalue);
				ICM_Credit_Index* New_CreditIndex = (ICM_Credit_Index*) CreditIndex->Clone() ;

				SetModel(ModelDef);
				
				New_CreditIndex->adoptDefCurve((ICM_DefaultCurve*)ModelDef->GetDefaultCurve()->Clone());
				// New_CreditIndex->SetForcedCurve((ICM_DefaultCurve*)ModelDef->GetDefaultCurve()->Clone());
				New_CreditIndex->SetForcedCurve( ModelDef->GetDefaultCurve() );
				
				Feeleg->SetCreditIndex((ICM_Credit_Index*)New_CreditIndex->Clone()) ;

				modifiedprice = ComputePrice(qCMPPRICE);

				SetModel(DefModel); //On reset le model initial		

				Feeleg->SetCreditIndex((ICM_Credit_Index*)CreditIndex->Clone()) ;

				if (ModelDef)
					delete ModelDef;
				ModelDef = NULL;

				if (New_CreditIndex)
					delete New_CreditIndex;
				New_CreditIndex = NULL;
				
		}
		break;
		case ICM_THETA_Type :
		{
			ARM_ZeroCurve* InitShortCurve = (ARM_ZeroCurve*) DefModel->GetZeroCurve();			
			ICM_DefaultCurve* DefCurve = (ICM_DefaultCurve*) DefModel->GetDefaultCurve();
			ARM_ZeroCurve* ShortCurve2 = (ARM_ZeroCurve*) DefCurve->GetZeroCurve();
			ARM_VolCurve* volcurve = (ARM_VolCurve*) DefModel->GetVolCurve();
			
			ARM_Date Asof = InitShortCurve->GetAsOfDate();
			ARM_Date NextBusDay = InitShortCurve->GetAsOfDate().NextBusinessDay(1) ;
			
			DefCurve->SetAsOfDate(NextBusDay);
			InitShortCurve->SetAsOfDate(NextBusDay) ;
			ShortCurve2->SetAsOfDate(NextBusDay);
			volcurve->SetAsOfDate(NextBusDay) ;
						
			// DefCurve->GenerateShiftCurve(0.0,ICMIRCURVE_TYPE);

			ICM_DefaultCurve* DefCurveIndex = (ICM_DefaultCurve*) CreditIndex->GetDefCurve();
			ICM_DefaultCurve* ForcedCurveIndex = (ICM_DefaultCurve*) CreditIndex->GetForcedCurve();
			ARM_ZeroCurve* ShortCurveIndex1 = (ARM_ZeroCurve*) DefCurveIndex->GetZeroCurve();
			ARM_ZeroCurve* ShortCurveIndex2 = (ARM_ZeroCurve*) ForcedCurveIndex->GetZeroCurve();
			ARM_VolCurve* VolCurveIndex = (ARM_VolCurve*)DefModel->GetVolCurve();

			DefCurveIndex->SetAsOfDate(NextBusDay);
			ShortCurveIndex1->SetAsOfDate(NextBusDay);
			VolCurveIndex->SetAsOfDate(NextBusDay);
			ForcedCurveIndex->SetAsOfDate(NextBusDay);
			ShortCurveIndex2->SetAsOfDate(NextBusDay);

			//ForcedCurveIndex->GenerateShiftCurve(0.0,ICMIRCURVE_TYPE);

			DefModel->SetStartDate(NextBusDay) ;
			
			ARM_Date NewExecutionDate = DefModel->GetStartDate();
			
			int DayCount = GetSecurity()->GetDayCount() ;
			int nbdays = DaysBetweenDates(DayCount,Asof, NextBusDay) ;
			
			ResetPricer() ;
			
			modifiedprice = ComputePrice(qCMPPRICE)/nbdays;
					
			DefCurve->SetAsOfDate(Asof);
			InitShortCurve->SetAsOfDate(Asof);
			ShortCurve2->SetAsOfDate(Asof);
			volcurve->SetAsOfDate(Asof);
			
			DefCurveIndex->SetAsOfDate(Asof);
			ShortCurveIndex1->SetAsOfDate(Asof);
			VolCurveIndex->SetAsOfDate(Asof);
			ForcedCurveIndex->SetAsOfDate(Asof);
			CreditIndex->GetForcedCurve()->GetZeroCurve()->SetAsOfDate(Asof);

			// DefCurve->GenerateShiftCurve(0.0,ICMIRCURVE_TYPE);
			
			// ForcedCurveIndex->GenerateShiftCurve(0.0,ICMIRCURVE_TYPE);

			DefModel->SetStartDate(Asof) ;

			ResetPricer() ;
			
			initialprice =initialprice/nbdays ;
		}
		break;
		default :
			result = -99999999.0;
		}

		ResetPricer();

		if (!result)
			result = modifiedprice - initialprice;

	}
 

	return (result);
}

double ICM_Pricer_CapFloorCMCDS::Accrued()
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
}
double ICM_Pricer_CapFloorCMCDS::Compute_Fwd_Spread(const ARM_Date&,const ARM_Date&, double& dur)
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
}
double ICM_Pricer_CapFloorCMCDS::ComputeDuration() 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
}
double ICM_Pricer_CapFloorCMCDS::ComputeImpliedVol(const double&)
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
}
