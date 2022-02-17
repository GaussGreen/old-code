#include "firsttoinc.h"
#include "ICMKernel/util/icm_utils.h"
#include "ICMKernel/pricer/icm_basepricer.h"
#include "ICMKernel/pricer/icm_pricer_cds.h"
#include "ICMKernel/inst/icm_credit_index.h"
#include "ICMKernel/mod/modelmulticurves.h"
#include "ICMKernel\util\icm_rootfinder1D.h"


//	-----------------------------------------------------------------------------------------------------
ICM_BasePricer::ICM_BasePricer()
{
	Init();
}
//	-----------------------------------------------------------------------------------------------------
// virtual 
ICM_BasePricer::~ICM_BasePricer()
{
	if (itspPricer) delete itspPricer; 
}
//	-----------------------------------------------------------------------------------------------------
void 
ICM_BasePricer::Set(ICM_Credit_Index* sec,ICM_ModelMultiCurves*mod,const ICM_Parameters&params,const ARM_Date&asof)
{
	ICM_Pricer::Set(sec,mod,params,&asof); 
}
//	-----------------------------------------------------------------------------------------------------
// virtual 
void 
ICM_BasePricer::DoPriceVector(qVECTMETH measure)
{
	switch(measure)
	{
	case qINDEXBASE: 
		doPriceIndexBase(); break; 
	case qMARKETBASE :
		doPriceMarketBase(); break; 
	case qTRUEMARKETBASE :
	case qTRUEMARKETSPREAD: 
		doPriceTrueMarketBase(); break; 
	default:
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_BasePricer::DoPriceVector: unhandled measure "<<measure); 
	}
}

//-------------------------------------------------------------------------------------------------------
void
ICM_BasePricer::doPriceMarketBase()
{
	ICM_Credit_Index* index = dynamic_cast<ICM_Credit_Index*>(GetSecurity()) ;
	ICM_ModelMultiCurves* mod = dynamic_cast<ICM_ModelMultiCurves*>(GetModel()) ;
	// 
	int basketSize = index->GetNbCurves(); 
	const ICM_DefaultCurve* indexCurve = mod->GetDefaultCurve(index->GetCreditIndexName()); 
	const ARM_Vector* curveDates  = indexCurve->GetDates(); 
	const ARM_Vector* spreadValues = indexCurve->GetRates(); 
	int nbPoints = curveDates->GetSize(); 
	ARM_Vector bases(nbPoints-1,0.);
	ARM_Vector spreadIntrinseques(nbPoints-1,0.);
	double PV_CDSs = 0.;
	//FILE * pFile= NULL;
	//pFile = fopen("C:\\temp\\testMarketBase.txt","w");

	for(int j=1;j<nbPoints;j++) 
	{
		// creation CDS:
		double Impliedspread = indexCurve->ImpliedSpreadInterpol(curveDates->Elt(j)); // en BP
		Impliedspread /=10000; // en valeur
		//fprintf(pFile,"\t Impliedspread for piller no %i : %3.10lf\n", j, Impliedspread*10000);
		std::auto_ptr<ICM_Cds> pCds((ICM_Cds*) indexCurve->stdCDS(this->GetAsOfDate(), 
																curveDates->Elt(j),Impliedspread,1,true));	
		
		PV_CDSs = 0.;
		// creation FlatCurve
		std::auto_ptr<ICM_DefaultCurve> pFlatCurve(indexCurve->createFlatCurve(curveDates->Elt(j)));
		double FlatDuration = pFlatCurve->RiskyDuration(curveDates->Elt(j));
		//fprintf(pFile,"\tFlatDuartion for piller no %i : %3.10lf\n", j, FlatDuration);
		for(int i=0;i<basketSize;i++) 
		{
			const ICM_DefaultCurve* pDefCurve = mod->GetDefaultCurve(index->GetLabel(i));			
			// creation du model :
			ICM_DefaultCurveModel DefModel(pDefCurve,mod->GetZeroCurve(),NULL, false);
			//Pricing
			ICM_Pricer_Cds PricerCDS;
			ICM_Parameters aICM_Parameters;
			PricerCDS.Set(pCds.get(), &DefModel,aICM_Parameters, GetAsOfDate());
			double pvCds = PricerCDS.ComputePrice(qCMPPRICE);
			PV_CDSs += pvCds;
			//DefModel.View("",pFile);
			//pCds->View("",pFile);
			//fprintf(pFile,"\t \tPV CDS for piller no %i, Name no %i : %3.10lf\n", j, i,  10000*pvCds);
			//fprintf(pFile,"\t \tFeeLeg CDS for piller no %i, : %lf\n", i, PricerCDS.ComputePrice(qCMPFEELEGPV) );
			//fprintf(pFile,"\t \tDEfLeg CDS for piller no %i, : %lf\n", i, PricerCDS.ComputePrice(qCMPDEFLEGPV) );	
		}
		PV_CDSs /=basketSize;
		//fprintf(pFile,"\t PV CDS for average no %i,  : %3.10lf\n", j,  10000*PV_CDSs);
		spreadIntrinseques.Elt(j-1) = - PV_CDSs/FlatDuration + Impliedspread;
		//fprintf(pFile,"\t spreadIntrinseques  no %i,  : %3.10lf\n", j,  10000*spreadIntrinseques.Elt(j-1));
		bases.Elt(j-1)				= spreadIntrinseques.Elt(j-1)/Impliedspread;
		//fprintf(pFile,"\t Base  no %i,  : %3.10lf\n", j,  10000*bases.Elt(j-1) );
	}

	adoptVectorValue(qMARKETBASE,dyn_clone(&bases)); 
	//if(pFile) fclose(pFile);
}
//--------------------------------------------------------------------------------------------------------
void ICM_BasePricer::doPriceTrueMarketBase()
{
	ICM_Credit_Index* index = dynamic_cast<ICM_Credit_Index*>(GetSecurity()) ;
	ICM_ModelMultiCurves* mod = dynamic_cast<ICM_ModelMultiCurves*>(GetModel()) ;
	// 
	int basketSize = index->GetNbCurves(); 
	const ICM_DefaultCurve* indexCurve = mod->GetDefaultCurve(index->GetCreditIndexName()); 
	const ARM_Vector* curveDates  = indexCurve->GetDates(); 
	const ARM_Vector* curveYT  = indexCurve->GetYearTerms();
	const ARM_Vector* spreadValues = indexCurve->GetRates(); 
	int nbPlots = curveDates->GetSize(); 
	ICM_Parameters aICM_Parameters;
	ARM_Vector bases(nbPlots-1,0.);
	ARM_Vector spreads(nbPlots-1,0.);
	double result = 1.;
	ARM_Vector spreadIntrinseques(nbPlots-1,0.);
	double PV_CDSs = 0.;
	//FILE * pFile= NULL;
	//pFile = fopen("C:\\temp\\testTrueMarketBase.txt","w");
	// cpt FirstCptEffDate -------------------------------------------------------
	std::auto_ptr<ARM_Vector> SchedSubStart;
	ARM_Date tmpStartDate = this->GetAsOfDate();
	ARM_Date FirstCouponEffDate;
	const int stdFreq = indexCurve->GetSTDCDS_Frequency();
	string ccy = indexCurve->GetCurrency();
	tmpStartDate.AddPeriod(-1*stdFreq,(char*)(ccy.c_str()));

		SchedSubStart = std::auto_ptr<ARM_Vector>(CptStartDates(tmpStartDate,(ARM_Date)curveDates->Elt(1),stdFreq,/*fwdrule ?? */ 2,
									K_SHORTSTART,indexCurve->GetSTDCDS_Adjusted() /*? intrule*/,(char*) ccy.c_str(),/*adjFirstdate*/1) );
		if(	SchedSubStart->size() >1)
			FirstCouponEffDate = (ARM_Date)(SchedSubStart->Elt(1));
		else { 
			ICMTHROW(ERR_INVALID_ARGUMENT,"No FirstCouponEffDate can be computed "); 
		}

		for(int j=1;j<nbPlots;j++) 
		{
			// creation CDS:
			bool isYT = index->GetIsYT();
			if (isYT){
				if (curveYT->Elt(j) != index->GetMaturities().Elt(j-1) ) ICMTHROW(ERR_INVALID_ARGUMENT," YT of Index doesn't correspond to YT of defcurveIndex, Use DateFormat "); 
			}else{
				if (curveDates->Elt(j) != index->GetMaturities().Elt(j-1) )	ICMTHROW(ERR_INVALID_ARGUMENT," Dates of Index doesn't correspond to Dates of defcurveIndex "); 
			}
			double couponSpread = index->GetRunning().Elt(j-1); // en BP
			couponSpread /=10000; // en valeur
		//	fprintf(pFile,"\t couponSpread for piller no %i : %3.10lf\n", j, couponSpread*10000);
			std::auto_ptr<ICM_Cds> pCds((ICM_Cds*) indexCurve->stdCDS(FirstCouponEffDate, 
																	curveDates->Elt(j),
																	couponSpread,10000,false));	
			
			PV_CDSs = 0.;
			for(int i=0;i<basketSize;i++) 
			{
				const ICM_DefaultCurve* pDefCurve = mod->GetDefaultCurve(index->GetLabel(i));			
				// creation du model :
				ICM_DefaultCurveModel DefModel(pDefCurve,mod->GetZeroCurve(),NULL, false);
				//Pricing
				ICM_Pricer_Cds PricerCDS;
				PricerCDS.Set(pCds.get(), &DefModel,aICM_Parameters, GetAsOfDate());
				double pvCds = PricerCDS.ComputePrice(qCMPPRICE);
				PV_CDSs += pvCds;
				//DefModel.View("",pFile);
				//pCds->View("",pFile);
				//fprintf(pFile,"\t \tPV CDS for piller no %i, Name no %i : %3.10lf\n", j, i,  pvCds);
			}
			PV_CDSs /=basketSize;
			//fprintf(pFile,"\t PV CDS for average no %i,  : %3.10lf\n", j,  PV_CDSs);
			itsPVoptimum = PV_CDSs;
			// find spread from "MTM Cds(FlatCurve) = PV_CDSs
			std::auto_ptr<ICM_Cds> pCdsIndex((ICM_Cds*) indexCurve->stdCDS(FirstCouponEffDate, 
																	curveDates->Elt(j),couponSpread,10000,false));	
			//fprintf(pFile," pCdsIndex : \n");
			//pCdsIndex->View("",pFile);
			// creation FlatCurve
			std::auto_ptr<ICM_DefaultCurve> pFlatCurve(indexCurve->createFlatCurve(curveDates->Elt(j)));
			double impliedSpread = pFlatCurve->ImpliedSpreadInterpol(curveDates->Elt(j));
			ICM_DefaultCurveModel DefModelIndex(pFlatCurve.get(),mod->GetZeroCurve(),NULL, true);
			//DefModelIndex.View("",pFile);
			if (itspPricer) delete itspPricer;
			itspPricer = new ICM_Pricer_Cds();
			((ICM_Pricer_Cds*)itspPricer)->Set(pCdsIndex.get(), &DefModelIndex,aICM_Parameters, GetAsOfDate());
			// ne pas avoir un spread Null sinon l'algo ne converge pas.
			double _inf = MAX( MAX(impliedSpread/10000-1.E-4,0.0), MAX(impliedSpread/10000-1.E-2, 0.0)); // en valeur
			//double _inf = impliedSpread/10000-1.E-4;
			double _sup = impliedSpread/10000+1.E-2;
			const qDEFCURVE_CALIB_ALGO CalibrationAlgo = indexCurve->GetCalibrationAlgo(); //qDEFCURVE_NEWTON;
			double spread = 0.;
			result = RootFinder1D(ff1::mem_call(&ICM_BasePricer::Evaluate,(*this))).ZeroBracketDecreasing(_inf,_sup);	
			// rate cannot be null because of calibration in reset rate. Pv( rate =0 , lambda =0) => 0
			if(( _inf <0.) || (_inf ==0.)) {
				_inf = 1.E-6;
			}
			if (!result) 
				ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_BasePricer::doPriceTrueMarketBase :  can't bracket spread"
					<< "from inf="<<_inf<<" sup="<<_sup); 
			
			switch(CalibrationAlgo)
				{														
				case qDEFCURVE_DICHO:
					spread = RootFinder1D(ff1::mem_call(&ICM_BasePricer::Evaluate,(*this))).Dichotomy(_inf,_sup,100,1.E-5,1.E-2,false); break;
				case qDEFCURVE_NEWTON:
					spread = RootFinder1D(ff1::mem_call(&ICM_BasePricer::Evaluate,(*this))).NewtonRaphsonWithBisection(_inf,_sup,100,0.00001,1.E-5,1.E-2); break;
				default:
					ICMTHROW(ERR_INVALID_ARGUMENT,"Can't find Spread for FlatCurve "<<CalibrationAlgo); 
			}
			//fprintf(pFile,"\t Spread CDS for average no %i,  : %3.10lf\n", j,  10000*spread);
			spreads.Elt(j-1) = spread*10000;
			bases.Elt(j-1)=  spread/impliedSpread*10000;
			//fprintf(pFile,"\t Base for average no %i,  : %3.10lf\n", j,  bases.Elt(j-1));
		}

		adoptVectorValue(qTRUEMARKETBASE,dyn_clone(&bases)); 
		adoptVectorValue(qTRUEMARKETSPREAD,dyn_clone(&spreads)); 
	//	if(pFile) fclose(pFile);
}
//	-----------------------------------------------------------------------------------------------------
// 
void
ICM_BasePricer::doPriceIndexBase()
{
	ICM_Credit_Index* index = dynamic_cast<ICM_Credit_Index*>(GetSecurity()) ;
	ICM_ModelMultiCurves* mod = dynamic_cast<ICM_ModelMultiCurves*>(GetModel()) ;
	// 
	int basketSize = index->GetNbCurves(); 
	const ICM_DefaultCurve* indexCurve = mod->GetDefaultCurve(index->GetCreditIndexName()); 
	const ARM_Vector* curveDates  = indexCurve->GetDates(); 
	const ARM_Vector* spreadValues = indexCurve->GetRates(); 
	int nbPoints = curveDates->GetSize(); 
	// 
	// ICM_QMatrix<double> spreads(basketSize,nbPoints),durations(basketSize,nbPoints); 
	ARM_Vector rates(nbPoints-1,0.); 
	ARM_Vector sumdurations(nbPoints-1,0.); 
	for(int i=0;i<basketSize;i++) 
	{
		const ICM_DefaultCurve* cdsCurve = mod->GetDefaultCurve(index->GetLabel(i)); 
		// starting from 1, since the point#0 is not 'real'
		for(int j=1;j<nbPoints;j++) 
		{
			// spreads(i,j) = mod->GetDefaultCurve(index->GetLabel(i))->ImpliedSpreadInterpol(curveDates->Elt(j)); 
			// durations(i,j) = mod->GetDefaultCurve(index->GetLabel(i))->RiskyDuration(curveDates->Elt(j)); 
			double spread = cdsCurve->ImpliedSpreadInterpol(curveDates->Elt(j)); 
			double duration = cdsCurve->RiskyDuration(curveDates->Elt(j)); 
			rates.Elt(j-1) += spread * duration; 
			sumdurations.Elt(j-1) += duration ;
		}
	}
	//	
	for(int j=1;j<nbPoints;j++) 
	{
		rates.Elt(j-1) /= ( sumdurations.Elt(j-1) * spreadValues->Elt(j) *10000. ); 
	}
	// rates.Elt(0)=rates.Elt(1); // handling point#0 

	adoptVectorValue(qINDEXBASE,dyn_clone(&rates)); 
}
//	-----------------------------------------------------------------------------------------------------
// virtual 
double ICM_BasePricer::ComputeSpread(const double& MtM ) 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}
//	-----------------------------------------------------------------------------------------------------
// virtual	
double ICM_BasePricer::ComputeImpliedVol(const double& Price) 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}//	-----------------------------------------------------------------------------------------------------
// virtual 
double ICM_BasePricer::Accrued() 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}//	-----------------------------------------------------------------------------------------------------
// virtual 
double ICM_BasePricer::FeeLegPV () 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}//	-----------------------------------------------------------------------------------------------------
// virtual 
double ICM_BasePricer::DefLegPV () 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}//	-----------------------------------------------------------------------------------------------------
// virtual 
double ICM_BasePricer::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi,
	const std::string& plot,
	const std::string& label,
	double epsilon, double epsvalueGamma ) 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}//	-----------------------------------------------------------------------------------------------------
// virtual 
double ICM_BasePricer::ComputeDuration()
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not implemented"); 
}

double ICM_BasePricer::Evaluate(const double& spread)
{
	ResetSpread( spread );
	itspPricer->ResetPricer();
	double npv = itspPricer->Price(qCMPPRICE);		
	return (npv-itsPVoptimum);
}

// FIXMEFRED: mig.vc8 (28/05/2007 10:39:22):cast
void ICM_BasePricer::ResetSpread(const double& spread)
{
	ICM_DefaultCurveModel* pDefModel = (ICM_DefaultCurveModel*)itspPricer->GetModel();
	std::auto_ptr<ICM_DefaultCurve> pDefCurve(dyn_clone(pDefModel->GetDefaultCurve()));
	ARM_Vector vSpreads(2);
	vSpreads.Elt(0) = spread;
	vSpreads.Elt(1) = spread;
	pDefCurve->setRates(vSpreads);
	pDefCurve->Calibrate();
	pDefModel->SetDefaultCurve(pDefCurve.get()); // flat clone = Y -> delete ensuite
	((ICM_Pricer_Cds*)itspPricer)->LoadMarketData(pDefModel,GetAsOfDate());
}