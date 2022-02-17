#include "ARMKernel\glob\firsttoinc.h"
#include "icm_pricer_tree_binomial_cds.h"
#include "ICMKernel/glob/icm_mktdatamng.h"
//#include "ICMKernel\util\icm_macro.h"
// #include "ICMKernel\inst\icm_cds.h"
#include "ICMKernel/crv/icm_defaultcurve.h"
#include "ICMKernel/inst/icm_cds.h"
#include "ICMKernel\util\icm_utils.h"

#include <set>
#include <limits>

//	------------------------------------------------------------------------
std::string ICM_Pricer_Tree_Binomial_Cds::PARAM_FORCEVOL("forceVol") ; 
std::string ICM_Pricer_Tree_Binomial_Cds::PARAM_USEDEFCURVEPOINTS("useDefCurvePoints");  
std::string ICM_Pricer_Tree_Binomial_Cds::PARAM_OPTIONDISCSTEPSIZE("discStepSize"); 
std::string ICM_Pricer_Tree_Binomial_Cds::PARAM_OPTIONDISCSTEPNUMBER("discStepNumber"); 
std::string ICM_Pricer_Tree_Binomial_Cds::PARAM_USECONSTANTPROBA("useConstantProba"); 
// std::string ICM_Pricer_Tree_Binomial_Cds::PARAM_USEOLDDEFLEGPRICING("useOldDefLegPricing"); 
std::string ICM_Pricer_Tree_Binomial_Cds::PARAM_DISCTOLERANCE("daysDiscTolerance"); 
std::string ICM_Pricer_Tree_Binomial_Cds::PARAM_DISCMAXITER("discMaxIter"); 
std::string ICM_Pricer_Tree_Binomial_Cds::PARAM_DISCVARIANCETOL("discVarianceTol"); 

//	------------------------------------------------------------------------
ICM_Pricer_Tree_Binomial_Cds::ICM_Pricer_Tree_Binomial_Cds() 
{ 
	Init();
}
//	------------------------------------------------------------------------
ICM_Pricer_Tree_Binomial_Cds::ICM_Pricer_Tree_Binomial_Cds(ARM_Security *sec, ARM_Object *mod,const ICM_Parameters& parameters,const ARM_Date&asof)
{
	Init();
	Set(sec, mod,parameters,asof);
}

//	------------------------------------------------------------------------
//	virtual 
ICM_Pricer_Tree_Binomial_Cds::~ICM_Pricer_Tree_Binomial_Cds()
{
}
//	------------------------------------------------------------------------
// virtual 
void ICM_Pricer_Tree_Binomial_Cds::MarketData(ARM_Security* sec,vector<string>& DMKT)
{
	int nbdmk=0;DMKT.resize(nbdmk+1);
	ICM_Cds* cds = (ICM_Cds*)sec;
	DMKT[nbdmk]= GetZeroCurveName(cds->GetCurrencyUnit()->GetCcyName(),GetAsOfDate());nbdmk++;DMKT.resize(nbdmk+1);
	DMKT[nbdmk]= GetDefaultCurveName(cds->GetSingleName(),GetAsOfDate());
	return;
}

//	------------------------------------------------------------------------
//virtual 
ARM_CLASS_NAME 
ICM_Pricer_Tree_Binomial_Cds::GetRootName()
{
	return ICM_PRICER; 
}
//	------------------------------------------------------------------------
void ICM_Pricer_Tree_Binomial_Cds::Init()
{
	SetName(ICM_BINOMIAL_TREE_PRICER_CDS);
	itsTree.clear(); 
}
//	------------------------------------------------------------------------
//	virtual 
void 
ICM_Pricer_Tree_Binomial_Cds::Reset(void)
{
	parent::Reset();
	itsTree.clear(); 
}
//	------------------------------------------------------------------------
void 
ICM_Pricer_Tree_Binomial_Cds::Set(ARM_Security *_sec, ARM_Object *mod,const ICM_Parameters&parameters,const ARM_Date&asof)
{
	if (!_sec) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Tree_Binomial_Cds::Set: no security"); 
	ICM_Cds& cds=dynamic_cast<ICM_Cds&>(*_sec); 
	ICM_Leg& feeLeg=*cds.GetFeeLeg(); 
	ICM_Security&sec=*feeLeg.GetCreditInfos(); 
	parent::Set(&cds,mod,parameters,&asof);

	ARM_Date asofDate = GetAsOfDate(); 
	ICM_DefaultCurve* defCurve = GetMktDataMng()->GetDefaultCurve(cds.GetSingleName(),asofDate); 
	if (!defCurve) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Tree_Binomial_Cds::Set: Can't fine Default Curve "<<cds.GetSingleName()); 
	double vol; 
	getParam(ICM_Pricer_Tree_Binomial_Cds::PARAM_FORCEVOL,vol); 
	// 
	
	
	//	-- CDS relevant dates are the feeleg pay dates
	ICMLOG("ICM_Pricer_Tree_Binomial_Cds::Set: Calibrating tree for "<<cds.GetSingleName());  
	const ARM_Vector&payDates = sec.GetPayDates() ; 
	

	unsigned long first; 
	unsigned long size=sec.countPaymentFlows(first,asofDate); 
	
	//	Create the SET with appropriate comparison scheme.
	double yfTol = std::numeric_limits<double>::epsilon(); 
	long daysTol=0; 
	getParam(PARAM_DISCTOLERANCE,daysTol,false); 
	if (daysTol!=0) yfTol=fabs(daysTol/365.); 
	tol_comp myComp(yfTol) ;
	std::set<double,tol_comp > times ( myComp ) ; 

	times.insert(0);

	//	The tree discretisation should match all CDS payment dates
	
	cds.GetFeeLeg()->ComputeYF(GetAsOfDate()) ;

	const ARM_Vector& YFAccstartDates = feeLeg.GetCreditInfos()->GetYFAccStartDates();
	const ARM_Vector& YFAccPayDates   = feeLeg.GetCreditInfos()->GetYFPayDates();
	
	double  yfMaturity = YFAccPayDates.Elt(size-1); 
	
	//for(unsigned int i=0;i<size;i++) 
	// 	times.insert( yfMaturity =CountYears(KACTUAL_365,asofDate,ARM_Date(payDates->Elt(first+i))) ) ;

	times.insert(YFAccstartDates.Elt(0));

	for(unsigned int i=0;i<size;i++) 
	 	times.insert(YFAccPayDates.Elt(i)) ;

		
	//	We optionnaly don't add the DefCurve points
	long useDefCurvePoints(1); // true , C convention
	getParam(PARAM_USEDEFCURVEPOINTS,useDefCurvePoints,false) ; 

	if (useDefCurvePoints) 
	{
		for(unsigned int i=0;i<defCurve->GetSize();i++) 
			times.insert(CountYears(KACTUAL_365,asofDate,
				ARM_Date(defCurve->GetDate(i)))
				); 
	}

	//	Adapt the discretization so that variance changes are minimal
	// 
	double  userStepSize ;
	long	userStepNumber ;
	double yfStep ;
	if (!getParam(PARAM_OPTIONDISCSTEPNUMBER,userStepNumber ,false)) userStepNumber=0; 
	if (!getParam(PARAM_OPTIONDISCSTEPSIZE,userStepSize ,false)) userStepSize=0; 
	yfStep = yfMaturity; 
	if (userStepNumber!=0 || userStepSize!=0)  
	{
		if (userStepNumber==0 && userStepSize>0) 
			yfStep = std::_cpp_min(userStepSize,yfStep) ; 
		else if (userStepNumber>0 && userStepSize==0) 
		{	yfStep=std::_cpp_min(yfMaturity/userStepNumber,yfStep) ; }
		else ICMTHROW(ERR_INVALID_ARGUMENT,"stepSize="<<userStepSize<<" stepNumber="<<userStepNumber) ; 
	}
	// now we consider that yfStep should be the min for all time steps. 
	std::set<double,tol_comp>::iterator it (times.begin()); 
	std::set<double,tol_comp>::iterator next_it(it); 
	++next_it ; 
	while ( next_it != times.end())
	{
		yfStep = std::_cpp_min(yfStep,*next_it -*it); 
		++it; ++next_it; 
	}
	double minDT = yfStep; 
	//	Here we iterate over yfStep in order to minimize variance changes
	//
	//		yfStep(n+1) = min(DT) / ( [min(DT)/yfStep] +1 )
	//
	//		criterion is 
	//			*	maxStep is attained (might be set to 0)
	//			*	varianceChange <= varianceTol 
	//		variance is given by vol*vol*yfStep
	double discVarianceTol (yfStep*vol*vol+1) ;
	long discMaxIter(0); 
	getParam(PARAM_DISCMAXITER,discMaxIter,false) ;
	getParam(PARAM_DISCVARIANCETOL,discVarianceTol,false) ;

	//	now we are looking for the it that will minimize the local variance, 
	//	here at constant volatility this is the one minimizing 
	long currentIter(0); 
	double varMin ; 
	double varMax ; 
	double varCriteria; 
	std::set<double,tol_comp> adaptativeTimes(times); 
	do 
	{
		if (currentIter!=0) yfStep = minDT/(int(minDT/yfStep) +1. ); 
		ICMLOG("yfStep="<<yfStep<<" at iteration "<<currentIter) ;
		double targetVariance = vol*vol*yfStep ; 
		it = times.begin(); next_it=it; 
		++next_it; 
		while (next_it != times.end())
		{
			targetVariance = std::_cpp_min(targetVariance,vol*vol*(*next_it- *it) ); 
			++it; ++next_it; 
		}
		//	now for each of the defined step, we should have Ni = [ localVariance / targetVariance ] 
		//	division of length DTi/[Ni]. 
		//	We have to work with a copy of the initial discretisation. 
		adaptativeTimes=times ; 
		it=times.begin(); next_it=it; 
		++next_it ;
		varMin = 100; 
		varMax = 0 ;
		while (next_it != times.end())
		{
			double DT = *next_it - *it; 
			double N = vol*vol*DT / targetVariance ;
			int Nint = N +0.5 ;		// this is the rouding to closest integer
			// int Nint=N; 
			//		hence	N=1.1 : we subdivise in 1 interval
			//				N=1.9 : we subdivise in 2 interval. 
			for(unsigned int j=1;j<Nint;j++) adaptativeTimes.insert(*it + j*DT/Nint) ;
			if (Nint>1) 
			{
				varMin = std::_cpp_min(varMin,vol*vol*DT/Nint); 
				varMax = std::_cpp_max(varMax,vol*vol*DT/Nint); 
			}
			else 
			{
				varMin = std::_cpp_min(varMin,vol*vol*DT); 
				varMax = std::_cpp_max(varMax,vol*vol*DT); 
			}
			++it; 
			++next_it; 
		}
		varCriteria=fabs(varMax-varMin) ; 
		ICMLOG("yfStep="<<yfStep<<" at iteration "<<currentIter<<" varCriteria="<<varCriteria<<" varMin="<<varMin<<" varMax="<<varMax) ;
		currentIter++; 
	} 
	while (currentIter<discMaxIter && gt(varCriteria,discVarianceTol) ) ;

	// creating the mesh
	ICM_QMatrix<double> theTimes(adaptativeTimes.size(),1); 
	it= adaptativeTimes.begin(); 
	i=0; 
	while (it!=adaptativeTimes.end()) 
	{ 
// 		ICMLOG(" Adding time yf="<<*it); 
		theTimes(i,0)=*it; 
		++it; 
		++i; 
	}


	itsTree.setTimes(theTimes); 

	long useConstantProba(0); 
	if (!getParam(PARAM_USECONSTANTPROBA,useConstantProba,false)) useConstantProba=0; 
	if (useConstantProba) 
		TreeServices::calibrate(itsTree,
			*defCurve,
			vol
			); 
	else 
		TreeServices::calibrate2(itsTree,
			*defCurve,
			vol
			); 

	// 
	// itsUseOldDefLegPricing=0; 
	// getParam(PARAM_USEOLDDEFLEGPRICING,itsUseOldDefLegPricing,false) ; 
}

//	------------------------------------------------------------------------
// virtual 
// double 
// ICM_Pricer_Tree_Binomial_Cds::Price(const int& mode, char* ExecutionDate ) 
// {
// 	return ICM_Pricer::Price(mode,ExecutionDate); 
// }
//	------------------------------------------------------------------------
// virtual 
double 
ICM_Pricer_Tree_Binomial_Cds::ComputePrice(qCMPMETH measure)
{
	return ICM_Pricer::ComputePrice(measure); 
}
//	------------------------------------------------------------------------
// virtual 
double ICM_Pricer_Tree_Binomial_Cds::ComputeSpread(const double& MtM )  
{
	if (GetSpreadFlg()) return GetSpread() ; 

	ICM_Cds& cds=dynamic_cast<ICM_Cds&>(* GetSecurity() ); 
	ICM_Leg& feeLeg=*cds.GetFeeLeg(); 
	
	ARM_Date useless; 

	double Spread = feeLeg.GetCreditSpread()* DefLegPV() / FeeLegPV() *100. ;

	SetSpread(Spread);
	
	return Spread;
		
}
//	------------------------------------------------------------------------
// virtual 
double ICM_Pricer_Tree_Binomial_Cds::Accrued() 
{
	// Here we use the execution date 

	//	--	Retrieve the fee leg
	ICM_Cds& cds=dynamic_cast<ICM_Cds&>(* GetSecurity() ); 
	ICM_Leg& feeLeg=*cds.GetFeeLeg(); 
	ICM_Security&sec=*feeLeg.GetCreditInfos(); 
	
	// 14514 int NumFlows = sec.GetNumFlows();
	const ARM_Vector& AccStartDates = sec.GetAccStartDates();
	const ARM_Vector& AccEndDates = sec.GetAccEndDates();
	const ARM_Vector& PayDates = sec.GetPayDates();
	int NumFlows = AccStartDates.size();

	double Accrued = 0.,Coupon =0.,DP=0.,SP=0.;
	bool AccruedFlag = false;

	ARM_ZeroCurve* discCurve= GetMktDataMng()->GetZeroCurve(
		cds.GetCurrencyUnit()->GetCcyName(), GetAsOfDate() );
	ICM_DefaultCurve * defCurve= GetMktDataMng()->GetDefaultCurve(cds.GetSingleName(),GetAsOfDate()) ;
	// this is copy/paste from the standard CDS pricer. 
	// 
	for (int i=0; i<NumFlows; i++)
	{
		if ((GetAsOfDate().GetJulian()>AccStartDates.Elt(i)) &&  (GetAsOfDate().GetJulian()>=AccEndDates.Elt(i)))
		{Coupon = 0.;}
//		JLA. encore un example assez instructif... 
//		on checke l'ExecutionDate (qui n'a pas de sens en général des les autres pricing) 
//		pour finalement utiliser l'AsOfDate ... 
//		else if ((ExecutionDate.GetJulian()>AccStartDates->Elt(i)) &&  (ExecutionDate.GetJulian()<AccEndDates->Elt(i)))
//		{Coupon = sec.AccruedCoupon(GetAsOfDate());	AccruedFlag = true;}		
		else if ((GetAsOfDate().GetJulian()>AccStartDates.Elt(i)) &&  (GetAsOfDate().GetJulian()<AccEndDates.Elt(i)))
		{Coupon = sec.AccruedCoupon(GetAsOfDate());	AccruedFlag = true;}

		/*
		if (AccruedFlag)
		{
		// Computes Accrued
		switch (cds.GetFeeLeg()->GetAccruedOnDefault())
		{
		case qACC_In_Risk : 
			DP = discCurve->DiscountPrice((ARM_Date)PayDates.Elt(i));
			SP = defCurve->SurvivalProba((ARM_Date)AccEndDates.Elt(i));
			Accrued = Coupon * DP * SP;
			break;
		case qNOS :
			Accrued = 0.;
			break;
		case qACCRUED_SETTLED : 
		default : 
			Accrued = Coupon;
		}
		break;
		}
		*/
	}

	// ICMTHROW(ERR_INVALID_MODEL,"Not Implemented") ; 
	SetAccrued(Accrued); 
	return Accrued; 
}

//	------------------------------------------------------------------------
//	virtual 
void 
ICM_Pricer_Tree_Binomial_Cds::View(char* id, FILE* ficOut)
{
	FILE* fOut;
	char  fOutName[200];

	if ( ficOut == NULL )
	{
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w"); 
	}
	else	fOut = ficOut;

	int size =0;

	fprintf(fOut, "\t\t\t ----------------- Binomial Tree CDS Pricer ----------------- \n\n");
	std::stringstream sstr ;
	TreeServices::dumpIntensities(sstr,itsTree); 
	sstr<<std::endl; 
	TreeServices::dumpProbas(sstr,itsTree); 
	sstr<<std::endl; 
	sstr<<itsTree; 
	fprintf(fOut,"%s\n",sstr.str().c_str()); 
	fprintf(fOut, "\n");


	parent::View(id,fOut); 
	if ( ficOut == NULL )fclose(fOut);
}

//	------------------------------------------------------------------------
// virtual 
double ICM_Pricer_Tree_Binomial_Cds::DefLegPV () 
{

	if (GetDefLegPriceFlg()) return GetDefLegPrice();

	//	--	Retrieve the def leg
	ICM_Cds& cds=dynamic_cast<ICM_Cds&>(* GetSecurity() ); 
	ICM_Leg& defLeg=*cds.GetDefLeg(); 
	ICM_Security&sec=*defLeg.GetCreditInfos(); 

	ARM_ZeroCurve* discCurve= GetMktDataMng()->GetZeroCurve(
		cds.GetCurrencyUnit()->GetCcyName(), GetAsOfDate() );	

	double losses = 
		(1.-GetMktDataMng()->GetDefaultCurve(cds.GetSingleName(),GetAsOfDate())->GetRecovery()) 
		* sec.GetNotionals().Elt(0) ; 
		
	double ret; 

	const ARM_Vector& payDates=sec.GetPayDates(); 
	const ARM_Vector& startDates=sec.GetAccStartDates(); 

	unsigned long first,size; 
	size=sec.countPaymentFlows(first,GetAsOfDate());

	long backwardFrom(itsTree.depth()-1); 
	long backwardTo(0); 

	ICM_QMatrix<double> flows_def(backwardFrom,2);
	
	ARM_Date EffectiveDate = startDates.Elt(first);
	double aux = CountYears(KACTUAL_365,GetAsOfDate(),EffectiveDate) ;

	for(unsigned int i=0;i<backwardFrom;i++)
	{
		if (itsTree.Time(i)<aux)
		{
			flows_def(i,0) = 0. ;
			flows_def(i,1) = 0. ;
		}
		else
		{
			flows_def(i,0)= 0.5 *(itsTree.Time(i)+itsTree.Time(i+1)) + cds.GetCreditLag()/365.;
			flows_def(i,1)=losses ;
		}
	}
	
	ICM_QMatrix<double> statesFrom(itsTree.slice(backwardFrom).size(),1); 
	ICM_QMatrix<double> ret2 ;
	ICM_QMatrix<double> ret2_greeks ;

	TreeServices::backwardprice_default3(
								itsTree,
								backwardFrom,
								backwardTo,
								flows_def,
								*discCurve,
								statesFrom, 
								ret2,
								ret2_greeks); 

	SetDefLegPrice(ret=ret2(0,0)); 

	SetDefLegPrice(ret);

	return ret; 
}


//----------------------------------------------------------------------
// Feeleg PV Calculation. integration steps are tree dates
//----------------------------------------------------------------------

double ICM_Pricer_Tree_Binomial_Cds::FeeLegPV ( ) 
{

	if (GetFeeLegPriceFlg()) return GetFeeLegPrice();

	//	--	Retrieve the fee leg
	ICM_Cds& cds=dynamic_cast<ICM_Cds&>(* GetSecurity() ); 
	ICM_Leg& feeLeg=*cds.GetFeeLeg(); 
	ICM_Security&sec=*feeLeg.GetCreditInfos(); 

	const ARM_Vector& payDates=sec.GetPayDates(); 
	const ARM_Vector& startDates=sec.GetAccStartDates(); 
	//	--	
	unsigned long first,size; 
	size=sec.countPaymentFlows(first,GetAsOfDate()); 
	
	ICM_QMatrix<double> flows(size,2); 
	
	for(unsigned int i=0;i<size;i++)
	{
		flows(i,0)=CountYears(KACTUAL_365,GetAsOfDate(),ARM_Date(payDates.Elt(first+i))) ;
		//flows(i,1)=sec.FullCoupon(startDates.Elt(first+i)); 
		flows(i,1)=sec.FullCoupon(first+i); 
	}

	ARM_ZeroCurve* discCurve= GetMktDataMng()->GetZeroCurve(
		cds.GetCurrencyUnit()->GetCcyName(), GetAsOfDate() );
	
	long backwardFrom(itsTree.depth()-1); 
	long backwardTo(0); 
	ICM_QMatrix<double> statesFrom(itsTree.slice(backwardFrom).size(),1); 
	ICM_QMatrix<double> ret2 ;
	ICM_QMatrix<double> ret2_greeks ;

	TreeServices::backwardprice(itsTree,
								backwardFrom,
								backwardTo,
								flows,
								*discCurve,
								statesFrom,
								ret2,
								ret2_greeks);

	ICM_QMatrix<double> statesFrom2(itsTree.slice(backwardFrom).size(),1); 
	ICM_QMatrix<double> ret3 ;
	double ret ;

	//FeeLeg Accrued Flows

	flows.Resize(backwardFrom+1,2);
	
	for(i=0;i<=backwardFrom;i++)
	{
		if ( gt(itsTree.Time(i),CountYears(KACTUAL_365,GetAsOfDate(),ARM_Date(startDates.Elt(first)))))
		{
			int l  = 0, i0 = 0.;
			
			while (l<size && lt(CountYears(KACTUAL_365,GetAsOfDate(),ARM_Date(startDates.Elt(first+l))),itsTree.Time(i)))
				l++;
			
			double EffectiveDate = CountYears(KACTUAL_365,GetAsOfDate(),ARM_Date(startDates.Elt(first+l-1))) ; // la bonne start date
			double PayDate = CountYears(KACTUAL_365,GetAsOfDate(),ARM_Date(payDates.Elt(first+l-1))) ;
			
			if(fabs(EffectiveDate-itsTree.Time(i))*365.<=4.00001 && l>1) // le 4.00001 est à revoir!!!!
			{
				EffectiveDate = CountYears(KACTUAL_365,GetAsOfDate(),ARM_Date(startDates.Elt(first+l-2))) ; // la bonne start date
				PayDate = CountYears(KACTUAL_365,GetAsOfDate(),ARM_Date(payDates.Elt(first+l-2))) ; // La vraie PayDate
			}

			while (i0<=backwardFrom && gt(PayDate,itsTree.Time(i0)))
				i0++;

			if (i0-1>= 0)
			{
				if ( fabs(PayDate - itsTree.Time(i0-1))<= fabs(PayDate - itsTree.Time(i0)))
					i0 = i0-1 ;
			}

			double aux = (double) ((double) GetAsOfDate().GetJulian() + ((double)EffectiveDate)*365.) ;
			ARM_Date AccruedStartDate =  ARM_Date(aux);
			// we can't apply + to ARM_Date... 
			// ARM_Date AccruedEndDate = ARM_Date(GetAsOfDate() + ((double)(std::_cpp_min(itsTree.Time(i),PayDate)))*365.) ;
			ARM_Date AccruedEndDate( GetAsOfDate()) ; 
			AccruedEndDate.AddDays(std::_cpp_min(itsTree.Time(i),PayDate)*365.) ;
			
			
			flows(i,0) = itsTree.Time(i0) ; // Le noeud qui se rapproche le plus de la vraie PayDate
			
			if (fabs(itsTree.Time(i)-PayDate)*365. <= 2.) // le 2 est à revoir!!
				flows(i,1) = 0. ;
			else
				flows(i,1)= sec.DiffAccrued(AccruedStartDate, AccruedEndDate);
		}
		else
		{
			flows(i,0) = 0. ;
			flows(i,1) = 0. ;
		}
	}


	ICM_QMatrix<double> ret3_greeks ;

	TreeServices::backwardprice_default3(
										itsTree,
										backwardFrom,
										backwardTo,
										flows,
										*discCurve,
										statesFrom2,
										ret3,
										ret3_greeks); 

	SetFeeLegPrice(ret = ret2(0,0)+ret3(0,0)); 

	SetFeeLegPrice(ret);

	return ret; 
}


// ------------------------------------------------------------
// Compute Sensitivity Method
// ------------------------------------------------------------
//	virtual protected
double 
ICM_Pricer_Tree_Binomial_Cds::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
												 const std::string& plot,
												 const std::string& label,
												 double epsvalue , double epsilonGamma // useless
												 )
{
    double sensitivity =0.;

	ICM_MktDataMng* MktDataMng = (ICM_MktDataMng*) GetMktDataMng();
		
	ARM_Date ExecutionDate = GetModel()->GetStartDate();

	double result = 0.0;
	double initialprice = 0.0;
	double modifiedprice = 0.0;

     
    {

		if (GetInitialPriceFlg())
			initialprice = GetInitialPrice();
		else
			initialprice = Price(qCMPPRICE);

		StdBinomialTree OldTree = itsTree ;

		//ResetPricer();
		ResetRootPricer() ;

		switch (typesensi)
		{
			case ICMRECOVERY_TYPE :
			case ICMIRCURVE_TYPE :
			case ICMSPREAD_TYPE :
			case ICMSPRELSHIFT_TYPE :
			case ICM_DTR_TYPE :
			{
				vector<string>* Tenors = new vector<string>;

				// if (strcmp(plot,"NONE"))
				if (plot!="NONE")
				{
					Tenors->resize(1) ;
					(*Tenors)[0] = plot;
				}
				
				ICM_BASIS_BUMP Parameters(*Tenors,label,epsvalue) ;

				ICM_BUMP<ICM_BASIS_BUMP>* bump = new ICM_BUMP<ICM_BASIS_BUMP>(typesensi, epsvalue) ;
				bump->Push(Parameters);

				ICM_MktDataMng* NewMktDataMng = MktDataMng->GenerateShiftMng(*bump);

				SetMktDataMng(NewMktDataMng);

				double vol; 
				getParam(ICM_Pricer_Tree_Binomial_Cds::PARAM_FORCEVOL,vol); 
				
				ICM_Cds& cds=dynamic_cast<ICM_Cds&>(*GetSecurity());
				ICM_DefaultCurve* defCurve = NewMktDataMng->GetDefaultCurve(cds.GetSingleName(),GetAsOfDate()); 

				TreeServices::calibrate2(itsTree,
										*defCurve,
										vol);

				modifiedprice = ComputePrice(qCMPPRICE);

				if (NewMktDataMng)
					delete NewMktDataMng;
				NewMktDataMng = NULL;

				if (Tenors)
					delete Tenors;
				Tenors = NULL;

				if (bump)
					delete bump;
				bump = NULL;
				
				SetMktDataMng(MktDataMng); //On reset le model initial

				itsTree=OldTree; 

			}
			break;
			default :
			result = -99999999.0;
		}

	ResetRootPricer();
	ComputePrice(qCMPPRICE);

	if (!result)
		result=modifiedprice - initialprice;


	}
   

	return (result);

}

double ICM_Pricer_Tree_Binomial_Cds::ComputeDuration(void) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
double ICM_Pricer_Tree_Binomial_Cds::ComputeImpliedVol(const double& Price) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
double ICM_Pricer_Tree_Binomial_Cds::Compute_Fwd_Spread(const ARM_Date &Mty, const ARM_Date & CDS_ExpiryDate, double& dur)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
void ICM_Pricer_Tree_Binomial_Cds::SetModel(ARM_Model* model)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
void ICM_Pricer_Tree_Binomial_Cds::computelossunit() 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
// 17783 void ICM_Pricer_Tree_Binomial_Cds::PerturbDefaultCurves() 
// 17783 	{
// 17783 		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
// 17783 	}


