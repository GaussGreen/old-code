#include "ARMKernel/glob/firsttoinc.h" 
#include "ICMKernel\pricer\ICM_Pricer_Tree_Binomial_Cds.h"
#include "ICMKernel\inst\icm_cds.h"
#include "ICMKernel/inst/icm_spreadoption.h"
#include "ICMKernel/mod/icm_fwdzerocurve.h"
#include "icm_pricer_tree_binomial_spreadoption.h"
#include "ICMKernel\pricer\icm_pricer_cds.h"

#include "ICMKernel/glob/icm_mktdatamng.h"
#include "ICMKernel/mod/icm_defcurvemodel.h"


#include <set>
#include <algorithm>



//	------------------------------------------------------------------------
std::string ICM_Pricer_Tree_Binomial_SpreadOption::PARAM_FORCEVOL("forceVol"); 
std::string ICM_Pricer_Tree_Binomial_SpreadOption::PARAM_OPTIONDISCSTEPSIZE("discStepSize"); 
std::string ICM_Pricer_Tree_Binomial_SpreadOption::PARAM_OPTIONDISCSTEPNUMBER("discStepNumber"); 
std::string ICM_Pricer_Tree_Binomial_SpreadOption::PARAM_USECONSTANTPROBA("useConstantProba"); 
std::string ICM_Pricer_Tree_Binomial_SpreadOption::PARAM_DISCTOLERANCE("daysDiscTolerance"); 
std::string ICM_Pricer_Tree_Binomial_SpreadOption::PARAM_DISCMAXITER("discMaxIter"); 
std::string ICM_Pricer_Tree_Binomial_SpreadOption::PARAM_DISCVARIANCETOL("discVarianceTol"); 
//	------------------------------------------------------------------------
ICM_Pricer_Tree_Binomial_SpreadOption::ICM_Pricer_Tree_Binomial_SpreadOption() 
{ 
	Init();	
}
//	------------------------------------------------------------------------
ICM_Pricer_Tree_Binomial_SpreadOption::ICM_Pricer_Tree_Binomial_SpreadOption(ARM_Security *sec, ARM_Object *mod,const ICM_Parameters& parameters,const ARM_Date&asof)
{
	Init();

	itsDeltaFlg = false;
	itsGammaFlg = false;
	itsVegaFlg = false ;
	itsThetaFlg = false ;
	itsRhoFlg  = false ;

	Set(sec, mod,parameters,asof);
}

//	------------------------------------------------------------------------
//	virtual 
ICM_Pricer_Tree_Binomial_SpreadOption::~ICM_Pricer_Tree_Binomial_SpreadOption()
{}
//	------------------------------------------------------------------------
//virtual 
ARM_CLASS_NAME 
ICM_Pricer_Tree_Binomial_SpreadOption::GetRootName()
{
	return ICM_PRICER; 
}
//	------------------------------------------------------------------------
void ICM_Pricer_Tree_Binomial_SpreadOption::Init()
{
	SetName(ICM_BINOMIAL_TREE_PRICER_SPREADOPTION);
	itsTree.clear(); 
	itsCdsReg.clear(); 
}
//	------------------------------------------------------------------------
//	virtual 
void 
ICM_Pricer_Tree_Binomial_SpreadOption::Reset(void)
{
	parent::Reset();
	itsTree.clear(); 
	itsCdsReg.clear() ;
	
	itsDeltaFlg = false;
	itsGammaFlg = false;
	itsVegaFlg = false ;
	itsThetaFlg = false ;
	itsRhoFlg  = false ;
}
//	------------------------------------------------------------------------
void 
ICM_Pricer_Tree_Binomial_SpreadOption::Set(ARM_Security *_sec, ARM_Object *mod,const ICM_Parameters&parameters,const ARM_Date&asof)
{
	if (!_sec) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Tree_Binomial_SpreadOption::Set: no security"); 
	ICM_SpreadOption& item=dynamic_cast<ICM_SpreadOption&>(*_sec); 

	//	--	Restrictions : the pricer will fail if not configured 
	//		with an european option KO not accelerated

	
	// if (item.IsAccelerated()) 
	// 	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Tree_Binomial_SpreadOption::Set: can't price ACCELERATED product"); 
	// if (!item.IsKO()) 
	// 	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Tree_Binomial_SpreadOption::Set: can't price Non KO product"); 
	// if (item.getExerciseStyle()!=K_EUROPEAN) 
	// 	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Tree_Binomial_SpreadOption::Set: can't price Non EUROPEAN product"); 
	if (item.getExerciseDates().GetNumLines()!=1) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Tree_Binomial_SpreadOption::Set: exerciseDate should be 1 line length"); 
	if (item.getExpiryDate()<asof) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Tree_Binomial_SpreadOption::Set: can't price after expiryDate"); 

	parent::Set(&item,mod,parameters,&asof);
	itsCdsReg.clear(); 

	// We construct the pricing tree. 
	ARM_Date asofDate = GetAsOfDate(); 

	double vol; 
	getParam(ICM_Pricer_Tree_Binomial_SpreadOption::PARAM_FORCEVOL,vol); 

	
	// 

	const ARM_Vector& undPayDates = item.getSecurity().GetFeeLeg().GetCreditInfosRef().GetPayDates(); 
	unsigned long first; 
	unsigned long size=item.getSecurity().GetFeeLeg().GetCreditInfosRef().countPaymentFlows(first,asofDate); 
	
	//	Create the SET with appropriate comparison scheme.
	double yfTol = std::numeric_limits<double>::epsilon(); 
	long daysTol=0; 
	getParam(PARAM_DISCTOLERANCE,daysTol,false); 
	if (daysTol!=0) yfTol=fabs(daysTol/365.); 
	tol_comp myComp(yfTol) ;
	std::set<double,tol_comp > times ( myComp ) ; 
	
	// std::set<double, lt_comp<double> > times; 
	//	first step - include in the discretisation all the payment dates
	times.insert(0);  
	double yfEndCDS;

	//for(unsigned int i=0;i<size;i++) 
	//	times.insert( yfEndCDS=CountYears(KACTUAL_365,asofDate,ARM_Date(undPayDates.Elt(first+i))) ) ;

	item.getSecurity().GetFeeLeg().ComputeYF(asofDate) ;
	const ARM_Vector& YFAccPayDates   = item.getSecurity().GetFeeLeg().GetCreditInfos()->GetYFPayDates();
	
	//for(unsigned int i=0;i<size;i++) 
	//	times.insert( yfEndCDS=CountYears(KACTUAL_365,asofDate,ARM_Date(undPayDates.Elt(first+i))) ) ;

	for(unsigned int i=0;i<size;i++) 
	 	times.insert(YFAccPayDates.Elt(i)) ;

	// second step - include the option maturity 
	double yfMaturity = CountYears(KACTUAL_365,asofDate,item.getExpiryDate()) ;
	times.insert(yfMaturity);

	// third step - include all the Roll dates (4 std maturity dates for cds per year)
	//yfEndCDS = CountYears(KACTUAL_365,asofDate,ARM_Date(undPayDates.Elt(size-1))) ;
	yfEndCDS = YFAccPayDates.Elt(size-1) ;

	ARM_Date RollDate ;
	int nb = 0 ;
	
	char TEMP[4];
	sprintf(TEMP,"%iM",nb);
	RollDate = AddPeriod(asofDate,TEMP, string(_sec->GetCurrencyUnit()->GetCcyName()),/*ARM_DEFAULT_COUNTRY*/true,qCredit_Adjust20);
	double yfRollDate = CountYears(KACTUAL_365,asofDate,RollDate) ;
	
	while (yfRollDate <= yfEndCDS)
	{
		times.insert(yfRollDate);
		nb++;
		sprintf(TEMP,"%iM",nb*3);
		RollDate = AddPeriod(asofDate,TEMP,string(_sec->GetCurrencyUnit()->GetCcyName()),/*ARM_DEFAULT_COUNTRY*/true,qCredit_Adjust20);
		ARM_Date aux = RollDate ;
		if(aux.NextBusinessDay(1).PreviousBusinessDay(1) != RollDate)
			RollDate.NextBusinessDay(1) ;		
		yfRollDate = CountYears(KACTUAL_365,asofDate,RollDate) ;
	}
	
	//	Adapt the discretization so that variance changes are minimal
	// 
	double  userStepSize ;
	long	userStepNumber ;
	double yfStep ;
	if (!getParam(PARAM_OPTIONDISCSTEPNUMBER,userStepNumber ,false)) userStepNumber=0; 
	if (!getParam(PARAM_OPTIONDISCSTEPSIZE,userStepSize ,false)) userStepSize=0; 
	yfStep = yfEndCDS; 
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
	double discVarianceTol (yfStep*vol*vol+1) ;
	long discMaxIter(0); 
	getParam(PARAM_DISCMAXITER,discMaxIter,false) ;
	getParam(PARAM_DISCVARIANCETOL,discVarianceTol,false) ;

	//	now we are looking for the it that will minimize the local variance, 
	//	here at constant volatility this is the one minimizing 
	long currentIter(0); 
	double varMin ; 
	double varMax ; 
	double varCriteria=1 ; 
	std::set<double,tol_comp> adaptativeTimes(times); 
	do 
	{
		if (currentIter!=0) yfStep = minDT/(int(minDT/yfStep) +1.); 
		//ICMLOG("yfStep="<<yfStep<<" at iteration "<<currentIter) ;
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
		adaptativeTimes = times; 
		it=times.begin(); next_it=it; 
		++next_it ;
		varMin = 100.; 
		varMax = 0 ; 
		while (next_it != times.end())
		{
			double DT = *next_it - *it; 
			double N = vol*vol*DT / targetVariance ;
			int Nint = N +0.5 ;		// this is the rouding to closest integer
			// int Nint = N   ;		// this is the rouding to closest integer
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
		varCriteria=fabs(varMax-varMin); 
		//ICMLOG("yfStep="<<yfStep<<" at iteration "<<currentIter<<" varCriteria="<<varCriteria<<" varMin="<<varMin<<" varMax="<<varMax) ;
		currentIter++; 
	} 
	while (currentIter<discMaxIter && gt(varCriteria,discVarianceTol) ) ;
	// ICMLOG("Final yfStep="<<yfStep<<" varMin="<<varMin<<" varMax="<<varMax) ;

	ICM_QMatrix<double> theTimes(adaptativeTimes.size(),1); 
	it =adaptativeTimes.begin(); 
	i=0; 
	while (it!=adaptativeTimes.end()) 
	{ 
		// ICMLOG(" Adding time yf="<<*it); 
		theTimes(i,0)=*it; 
		++it; 
		++i; 
	}

	
	itsTree.setTimes(theTimes); 

	long useConstantProba(0); 
	if (!getParam(PARAM_USECONSTANTPROBA,useConstantProba,false)) useConstantProba=0; 
	if (useConstantProba) TreeServices::calibrate(itsTree,getDefCurve(), vol); 
	else TreeServices::calibrate2(itsTree,getDefCurve(),vol); 
		
	// 
	// itsUseOldDefLegPriging=0; 
	// getParam(PARAM_USEOLDDEFLEGPRICING,itsUseOldDefLegPriging,false) ; 
}

static std::string getProductTypeKey(bool isAcc,bool isKo,int exStyle,qUnderlying_Maturity_Style undMatuStyle)
{
	std::stringstream sstr ;
	sstr<<"["<<ICM_SpreadOption::TYPE()<<":ICM_SpreadOption]" ;
	sstr<<"["<<ICM_SpreadOption::ACCELERATED()<<":"<<isAcc<<"]" ;
	sstr<<"["<<ICM_SpreadOption::EXSTYLE()<<":"<<exStyle<<"]" ;
	sstr<<"["<<ICM_SpreadOption::UNDMATUSTYLE()<<":"<<ICM_EnumsCnv::toString(undMatuStyle)<<"]" ;
	sstr<<"["<<ICM_SpreadOption::KO()<<":"<<isKo<<"]" ;
	return sstr.str(); 
}
//	------------------------------------------------------------------------
// virtual 
double 
ICM_Pricer_Tree_Binomial_SpreadOption::ComputePrice(qCMPMETH measure )
{
	//	
	//		Temporary:  
	//			we register here the pricing function according to "extended product type" 
	//			delivererd by the product . 
	// 
	typedef double (ICM_Pricer_Tree_Binomial_SpreadOption::*calcfunc_t )() ;
	static std::map<std::string,calcfunc_t> methodRegistry ;
	if (methodRegistry.size()==0) 
	{
		ICMLOG("ICM_Pricer_Tree_Binomial_SpreadOption::ComputePrice: setting up the registry"); 
		//	register the EUR pricing cases
		methodRegistry[getProductTypeKey(false/*acc*/,true/*ko*/,K_EUROPEAN,qConstantMaturity)]= &ICM_Pricer_Tree_Binomial_SpreadOption::EuropeanPrice ; 
		methodRegistry[getProductTypeKey(true/*acc*/,false/*ko*/,K_EUROPEAN,qConstantMaturity)]= &ICM_Pricer_Tree_Binomial_SpreadOption::EuropeanPrice ; 
		methodRegistry[getProductTypeKey(false/*acc*/,false/*ko*/,K_EUROPEAN,qConstantMaturity)]= &ICM_Pricer_Tree_Binomial_SpreadOption::EuropeanPrice ; 
		//	register the BER pricing cases
		methodRegistry[getProductTypeKey(false/*acc*/,false/*ko*/,K_BERMUDAN,qConstantMaturity)]= &ICM_Pricer_Tree_Binomial_SpreadOption::AmericanPrice; 
		methodRegistry[getProductTypeKey(true/*acc*/,false/*ko*/,K_BERMUDAN,qConstantMaturity)]= &ICM_Pricer_Tree_Binomial_SpreadOption::AmericanPrice ; 
		methodRegistry[getProductTypeKey(false/*acc*/,true/*ko*/,K_BERMUDAN,qConstantMaturity)]= &ICM_Pricer_Tree_Binomial_SpreadOption::AmericanPrice ; // for puts
		methodRegistry[getProductTypeKey(false/*acc*/,false/*ko*/,K_BERMUDAN,qResidualMaturity)]= &ICM_Pricer_Tree_Binomial_SpreadOption::AmericanPrice ; 
		methodRegistry[getProductTypeKey(true/*acc*/,false/*ko*/,K_BERMUDAN,qResidualMaturity)]= &ICM_Pricer_Tree_Binomial_SpreadOption::AmericanPrice ; 
		methodRegistry[getProductTypeKey(false/*acc*/,true/*ko*/,K_BERMUDAN,qResidualMaturity)]= &ICM_Pricer_Tree_Binomial_SpreadOption::AmericanPrice ; 
		// register the AME pricing cases
		methodRegistry[getProductTypeKey(true/*acc*/,false/*ko*/,K_AMERICAN,qConstantMaturity)]= &ICM_Pricer_Tree_Binomial_SpreadOption::AmericanPrice; 
		methodRegistry[getProductTypeKey(false/*acc*/,true/*ko*/,K_AMERICAN,qConstantMaturity)]= &ICM_Pricer_Tree_Binomial_SpreadOption::AmericanPrice;  // for puts
		methodRegistry[getProductTypeKey(true/*acc*/,false/*ko*/,K_AMERICAN,qResidualMaturity)]= &ICM_Pricer_Tree_Binomial_SpreadOption::AmericanPrice ;
		methodRegistry[getProductTypeKey(false/*acc*/,true/*ko*/,K_AMERICAN,qResidualMaturity)]= &ICM_Pricer_Tree_Binomial_SpreadOption::AmericanPrice ; // for puts	
		ICMLOG("ICM_Pricer_Tree_Binomial_SpreadOption::ComputePrice: done."); 
	}
	//
	//
	//
	// temoprary here... 
	switch (measure) 
	{
		case qCMPPRICE:
		case qCMPPREMIUM:
		{
			std::string prdtype = getSpreadOption().getProductTypeKey(); 
			ICMLOG("ICM_Pricer_Tree_Binomial_SpreadOption::ComputePrice: now pricing "<<prdtype); 
			std::map<std::string,calcfunc_t>::iterator it = methodRegistry.find(prdtype); 
			if (it==methodRegistry.end()) 
				ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Tree_Binomial_SpreadOption::ComputePrice: unsupported "<<prdtype); 
			calcfunc_t  calcFunction = it->second ; 
			if(!GetPriceFlg())	(this->*calcFunction)() ; 
			return GetPrice() ;
		}
		break;
		case qCMPFWDSPREAD:
		{
			if(!GetSpreadFlg())	this->FwdSpreadPrice() ; 
			return   GetSpread() ;
		}
		break ;
	default:
		ICMTHROW(ERR_INVALID_ARGUMENT,"Can't price for pricing mode "<<measure); 
	} ;
}

//	------------------------------------------------------------------------
//	virtual 
void 
ICM_Pricer_Tree_Binomial_SpreadOption::View(char* id, FILE* ficOut)
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

	fprintf(fOut, "\t\t\t ----------------- Binomial Tree SpreadOption Pricer ----------------- \n\n");
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


double 
ICM_Pricer_Tree_Binomial_SpreadOption::GenericLegPrice(const int& mode,double* Data)
{
	double Result = 0. ;

	ICM_SpreadOption*  item=dynamic_cast<ICM_SpreadOption*>(GetSecurity()); 
	if (!item)
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Tree_Binomial_SpreadOption: No Spreadoption defined. "); 

	//	Step 1	- We price PL and DL of the underlying CDS
	//
	//		This is done from CDS last payment date to optionMaturity date
	//		With a forward discount curve
	//		We get for DL and PL a vector of output (one value per tree node)
	//	
	unsigned long posMaturity(0); 
	if (!itsTree.findpos( CountYears(KACTUAL_365,GetAsOfDate(),item->getExpiryDate()) ,posMaturity ) ) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Cant find posMaturity in the tree"); 
	
	//	--	Retrieve the fee leg

	const ARM_Vector& payDates=   item->getSecurity().GetFeeLeg().GetCreditInfos()->GetPayDates(); 
	const ARM_Vector& startDates=  item->getSecurity().GetFeeLeg().GetCreditInfos()->GetAccStartDates(); 
	// 14514 double notional = item->getSecurity().GetDefLeg().GetCreditInfosRef().GetInitialNotional() ;
	double notional = item->getSecurity().GetDefLeg().GetCreditInfosRef().GetNotionals().Elt(0) ;
	double losses = ( 1.-GetMktDataMng()->GetDefaultCurve(item->getSecurity().GetSingleName(),GetAsOfDate())->GetRecovery()) * notional ; 
	
	unsigned long posEndCDS(0); 
	if (!itsTree.findpos(CountYears(KACTUAL_365,GetAsOfDate(), payDates.Elt(payDates.GetNumLines()-1)),posEndCDS)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Cant find posEndCDS in the tree"); 
	//	--	
	unsigned long first,size; 
	size=  item->getSecurity().GetFeeLeg().GetCreditInfos()->countPaymentFlows(first,GetAsOfDate());
	ICM_QMatrix<double> flows(size,2); 
	
	const ARM_Vector& YFAccstartDates = item->getSecurity().GetFeeLeg().GetCreditInfos()->GetYFAccStartDates();
	const ARM_Vector& YFAccpayDates   = item->getSecurity().GetFeeLeg().GetCreditInfos()->GetYFPayDates();

	//FeeLeg Flows

	for(unsigned int i=0;i<size;i++)
	{
		flows(i,0)= YFAccpayDates.Elt(i);
		//flows(i,1)= item->getSecurity().GetFeeLeg().GetCreditInfos()->FullCoupon(startDates.Elt(first+i)); 
		flows(i,1)= item->getSecurity().GetFeeLeg().GetCreditInfos()->FullCoupon(first+i); 
	}

	ARM_ZeroCurve* discCurve= GetMktDataMng()->GetZeroCurve(
		item->GetCurrencyUnit()->GetCcyName(), GetAsOfDate() );

	ICM_QMatrix<double> statesFrom(itsTree.slice(posEndCDS).size(),1) ;
	ICM_QMatrix<double> feelegFwd,feelegFwd0 ;
	ICM_QMatrix<double> ret_greeks ; 
	
	ICM_ZeroCurveForward fwdCurve(*discCurve,item->getExpiryDate());
	long backwardFrom(itsTree.depth()-1);
	
	long backwardTo(0);
	unsigned int j = 0;
	double aux = 0. ;
	
	ICM_QMatrix<double> statesFrom_def(itsTree.slice(backwardFrom).size(),1);

	TreeServices::backwardprice(itsTree,posEndCDS,posMaturity,flows,*discCurve,statesFrom_def,feelegFwd,ret_greeks); 

	flows.Resize(backwardFrom+1,2);
	ICM_QMatrix<double> flows_def(backwardFrom,2);
	ICM_QMatrix<double> indexes(2,1);
	
	for(i=0;i<=backwardFrom;i++)
	{
		//FeeLeg Accrued Flows

		if (i>posMaturity)
		{
			double EffectiveDate = itsTree.Time(posMaturity);
			double PayDate = 0.;				

			int l  = 0, i0 = 0;;
			
			while (l<size && lt(YFAccstartDates.Elt(l),itsTree.Time(i)))
				l++;
			
			EffectiveDate = YFAccstartDates.Elt(l-1);				// la bonne start date
			
			PayDate = YFAccpayDates.Elt(l-1);
			
			if(fabs(EffectiveDate-itsTree.Time(i))*365.<=4.00001 && l>1) // le 4.00001 est à revoir!!!!
			{
				EffectiveDate = YFAccstartDates.Elt(l-2) ;			// la bonne start date
				PayDate = YFAccpayDates.Elt(l-2) ;					// La vraie PayDate
			}

			TreeServices::dichotomic_research(itsTree,PayDate ,0 , backwardFrom, indexes) ;
			i0 = indexes(1,0) ;

			if (i0-1>= 0)
			{
				if ( fabs(PayDate - itsTree.Time(i0-1))<= fabs(PayDate - itsTree.Time(i0)))
					i0 = i0-1 ;
			}

			flows(i,0) = itsTree.Time(i0) ; // Le noeud qui se rapproche le plus de la vraie PayDate

			if (fabs(itsTree.Time(i)-PayDate)*365. <= 2.) // le 2 est à revoir!!
				flows(i,1) = 0. ;
			else
			{
				aux = (double) ((double) GetAsOfDate().GetJulian() + ((double)(std::_cpp_max((double)itsTree.Time(posMaturity),(double)EffectiveDate)))*365.) ;
				
				ARM_Date AccruedStartDate =  ARM_Date(aux);
				// ARM_Date AccruedEndDate = ARM_Date(GetAsOfDate() + ((double) (std::_cpp_min(itsTree.Time(i),PayDate)))*365.) ;
				ARM_Date AccruedEndDate (GetAsOfDate() ) ;
				AccruedEndDate.AddDays( std::_cpp_min(itsTree.Time(i),PayDate)*365. )  ;
			
				flows(i,1)= item->getSecurity().GetFeeLeg().GetCreditInfos()->DiffAccrued(AccruedStartDate, AccruedEndDate);
			}
		}
		else
		{
			flows(i,0) = 0. ;
			flows(i,1) = 0. ;
		}
		
		//DefLeg flows
		
		//flows_def(i,0) représente les payment Dates du flux du Node i
		//flows_def(i,1) représente le flux au Node i, soit La Loss à cette date

		if (i<backwardFrom)
		{
			aux = YFAccstartDates.Elt(0);

			if (itsTree.Time(i)<aux)
			{
				flows_def(i,0) = 0. ;
				flows_def(i,1) = 0. ;
			}
			else
			{	
				flows_def(i,0)= 0.5 *(itsTree.Time(i)+itsTree.Time(i+1)) + item->getSecurity().GetCreditLag()/365.;
				flows_def(i,1)=losses ;
			}
		}
	}

	ICM_QMatrix<double> statesFrom2(itsTree.slice(backwardFrom).size(),1); 
	ICM_QMatrix<double> AccruedFL ;
	
	TreeServices::backwardprice_default3(itsTree, backwardFrom,posMaturity,flows,*discCurve,statesFrom2,AccruedFL,ICM_QMatrix<double>(6,1));
	
	for(i=0;i<feelegFwd.Getnbrows();i++) 
		feelegFwd(i,0) = feelegFwd(i,0) + AccruedFL(i,0);


	ICM_QMatrix<double> deflegFwd,deflegFwd0 ; 

	TreeServices::backwardprice_default3(itsTree,posEndCDS,posMaturity,flows_def,*discCurve,statesFrom_def,deflegFwd,ICM_QMatrix<double>(6,1)); 
	
	//Discountage
		
	TreeServices::backwardprice(itsTree,posMaturity,0,ICM_QMatrix<double>(1,2),*discCurve,feelegFwd,feelegFwd0,ret_greeks); 
		
	TreeServices::backwardprice(itsTree,posMaturity,0,ICM_QMatrix<double>(1,2),*discCurve,deflegFwd,deflegFwd0,ret_greeks); 

	//	Step 2	- KO European Option pricing is a backward pricing 
	//
	//		CALL	: We assign to the end value the max(fwd-K,0)*PV01fwd for each node
	//		PUT		: We assign to the end value the max(K-fwd,0)*PV01fwd for each node 
	//

	ICM_QMatrix<double> payoff(feelegFwd.Getnbrows(),1); 
	
	for(i=0;i<payoff.Getnbrows();i++) 
	{
		double spread,pv01; 
		if (eq(feelegFwd(i,0),0.)) 
		{
			spread=pv01=0; 
		} 
		else 
		{
			spread=deflegFwd(i,0)/feelegFwd(i,0)*100.*item->getSecurity().GetFeeLeg().GetCreditSpread(); 
			pv01 = 100. * feelegFwd(i,0) / (notional *item->getSecurity().GetFeeLeg().GetCreditSpread()); 
		}
// 		ICMLOG ("noeud "<<i<<" spread="<<spread<<" pv01="<<pv01); 

		payoff(i,0) =  std::_cpp_max((item->IsCall()?1:-1)*(spread-item->getStrike()),0.) * pv01 ;
	}
; 
	ICM_QMatrix<double> optionPrice; 

	TreeServices::backwardprice(itsTree,posMaturity,0,ICM_QMatrix<double>(1,2),*discCurve,payoff,optionPrice,ret_greeks); 

	Result = optionPrice(0,0) ; 

	 // If No KO

	if ((Data[0] == 0 || Data[0] == 1) )
	{

		double OptionPrice_NoKO = 0.;
		
		std::auto_ptr<ICM_Cds> Shortcds(
			(ICM_Cds*) GetMktDataMng()->GetDefaultCurve(item->getSecurity().GetSingleName(),GetAsOfDate())->stdCDS(GetAsOfDate(),item->getExpiryDate())
			) ;

		const ARM_Vector& startDates2= Shortcds->GetFeeLeg()->GetCreditInfos()->GetAccStartDates();

		unsigned long first2,size2; 
		
		Shortcds->GetFeeLeg()->ComputeYF(GetAsOfDate()) ;
		const ARM_Vector& YFstartDates2 = Shortcds->GetFeeLeg()->GetCreditInfos()->GetYFAccStartDates();

		size2 = Shortcds->GetFeeLeg()->GetCreditInfos()->countPaymentFlows(first2,GetAsOfDate());

		//ARM_Date EffectiveDate = startDates2.Elt(first);

		ICM_QMatrix<double> flows_def2(backwardFrom,2);
		
		
		for(i=0;i<backwardFrom;i++)
		{
			aux = YFstartDates2.Elt(0) ;
			
			if (itsTree.Time(i)<aux)
			{
				flows_def2(i,0) = 0.;
				flows_def2(i,1) = 0.;
			}
			else
			{
				if (Data[0] == 0) // No KO + with Acceleration
					flows_def2(i,0)= itsTree.Time(i+1);
				if (Data[0] == 1) // No KO + without Acceleration
					flows_def2(i,0)= itsTree.Time(posMaturity);
				
				if(item->IsCall())
					flows_def2(i,1)=losses ;
			}
		}

		
		TreeServices::backwardprice_default3(itsTree,posMaturity,0,flows_def2,*discCurve,statesFrom_def,deflegFwd,ICM_QMatrix<double>(6,1));

		OptionPrice_NoKO = optionPrice(0,0) + deflegFwd(0,0)*10000./notional ;

		Result = OptionPrice_NoKO ;
	}

	// End of No KO additional value Pricing

	if(Data[0] == -2)
	{
		double spreadFwd; 
		
		if (eq(feelegFwd0(0,0),0.))
			spreadFwd=0; 
		else
			spreadFwd=(deflegFwd0(0,0)/feelegFwd0(0,0))*100.*item->getSecurity().GetFeeLeg().GetCreditSpread();
		
		Result = spreadFwd ;
	}
	

	ICM_QMatrix<double> Times_AM(posMaturity+1,1,0.);
	ICM_QMatrix<double> flowsDef_AM(posMaturity+1,2,0.);
	StdBinomialTree flows_AM ;
	ICM_QMatrix<double> RET(posMaturity+1,posMaturity+1);
	ICM_QMatrix<double> statesFrom_AM(itsTree.slice(posMaturity).size(),1) ; 
	ICM_QMatrix<double> optionPrice2;
	

	if(Data[0] == 2) // American Options
	{
		CpteSpreadsAndPv01s(RET, 1.); // Constant Maturity Spreads and PV01s

		//Option Américaine (forcément avec KO et avec Accélération)
		
		flows_AM.setTimes(itsTree.Time());

		//flows_AM := valeurs d'exercices stockées en h	

		//Times_AM(i,0) := date d'exercice juste après le noeud i
		//flowsDef_AM(i,0 := Paiement en cas de défaut au noeud i
		//flowsDef_AM(i,1) := date de paiement en cas de défaut entre i et i+1

		for(i=0;i<=posMaturity; i++) 
		{
			StdBinomialTree::slice_t aSlice=flows_AM.slice(i); 
			StdBinomialTree::slice_t::iterator it = aSlice.begin(); 
			StdBinomialTree::node_data_t& baseData= it->data(); 
			
			for (j=1;j<aSlice.size(); j++)
			{	
				double aux = RET.Getvalue(i,j);	
				aSlice.data(j).h = aux;	
				//ICMLOG("Exercice value at Node"<<i<<","<<j<<" is "<<aux) ;
			}

			Times_AM(i,0)=itsTree.Time(i);
			
			if(item->IsCall())
				flowsDef_AM(i,0)=losses * 10000./notional;
						
			// No KO + with Acceleration
			flowsDef_AM(i,1)=itsTree.Time(i+1);		
		}

		TreeServices::backwardprice_BN(itsTree,
									   posMaturity,
									   0,
									   flows_AM,
									   Times_AM,
									   flowsDef_AM,
									   *discCurve,
									   statesFrom_AM,
									   optionPrice2,
									   ICM_QMatrix<double>(6,1));
		Result = optionPrice2(0,0) ;
	}

	double ExerciceFrequency = Data[1] ;

	if(Data[0] == 3 || Data[0] == 4) // Bermuda Options
	{
		//Option Bermuda (forcément avec KO mais avec ou sans Accélération )

		CpteSpreadsAndPv01s(RET, 1.); // Constant Maturity Spreads and PV01s

		flows_AM.setTimes(itsTree.Time());

		//flows_AM := valeurs d'exercices stockées en h	

		//Times_AM(i,0) := date d'exercice juste après le noeud i
		//flowsDef_AM(i,0 := Paiement en cas de défaut au noeud i
		//flowsDef_AM(i,1) := date de paiement en cas de défaut entre i et i+1		

		for(i=0;i<=posMaturity; i++) 
		{
			StdBinomialTree::slice_t aSlice=flows_AM.slice(i); 
			StdBinomialTree::slice_t::iterator it = aSlice.begin(); 
			StdBinomialTree::node_data_t& baseData= it->data(); 
			
			for (j=1;j<aSlice.size(); j++)
			{	
				double aux = RET.Getvalue(i,j);	
				aSlice.data(j).h = aux;	
				//ICMLOG("Exercice value at Node"<<i<<","<<j<<" is "<<aux) ;
			}

			j = 1;

			if(i==0)
				j = 1;
			else
			{
				while(itsTree.Time(i-1)>j*(1./ExerciceFrequency))
					j++;
			}

			Times_AM(i,0)= (1./ExerciceFrequency) *j;
			
			if(item->IsCall())
				flowsDef_AM(i,0)=losses * 10000./notional ;
			if (Data[0] == 3) // No KO + with Acceleration
				flowsDef_AM(i,1)=itsTree.Time(i+1);
			if (Data[0] == 4) // No KO + without Acceleration
				flowsDef_AM(i,1) = (1./ExerciceFrequency) *j;


		}

		TreeServices::backwardprice_BN(itsTree,
									   posMaturity,
									   0,
									   flows_AM,
									   Times_AM,
									   flowsDef_AM,
									   *discCurve,
									   statesFrom_AM,
									   optionPrice2,
									   ICM_QMatrix<double>(6,1));
		Result = optionPrice2(0,0) ;
	}

	if(Data[0] == 5) // American Option on residual CDS
	{
		CpteSpreadsAndPv01s(RET, 0.); // Residual Maturity Spreads and PV01s
		flows_AM.setTimes(itsTree.Time());

		//flows_AM := valeurs d'exercices stockées en h	

		//Times_AM(i,0) := date d'exercice juste après le noeud i
		//flowsDef_AM(i,0 := Paiement en cas de défaut au noeud i
		//flowsDef_AM(i,1) := date de paiement en cas de défaut entre i et i+1

		for(i=0;i<=posMaturity; i++) 
		{
			StdBinomialTree::slice_t aSlice=flows_AM.slice(i); 
			StdBinomialTree::slice_t::iterator it = aSlice.begin(); 
			StdBinomialTree::node_data_t& baseData= it->data(); 
			
			for (j=1;j<aSlice.size(); j++)
			{	
				double aux = RET.Getvalue(i,j);	
				aSlice.data(j).h = aux;	
				//ICMLOG("Exercice value at Node"<<i<<","<<j<<" is "<<aux) ;
			}

			Times_AM(i,0)=itsTree.Time(i);
			
			if(item->IsCall())
				flowsDef_AM(i,0)=losses * 10000./notional ;
			
			// No KO + with Acceleration
			flowsDef_AM(i,1)=itsTree.Time(i+1);
		}
	
		TreeServices::backwardprice_BN(itsTree,
									   posMaturity,
									   0,
									   flows_AM,
									   Times_AM,
									   flowsDef_AM,
									   *discCurve,
									   statesFrom_AM,
									   optionPrice2,
									   ICM_QMatrix<double>(6,1));
		Result = optionPrice2(0,0) ;
	}

	if(Data[0] == 6 ||Data[0] == 7) // Bermuda Options on Residual CDS
	{
		//Option Bermuda (forcément avec KO mais avec ou sans Accélération )

		CpteSpreadsAndPv01s(RET, 0.); // Residual Maturity Spreads and PV01s

		flows_AM.setTimes(itsTree.Time());

		//flows_AM := valeurs d'exercices stockées en h	

		//Times_AM(i,0) := date d'exercice juste après le noeud i
		//flowsDef_AM(i,0 := Paiement en cas de défaut au noeud i
		//flowsDef_AM(i,1) := date de paiement en cas de défaut entre i et i+1

		for(i=0;i<=posMaturity; i++) 
		{
			StdBinomialTree::slice_t aSlice=flows_AM.slice(i); 
			StdBinomialTree::slice_t::iterator it = aSlice.begin(); 
			StdBinomialTree::node_data_t& baseData= it->data(); 
			
			for (j=1;j<aSlice.size(); j++)
			{	
				double aux = RET.Getvalue(i,j);	
				aSlice.data(j).h = aux;	
				//ICMLOG("Exercice value at Node"<<i<<","<<j<<" is "<<aux) ;
			}

			j = 1;

			while(itsTree.Time(i)>j*(1./ExerciceFrequency))
				j++;

			Times_AM(i,0)= (1./ExerciceFrequency) *j;
			
			if(item->IsCall())
				flowsDef_AM(i,0)=losses * 10000./notional ;
			if (Data[0] == 6) // No KO + with Acceleration
				flowsDef_AM(i,1)=itsTree.Time(i+1);
			if (Data[0] == 7) // No KO + without Acceleration
				flowsDef_AM(i,1) = (1./ExerciceFrequency) *j;
		}

		TreeServices::backwardprice_BN(itsTree,
									   posMaturity,
									   0,
									   flows_AM,
									   Times_AM,
									   flowsDef_AM,
									   *discCurve,
									   statesFrom_AM,
									   optionPrice2,
									   ICM_QMatrix<double>(6,1));
		Result = optionPrice2(0,0) ;
	}
	
	return Result ;
}

void
ICM_Pricer_Tree_Binomial_SpreadOption::CpteSpreadsAndPv01s(ICM_QMatrix<double>&payoff, double mode)
{
	ICM_SpreadOption*  item=dynamic_cast<ICM_SpreadOption*>(GetSecurity()); 
	
	if (!item)
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Tree_Binomial_SpreadOption: No Spreadoption defined. "); 
	
	unsigned long SizeTree = itsTree.depth();

	ICM_Cds Sec = item->getSecurity();	

	const ARM_Vector& YFAccStartDates = item->getSecurity().GetFeeLeg().GetCreditInfos()->GetYFAccStartDates();
	const ARM_Vector& YFAccPayDates   = item->getSecurity().GetFeeLeg().GetCreditInfos()->GetYFPayDates();

	const ARM_Vector& PayDates=   Sec.GetFeeLeg()->GetCreditInfos()->GetPayDates(); 
	
	unsigned long PosMaturity(0); 
		if (!itsTree.findpos( CountYears(KACTUAL_365,GetAsOfDate(),item->getExpiryDate()) ,PosMaturity ) ) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"Cant find posMaturity in the tree"); 

	unsigned long PosEndCDS(0); 
		if (!itsTree.findpos(CountYears(KACTUAL_365,GetAsOfDate(), PayDates.Elt(PayDates.GetNumLines()-1)),PosEndCDS)) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"Cant find posEndCDS in the tree"); 

	for (unsigned j = PosMaturity; j > 0 ; j--)
	{
		ARM_Date Maturity;

		if (mode == 1.0)
		{
			Maturity = ARM_Date(GetAsOfDate()); 
			Maturity.AddDays(itsTree.Time(j)*365).AddDays((itsTree.Time(PosEndCDS)-itsTree.Time(PosMaturity))*365);
		}
		else
		{
			Maturity = ARM_Date(PayDates.Elt(PayDates.GetNumLines()-1));
		}

		ARM_Date startDate(	GetAsOfDate() ) ;	
		startDate.AddDays(itsTree.Time(j)*365); 
		std::auto_ptr<ICM_Cds> Shortcds((ICM_Cds*) GetMktDataMng()->GetDefaultCurve(item->getSecurity().GetSingleName(),GetAsOfDate())->stdCDS(
							startDate, 
							Maturity)) ;
		
		Shortcds->SetCreditSpread(100.*Sec.GetFeeLeg()->GetCreditSpread()) ;
		Shortcds->GetFeeLeg()->GetCreditInfos()->SetAccrualBasis(Sec.GetAccrualBasis());
//		Shortcds->GetFeeLeg()->GetCreditInfos()->SetBasis(Sec.GetBasis());
								
		//	--	Retrieve the fee leg

		const ARM_Vector& payDates = Shortcds->GetFeeLeg()->GetCreditInfos()->GetPayDates(); 
		const ARM_Vector& startDates=  Shortcds->GetFeeLeg()->GetCreditInfos()->GetAccStartDates(); 

		Shortcds->GetFeeLeg()->ComputeYF(GetAsOfDate()) ;

		const ARM_Vector& YFAccstartDates =  Shortcds->GetFeeLeg()->GetCreditInfos()->GetYFAccStartDates();
		const ARM_Vector& YFAccpayDates   =  Shortcds->GetFeeLeg()->GetCreditInfos()->GetYFPayDates();

		// 14514 double notional = Sec.GetDefLeg()->GetCreditInfosRef().GetInitialNotional() ;
		double notional = Sec.GetDefLeg()->GetCreditInfosRef().GetNotionals().Elt(0) ;
		double losses = ( 1.-GetMktDataMng()->GetDefaultCurve(Sec.GetSingleName(),GetAsOfDate())->GetRecovery())* notional ; 

		unsigned long first,size; 
		long backwardFrom(itsTree.depth()-1);
		long backwardTo(0);

		size=  Shortcds->GetFeeLeg()->GetCreditInfos()->countPaymentFlows(first,GetAsOfDate());
		
		unsigned int j0 = 0;
		
		double aux = YFAccpayDates.Elt(size-1) ;

		ICM_QMatrix<double> indexes(2,1);

		TreeServices::dichotomic_research(itsTree,aux ,0 , backwardFrom, indexes) ;
		j0 = indexes(1,0) ;
				
		if(fabs(aux - itsTree.Time(j0-1))*365 <= 1.){j0--;}

		if (j0-1>= 0)
		{
			if ( fabs(aux - itsTree.Time(j0-1))<= fabs(aux - itsTree.Time(j0)))
				j0=j0-1 ;
		}
				
		unsigned long posEndCDS(j0); //slice correspondant à la End Date du CDS
		unsigned long posMaturity(j); //slice correspondant à la startDate du CDS

		ICM_QMatrix<double> flows(posEndCDS+1,2); 
		
		ICM_QMatrix<double> flows_def(backwardFrom,2);

		
		for(unsigned int i=0;i<=posEndCDS;i++)
		{
			if (i>posMaturity)
			{	

				//FeeLeg
				
				//flows(i,0) représente les payment Dates du flux du Node i
				//flows(i,1) représente le flux au Node i, soit le coupon couru à cette date
				
				// EffectiveDate la start Date juste avant Node i , il y en a forcément une
				// PayDate la payment Date juste après Node i , il y en a forcément une

				double EffectiveDate = itsTree.Time(posMaturity);				
				double PayDate = 0.;				
				
				int l  = 0, i0 = 0;;
				
				while (l<size && lt(YFAccstartDates.Elt(first+l),itsTree.Time(i)))
					l++;
				
				EffectiveDate = YFAccstartDates.Elt(first+l-1) ;		// start date just strictly before slice i

				PayDate = YFAccpayDates.Elt(first+l-1) ;
				
				if(fabs(EffectiveDate-itsTree.Time(i))*365.<=4.00001 && l>1) // le 4.00001 est à revoir!!!!
				{
					EffectiveDate = YFAccstartDates.Elt(first+l-2) ; // la bonne start date

					PayDate = YFAccpayDates.Elt(first+l-2) ; // La vraie PayDate
				}

				TreeServices::dichotomic_research(itsTree,PayDate ,0 , backwardFrom, indexes) ;
				i0 = indexes(1,0) ;

				if (i0-1>= 0)
				{
					if ( fabs(PayDate - itsTree.Time(i0-1))<= fabs(PayDate - itsTree.Time(i0)))
						i0 = i0-1 ;
				}

				aux = (double) ((double) GetAsOfDate().GetJulian() + ((double)(std::_cpp_max((double)itsTree.Time(posMaturity),(double)EffectiveDate)))*365.) ;

				ARM_Date AccruedStartDate =  ARM_Date(aux);
				ARM_Date AccruedEndDate ( GetAsOfDate() ) ; 
				AccruedEndDate .AddDays( std::_cpp_min(itsTree.Time(i),PayDate)*365.) ;
				
				flows(i,0) = itsTree.Time(i0) ; // Le noeud qui se rapproche le plus de la vraie PayDate
				flows(i,1)= Shortcds->GetFeeLeg()->GetCreditInfos()->DiffAccrued(AccruedStartDate, AccruedEndDate)* notional/100.;
			}
			else
			{
				flows(i,0) = 0. ;
				flows(i,1)= 0. ;
			}

			//DefLeg Flows
			
			//flows_def(i,0) représente les payment Dates du flux du Node i
			//flows_def(i,1) représente le flux au Node i, soit La Loss à cette date
		
			if(i<backwardFrom)
			{
				//ARM_Date EffectiveDate = startDates.Elt(first);
				//aux = CountYears(KACTUAL_365,GetAsOfDate(),EffectiveDate) ;
				
				aux = YFAccstartDates.Elt(first) ;
				
				if (itsTree.Time(i)<aux)
				{
					flows_def(i,0) = 0. ;
					flows_def(i,1) = 0. ;
				}
				else
				{
					flows_def(i,0)= 0.5 *(itsTree.Time(i)+itsTree.Time(i+1)) + item->getSecurity().GetCreditLag()/365.;
					flows_def(i,1)=losses ;
				}
			}
		}
			
		ARM_ZeroCurve* discCurve= GetMktDataMng()->GetZeroCurve(
			item->GetCurrencyUnit()->GetCcyName(), GetAsOfDate() );

		ICM_QMatrix<double> statesFrom(itsTree.slice(posEndCDS-(posMaturity-j)).size(),1) ;
		ICM_QMatrix<double> statesFrom_def(itsTree.slice(posEndCDS-(posMaturity-j)).size(),1); 
		ICM_QMatrix<double> feelegFwd,feelegFwd0 ; 
		ICM_QMatrix<double> deflegFwd,deflegFwd0 ;
		ICM_QMatrix<double> ret_greeks ;
		
		TreeServices::backwardprice(itsTree,posEndCDS-(posMaturity-j),j,flows,*discCurve,statesFrom,feelegFwd,ret_greeks);
		TreeServices::backwardprice_default3(itsTree,posEndCDS-(posMaturity-j),j,flows_def,*discCurve,statesFrom_def,deflegFwd,ICM_QMatrix<double>(6,1)); 
		
		// Partie test
		TreeServices::backwardprice(itsTree,j,0,ICM_QMatrix<double>(1,2),*discCurve,feelegFwd,feelegFwd0,ret_greeks); 
		TreeServices::backwardprice(itsTree,j,0,ICM_QMatrix<double>(1,2),*discCurve,deflegFwd,deflegFwd0,ret_greeks); 

		double spreadFwd; 
		if (eq(feelegFwd0(0,0),0.))
			spreadFwd=0; 
		else
			spreadFwd=(deflegFwd0(0,0)/feelegFwd0(0,0))*100.*item->getSecurity().GetFeeLeg().GetCreditSpread(); 
		
		ICM_DefaultCurve* DefCurve = (ICM_DefaultCurve*) GetMktDataMng()->GetDefaultCurve(item->getSecurity().GetSingleName(),GetAsOfDate());
		ARM_ZeroCurve* InitShortCurve = (ARM_ZeroCurve*) DefCurve->GetZeroCurve();
		
		ICM_DefaultCurveModel* DefModel = new ICM_DefaultCurveModel(DefCurve, InitShortCurve);

		// ICM_Pricer_Cds* PricerCds = new ICM_Pricer_Cds(Shortcds.get(),DefModel);
		ICM_Pricer_Cds  PricerCds ; PricerCds.Set(Shortcds.get(),DefModel,ICM_Parameters(),GetAsOfDate());

		PricerCds.SetFaster(true);
		ARM_Date ExecutionDate = DefModel->GetStartDate();
			
		// double DL = PricerCds.DefLegPV() ;
		// double PL = PricerCds.FeeLegPV() ;
		double DL = PricerCds.Price(qCMPDEFLEGPV) ;
		double PL = PricerCds.Price(qCMPFEELEGPV) ;
		double Spread = PricerCds.ComputeSpread(0.) ;

		// if (PricerCds)
		// 	delete PricerCds;
		// PricerCds = NULL;

		if (DefModel)
			delete DefModel;
		DefModel = NULL;

		// fin de la partie tests
		
		for(i=0;i<feelegFwd.Getnbrows();i++) 
		{
			double spread,pv01; 
			if (eq(feelegFwd(i,0),0.)) 
				spread=pv01=0; 
			else 
			{
				spread=deflegFwd(i,0)/feelegFwd(i,0)*100.*Sec.GetFeeLeg()->GetCreditSpread(); 
				pv01 = 100. * feelegFwd(i,0) / (notional * Sec.GetFeeLeg()->GetCreditSpread());
			}
			payoff(j,i) =  std::_cpp_max((item->IsCall()?1:-1)*(spread-item->getStrike()),0.) * pv01 ;
		}
	}
}

// ------------------------------------------------------------
// Compute Sensitivity Method
// ------------------------------------------------------------
// virtual protected 
double 
ICM_Pricer_Tree_Binomial_SpreadOption::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
														  const std::string& plot,
														  const std::string&label, 
														  double  epsvalue,  double epsilonGamma //useless 
														  )
{
    double sensitivity =0.;

	ICM_MktDataMng* MktDataMng = (ICM_MktDataMng*) GetMktDataMng();
		
	ARM_Date ExecutionDate = GetModel()->GetStartDate();

	double result = 0.0;
	double initialprice = 0.0;
	double modifiedprice = 0.0;
	ICM_SpreadOption*  item=dynamic_cast<ICM_SpreadOption*>(GetSecurity());


    {
		if (GetInitialPriceFlg())
			initialprice = GetInitialPrice();
		else
			initialprice = Price(qCMPPRICE);

		StdBinomialTree OldTree = itsTree ;
		double delta, gamma, vega, theta, rho, spreadfwd ;
		bool bool1 = false, bool2 = false ;

		delta = GetDelta();
		gamma = GetGamma();
		theta = GetTheta();
		spreadfwd =  GetSpread() ;

		if(GetRhoFlg())
		{
			rho = GetRho();
			bool1 = true ;
		}
		
		if(GetVegaFlg())
		{
			vega = GetVega() ;
			bool2 = true ;
		}

		ResetRootPricer() ;

		switch (typesensi)
		{
			case ICMRECOVERY_TYPE :
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
				getParam(ICM_Pricer_Tree_Binomial_SpreadOption::PARAM_FORCEVOL,vol); 
						
				ICM_DefaultCurve* defCurve = (ICM_DefaultCurve*) NewMktDataMng->GetDefaultCurve(item->getSecurity().GetSingleName(),GetAsOfDate());

				TreeServices::calibrate2(itsTree,
										*defCurve,
										vol);

				modifiedprice = ComputePrice(qCMPPRICE);

				
				// On reset les valeurs initiales

				SetMktDataMng(MktDataMng); //On reset le model initial
				itsTree=OldTree;
				SetInitialPrice(initialprice);
				SetPrice(initialprice) ;
				SetSpread(spreadfwd);
				SetDelta(delta) ;
				SetGamma(gamma);
				SetTheta(theta);
				
				if(bool2)
					SetVega(vega);

				if(bool1)
					SetRho(rho);
				
				if (NewMktDataMng)
					delete NewMktDataMng;
				NewMktDataMng = NULL;

				if (Tenors)
					delete Tenors;
				Tenors = NULL;

				if (bump)
					delete bump;
				bump = NULL;

			}
			break;
			case ICMIRCURVE_TYPE :
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
				getParam(ICM_Pricer_Tree_Binomial_SpreadOption::PARAM_FORCEVOL,vol); 
						
				ICM_DefaultCurve* defCurve = (ICM_DefaultCurve*) NewMktDataMng->GetDefaultCurve(item->getSecurity().GetSingleName(),GetAsOfDate());

				TreeServices::calibrate2(itsTree,
										*defCurve,
										vol);

				modifiedprice = ComputePrice(qCMPPRICE);
				
				SetRho(modifiedprice - initialprice) ;

				// On reset les valeurs initiales

				SetMktDataMng(MktDataMng); //On reset le model initial
				itsTree=OldTree;
				SetInitialPrice(initialprice);
				SetPrice(initialprice) ;
				SetSpread(spreadfwd);
				SetDelta(delta) ;
				SetGamma(gamma);
				SetTheta(theta);

				if(bool2)
					SetVega(vega);

				if (NewMktDataMng)
					delete NewMktDataMng;
				NewMktDataMng = NULL;

				if (Tenors)
					delete Tenors;
				Tenors = NULL;

				if (bump)
					delete bump;
				bump = NULL;

				
			}
			break;
			case ICM_GREEK_VEGA_TYPE :
			{
				double vol;
				getParam(ICM_Pricer_Tree_Binomial_SpreadOption::PARAM_FORCEVOL,vol); 
				vol += epsvalue;

				ICM_SpreadOption*  item=dynamic_cast<ICM_SpreadOption*>(GetSecurity());
				ICM_DefaultCurve* defCurve = (ICM_DefaultCurve*) MktDataMng->GetDefaultCurve(item->getSecurity().GetSingleName(),GetAsOfDate());
				
				TreeServices::calibrate2(itsTree,
										*defCurve,
										vol);

				modifiedprice = ComputePrice(qCMPPRICE);

				SetVega(modifiedprice - initialprice) ;
				
				// On reset les valeurs initiales

				SetMktDataMng(MktDataMng); //On reset le model initial
				itsTree=OldTree; 
				SetInitialPrice(initialprice);
				SetPrice(initialprice) ;
				SetSpread(spreadfwd);
				SetDelta(delta) ;
				SetGamma(gamma);
				SetTheta(theta);
				
				if(bool1)
					SetRho(rho);

			}
			break;
			case ICM_GREEK_DELTA_TYPE :
			{
				double cdsspread = 100.*item->getSecurity().GetFeeLeg().GetCreditSpread();
				// 14514 double notional = item->getSecurity().GetDefLeg().GetCreditInfosRef().GetInitialNotional() ;
				double notional = item->getSecurity().GetDefLeg().GetCreditInfosRef().GetNotionals().Elt(0) ;
				
				const ICM_Cds& cds = dynamic_cast<const ICM_Cds&>(item->getSecurity()); 
				int nb = CountYears(KACTUAL_365,GetAsOfDate(),cds.GetEndDateNA()) ; // ((ICM_Cds&)(item->getSecurity())).GetEndDate());		
				char TEMP[4];
				sprintf(TEMP,"%iY",nb);
				
				//double OptionSensi = ComputeSensitivity(ICMSPREAD_TYPE,TEMP,label, epsvalue) ;
				double OptionSensi = ComputeSensitivity(ICMSPREAD_TYPE,plot,label, epsvalue) ;

				OptionSensi *= notional / 10000.;

				ICM_MktDataMng* DefModel = GetMktDataMng();
				
				//std::auto_ptr<ICM_Cds> Shortcds(
				//					(ICM_Cds*) GetMktDataMng()->GetDefaultCurve(item->getSecurity().GetSingleName(),GetAsOfDate())->stdCDS(
				//					GetAsOfDate(),
				//					((ICM_Cds&)(item->getSecurity())).GetEndDate())) ;

				
				ARM_Date EndDate;
				
				// if (!strcmp(plot,"NONE"))  //Detection d'un parallel Shift
				if (plot=="NONE") 
					EndDate = cds.GetEndDateNA() ; // ((ICM_Cds&)(item->getSecurity())).GetEndDate() ;
				else
					EndDate = AddPeriod(GetAsOfDate(),
											plot,string(item->getSecurity().GetCurrencyUnit()->GetCcyName()),
											/*ARM_DEFAULT_COUNTRY*/true,qCredit_Adjust20);
				
				std::auto_ptr<ICM_Cds> Shortcds(
									(ICM_Cds*) GetMktDataMng()->GetDefaultCurve(item->getSecurity().GetSingleName(),GetAsOfDate())->stdCDS(
									GetAsOfDate(),
									EndDate)) ;
				
				const string& name = item->getSecurity().GetSingleName();

				Shortcds->SetSingleName(name);

				ARM_ReferenceValue refvalfix(notional , 1 /* price */, 0 /* K_CONSTANT */);
				
				(Shortcds->GetFeeLeg())->SetCreditSpread(cdsspread);
				(Shortcds->GetFeeLeg())->SetAmount(&refvalfix,100.0);
				(Shortcds->GetDefLeg())->SetAmount(&refvalfix,100.0);
				
				// ICM_Parameters*  Params =  (ICM_Parameters*) GetParameters() ;

				ICM_Pricer_Tree_Binomial_Cds* PricerCds = new ICM_Pricer_Tree_Binomial_Cds(
																						Shortcds.get(),
																						DefModel,
																						GetParameters(),
																						GetAsOfDate());

				PricerCds->SetFaster(true);
								
				//double CdsSensi = - PricerCds->ComputeSensitivity(ICMSPREAD_TYPE,TEMP,label,epsvalue);
				double CdsSensi = - PricerCds->Hedge(ICMSPREAD_TYPE,plot,label,epsvalue);

				if (PricerCds)
					delete PricerCds;
				PricerCds = NULL;

				modifiedprice = (OptionSensi/CdsSensi) + initialprice;

				SetDelta(OptionSensi/CdsSensi) ;

				// On reset les valeurs initiales

				SetMktDataMng(MktDataMng); //On reset le model initial
				itsTree=OldTree;
				SetInitialPrice(initialprice);
				SetPrice(initialprice) ;
				SetSpread(spreadfwd);	
				SetGamma(gamma);
				SetTheta(theta);
				
				if(bool1)
					SetRho(rho);
				if(bool2)
					SetVega(vega);

			}
			break;			
			case ICM_GREEK_GAMMA_TYPE :
			{
				vector<string>* Tenors = new vector<string>;
				
				ICM_BASIS_BUMP Parameters(*Tenors,label,0.02) ;

				ICM_BUMP<ICM_BASIS_BUMP>* bump = new ICM_BUMP<ICM_BASIS_BUMP>(ICMSPREAD_TYPE, 0.02) ;
				bump->Push(Parameters);

				ICM_MktDataMng* NewMktDataMng = MktDataMng->GenerateShiftMng(*bump);

				SetMktDataMng(NewMktDataMng);

				double vol; 
				getParam(ICM_Pricer_Tree_Binomial_SpreadOption::PARAM_FORCEVOL,vol); 
						
				ICM_DefaultCurve* defCurve = (ICM_DefaultCurve*) NewMktDataMng->GetDefaultCurve(item->getSecurity().GetSingleName(),GetAsOfDate());

				TreeServices::calibrate2(itsTree,
										*defCurve,
										vol);

				SetInitialPriceFlg(false);
				double delta1 = ComputeSensitivity(ICM_GREEK_DELTA_TYPE,"NONE",label, epsvalue) ;

				ResetRootPricer();

				if (NewMktDataMng)
					delete NewMktDataMng;
				NewMktDataMng = NULL;

				if (bump)
					delete bump;
				bump = NULL;

				ICM_BASIS_BUMP Parameters2(*Tenors,label,-0.02) ;

				bump = new ICM_BUMP<ICM_BASIS_BUMP>(ICMSPREAD_TYPE, -0.02) ;
				bump->Push(Parameters2);

				NewMktDataMng = MktDataMng->GenerateShiftMng(*bump);

				SetMktDataMng(NewMktDataMng);

				defCurve = (ICM_DefaultCurve*) NewMktDataMng->GetDefaultCurve(item->getSecurity().GetSingleName(),GetAsOfDate());

				TreeServices::calibrate2(itsTree,
										*defCurve,
										vol);

				SetInitialPriceFlg(false);
				double delta2 = ComputeSensitivity(ICM_GREEK_DELTA_TYPE,"NONE",label, epsvalue) ;

				// On reset les valeurs initiales
				
				SetMktDataMng(MktDataMng); //On reset le model initial
				SetInitialPrice(initialprice);
				SetPrice(initialprice) ;
				SetSpread(spreadfwd);
				itsTree=OldTree;
				SetDelta(delta) ;
				SetTheta(theta);
				
				if(bool1)
					SetRho(rho);
				if(bool2)
					SetVega(vega);
				
				double aux1 = (delta1-delta2)/4.;

				modifiedprice = aux1 + initialprice;

				SetGamma(aux1);

				if (NewMktDataMng)
					delete NewMktDataMng;
				NewMktDataMng = NULL;

				if (Tenors)
					delete Tenors;
				Tenors = NULL;
				
				if (bump)
					delete bump;
				bump = NULL;
			}
			break;
			default :
			result = -99999999.0;
		}

	if (!result)
		result=modifiedprice - initialprice;
	}
    

	return (result);

}
//	----------------------------------------------------------------------------------------------
ICM_DefaultCurve&	
ICM_Pricer_Tree_Binomial_SpreadOption::getDefCurve() 
{
	ICM_DefaultCurve* item  = dynamic_cast<ICM_DefaultCurve*> 
		( GetMktDataMng()->GetDefaultCurve(getSpreadOption().getSecurity().GetSingleName(),GetAsOfDate()) ) ;
	if (!item) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Tree_Binomial_SpreadOption::getDefCurve: can't get "<<
		getSpreadOption().getSecurity().GetSingleName() << " " << GetAsOfDate() ); 
	return *item; 
}
//	----------------------------------------------------------------------------------------------
ARM_ZeroCurve&		
ICM_Pricer_Tree_Binomial_SpreadOption::getDiscountCurve()
{
	ARM_ZeroCurve* item = dynamic_cast<ARM_ZeroCurve*>
		( GetMktDataMng()->GetZeroCurve(
		getSpreadOption().GetCurrencyUnit()->GetCcyName(), GetAsOfDate() )) ; 
	if (!item) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Tree_Binomial_SpreadOption::getDiscountCurve: can't get "<<
			item->GetCurrencyUnit()->GetCcyName() << " " << GetAsOfDate()
		); 
	return *item; 
}
//	----------------------------------------------------------------------------------------------
ICM_SpreadOption&	
ICM_Pricer_Tree_Binomial_SpreadOption::getSpreadOption()
{
	ICM_SpreadOption*  item=dynamic_cast<ICM_SpreadOption*>(GetSecurity()); 
	if (!item) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Tree_Binomial_SpreadOption::getSpreadOption: no spreadoption defined"); 
	return *item; 
}
//	----------------------------------------------------------------------------------------------
ICM_Cds&	
ICM_Pricer_Tree_Binomial_SpreadOption::getStdCDS(const ARM_Date& startDate,const ARM_Date&endDate)
{
	ICM_Cds* item = itsCdsReg.get(startDate) ;
	if (!item) 
	{
		item = dynamic_cast<ICM_Cds*> ( getDefCurve().stdCDS(startDate,endDate) ) ;
		if (!item) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Tree_Binomial_SpreadOption::getStdCDS; cant create stdCDS"); 
		itsCdsReg.adopt(startDate,item); 
		item->SetCreditSpread(100.* getSpreadOption().getSecurity().GetFeeLeg().GetCreditSpread()) ;
		item->GetFeeLeg()->GetCreditInfos()->SetAccrualBasis(getSpreadOption().getSecurity().GetAccrualBasis());
//		item->GetFeeLeg()->GetCreditInfos()->SetBasis(getSpreadOption().getSecurity().GetBasis());
	}
	return *item; 	
}
//	----------------------------------------------------------------------------------------------
void
ICM_Pricer_Tree_Binomial_SpreadOption::CpteSpreadsAndPv01s(StdBinomialTree& aTree)
{
	ICM_SpreadOption*  item = & getSpreadOption(); 
	aTree.setTimes(itsTree.Time());	
	unsigned long SizeTree = itsTree.depth();

	const ARM_Vector& PayDates=  getSpreadOption().getSecurity().GetFeeLeg().GetCreditInfos()->GetPayDates(); 
	
	unsigned long PosMaturity(0); 
		if (!itsTree.findpos( CountYears(KACTUAL_365,GetAsOfDate(),item->getExpiryDate()) ,PosMaturity ) ) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"Cant find posMaturity in the tree"); 

	unsigned long PosEndCDS(0); 
		if (!itsTree.findpos(CountYears(KACTUAL_365,GetAsOfDate(), PayDates.Elt(PayDates.GetNumLines()-1)),PosEndCDS)) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"Cant find posEndCDS in the tree"); 

	/*int k0 = 0;
	ICM_QMatrix<double> ret_greeks_PL, ret_greeks_DL, ret_greeks ;
	ret_greeks.Resize(6,1) ;*/

	for (int j = PosMaturity; j >= 0 ; j--)
	{
		ARM_Date jDate ( GetAsOfDate() ); 
		jDate.AddDays(itsTree.Time(j)*365).PreviousBusinessDay(1); 
		ARM_Date Maturity;
		
		if(j==PosMaturity)
			Maturity = Maturity = ARM_Date(PayDates.Elt(PayDates.GetNumLines()-1));
		else
		{
			if (item->underlyingMaturityStyle() == qResidualMaturity) 
				Maturity = ARM_Date(PayDates.Elt(PayDates.GetNumLines()-1));
			else if (item->underlyingMaturityStyle() == qConstantMaturity) 
			{
				int nb = (itsTree.Time(PosEndCDS)-itsTree.Time(PosMaturity))*12 ;		
				char TEMP[4];
				sprintf(TEMP,"%iM",nb);
				Maturity = AddPeriod(jDate.GetJulian(),TEMP,string(item->getSecurity().GetCurrencyUnit()->GetCcyName()),
					/*ARM_DEFAULT_COUNTRY*/true,qCredit_Adjust20);
			}
			else 
				ICMTHROW(ERR_INVALID_ARGUMENT,"Unknonwn Maturity Style "<<item->underlyingMaturityStyle() ); 
		}

		ICM_Cds& Shortcds = getStdCDS(jDate,Maturity) ;			// Underlying CDS priced on the lattice
								
		//	Look for the slices corresponding to the start and end dates of this CDS
		
		const ARM_Vector& payDates   =  Shortcds.GetFeeLeg()->GetCreditInfos()->GetPayDates(); 
		const ARM_Vector& startDates =  Shortcds.GetFeeLeg()->GetCreditInfos()->GetAccStartDates(); 

		Shortcds.GetFeeLeg()->ComputeYF(GetAsOfDate()) ;

		const ARM_Vector& YFAccstartDates =  Shortcds.GetFeeLeg()->GetCreditInfos()->GetYFAccStartDates();
		const ARM_Vector& YFAccpayDates   =  Shortcds.GetFeeLeg()->GetCreditInfos()->GetYFPayDates();

		// 14514 double notional = getSpreadOption().getSecurity().GetDefLeg().GetCreditInfosRef().GetInitialNotional() ;
		double notional = getSpreadOption().getSecurity().GetDefLeg().GetCreditInfosRef().GetNotionals().Elt(0) ;
		double lossesAmount =(1.0 - getDefCurve().GetRecovery()) * notional; 

		unsigned long first,size; 
		double aux = 0.;

		long backwardFrom(itsTree.depth()-1);
		long backwardTo(0);

		size = Shortcds.GetFeeLeg()->GetCreditInfos()->countPaymentFlows(first,GetAsOfDate());
		
		unsigned int j0 = 0;
		aux = YFAccpayDates.Elt(size-1);
		
		ICM_QMatrix<double> indexes(2,1);

		TreeServices::dichotomic_research(itsTree,aux,0 , backwardFrom, indexes) ;
		j0 = indexes(1,0) ;

		unsigned long posMaturity(j); //slice correspondant à la startDate du CDS
		unsigned long posEndCDS(j0);  //slice correspondant à la End Date du CDS
		
		ICM_QMatrix<double> flows(posEndCDS+1,2);

		ICM_QMatrix<double> flows_def(backwardFrom,2);
		
		//--	Retrieve the fee leg

		for(unsigned int i=0;i<=posEndCDS;i++)
		{
			//FeeLeg Flows (with Accrued if any)

			//flows(i,0) représente les payment Dates du flux du Node i
			//flows(i,1) représente le flux au Node i, soit le coupon couru à cette date
			
			if (i>posMaturity)
			{	
				// EffectiveDate : la start Date juste avant Node i , il y en a forcément une
				// PayDate		 : la payment Date juste après Node i , il y en a forcément une

				double EffectiveDate = itsTree.Time(posMaturity);				
				double PayDate = 0.;				
				
				int l  = 0 ,j0 = posEndCDS, i0 = 0;
								
				while (l<size && lt(YFAccstartDates.Elt(first+l),itsTree.Time(i)))
					l++;
		
				EffectiveDate = YFAccstartDates.Elt(first+l-1) ;		// start date just strictly before slice i
				PayDate = YFAccpayDates.Elt(first+l-1);					// Pay Date just after or on slice i
				
				TreeServices::dichotomic_research(itsTree, PayDate, 0, backwardFrom, indexes) ;
				i0 = indexes(1,0) ;

				aux = GetAsOfDate().GetJulian() + 365.*std::_cpp_max(itsTree.Time(posMaturity),EffectiveDate) ;

				ARM_Date AccruedStartDate =  ARM_Date(aux);
				ARM_Date AccruedEndDate = ARM_Date(GetAsOfDate() ); 
				AccruedEndDate.AddDays( std::_cpp_min(itsTree.Time(i),PayDate)*365.) ;
				
				flows(i,0) = itsTree.Time(i0) ;							// the closest node to the real payment date
				flows(i,1)= Shortcds.GetFeeLeg()->GetCreditInfos()->DiffAccrued(
																AccruedStartDate, AccruedEndDate)* notional/100.;
			}
		
			else
			{
				flows(i,0) = 0. ;
				flows(i,1)= 0. ;
			}

			//DefLeg flows
			
			//flows_def(i,0) représente les payment Dates du flux du Node i
			//flows_def(i,1) représente le flux au Node i, soit La Loss à cette date
			
			if(i<backwardFrom)
			{
				aux = YFAccstartDates.Elt(first) ;
			
				if (itsTree.Time(i)<aux)
				{
					flows_def(i,0) = 0. ;
					flows_def(i,1) = 0. ;
				}
				else
				{
					flows_def(i,0)= 0.5 *(itsTree.Time(i)+itsTree.Time(i+1)) + item->getSecurity().GetCreditLag()/365.;
					flows_def(i,1)=lossesAmount ;
				}
			}
		}
					
		ARM_ZeroCurve* discCurve= GetMktDataMng()->GetZeroCurve(
			item->GetCurrencyUnit()->GetCcyName(), GetAsOfDate() );

		ICM_QMatrix<double> statesFrom(itsTree.slice(posEndCDS-(posMaturity-j)).size(),1) ;
		ICM_QMatrix<double> statesFrom_def(itsTree.slice(posEndCDS-(posMaturity-j)).size(),1); 
		ICM_QMatrix<double> feelegFwd , feelegFwd0 ; 
		ICM_QMatrix<double> deflegFwd ,deflegFwd0 ; 

		TreeServices::backwardprice(itsTree,posEndCDS-(posMaturity-j),j,flows,*discCurve,statesFrom,feelegFwd,ICM_QMatrix<double>(6,1));
		TreeServices::backwardprice_default3(itsTree,posEndCDS-(posMaturity-j),j,flows_def,*discCurve,statesFrom_def,deflegFwd,ICM_QMatrix<double>(6,1));

		/*if (j <= 2)
		{
			for (int k = 0 ; k <= j ; k ++)
			{
				ret_greeks(k0,0) = deflegFwd(j-k,0) - feelegFwd(j-k,0) ;
				k0++ ;
			}
		}*/
		/*if (j == PosMaturity)
		{
			TreeServices::backwardprice(itsTree,j,0,ICM_QMatrix<double>(1,2),*discCurve,feelegFwd,feelegFwd0,ret_greeks_PL); 
			TreeServices::backwardprice(itsTree,j,0,ICM_QMatrix<double>(1,2),*discCurve,deflegFwd,deflegFwd0,ret_greeks_DL); 

			for (int i0 = 0 ; i0 < 6 ; i0 ++)
				ret_greeks(i0,0) = ret_greeks_DL(i0,0) - ret_greeks_PL(i0,0) ;
		}*/
		
/******************
		// Partie test
		ICM_QMatrix<double> feelegFwd0 ; 
		ICM_QMatrix<double> deflegFwd0 ; 
		TreeServices::backwardprice(itsTree,j,0,ICM_QMatrix<double>(1,2),*discCurve,feelegFwd,feelegFwd0); 
		TreeServices::backwardprice(itsTree,j,0,ICM_QMatrix<double>(1,2),*discCurve,deflegFwd,deflegFwd0); 

		double spreadFwd; 
		if (eq(feelegFwd0(0,0),0.))
			spreadFwd=0; 
		else
			spreadFwd=(deflegFwd0(0,0)/feelegFwd0(0,0))*100.*item->getSecurity().GetFeeLeg().GetCreditSpread(); 
		
		// ICM_DefaultCurve* DefCurve = (ICM_DefaultCurve*) GetMktDataMng()->GetDefaultCurve(item->getSecurity().GetSingleName(),GetAsOfDate());
		ICM_DefaultCurve& DefCurve = getDefCurve(); // (ICM_DefaultCurve*) GetMktDataMng()->GetDefaultCurve(item->getSecurity().GetSingleName(),GetAsOfDate());
		ARM_ZeroCurve* InitShortCurve =  DefCurve.GetZeroCurve();
		
		ICM_DefaultCurveModel DefModel(&DefCurve, InitShortCurve);

		ICM_Pricer_Cds PricerCds (&Shortcds ,&DefModel);

		PricerCds.SetFaster(true);
		ARM_Date ExecutionDate = DefModel.GetStartDate();
			
		double DL = PricerCds.DefLegPV(ExecutionDate) ;
		double PL = PricerCds.FeeLegPV(ExecutionDate) ;
		double Spread = PricerCds.ComputeSpread(0.) ;

		// fin de la partie tests

****************/ 
		
		for(i=0;i<feelegFwd.Getnbrows();i++) 
		{
			double spread,pv01; 
			if (eq(feelegFwd(i,0),0.)) 
				spread=pv01=0; 
			else 
			{
				spread=deflegFwd(i,0)/feelegFwd(i,0)*100.*getSpreadOption().getSecurity().GetFeeLeg().GetCreditSpread(); 
				pv01 = 100. * feelegFwd(i,0) / (notional * getSpreadOption().getSecurity().GetFeeLeg().GetCreditSpread()); 
			}
			aTree.node(j,i).data().h = std::_cpp_max((item->IsCall()?1:-1)*(spread-item->getStrike()),0.) * pv01 ;
		}
	}
}
			


//	------------------------------------------------------------------------------------
double	
ICM_Pricer_Tree_Binomial_SpreadOption::EuropeanPrice()
{

	double Result = 0. ;
	ICM_SpreadOption&  item= getSpreadOption(); 

	//	Step 1	- We price PL and DL of the underlying CDS
	//
	//		This is done from CDS last payment date to optionMaturity date
	//		With a forward discount curve
	//		We get for DL and PL a vector of output (one value per tree node)
	//	
	unsigned long posMaturity(0); 
	if (!itsTree.findpos( CountYears(KACTUAL_365,GetAsOfDate(),item.getExpiryDate()) ,posMaturity ) ) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Cant find posMaturity in the tree"); 
	
	const ARM_Vector& payDates   =  item.getSecurity().GetFeeLeg().GetCreditInfos()->GetPayDates(); 
	const ARM_Vector& startDates =  item.getSecurity().GetFeeLeg().GetCreditInfos()->GetAccStartDates(); 

	const ARM_Vector& YFAccstartDates = item.getSecurity().GetFeeLeg().GetCreditInfos()->GetYFAccStartDates();
	const ARM_Vector& YFAccpayDates   = item.getSecurity().GetFeeLeg().GetCreditInfos()->GetYFPayDates();

	// 14514 double notional = item.getSecurity().GetDefLeg().GetCreditInfosRef().GetInitialNotional() ;
	double notional = item.getSecurity().GetDefLeg().GetCreditInfosRef().GetNotionals().Elt(0) ;
	double lossesAmount = (1.0- getDefCurve().GetRecovery()) * notional ;
	
	unsigned long posEndCDS(0); 
	if (!itsTree.findpos(CountYears(KACTUAL_365,GetAsOfDate(), payDates.Elt(payDates.GetNumLines()-1)),posEndCDS)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Cant find posEndCDS in the tree"); 
	
	unsigned long first,size; 
	size=  item.getSecurity().GetFeeLeg().GetCreditInfos()->countPaymentFlows(first,GetAsOfDate());
	ICM_QMatrix<double> flows(size,2); 
	
	//FeeLeg Flows

	for(unsigned int i=0;i<size;i++)
	{
		flows(i,0) = YFAccpayDates.Elt(first+i);
		//flows(i,1) = item.getSecurity().GetFeeLeg().GetCreditInfos()->FullCoupon(startDates.Elt(first+i)); 
		flows(i,1) = item.getSecurity().GetFeeLeg().GetCreditInfos()->FullCoupon(first+i); 
	}

	ARM_ZeroCurve& discCurve = getDiscountCurve() ; 
	
	ICM_QMatrix<double> statesFrom(itsTree.slice(posEndCDS).size(),1) ;
	ICM_QMatrix<double> feelegFwd,feelegFwd0 ; 
	ICM_QMatrix<double> ret_greeks_PL, ret_greeks_DL,ret_greeks_NPV, ret_greeks_Price, ret_greeks_defleg; 
	
	ICM_ZeroCurveForward fwdCurve(discCurve,item.getExpiryDate());
	long backwardFrom(itsTree.depth()-1);
	
	long backwardTo(0);
	unsigned int j = 0;
	double aux = 0. ;
	
	ICM_QMatrix<double> statesFrom_def(itsTree.slice(backwardFrom).size(),1);

	TreeServices::backwardprice(itsTree,posEndCDS,posMaturity,flows,discCurve,statesFrom_def,feelegFwd,ICM_QMatrix<double>(6,1)); 

	flows.Resize(backwardFrom+1,2);
		
	ICM_QMatrix<double> flows_def(backwardFrom,2);
	
	for(i=0;i<=backwardFrom;i++)
	{
		//FeeLeg Accrued Flows

		if (i>posMaturity)
		{
			double EffectiveDate = itsTree.Time(posMaturity);				
			double PayDate = 0.;				

			int l  = 0, i0 = 0;;
			
			while (l<size && lt(YFAccstartDates.Elt(first+l),itsTree.Time(i)))
				l++;
			
			EffectiveDate = YFAccstartDates.Elt(first+l-1) ;	// start date just strictly before slice i
			PayDate = YFAccpayDates.Elt(first+l-1) ;			// Pay Date just after or on slice i
			
			ICM_QMatrix<double> indexes(2,1);

			TreeServices::dichotomic_research(itsTree,PayDate ,0 , backwardFrom, indexes) ;
			i0 = indexes(1,0) ;

			flows(i,0) = itsTree.Time(i0) ;						// the closest node to the real payment date
			
			if (fabs(itsTree.Time(i)-PayDate)*365. < 1.)		// The Slice is a Payment Date, flows already taken into account
				flows(i,1) = 0. ;
			else
			{
				aux = GetAsOfDate().GetJulian() + 365.*std::_cpp_max(itsTree.Time(posMaturity),EffectiveDate)  ;
				
				ARM_Date AccruedStartDate =  ARM_Date(aux);
				ARM_Date AccruedEndDate = ARM_Date(GetAsOfDate() ) ;
				AccruedEndDate .AddDays( std::_cpp_min(itsTree.Time(i),PayDate)*365.) ;
				
				flows(i,1)= item.getSecurity().GetFeeLeg().GetCreditInfos()->DiffAccrued(AccruedStartDate, AccruedEndDate);
			}
		}
		else
		{
			flows(i,0) = 0. ;
			flows(i,1) = 0. ;
		}

		//DefLeg flows
		
		//flows_def(i,0) représente les payment Dates du flux du Node i
		//flows_def(i,1) représente le flux au Node i, soit La Loss à cette date

		if (i<backwardFrom)
		{
			aux = YFAccstartDates.Elt(first) ;
			
			if (itsTree.Time(i)<aux)
			{
				flows_def(i,0) = 0. ;
				flows_def(i,1) = 0. ;
			}
			else
			{	
				flows_def(i,0)= 0.5 *(itsTree.Time(i)+itsTree.Time(i+1)) + item.getSecurity().GetCreditLag()/365.;
				flows_def(i,1)=lossesAmount ;
			}
		}
	}

	ICM_QMatrix<double> statesFrom2(itsTree.slice(backwardFrom).size(),1); 
	ICM_QMatrix<double> AccruedFL ;
	
	TreeServices::backwardprice_default3(itsTree, backwardFrom,posMaturity,flows,discCurve,statesFrom2,AccruedFL,ICM_QMatrix<double>(6,1));
	
	for(i=0;i<feelegFwd.Getnbrows();i++) 
		feelegFwd(i,0) = feelegFwd(i,0) + AccruedFL(i,0);

	ICM_QMatrix<double> deflegFwd,deflegFwd0 ; 

	TreeServices::backwardprice_default3(itsTree,posEndCDS,posMaturity,flows_def,discCurve,statesFrom_def,deflegFwd,ICM_QMatrix<double>(6,1)); 
	
	//Discountage
		
	TreeServices::backwardprice(itsTree,posMaturity,0,ICM_QMatrix<double>(1,2),discCurve,feelegFwd,feelegFwd0,ICM_QMatrix<double>(6,1));
	TreeServices::backwardprice(itsTree,posMaturity,0,ICM_QMatrix<double>(1,2),discCurve,deflegFwd,deflegFwd0,ICM_QMatrix<double>(6,1)); 

	
	//	Step 2	- KO European Option pricing is a backward pricing 
	//
	//		CALL	: We assign to the end value the max(fwd-K,0)*PV01fwd for each node
	//		PUT		: We assign to the end value the max(K-fwd,0)*PV01fwd for each node 
	//

	ICM_QMatrix<double> payoff(feelegFwd.Getnbrows(),1); 
	
	for(i=0;i<payoff.Getnbrows();i++) 
	{
		double spread,pv01; 
		if (eq(feelegFwd(i,0),0.)) 
		{
			spread=pv01=0; 
		} 
		else 
		{
			spread=deflegFwd(i,0)/feelegFwd(i,0)*100.*item.getSecurity().GetFeeLeg().GetCreditSpread(); 
			pv01 = 100. * feelegFwd(i,0) / (notional * item.getSecurity().GetFeeLeg().GetCreditSpread()); 
		}
// 		ICMLOG ("noeud "<<i<<" spread="<<spread<<" pv01="<<pv01); 

		payoff(i,0) =  std::_cpp_max((item.IsCall()?1:-1)*(spread-item.getStrike()),0.) * pv01 ;
	}

	ICM_QMatrix<double> optionPrice; 

	TreeServices::backwardprice(itsTree,posMaturity,0,ICM_QMatrix<double>(1,2),discCurve,payoff,optionPrice,ret_greeks_Price); 

	// this is the EUR with KO price ...

	Result = optionPrice(0,0) ; 

	double spreadFwd; 
		
	if (eq(feelegFwd0(0,0),0.))
		spreadFwd=0; 
	else
		spreadFwd=(deflegFwd0(0,0)/feelegFwd0(0,0))*100.*item.getSecurity().GetFeeLeg().GetCreditSpread();
	
	SetSpread(spreadFwd) ; 

	// If No KO

	if (item.IsKO()==false) 
	{
		double OptionPrice_NoKO = 0.;
		
		ARM_Date CdsStartDate = GetAsOfDate() ;
	
		CdsStartDate = CdsStartDate.AddDays(2) ;

		ICM_Cds* Shortcds = & getStdCDS(CdsStartDate,item.getExpiryDate()); 
		
		unsigned long first2,size2; 
		
		size2 = Shortcds->GetFeeLeg()->GetCreditInfos()->countPaymentFlows(first2,GetAsOfDate());

		Shortcds->GetFeeLeg()->ComputeYF(GetAsOfDate()) ;
		const ARM_Vector& YFstartDates2 = Shortcds->GetFeeLeg()->GetCreditInfos()->GetYFAccStartDates();

		ICM_QMatrix<double> flows_def2(backwardFrom,2);
			
		aux = YFstartDates2.Elt(0) ;

		for(i=0;i<=posMaturity;i++)
		{
			if (itsTree.Time(i)<aux)
			{
				flows_def2(i,0) = 0.;
				flows_def2(i,1) = 0.;
			}
			else
			{
				if (item.IsAccelerated())								//   No KO + with Acceleration
					flows_def2(i,0)= itsTree.Time(i+1) ;
									//+ Shortcds->GetCreditLag()/365.;
				else													// No KO + without Acceleration
					flows_def2(i,0)= itsTree.Time(posMaturity);
				
				if(item.IsCall())
					flows_def2(i,1)=lossesAmount ;
			}
		}
		TreeServices::backwardprice_default3(itsTree,posMaturity,0,flows_def2,discCurve,statesFrom_def,deflegFwd, ret_greeks_defleg);

		OptionPrice_NoKO = optionPrice(0,0) + deflegFwd(0,0)*10000./notional ;

		Result = OptionPrice_NoKO ;

		for (int i0 = 0 ; i0 < 6 ; i0  ++)
			ret_greeks_Price(i0,0) = ret_greeks_Price(i0,0) +  ret_greeks_defleg(i0,0)*10000./notional;
 	}

	double Theta = 0. ;

	/*ICM_Pricer_Tree_Binomial_SpreadOption::ComputeGreeksData(ret_greeks_NPV) ;
	
	double Delta = 0., Gamma = 0., Theta = 0. ;
		
	Delta = (notional/10000) *(ret_greeks_Price(4,0)-ret_greeks_Price(3,0))/(ret_greeks_NPV(4,0)-ret_greeks_NPV(3,0));

	//double PV01_CT, PV01_LT, PV01_Fwd ;

	//PV01_CT = getDefCurve().RiskyPV01(GetAsOfDate(), item.getExpiryDate()) ;
	//PV01_LT = getDefCurve().RiskyPV01(GetAsOfDate(), ((ICM_Cds&)(item.getSecurity())).GetEndDate()) ;
	//PV01_Fwd = PV01_LT - PV01_CT ;
	
	//if (item.IsKO()==true)
	//	Delta *= PV01_LT/PV01_Fwd ;
	//else
	//	Delta = (Delta * PV01_LT/PV01_Fwd)- PV01_CT/PV01_Fwd ;
	
	double Delta_Up = 0. , Delta_dpwn = 0. ;
	Delta_Up = (notional/10000) * (ret_greeks_Price(2,0)-ret_greeks_Price(1,0))/(ret_greeks_NPV(2,0)-ret_greeks_NPV(1,0));
	Delta_dpwn = (notional/10000) *(ret_greeks_Price(1,0)-ret_greeks_Price(0,0))/(ret_greeks_NPV(1,0)-ret_greeks_NPV(0,0));
	
	Gamma = (Delta_Up - Delta_dpwn)/4. ;*/

	Theta = (ret_greeks_Price(1,0)-ret_greeks_Price(5,0))/(365. * itsTree.Time(2)) ;
	
	/*SetDelta(Delta) ;
	SetGamma(Gamma)	;*/
	SetTheta(Theta) ;
	
	SetPrice(Result); 
	
	return Result;
}
//	------------------------------------------------------------------------------------
double	
ICM_Pricer_Tree_Binomial_SpreadOption::AmericanPrice()
{

	double Result = 0. ;

	ICM_SpreadOption&item=getSpreadOption(); 

	// ICM_SpreadOption*  item=dynamic_cast<ICM_SpreadOption*>(GetSecurity()); 
	// if (!item)
	// 	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Tree_Binomial_SpreadOption: No Spreadoption defined. "); 

	//	Step 1	- We price PL and DL of the underlying CDS
	//
	//		This is done from CDS last payment date to optionMaturity date
	//		With a forward discount curve
	//		We get for DL and PL a vector of output (one value per tree node)
	//	
	unsigned long posMaturity(0); 
	if (!itsTree.findpos( CountYears(KACTUAL_365,GetAsOfDate(),item.getExpiryDate()) ,posMaturity ) ) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Cant find posMaturity in the tree"); 
	
	// 14514 double notional = item.getSecurity().GetDefLeg().GetCreditInfosRef().GetInitialNotional() ;
	double notional = item.getSecurity().GetDefLeg().GetCreditInfosRef().GetNotionals().Elt(0) ;
	double lossesAmount = (1.0 - getDefCurve().GetRecovery()) * notional ; 

	ICM_QMatrix<double> Times_AM(posMaturity+1,1,0.);
	ICM_QMatrix<double> flowsDef_AM(posMaturity+1,2,0.);
	StdBinomialTree flows_AM ;
	ICM_QMatrix<double> RET(posMaturity+1,posMaturity+1);
	ICM_QMatrix<double> statesFrom_AM(itsTree.slice(posMaturity).size(),1) ; 
	ICM_QMatrix<double> optionPrice;
	ICM_QMatrix<double> ret_greeks_NPV;
	unsigned int i,j ;	
	
	CpteSpreadsAndPv01s(flows_AM); 

	//	Populate the Exercise Dates
	
	if (item.getExerciseStyle()==K_AMERICAN) 
		for(i=0;i<=posMaturity;i++)
			Times_AM(i,0)=itsTree.Time(i);

	else if (item.getExerciseStyle()==K_BERMUDAN) 
	{
		ARM_Date ExerciceDate ;
		int nb = 12/item.getExerciseFrequency() ;
				
		char TEMP[4];
		sprintf(TEMP,"%iM",nb);
		ExerciceDate = AddPeriod(GetAsOfDate().GetJulian(),TEMP,
			string(item.getSecurity().GetCurrencyUnit()->GetCcyName()),/*ARM_DEFAULT_COUNTRY,*/true,qCredit_Adjust20);
		double yfExerciceDate = CountYears(KACTUAL_365,GetAsOfDate(),ExerciceDate) ;

		if (yfExerciceDate > itsTree.Time(posMaturity))
		{
			for(i=0;i<=posMaturity;i++)
				Times_AM(i,0)= itsTree.Time(posMaturity);
		}
		else
		{	
			i=0; j=2;
			while (yfExerciceDate <= itsTree.Time(posMaturity))
			{
				while(i<=posMaturity && yfExerciceDate >= itsTree.Time(i))
				{
					Times_AM(i,0) = yfExerciceDate ;
					i++ ;
				}
				
				sprintf(TEMP,"%iM",nb*j);
				ExerciceDate = AddPeriod(GetAsOfDate().GetJulian(),TEMP,
					string(item.getSecurity().GetCurrencyUnit()->GetCcyName()),/*ARM_DEFAULT_COUNTRY*/true,qCredit_Adjust20);
				ARM_Date aux = ExerciceDate ;

				if(aux.NextBusinessDay(1).PreviousBusinessDay(1) != ExerciceDate)
					ExerciceDate.NextBusinessDay(1) ;		
				yfExerciceDate = CountYears(KACTUAL_365,GetAsOfDate(),ExerciceDate) ;
				j++ ;
			}
		}
	}

	//	should be handled at an upper level 
	//	else ICMTHROW("ICM_Pricer_Tree_Binomial_SpreadOption::AmericanPrice: unsupported exerciseStyle "<< item.getExerciseStyle()) ; 

	//	Populate the flows_Def
	
	if (item.IsCall()) 
	{
		for(i=0;i<=posMaturity;i++) 
			flowsDef_AM(i,0)= lossesAmount * 10000./notional ; 
	}

	if (item.IsAccelerated()) 
	{
		for(i=0;i<=posMaturity;i++) 
			flowsDef_AM(i,1)= itsTree.Time(i+1); 
	}
	else 
	{
		for(i=0;i<=posMaturity;i++) 
			flowsDef_AM(i,1)= Times_AM(i,0); 
	}
	ICM_QMatrix<double> ret_greeks_Price ;

	TreeServices::backwardprice_BN(itsTree,
								   posMaturity, 0,
								   flows_AM, Times_AM, flowsDef_AM,
								   getDiscountCurve(),
								   statesFrom_AM,
								   optionPrice,
								   ret_greeks_Price);
	Result = optionPrice(0,0) ;

	//ICM_Pricer_Tree_Binomial_SpreadOption::ComputeGreeksData(ret_greeks_NPV) ;

	double Delta = 0., Gamma = 0., Theta = 0. ;
	
	//Delta = (notional/10000) *(ret_greeks_Price(4,0)-ret_greeks_Price(3,0))/(ret_greeks_NPV(4,0)-ret_greeks_NPV(3,0));

	/*double PV01_CT, PV01_LT, PV01_Fwd ;

	PV01_CT = getDefCurve().RiskyPV01(GetAsOfDate(), item.getExpiryDate()) ;
	PV01_LT = getDefCurve().RiskyPV01(GetAsOfDate(), ((ICM_Cds&)(item.getSecurity())).GetEndDate()) ;
	PV01_Fwd = PV01_LT - PV01_CT ;
	
	if (item.IsKO()==true)
		Delta *= PV01_LT/PV01_Fwd ;
	else
		Delta = (Delta * PV01_LT/PV01_Fwd)- PV01_CT/PV01_Fwd ;*/
	
	/*double Delta_Up = 0. , Delta_dpwn = 0. ;
	Delta_Up = (notional/10000) * (ret_greeks_Price(2,0)-ret_greeks_Price(1,0))/(ret_greeks_NPV(2,0)-ret_greeks_NPV(1,0));
	Delta_dpwn = (notional/10000) *(ret_greeks_Price(1,0)-ret_greeks_Price(0,0))/(ret_greeks_NPV(1,0)-ret_greeks_NPV(0,0));
	
	Gamma = (Delta_Up - Delta_dpwn)/4. ;*/
	Theta = (ret_greeks_Price(1,0)-ret_greeks_Price(5,0))/(365. * itsTree.Time(2)) ;
	
	/*SetDelta(Delta) ;
	SetGamma(Gamma)	;*/
	SetTheta(Theta) ;
	
	SetPrice(Result); 
	
	return Result; 
}
//	------------------------------------------------------------------------------------
double	
ICM_Pricer_Tree_Binomial_SpreadOption::FwdSpreadPrice()
{

	double Result = 0. ;
	ICM_SpreadOption&  item= getSpreadOption(); 

	//	Step 1	- We price PL and DL of the underlying CDS
	//
	//		This is done from CDS last payment date to optionMaturity date
	//		With a forward discount curve
	//		We get for DL and PL a vector of output (one value per tree node)
	//	
	unsigned long posMaturity(0); 
	if (!itsTree.findpos( CountYears(KACTUAL_365,GetAsOfDate(),item.getExpiryDate()) ,posMaturity ) ) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Cant find posMaturity in the tree"); 
		
	const ARM_Vector& payDates= item.getSecurity().GetFeeLeg().GetCreditInfos()->GetPayDates(); 
	const ARM_Vector& startDates=  item.getSecurity().GetFeeLeg().GetCreditInfos()->GetAccStartDates(); 

	const ARM_Vector& YFAccstartDates = item.getSecurity().GetFeeLeg().GetCreditInfos()->GetYFAccStartDates();
	const ARM_Vector& YFAccpayDates   = item.getSecurity().GetFeeLeg().GetCreditInfos()->GetYFPayDates();
	
	// 14514 double notional = item.getSecurity().GetDefLeg().GetCreditInfosRef().GetInitialNotional() ;
	double notional = item.getSecurity().GetDefLeg().GetCreditInfosRef().GetNotionals().Elt(0) ;
	double lossesAmount = (1.0 - getDefCurve().GetRecovery()) * notional ; 
	
	unsigned long posEndCDS(0); 
	if (!itsTree.findpos(CountYears(KACTUAL_365,GetAsOfDate(), payDates.Elt(payDates.GetNumLines()-1)),posEndCDS)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Cant find posEndCDS in the tree"); 
	
	unsigned long first,size; 
	size=  item.getSecurity().GetFeeLeg().GetCreditInfos()->countPaymentFlows(first,GetAsOfDate());
	ICM_QMatrix<double> flows(size,2); 
	
	//FeeLeg Flows

	for(unsigned int i=0;i<size;i++)
	{
		flows(i,0)=YFAccpayDates.Elt(first+i) ;
		//flows(i,1)= item.getSecurity().GetFeeLeg().GetCreditInfos()->FullCoupon(startDates.Elt(first+i)); 
		flows(i,1)= item.getSecurity().GetFeeLeg().GetCreditInfos()->FullCoupon(first+i); 
	}

	ARM_ZeroCurve& discCurve = getDiscountCurve() ; 
	
	ICM_QMatrix<double> statesFrom(itsTree.slice(posEndCDS).size(),1) ;
	ICM_QMatrix<double> feelegFwd,feelegFwd0 ; 
		
	ICM_ZeroCurveForward fwdCurve(discCurve,item.getExpiryDate());
	long backwardFrom(itsTree.depth()-1);
	
	long backwardTo(0);
	unsigned int j = 0;
	double aux = 0. ;
	
	ICM_QMatrix<double> statesFrom_def(itsTree.slice(backwardFrom).size(),1);

	TreeServices::backwardprice(itsTree,posEndCDS,posMaturity,flows,discCurve,statesFrom_def,feelegFwd,ICM_QMatrix<double>(6,1)); 
	ICM_QMatrix<double> flows_def(backwardFrom,2);

	flows.Resize(backwardFrom+1,2);
	
	for(i=0;i<=backwardFrom;i++)
	{
		//FeeLeg Accrued Flows

		if (i>posMaturity)
		{
			double EffectiveDate = itsTree.Time(posMaturity);
			double PayDate = 0.;				

			int l  = 0, i0 = 0;;

			while (l<size && lt(YFAccstartDates.Elt(first+l),itsTree.Time(i)))
				l++;
			
			EffectiveDate = YFAccstartDates.Elt(first+l-1) ;	// start date just strictly before slice i
			PayDate = YFAccpayDates.Elt(first+l-1);				// Pay Date just after or on slice i
		
			ICM_QMatrix<double> indexes(2,1);
			
			TreeServices::dichotomic_research(itsTree,PayDate ,0 , backwardFrom, indexes) ;
			i0 = indexes(1,0) ;

			flows(i,0) = itsTree.Time(i0) ;						// the closest node to the real payment date

			if (fabs(itsTree.Time(i)-PayDate)*365. < 1.)		// The Slice is a Payment Date, flows already taken into account
				flows(i,1) = 0. ;
			else
			{
				aux = GetAsOfDate().GetJulian() + 365.*std::_cpp_max(itsTree.Time(posMaturity),EffectiveDate)  ;
				
				ARM_Date AccruedStartDate =  ARM_Date(aux);
				ARM_Date AccruedEndDate = ARM_Date(GetAsOfDate() ); 
				AccruedEndDate .AddDays( std::_cpp_min(itsTree.Time(i),PayDate)*365.) ;
				
				flows(i,1)= item.getSecurity().GetFeeLeg().GetCreditInfos()->DiffAccrued(AccruedStartDate, AccruedEndDate);
			}
		}
		else
		{
			flows(i,0) = 0. ;
			flows(i,1) = 0. ;
		}

		//DefLeg Flows
 
		if (i <backwardFrom)
		{
			aux = YFAccstartDates.Elt(first) ;

			if (itsTree.Time(i)<aux)
			{
				flows_def(i,0) = 0. ;
				flows_def(i,1) = 0. ;
			}
			else
			{	
				flows_def(i,0)= 0.5 *(itsTree.Time(i)+itsTree.Time(i+1)) + item.getSecurity().GetCreditLag()/365.;
				flows_def(i,1)=lossesAmount ;
			}
		}
	}

	ICM_QMatrix<double> statesFrom2(itsTree.slice(backwardFrom).size(),1); 
	ICM_QMatrix<double> AccruedFL ;
	
	TreeServices::backwardprice_default3(itsTree, backwardFrom,posMaturity,flows,discCurve,statesFrom2,AccruedFL,ICM_QMatrix<double>(6,1));
	
	for(i=0;i<feelegFwd.Getnbrows();i++) 
		feelegFwd(i,0) = feelegFwd(i,0) + AccruedFL(i,0);

	ICM_QMatrix<double> deflegFwd,deflegFwd0 ; 

	TreeServices::backwardprice_default3(itsTree,posEndCDS,posMaturity,flows_def,discCurve,statesFrom_def,deflegFwd,ICM_QMatrix<double>(6,1));
	
	//Discountage
		
	TreeServices::backwardprice(itsTree,posMaturity,0,ICM_QMatrix<double>(1,2),discCurve,feelegFwd,feelegFwd0,ICM_QMatrix<double>(6,1));
	TreeServices::backwardprice(itsTree,posMaturity,0,ICM_QMatrix<double>(1,2),discCurve,deflegFwd,deflegFwd0,ICM_QMatrix<double>(6,1));

	if (eq(feelegFwd0(0,0),0.))
		Result=0; 
	else
		Result=(deflegFwd0(0,0)/feelegFwd0(0,0))*100.*item.getSecurity().GetFeeLeg().GetCreditSpread();
	
	SetSpread(Result) ; 
	return Result; 
}


// *************************************************************
// Get the Greeks
// *************************************************************

double
ICM_Pricer_Tree_Binomial_SpreadOption::ComputeGreeks (const int& type)
{
	qSENSITIVITY_TYPE Type = (qSENSITIVITY_TYPE) type;

	double Greek = 0. ;
	
	switch (Type)
	{
		case ICM_GREEK_DELTA_TYPE :
		{
			Greek = GetDelta();
		}
		break;
		case ICM_GREEK_GAMMA_TYPE :
		{
			Greek = GetGamma();
		}
		break;
		case ICM_GREEK_VEGA_TYPE :
		{
			Greek = GetVega();
		}
		break;
		case ICM_GREEK_THETA_TYPE :
		{
			Greek = GetTheta();
		}
		break;
		case ICM_GREEK_RHO_TYPE :
		{
			Greek = GetRho();
		}
		break;
		default :
			Greek = GetDelta();
	}

	return (Greek) ;
}

// *******************************************************************************************
// Get the Data for the Long Term CDS in order to compute the greeks on the tree
// *******************************************************************************************


/*void
ICM_Pricer_Tree_Binomial_SpreadOption::ComputeGreeksData(ICM_QMatrix<double>& ret)
{
	ICM_SpreadOption*  item=dynamic_cast<ICM_SpreadOption*>(GetSecurity());
	
	ARM_Date LongMty = ((ICM_Cds&)(item->getSecurity())).GetEndDate() ;
	
	std::auto_ptr<ICM_Cds> Longcds((ICM_Cds*) GetMktDataMng()->GetDefaultCurve(item->getSecurity().GetSingleName(),GetAsOfDate())->stdCDS(
							GetAsOfDate(),
							LongMty)) ;

	ICM_Cds Sec = item->getSecurity();
		
	Longcds->SetCreditSpread(100.*Sec.GetFeeLeg()->GetCreditSpread()) ;
	Longcds->GetFeeLeg()->GetCreditInfos()->SetAccrualBasis(Sec.GetAccrualBasis());
	Longcds->GetFeeLeg()->GetCreditInfos()->SetBasis(Sec.GetBasis());
	
	ret.Resize(6,1) ;
	
	//	Look for the slices corresponding to the start and end dates of this CDS
	
	ARM_Vector& payDates   = *Longcds->GetFeeLeg()->GetCreditInfos()->GetPayDates(); 
	ARM_Vector& startDates = *Longcds->GetFeeLeg()->GetCreditInfos()->GetAccStartDates(); 

	Longcds->GetFeeLeg()->ComputeYF(GetAsOfDate()) ;

	ARM_Vector& YFAccstartDates = * Longcds->GetFeeLeg()->GetCreditInfos()->GetYFAccStartDates();
	ARM_Vector& YFAccpayDates   = * Longcds->GetFeeLeg()->GetCreditInfos()->GetYFPayDates();

	double notional = getSpreadOption().getSecurity().GetDefLeg().GetCreditInfosRef().GetInitialNotional() ;
	double lossesAmount = getDefCurve().GetRecovery() * notional; 

	unsigned long first,size; 
	double aux = 0.;

	long backwardFrom(itsTree.depth()-1);
	long backwardTo(0);

	size = Longcds->GetFeeLeg()->GetCreditInfos()->countPaymentFlows(first,GetAsOfDate());
	
	unsigned int j0 = 0;
	aux = YFAccpayDates.Elt(size-1);
	
	ICM_QMatrix<double> indexes(2,1);
	
	TreeServices::dichotomic_research(itsTree,aux,0 , backwardFrom, indexes) ;
	j0 = indexes(1,0) ;

	unsigned long posMaturity(0); //slice correspondant à la startDate du CDS
	unsigned long posEndCDS(j0);  //slice correspondant à la End Date du CDS
	
	ICM_QMatrix<double> flows(posEndCDS+1,2);

	ICM_QMatrix<double> flows_def(backwardFrom,2);
	
	//--	Retrieve the fee leg

	for(unsigned int i=0;i<=posEndCDS;i++)
	{
		//FeeLeg Flows (with Accrued if any)

		//flows(i,0) représente les payment Dates du flux du Node i
		//flows(i,1) représente le flux au Node i, soit le coupon couru à cette date
		
		if (i>posMaturity)
		{	
			// EffectiveDate : la start Date juste avant Node i , il y en a forcément une
			// PayDate		 : la payment Date juste après Node i , il y en a forcément une

			double EffectiveDate = itsTree.Time(posMaturity);				
			double PayDate = 0.;				
			
			int l  = 0 ,j0 = posEndCDS, i0 = 0;
			
			while (l<size && lt(YFAccstartDates.Elt(first+l),itsTree.Time(i)))
				l++;
			
			EffectiveDate = YFAccstartDates.Elt(first+l-1) ;		// start date just strictly before slice i
			PayDate = YFAccpayDates.Elt(first+l-1);					// Pay Date just after or on slice i
			
			TreeServices::dichotomic_research(itsTree, PayDate, 0, backwardFrom, indexes) ;
			i0 = indexes(1,0) ;

			aux = GetAsOfDate().GetJulian() + 365.*std::_cpp_max(itsTree.Time(posMaturity),EffectiveDate) ;

			ARM_Date AccruedStartDate =  ARM_Date(aux);
			ARM_Date AccruedEndDate = ARM_Date(GetAsOfDate() + ((double) (std::_cpp_min(itsTree.Time(i),PayDate)))*365.) ;
			
			flows(i,0) = itsTree.Time(i0) ;							// the closest node to the real payment date
			flows(i,1)= Longcds->GetFeeLeg()->GetCreditInfos()->DiffAccrued(
																AccruedStartDate, AccruedEndDate)* notional/100.;
		}
		
		else
		{
			flows(i,0) = 0. ;
			flows(i,1)= 0. ;
		}

		//DefLeg flows
		
		//flows_def(i,0) représente les payment Dates du flux du Node i
		//flows_def(i,1) représente le flux au Node i, soit La Loss à cette date
		
		if(i<backwardFrom)
		{
			aux = YFAccstartDates.Elt(first) ;
			
			if (itsTree.Time(i)<aux)
			{
				flows_def(i,0) = 0. ;
				flows_def(i,1) = 0. ;
			}
			else
			{
				flows_def(i,0)= 0.5 *(itsTree.Time(i)+itsTree.Time(i+1)) + item->getSecurity().GetCreditLag()/365.;
				flows_def(i,1)=lossesAmount ;
			}
		}
	}
	
	ARM_ZeroCurve* discCurve= GetMktDataMng()->GetZeroCurve(
		item->GetCurrencyUnit()->GetCcyName(), GetAsOfDate() );

	ICM_QMatrix<double> statesFrom(itsTree.slice(posEndCDS-(posMaturity)).size(),1) ;
	ICM_QMatrix<double> statesFrom_def(itsTree.slice(posEndCDS-(posMaturity)).size(),1); 
	ICM_QMatrix<double> feelegFwd , feelegFwd0 ; 
	ICM_QMatrix<double> deflegFwd ,deflegFwd0 ; 
	ICM_QMatrix<double> ret_greeks_DL, ret_greeks_PL ;

	TreeServices::backwardprice(itsTree,posEndCDS,0,flows,*discCurve,statesFrom,feelegFwd,ret_greeks_PL);
	TreeServices::backwardprice_default3(itsTree,posEndCDS,0,flows_def,*discCurve,statesFrom_def,deflegFwd,ret_greeks_DL);

	for (int i0 = 0 ; i0 < 6 ; i0 ++)
		ret(i0,0) = ret_greeks_DL(i0,0) - ret_greeks_PL(i0,0) ;
			
}*/

/** 

void
ICM_Pricer_Tree_Binomial_SpreadOption::ComputeGreeksData(ICM_QMatrix<double>& ret)
{

	ICM_SpreadOption*  item=dynamic_cast<ICM_SpreadOption*>(GetSecurity());
	
	const ICM_Cds& cds = dynamic_cast<const ICM_Cds&>(item->getSecurity()); // cds.GetDefLeg().GetCreditInfosRef().GetAccEndDates().last()
	ARM_Date LongMty = cds.GetEndDateNA() ;// ((ICM_Cds&)(item->getSecurity())).GetEndDate() ;
	
	std::auto_ptr<ICM_Cds> Longcds((ICM_Cds*) GetMktDataMng()->GetDefaultCurve(item->getSecurity().GetSingleName(),GetAsOfDate())->stdCDS(
							GetAsOfDate(),
							LongMty)) ;

	ICM_Cds Sec = item->getSecurity();
		
	Longcds->SetCreditSpread(100.*Sec.GetFeeLeg()->GetCreditSpread()) ;
	Longcds->GetFeeLeg()->GetCreditInfos()->SetAccrualBasis(Sec.GetAccrualBasis());
	//Longcds->GetFeeLeg()->GetCreditInfos()->SetBasis(Sec.GetBasis());
	
	ret.Resize(6,1) ;
	

//	--	Retrieve the fee leg
	
	ICM_Leg& feeLeg=*Longcds->GetFeeLeg(); 
	ICM_Security&sec=*feeLeg.GetCreditInfos(); 

	const ARM_Vector& payDates=sec.GetPayDates(); 
	const ARM_Vector& startDates=sec.GetAccStartDates(); 

// 	14514 double notional = getSpreadOption().getSecurity().GetDefLeg().GetCreditInfosRef().GetInitialNotional() ;
	double notional = getSpreadOption().getSecurity().GetDefLeg().GetCreditInfosRef().GetNotionals().Elt(0) ;
	double losses = getDefCurve().GetRecovery() * notional; 

	//	--	
	unsigned long first,size; 
	size=sec.countPaymentFlows(first,GetAsOfDate()); 
	
	ICM_QMatrix<double> flows(size,2); 
	
	for(unsigned int i=0;i<size;i++)
	{
		flows(i,0)=CountYears(KACTUAL_365,GetAsOfDate(),ARM_Date(payDates.Elt(first+i))) ;
		//flows(i,1)=sec.FullCoupon(startDates.Elt(first+i))* notional/100; 
		flows(i,1)=sec.FullCoupon(first+i)* notional/100; 
	}

	ARM_ZeroCurve* discCurve= GetMktDataMng()->GetZeroCurve(
		Longcds->GetCurrencyUnit()->GetCcyName(), GetAsOfDate() );
	
	long backwardFrom(itsTree.depth()-1); 
	long backwardTo(0); 
	ICM_QMatrix<double> statesFrom(itsTree.slice(backwardFrom).size(),1); 
	ICM_QMatrix<double> ret1 ;
	ICM_QMatrix<double> ret_greeks_PL ;

	TreeServices::backwardprice(itsTree,
								backwardFrom,
								backwardTo,
								flows,
								*discCurve,
								statesFrom,
								ret1,
								ret_greeks_PL);

	ICM_QMatrix<double> statesFrom2(itsTree.slice(backwardFrom).size(),1); 
	ICM_QMatrix<double> ret2 ;

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
			ARM_Date AccruedEndDate = ARM_Date(GetAsOfDate() ); 
			AccruedEndDate .AddDays( std::_cpp_min(itsTree.Time(i),PayDate)*365.) ;
			
			flows(i,0) = itsTree.Time(i0) ; // Le noeud qui se rapproche le plus de la vraie PayDate
			
			if (fabs(itsTree.Time(i)-PayDate)*365. <= 2.) // le 2 est à revoir!!
				flows(i,1) = 0. ;
			else
				flows(i,1)= sec.DiffAccrued(AccruedStartDate, AccruedEndDate)*notional/100.;
		}
		else
		{
			flows(i,0) = 0. ;
			flows(i,1) = 0. ;
		}
	}


	ICM_QMatrix<double> ret_greeks_AccruedPL ;

	TreeServices::backwardprice_default3(
										itsTree,
										backwardFrom,
										backwardTo,
										flows,
										*discCurve,
										statesFrom2,
										ret2,
										ret_greeks_AccruedPL);

			
	ICM_QMatrix<double> flows_def(backwardFrom,2);
	
	ARM_Date EffectiveDate = startDates.Elt(first);
	double aux = CountYears(KACTUAL_365,GetAsOfDate(),EffectiveDate) ;

	for(i=0;i<backwardFrom;i++)
	{
		if (itsTree.Time(i)<aux)
		{
			flows_def(i,0) = 0. ;
			flows_def(i,1) = 0. ;
		}
		else
		{
			flows_def(i,0)= 0.5 *(itsTree.Time(i)+itsTree.Time(i+1)) + Longcds->GetCreditLag()/365.;
			flows_def(i,1)=losses ;
		}
	}
	
	ICM_QMatrix<double> statesFrom3(itsTree.slice(backwardFrom).size(),1); 
	ICM_QMatrix<double> ret3 ;
	ICM_QMatrix<double> ret_greeks_DL ;

	TreeServices::backwardprice_default3(
								itsTree,
								backwardFrom,
								backwardTo,
								flows_def,
								*discCurve,
								statesFrom3, 
								ret3,
								ret_greeks_DL);

	for (int i0 = 0 ; i0 < 6 ; i0 ++)
		ret(i0,0) = ret_greeks_DL(i0,0) - (ret_greeks_PL(i0,0)+ret_greeks_AccruedPL(i0,0));
}
**/ 
double ICM_Pricer_Tree_Binomial_SpreadOption::ComputeSpread(const double& MtM )  
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
double ICM_Pricer_Tree_Binomial_SpreadOption::Accrued() 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
double ICM_Pricer_Tree_Binomial_SpreadOption::FeeLegPV () 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
double ICM_Pricer_Tree_Binomial_SpreadOption::DefLegPV ()
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
double ICM_Pricer_Tree_Binomial_SpreadOption::ComputeDuration(void) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
double ICM_Pricer_Tree_Binomial_SpreadOption::ComputeImpliedVol(const double& Price) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}
double ICM_Pricer_Tree_Binomial_SpreadOption::Compute_Fwd_Spread(const ARM_Date &Mty, const ARM_Date &CDS_ExpiryDate, double& dur) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
	}

