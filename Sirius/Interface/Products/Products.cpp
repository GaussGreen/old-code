//	Products.cpp : Implementation of DLL Exports.
//
//	Note: Proxy/Stub Information
//        To build a separate proxy/stub DLL, 
//        run nmake -f Productsps.mk in the project directory.
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "resource.h"
#include <initguid.h>
#include "products.h"
#include "products_i.c"
#include "sirius_i.c"

#include "product_napoleon.h"
#include "product_cappedflooredcliquet.h"
#include "product_barriercliquet.h"
#include "product_generic.h"
#include "product_europeanoption.h"
#include "product_future.h"
#include "product_callableCliquetSwap.h"
#include "product_momentumcliquet.h"
#include "product_stock.h"
#include "MlEqPde.h"

#include "MlEqObjects.h"
#include "MlEqParameterList.h"
#include "threemoments.h"
#include "hpsort.h"
#include "product_himalayan.h"
#include "MonteCarlo.h"
#include "MlEqSwap.h"
#include "product_varianceswap.h"
#include "product_whale.h"
#include "product_basketcliquet.h"
#include "Heston.h"
#include "mleqshortmodels.h"
#include "localvol.h"
#include "PdeProducts.h"
#include "ran2.h"
#include "product_callableNote.h"
#include "product_autocallableNote.h"
#include "product_bermudeanNote.h"

CProductsComModule							_Module;

BEGIN_OBJECT_MAP(ObjectMap)
	OBJECT_ENTRY(CLSID_Napoleon, CNapoleon)
	OBJECT_ENTRY(CLSID_CappedFlooredCliquet, CCappedFlooredCliquet)
	OBJECT_ENTRY(CLSID_BarrierCliquet, CBarrierCliquet)
	OBJECT_ENTRY(CLSID_Generic, CGeneric)
	OBJECT_ENTRY(CLSID_Himalayan, CHimalayan)
	OBJECT_ENTRY(CLSID_EuropeanOption, CEuropeanOption)
	OBJECT_ENTRY(CLSID_Future, CFuture)
	OBJECT_ENTRY(CLSID_CallableCliquetSwap, CCallableCliquetSwap)
	OBJECT_ENTRY(CLSID_MomentumCliquet, CMomentumCliquet)
	OBJECT_ENTRY(CLSID_Stock, CStock)
	OBJECT_ENTRY(CLSID_varianceswap, CVarianceSwapPricer)
	OBJECT_ENTRY(CLSID_Whale, CWhale)
	OBJECT_ENTRY(CLSID_BasketCliquet, CBasketCliquet)
	OBJECT_ENTRY(CLSID_CallableNote, CCallableNote)
	OBJECT_ENTRY(CLSID_AutocallableNote, CAutocallableNote)
	OBJECT_ENTRY(CLSID_BermudeanNote, CBermudeanNote)
END_OBJECT_MAP()


/////////////////////////////////////////////////////////////////////////////
//	Dll[...] standard functions
//
//  Used to determine whether the DLL can be unloaded by OLE
//
STDAPI DllCanUnloadNow(void)
{
    return (_Module.GetLockCount()==0) ? S_OK : S_FALSE;
}
//
//	Returns a class factory to create an object of the requested type
//
STDAPI DllGetClassObject(REFCLSID rclsid, REFIID riid, LPVOID* ppv)
{
    return _Module.GetClassObject(rclsid, riid, ppv);
}
//
//	DLL Entry Point
//
extern "C" BOOL WINAPI DllMain(HINSTANCE hInstance, DWORD dwReason, LPVOID /*lpReserved*/)
{
    if (dwReason == DLL_PROCESS_ATTACH){
		_Module.Init(ObjectMap, hInstance, &LIBID_Products);
        DisableThreadLibraryCalls(hInstance);		
    } else if (dwReason == DLL_PROCESS_DETACH){        
		_Module.Term();
	}
    return TRUE;
}
//
//	Adds entries to the system registry - regsvr32.exe calls this function
//
STDAPI DllRegisterServer(void)
{
    // registers object, typelib and all interfaces in typelib
    return _Module.RegisterServer(TRUE);
}
//
//	Removes entries from the system registry - regsvr32.exe calls this function
//
STDAPI DllUnregisterServer(void)
{
    return _Module.UnregisterServer(TRUE);
}


/////////////////////////////////////////////////////////////////////////////
//	Emitted functions
//
//	all these functions are emitted to SiriusXL.xla
//


/****************************************************************
**	 
**	Routine: MLBlackScholes
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLBlackScholes(VARIANT Forward, VARIANT Strike, VARIANT Volatility, VARIANT Maturity, VARIANT DiscountFactor, VARIANT CallPut, VARIANT BlendOpt,VARIANT VolVolOpt,VARIANT nstdevOpt,VARIANT npointsOpt, VARIANT* pResult)
{
	begin_function
	map_parameter(Forward, double, fForward);
	map_parameter(Strike, double, fStrike);
	map_parameter(Volatility, double, fVolatility);
	map_parameter(Maturity, double, fMaturity);
	map_parameter(DiscountFactor, double, fDiscountFactor);
	map_parameter(CallPut, long, nCallPut);
	map_optional_parameter(BlendOpt, double, fBlend, 0.0);		
	map_optional_parameter(VolVolOpt, double, fvolvol, 0.0);		
	map_optional_parameter(nstdevOpt, long, nstdev, 4);		
	map_optional_parameter(npointsOpt, long, npoints, 12);		


	double res;
	if ( fabs(fvolvol) < 1e-3 || npoints == 0 ){
		res = Bs(fForward, fVolatility, fMaturity, fStrike, fDiscountFactor, nCallPut, fBlend);
	}
	else{
		res = Bs(fForward, fVolatility, fMaturity, fStrike, fDiscountFactor, nCallPut, fBlend,fvolvol,nstdev,npoints);
	}

	// return the result
	CParameterMap pm;
	pm.SetValue(res);
	return pm.GetValue(pResult);	
	end_function
}


HRESULT __stdcall MLBlackScholesGreeks(VARIANT Forward, VARIANT Strike, VARIANT Volatility, VARIANT Maturity, VARIANT DiscountFactor, VARIANT CallPut, VARIANT BlendOpt,VARIANT VolVolOpt,VARIANT nstdevOpt,VARIANT npointsOpt, VARIANT* pResult)
{
	begin_function
	map_parameter(Forward, double, fForward);
	map_parameter(Strike, double, fStrike);
	map_parameter(Volatility, double, fVolatility);
	map_parameter(Maturity, double, fMaturity);
	map_parameter(DiscountFactor, double, fDiscountFactor);
	map_parameter(CallPut, long, nCallPut);
	map_optional_parameter(BlendOpt, double, fBlend, 0.0);		
	map_optional_parameter(VolVolOpt, double, fvolvol, 0.0);		
	map_optional_parameter(nstdevOpt, long, nstdev, 4);		
	map_optional_parameter(npointsOpt, long, npoints, 12);		

	double eps=0.001;

	CVector res(4);

	res[0] = Bs(fForward, fVolatility, fMaturity, fStrike, fDiscountFactor, nCallPut, fBlend);

	res[1]=1/eps*(Bs(fForward, fVolatility, fMaturity, fStrike, fDiscountFactor, nCallPut, fBlend)-Bs(fForward + eps , fVolatility, fMaturity, fStrike, fDiscountFactor, nCallPut, fBlend));

	res[2]=1/eps/eps*(Bs(fForward-eps, fVolatility, fMaturity, fStrike, fDiscountFactor, nCallPut, fBlend)+Bs(fForward+eps, fVolatility, fMaturity, fStrike, fDiscountFactor, nCallPut, fBlend)-2*Bs(fForward, fVolatility, fMaturity, fStrike, fDiscountFactor, nCallPut, fBlend));

	res[3]=1/eps*(Bs(fForward, fVolatility+eps, fMaturity, fStrike, fDiscountFactor, nCallPut, fBlend)-Bs(fForward , fVolatility, fMaturity, fStrike, fDiscountFactor, nCallPut, fBlend));

	// return the result
	CParameterMap pm;
	pm.SetValue(res);
	return pm.GetValue(pResult);	
	end_function
}


/****************************************************************
**	 
**	Routine: MLBlackScholesKnockOutDigital
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLBlackScholesKnockOutDigital(VARIANT Forward,VARIANT Spot, VARIANT Volatility, VARIANT Maturity, VARIANT BarrierUp,VARIANT BarrierDown,VARIANT Rebate, VARIANT DiscountFactor, VARIANT CallPut, VARIANT BlendOpt,VARIANT VolVolOpt,VARIANT nstdevOpt,VARIANT npointsOpt, VARIANT* pResult)
{
	begin_function
	map_parameter(Forward, double, fForward);
	map_parameter(Spot, double, fSpot);
	map_parameter(BarrierUp, double, fbarrierUp);
	map_parameter(BarrierDown, double, fbarrierDown);
	map_parameter(Rebate, double, frebate);

	map_parameter(Volatility, double, fVolatility);
	map_parameter(Maturity, double, fMaturity);
	map_parameter(DiscountFactor, double, fDiscountFactor);
	map_parameter(CallPut, long, nCallPut);
	map_optional_parameter(BlendOpt, double, fBlend, 0.0);		
	map_optional_parameter(VolVolOpt, double, fvolvol, 0.0);		
	map_optional_parameter(nstdevOpt, long, nstdev, 4);		
	map_optional_parameter(npointsOpt, long, npoints, 12);		

	double price;


	double_knock_out_rebate(price,fMaturity,fSpot,fForward,
						  fDiscountFactor,fVolatility,fBlend,fbarrierUp,
						  fbarrierDown,frebate,fvolvol,nstdev,npoints);




	// return the result
	CParameterMap pm;
	pm.SetValue(price);
	return pm.GetValue(pResult);	
	end_function
}

/****************************************************************
**	 
**	Routine: MLKnockOutDigital
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLKnockOutDigital(VARIANT DownBarrier,VARIANT Spot, VARIANT Maturity,
									VARIANT Rate,VARIANT DvidentYield,
									VARIANT SpotBarrierVols,VARIANT SpotBarrierVolSlope,
									VARIANT FwdBarrierVols,VARIANT FwdBarrierVolsSlope,
									VARIANT PayAtTheEnd,
									VARIANT npoints,VARIANT nplotpoints,VARIANT plotCummulativeHitProb,
									VARIANT* pResult)
{
	begin_function

	map_parameter(DownBarrier, double, fDownBarrier);
	map_parameter(Spot, double, fSpot);
	map_parameter(Maturity, double, fMaturity);
	map_parameter(Rate, double, fRate);
	map_parameter(DvidentYield, double, fYield);
	map_parameter(SpotBarrierVols, double, fSpotBarrierVols);
	map_parameter(SpotBarrierVolSlope, double, fSpotBarrierVolSlope);
	map_parameter(FwdBarrierVols, double, fFwdBarrierVols);
	map_parameter(FwdBarrierVolsSlope, double, fFwdBarrierVolsSlope);
	map_parameter(npoints, long, nnpoints);
	map_parameter(nplotpoints, long, nnplotpoints);
	map_parameter(plotCummulativeHitProb, bool, bplotCummulativeHitProb);
	map_parameter(PayAtTheEnd, bool, bPayAtTheEnd);

	

	double price;
	CMatrix AmPrices;

	int nmethod = 3;
	AmericanDownInPut(price,AmPrices,nmethod,fDownBarrier,
					   fSpot,fMaturity,fRate,fYield,
					   fSpotBarrierVols,fSpotBarrierVolSlope,
					   fFwdBarrierVols,fFwdBarrierVolsSlope,
					   bPayAtTheEnd,
					   nnpoints,nnplotpoints,bplotCummulativeHitProb);


	CMatrix result;
	if ( AmPrices.rows() > 0 )
	{
		result.resize(AmPrices.rows()+1,2);
		result[0][0] = price;
		for ( int i = 0 ; i < AmPrices.rows(); i++ )
		{
			for ( int j = 0 ; j < AmPrices.cols(); j++ ){
				result[i+1][j] = AmPrices[i][j];
			}
		}
	}
	else
	{
		result.resize(1,1);
		result[0][0] = price;
	}


	// return the result
	CParameterMap pm;
	pm.SetValue(result);
	return pm.GetValue(pResult);	
	end_function
}




/****************************************************************
**	 
**	Routine: MLCholeskyDecomposition
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLCholeskyDecomposition(VARIANT CorrelationMatrix, VARIANT* pResult)
{
	begin_function
	map_parameter(CorrelationMatrix, CMatrix, correl);

	CMatrix cholesky;
	createCholesky(cholesky,correl);

	// return the result
	CParameterMap pm;
	pm.SetValue(cholesky);
	return pm.GetValue(pResult);

	end_function
}



/****************************************************************
**	 
**	Routine: testMultiAssetMC
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLtestMultiAssetMC(VARIANT assets, VARIANT mcDates, VARIANT modelParameters,VARIANT dummy,  VARIANT* pResult)
{

	begin_function
	
	map_object_vector_parameter(assets, Asset, vAssets);
	map_parameter(modelParameters, CVector, modelparams);
	map_parameter(mcDates, std::vector< long >, dates);

	CMatrix calibTimeIndex;
	CMatrix calibVols;
	vector<vector<MlEqStrikeHandle> > pCalibStrikes;

	int nassets = vAssets.size();
	GVector<int> rateIsStochastic(nassets);

	vector < MlEqStochBetaVolatilityStructureHandle > pbetaVolStruct(nassets);
	
	for ( int i = 0; i < nassets; i++ ){
		pbetaVolStruct[i] = dynamic_cast<MlEqStochBetaVolatilityStructure*>(&*vAssets[i]->GetVolatilityStructure());
	}

	MlEqConstDateHandle hDate = pbetaVolStruct[0]->getDateToDouble();
			
	CMultiAssetForwardSkewMC multiMc(hDate);
	multiMc.Initialize(vAssets, dates, pbetaVolStruct, rateIsStochastic, modelparams, calibTimeIndex, calibVols, pCalibStrikes);
	multiMc.generatePaths();

	end_function
}	


/****************************************************************
**	 
**	Routine: MLFitPolynomial
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/



HRESULT __stdcall MLFitPolynomial(VARIANT InitialGuess, VARIANT xData, VARIANT yData,VARIANT xValBounds, VARIANT FinalTolerance, VARIANT StoppingTolerance, VARIANT* pResult)
{		
	double								fInitialTolerance;
	CMatrix								objBounds;
	CMatrix								NonZeroPartialDerivativesSpecification;
	int									outputFlag = 0;
	fitPolynomial						fP;
	int									returnCode;
	int									maximizeFlag = 0;

	// map input parameters
	begin_function
	map_parameter(InitialGuess, CVector, vectorInitialGuess);
	map_parameter(xData, CVector, vxData);
	map_parameter(yData, CVector, vyData);
	map_parameter(xValBounds, CMatrix, matrix_xValBounds);
	map_parameter(FinalTolerance, double, fFinalTolerance);
	map_parameter(StoppingTolerance, double, fStoppingTolerance);
	
	if ( vxData.getsize() != vyData.getsize() ){
		throw("xData size must equal yData size");
	}

	CMatrix matrixData(vxData.getsize(),2);

	for ( int i = 0 ; i < vxData.getsize(); i++ ){

		matrixData[i][0] = vxData[i];
		matrixData[i][1] = vyData[i];
	}

	// do the work
	fInitialTolerance = 2.0 * fFinalTolerance;
	fP.initialize(vectorInitialGuess, 
				  matrixData, 
				  matrix_xValBounds, 
				  objBounds, 
				  fInitialTolerance, 
				  fFinalTolerance, 
				  fStoppingTolerance, 
				  NonZeroPartialDerivativesSpecification,
				  outputFlag);

	CVector results(vectorInitialGuess);
	fP.solve(results, maximizeFlag, returnCode);
	
	// return the result
	CParameterMap pm;
	pm.SetValue(results);
	return pm.GetValue(pResult);
	end_function
}


/****************************************************************
**	 
**	Routine: MLGetDiscountFactor
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/



HRESULT __stdcall MLGetDiscountFactor(VARIANT ZeroCurve, VARIANT DiscountDate ,VARIANT *pResult)
{			
	double								fDiscountFactor;

	// map the inputs
	begin_function		
	map_com_object_parameter(ZeroCurve, ZeroCurve, spZeroCurve)
	map_parameter(DiscountDate, double, fDate)

	// do the work
	if (spZeroCurve->GetDiscountFactor(fDate, &fDiscountFactor)) propagate_error;
	
	// return the result
	return CComVariant(fDiscountFactor).Detach(pResult);
	end_function
}

/****************************************************************
**	 
**	Routine: MLGetDiscountFactor
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/



HRESULT __stdcall MLGetForwardDiscountFactor(VARIANT ZeroCurve, VARIANT StartDate, VARIANT DiscountDate ,VARIANT *pResult)
{			
	double								fDiscountFactor1;
    double								fDiscountFactor2;


	// map the inputs
	begin_function		
	map_com_object_parameter(ZeroCurve, ZeroCurve, spZeroCurve)
	map_parameter(DiscountDate, double, fDate)
    map_parameter(StartDate, double, fStartDate)

	if ( fDate < fStartDate ){
		throw("DiscountDate must be greater or equal StartDate in MLGetForwardDiscountFactor");
	}

	// do the work
	if (spZeroCurve->GetDiscountFactor(fDate, &fDiscountFactor1)) propagate_error;
	if (spZeroCurve->GetDiscountFactor(fStartDate, &fDiscountFactor2)) propagate_error;
	
	// return the result
	return CComVariant(fDiscountFactor1/fDiscountFactor2).Detach(pResult);
	end_function
}



/****************************************************************
**	 
**	Routine: MLGetVolatility
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLGetVolatility(VARIANT VolatilityStructureHandle, VARIANT StrikeHandle, VARIANT VolDate, VARIANT BidOfferOpt, VARIANT* pResult)
{
	int									nStrikes;
	vector<MlEqStrikeHandle>			m_pStrikes;

	// map the inputs
	begin_function	
	map_object_parameter(VolatilityStructureHandle, VolatilityStructure, hVolatilityStructure);
	map_object_parameter(StrikeHandle, Strikes, hStrike);
	map_parameter(VolDate, double, fVolDate);
	map_optional_enum_parameter(BidOfferOpt, BidAskMidEnum, nBidOffer, Middle);
		
	// do the work
	m_pStrikes = hStrike->m_hStrikes;
	nStrikes = hStrike->m_hStrikes.size();
	CVector vector(nStrikes);
	for (int i = 0; i < nStrikes; i++){
		vector[i] = hVolatilityStructure->getVol(*(hStrike->m_hStrikes[i]), fVolDate, nBidOffer);
	}
	
	// return the result
	CParameterMap pm;
	pm.SetValue(vector);
	return pm.GetValue(pResult);	
	end_function
}


/****************************************************************
**	 
**	Routine: MLGetFutureVolatility
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLGetFutureVolatility(VARIANT VolatilityStructureHandle, VARIANT StrikeHandle, VARIANT VolStartDate, VARIANT VolEndDate, VARIANT BidOfferOpt, VARIANT* pResult)
{
	int									nStrikes;
	vector<MlEqStrikeHandle>			m_pStrikes;
				
	// map the inputs
	begin_function	
	map_object_parameter(VolatilityStructureHandle, VolatilityStructure, hVolatilityStructure);
	map_object_parameter(StrikeHandle, Strikes, hStrike);
	map_parameter(VolStartDate, double, sVolDate);
	map_parameter(VolEndDate, double, eVolDate);
	map_optional_enum_parameter(BidOfferOpt, BidAskMidEnum, nBidOffer, Middle);
		
	// do the work
	m_pStrikes = hStrike->m_hStrikes;
	nStrikes = hStrike->m_hStrikes.size();
	CVector vector(nStrikes);
	for (int i = 0; i < nStrikes; i++){
		vector[i] = hVolatilityStructure->getFutureVol(*(hStrike->m_hStrikes[i]), sVolDate,eVolDate, nBidOffer);
	}
	
	// return the result
	CParameterMap pm;
	pm.SetValue(vector);
	return pm.GetValue(pResult);	
	end_function
}

/****************************************************************
**	 
**	Routine: MLCalibLocalVol
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLCalibrateEffectiveLocalVol(VARIANT asset,VARIANT startDate,VARIANT endDate,VARIANT calibStrikes,VARIANT calibSpots,VARIANT targetVols,VARIANT initialGuessOpt,VARIANT nGauss, VARIANT calcStrikesOpt,VARIANT useParametricLVOpt,VARIANT tanhWingInfoOpt,VARIANT* pResult)
{	
	
	vector< MlEqStrikeHandle >			strikes;
	vector< MlEqStrikeHandle >			spotStrikes;
	vector< MlEqStrikeHandle >			cStrikes;
	
	begin_function	
		
	map_object_parameter(asset, Asset, hasset);
	map_object_parameter(startDate, Date, startD);
	map_parameter(endDate, long,enddate );
	map_parameter(targetVols, CVector,target );
	map_parameter(nGauss, long,ngauss );
	map_object_parameter(calibStrikes, Strikes, hStrike);
	map_object_parameter(calibSpots, Strikes, hSpots);
	map_optional_parameter(initialGuessOpt, CVector, fguess,CVector());
	map_object_parameter(calcStrikesOpt, Strikes, calcStrikes);
	map_optional_parameter(useParametricLVOpt, bool, usePLV,false);
	map_optional_parameter(tanhWingInfoOpt, CVector, tanhInfo,CVector());
	
	strikes = hStrike->m_hStrikes;
	int nstrikes = strikes.size();
	GVector< MlEqStrikeHandle > calibStrikes(nstrikes);
	for ( int i = 0 ; i < nstrikes; i++ ){
		calibStrikes[i] = strikes[i];
	}
	
	
	spotStrikes = hSpots->m_hStrikes;
	nstrikes = spotStrikes.size();
	GVector< MlEqStrikeHandle > calibSpots(nstrikes);
	for ( int i = 0 ; i < nstrikes; i++ ){
		calibSpots[i] = spotStrikes[i];
	}
	
	cStrikes = calcStrikes->m_hStrikes;
	calibrateEffectiveLocalVol effLv;
		
	double finalTol		= 1e-6;
	double stoppingTol	= 1e-6;

	effLv.init(*hasset,startD,enddate,
			   calibStrikes,calibSpots,target,fguess,ngauss,tanhInfo,finalTol,stoppingTol);

	
	CMatrix res;
	
	MlEqStrike stk;
	MlEqStrike::convertStrikes(stk,*cStrikes[0]);
	
	if ( stk.m_strike < 0 )
	{
		if ( !usePLV )
		{
//			effLv.calibrate(res);
		}
		else
		{
			effLv.calibrateLV(res,fguess);
		}

	}
	else
	{
		GVector< MlEqStrikeHandle > stk(cStrikes.size());
		for ( int i = 0 ; i < stk.getsize(); i++ ){
			stk[i] = cStrikes[i];
		}
		
		CVector result;
		res.resize(cStrikes.size(),1);
		effLv.calcOptions(result,fguess,stk);
		
		for ( int i = 0 ; i < stk.getsize(); i++ ){
			res[i][0] = result[i];
		}
	}
	
	CParameterMap pm;
	pm.SetValue(res);
	return pm.GetValue(pResult);
	
	
	end_function
}	


/****************************************************************
**	 
**	Routine: MLCalibLocalVol
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLCalibrateEffectiveLocalVolFull(VARIANT asset,VARIANT startDate,VARIANT endDates,VARIANT quadraticInterpolators,VARIANT calibStrikes,VARIANT calibSpots,VARIANT initialGuessOpt,VARIANT nGauss, VARIANT usePrevResultAsGuess,VARIANT calcStrikesOpt,VARIANT idateOpt,
												   VARIANT usePdeOpt,VARIANT lowSpotOpt,VARIANT highSpotOpt,VARIANT nxOpt,VARIANT ntOpt,
												   VARIANT* pResult)
{	

  
	begin_function	
		
	map_object_parameter(asset, Asset, hasset);
	map_object_vector_parameter(quadraticInterpolators,Interpolator,hinterpolators);
	map_object_parameter(startDate, Date, startD);
	map_parameter(endDates, GVector< long > ,enddate );
	map_optional_parameter(initialGuessOpt,CMatrix, fguess,CMatrix());
	map_optional_parameter(nGauss,long, ngauss,100);
	map_strikes_vector_parameter(calibStrikes, hStrikes);
	map_strikes_vector_parameter(calibSpots, hCStrikes);
	map_optional_parameter(idateOpt, long  ,idate,0 );
	map_optional_object_parameter(calcStrikesOpt, Strikes, hCalcStrike);
	map_parameter(usePrevResultAsGuess, bool ,usePrev );

	map_optional_parameter(usePdeOpt, bool ,usePde,false );
	map_optional_parameter(lowSpotOpt, double ,lowSpot,-1e99 );
	map_optional_parameter(highSpotOpt, double ,highSpot,-1e99 );
	map_optional_parameter(nxOpt, double ,nx,-1e99 );
	map_optional_parameter(ntOpt, double ,nt,-1e99 );

	if ( usePde )
	{
		if ( fabs(lowSpot+1e99) < 1e-11 || fabs(highSpot+1e99)< 1e-11 ||
			 fabs(nx+1e99) < 1e-11   || fabs(nt+1e99) < 1e-11 )
		{
			throw("pde set up not complete");
		}
	}


	int nStrikes ;

	if (  hCalcStrike == NULL ){
		nStrikes = 0;
	}
	else{
		nStrikes = hCalcStrike->m_hStrikes.size();

	}

	GVector< MlEqStrikeHandle > calcStrikes(nStrikes);
	for ( int i = 0 ; i < nStrikes; i++ ){
		calcStrikes[i] = hCalcStrike->m_hStrikes[i];
	}

	
	calibrateEffectiveLocalVol effLv;

	GVector <MlEqAnalyticCurveWithTanhWingEdgeHandle > quadraticInterpolator(hinterpolators.size());

	for ( int i = 0 ; i < quadraticInterpolator.getsize(); i++ )
	{
		quadraticInterpolator[i] = dynamic_cast<MlEqAnalyticCurveWithTanhWingEdge*>(&*hinterpolators[i]);		
		if (!quadraticInterpolator[i]) throw "Bad horsey interpolator type!";
	}
	

	int n = hStrikes.size();
	GVector < GVector  < MlEqStrikeHandle > > calibStk;
	calibStk.resize(hStrikes.size());
	for ( int i = 0 ; i < n; i++ )
	{
		calibStk[i].resize(hStrikes[i].size());
		for ( int j = 0 ; j < hStrikes[i].size(); j++ ){
			calibStk[i][j] = hStrikes[i][j];
		}
	}

	n = hCStrikes.size();
	GVector < GVector  < MlEqStrikeHandle > > calibSpots;
	calibSpots.resize(hCStrikes.size());
	for ( int i = 0 ; i < n; i++ )
	{
		calibSpots[i].resize(hCStrikes[i].size());
		for ( int j = 0 ; j < hCStrikes[i].size(); j++ ){
			calibSpots[i][j] = hCStrikes[i][j];
		}
	}


	CMatrix result;

	if ( !nStrikes )
	{
		if ( !usePde )
		{
			calibrateEffectiveLocalVolFull(result,*hasset,quadraticInterpolator,startD,
									   enddate,fguess,ngauss,usePrev,
									   calibStk,calibSpots);
		}
		else
		{
			 calibrateEffectiveLocalVolGridFull(result,*hasset, quadraticInterpolator,
									startD,enddate,fguess,usePrev,
									calibStk,calibSpots,
									lowSpot,highSpot,nx,nt);

		}
	}
	else
	{

		calibrateEffectiveLocalVol effLv;

		if ( idate > quadraticInterpolator.getsize()-1 ){
			throw("decrease idate or U're in trouble");
		}

		CVector targetvols(calibStk[0].getsize());

		if ( idate < 0 ){
			throw("idate must be positive");
		}

		startD->PutDate(enddate[idate]);

		CVector guess;

		effLv.init(*hasset,quadraticInterpolator[idate],startD,enddate[idate+1],
				   calibStk[0],calibSpots[0],targetvols,guess,ngauss);


		
		result.resize(calcStrikes.getsize(),1);

		CVector res(calcStrikes.getsize());

		const CVector& coeff = quadraticInterpolator[idate]->getCoeff();
		guess = coeff;
		effLv.calcOptions(res,guess,calcStrikes);
		
		for ( int i = 0 ; i < calcStrikes.getsize(); i++ ){
			result[i][0] = res[i];
		}


	}


	

	
	CParameterMap pm;
	pm.SetValue(result);
	return pm.GetValue(pResult);
	
	
	end_function
}	


/****************************************************************
**	 
**	Routine: MLStrikeInTermsOfStandartDeviations
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLStrikeInTermsOfStandartDeviations(VARIANT VolatilityStructureHandle, VARIANT StrikeHandle, VARIANT VolDate, VARIANT Forward, VARIANT BidOfferOpt, VARIANT* pResult)
{		
	int									nStrikes;
	vector<MlEqStrikeHandle>			m_pStrikes;
	
	// map the inputs
	begin_function
	map_object_parameter(VolatilityStructureHandle, VolatilityStructure, hVolatilityStructure);
	map_object_parameter(StrikeHandle, Strikes, hStrike);
	map_parameter(VolDate, double, fVolDate);
	map_optional_enum_parameter(BidOfferOpt, BidAskMidEnum, nBidOffer, Middle);
	map_parameter(Forward, double, fForward);
		
	// do the work
	m_pStrikes = hStrike->m_hStrikes;
	nStrikes = hStrike->m_hStrikes.size();
	CVector vector(nStrikes);
	for (int i = 0; i < nStrikes; i++){
		vector[i] = hVolatilityStructure->StrikeInTermsOfStandartDeviations(*(hStrike->m_hStrikes[i]), fVolDate, fForward, nBidOffer);
	}
	
	// return the result
	CParameterMap pm;
	pm.SetValue(vector);
	return pm.GetValue(pResult);	
	end_function
}



/****************************************************************
**	 
**	Routine: MLGetDensity
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/



HRESULT __stdcall MLGetDensity(VARIANT VolatilityStructureHandle, VARIANT StrikeHandle, VARIANT VolDate, VARIANT Forward, VARIANT BidOfferOpt,VARIANT CumulativeFlag, VARIANT* pResult)
{		
	int									nStrikes;
	vector<MlEqStrikeHandle>			m_hStrikes;
	
	// map the inputs
	begin_function
	map_object_parameter(VolatilityStructureHandle, VolatilityStructure, hVolatilityStructure);
	map_object_parameter(StrikeHandle, Strikes, hStrike);
	map_parameter(VolDate, double, fVolDate);
	map_optional_enum_parameter(BidOfferOpt, BidAskMidEnum, nBidOffer, Middle);
	map_parameter(Forward, double, fForward);
	map_optional_parameter(CumulativeFlag, long, ncummulative, 0);
	
	// do the work
	m_hStrikes = hStrike->m_hStrikes;
	nStrikes = hStrike->m_hStrikes.size();
	CVector vector(nStrikes);
	for (int i = 0; i < nStrikes; i++){
		vector[i] = hVolatilityStructure->getDensity(*(hStrike->m_hStrikes[i]), fVolDate, fForward, nBidOffer, ncummulative);
	}
	
	// return the result
	CParameterMap pm;
	pm.SetValue(vector);
	return pm.GetValue(pResult);	
	end_function
}



/****************************************************************
**	 
**	Routine: MLBSImpliedVolatility
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/



HRESULT __stdcall MLBSImpliedVolatility(VARIANT OptionPrice, VARIANT Forward, VARIANT Strike, VARIANT Maturity, VARIANT DiscountFactor, VARIANT CallPut, VARIANT BlendOpt, VARIANT AccuracyOpt, VARIANT* pResult)
{	
	int									rootfind_flag = 8;
	double								lower_x = 0.0;	
	double								upper_x	= 1.5;	

	begin_function
	map_parameter(OptionPrice, double, fOptionPrice);
	map_parameter(Forward, double, fForward);
	map_parameter(Strike, double, fStrike);
	map_parameter(Maturity, double, fMaturity);
	map_parameter(DiscountFactor, double, fDiscountFactor);
	map_parameter(CallPut, long, nCallPut);
	map_optional_parameter(BlendOpt, double, fBlend, 0.0);
	map_optional_parameter(AccuracyOpt, double, fAccuracy, 0.00001);

	double res = MlEqBSImpliedVol(fOptionPrice,	
								  fForward,			
								  fMaturity,		
								  fStrike,			
								  fDiscountFactor,
								  nCallPut,			
								  rootfind_flag,		
								  fAccuracy,		
								  lower_x,		
								  upper_x);

	return CComVariant(res).Detach(pResult);
	end_function
}




/****************************************************************
**	 
**	Routine: MLHermitSkew
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/



HRESULT __stdcall MLHermitSkew( VARIANT Forward, VARIANT Maturity,VARIANT Strikes , VARIANT HermitCoeff,VARIANT ElasticityOpt,VARIANT VolVolOpt,VARIANT NGaussOpt, VARIANT* pResult)
{	

	vector< MlEqStrikeHandle >			strikes;

	begin_function

	map_parameter(Forward, double, fForward);
	map_parameter(Maturity, double, fMaturity);
	map_parameter(HermitCoeff, CVector, fHermitCoeff);
	map_object_parameter(Strikes, Strikes, hStrike);

	map_optional_parameter(VolVolOpt, double, volvol, 0.0);
	map_optional_parameter(ElasticityOpt, double, elasticity, 0.0);
	map_optional_parameter(NGaussOpt, long, ngauss, 0);


	strikes = hStrike->m_hStrikes;
	int nstrikes = strikes.size();

	MlEqStrike stk;
	CVector cstrikes(nstrikes);

	for ( int i= 0 ; i < nstrikes; i++ )
	{
		MlEqStrike::convertStrikes(stk,*strikes[i]);
		cstrikes[i] = stk.m_strike;
	}

	CVector result(nstrikes);

	int returnVolFlag = 1;

	int ntimesteps		= 5; 
	int npaths			= 10;
	int nspacialPoints	= 100;
	int nstdev			= 6;

//	HermiteOptions(result, cstrikes,fHermitCoeff,fMaturity,fForward,returnVolFlag,elasticity,volvol,ngauss);
	HermiteOptionsNew(result,cstrikes,fHermitCoeff,fMaturity,fForward,returnVolFlag,
		ntimesteps,npaths,nspacialPoints,nstdev);

	// return the result
	CParameterMap pm;
	pm.SetValue(result);
	return pm.GetValue(pResult);	
	end_function

}	
	
	

/****************************************************************
**	 
**	Routine: MLHermitSkew
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/	
	
	
HRESULT __stdcall MLCalibHermiteHermitSkew( VARIANT Forward, VARIANT Maturity,VARIANT CalibStrikes , VARIANT HermitCoeff,VARIANT calibVols,
										    VARIANT DateToDouble,VARIANT ntimeSteps,VARIANT npaths,VARIANT nx,VARIANT nstdev,VARIANT VarSwapVolOpt, VARIANT* pResult)
{	
	
	vector< MlEqStrikeHandle >			strikes;
	

	begin_function

	map_parameter(Forward, double, fForward);
	map_parameter(Maturity, long, fMaturity);
	map_parameter(HermitCoeff, CVector, fHermitCoeff);
	map_object_parameter(CalibStrikes, Strikes, hStrike);
	map_object_parameter(DateToDouble, Date, dateToDouble);
	map_parameter(calibVols, CVector, hvols);

	map_parameter(ntimeSteps, long, nt);
	map_parameter(npaths, long, mnpaths);
	map_parameter(nx, long, mnx);
	map_parameter(nstdev, long, mnstdev);
	map_optional_parameter(VarSwapVolOpt, CVector, fvarvol,CVector());


	strikes = hStrike->m_hStrikes;
	GVector< MlEqStrikeHandle > CalibStrikes(strikes.size());

	for ( int i = 0 ; i < strikes.size(); i++ ){
		CalibStrikes[i] = strikes[i];
	}


	
	calibrateHermiteOptions cOptions;


	cOptions.initialize(CalibStrikes,
						hvols,fForward,
										  dateToDouble,
										  fMaturity,
										  fvarvol,
										  nt,
										  mnpaths,
										  mnx,
										  mnstdev,
										  fHermitCoeff);


	CVector result;

	cOptions.calibrate(result);



	// return the result
	CParameterMap pm;
	pm.SetValue(result);
	return pm.GetValue(pResult);	
	end_function

}


/****************************************************************
**	 
**	Routine: MLHermitSkew
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/



HRESULT __stdcall MLHermitVolOptions( VARIANT Forward, VARIANT Maturity,VARIANT Strikes , VARIANT HermitCoeff,VARIANT ElasticityOpt,VARIANT VolVolOpt,VARIANT NGaussOpt, VARIANT* pResult)
{
	


	vector< MlEqStrikeHandle >			strikes;

	begin_function

	map_parameter(Forward, double, fForward);
	map_parameter(Maturity, double, fMaturity);
	map_parameter(HermitCoeff, CVector, fHermitCoeff);
	map_object_parameter(Strikes, Strikes, hStrike);

	map_optional_parameter(VolVolOpt, double, volvol, 0.0);
	map_optional_parameter(ElasticityOpt, double, elasticity, 0.0);
	map_optional_parameter(NGaussOpt, long, ngauss, 0);

	strikes = hStrike->m_hStrikes;
	int nstrikes = strikes.size();

	MlEqStrike stk;
	CVector cstrikes(nstrikes);

	for ( int i= 0 ; i < nstrikes; i++ )
	{
		MlEqStrike::convertStrikes(stk,*strikes[i]);
		cstrikes[i] = stk.m_strike;
	}

	CMatrix result;

	int returnVolFlag = 1;

	int ntimesteps		= 5; 
	int npaths			= 10;
	int nspacialPoints	= 100;
	int nstdev			= 6;

	VolOptions(result,cstrikes,fHermitCoeff,fMaturity,fForward,returnVolFlag,
		ntimesteps,npaths,nspacialPoints,nstdev);


	// return the result
	CParameterMap pm;
	pm.SetValue(result);
	return pm.GetValue(pResult);	
	end_function

}



/****************************************************************
**	 
**	Routine: MLHermitSkew
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/



HRESULT __stdcall MLAdjustHermitSkew( VARIANT Forward, VARIANT Maturity,VARIANT Strikes , VARIANT HermitCoeff,VARIANT TargetSlope, VARIANT* pResult)
{	

	vector< MlEqStrikeHandle >			strikes;

	begin_function

	map_parameter(Forward, double, fForward);
	map_parameter(Maturity, double, fMaturity);
	map_parameter(HermitCoeff, CVector, fHermitCoeff);
	map_object_parameter(Strikes, Strikes, hStrike);
	map_parameter(TargetSlope, double, fTargetSlope);

	strikes = hStrike->m_hStrikes;
	int nstrikes = strikes.size();

	MlEqStrike stk;
	CVector cstrikes(nstrikes);

	for ( int i= 0 ; i < nstrikes; i++ )
	{
		MlEqStrike::convertStrikes(stk,*strikes[i]);
		cstrikes[i] = stk.m_strike;
	}


	double eps = 0.05;
	hermiteSkew hSk;
	hSk.initialize(fMaturity,fForward,fHermitCoeff);
	hSk.adjustSlope(fTargetSlope,fForward,eps);


	// return the result
	CParameterMap pm;
	pm.SetValue(hSk.m_hermitCoeff);
	return pm.GetValue(pResult);	
	end_function

}


/****************************************************************
**	 
**	Routine: MLCalculateOptions
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLCalculateOptions(VARIANT ForwardSkewHandle, VARIANT Strikes, VARIANT startDateIndex,VARIANT endDateIndex, VARIANT ImpliedVolFlag,VARIANT Asset, VARIANT* pResult)
{		
	vector<vector<MlEqStrikeHandle> >	pStrikes;
	vector< MlEqStrikeHandle >			strikes;

	int									i;	

	// map the inputs
	begin_function			

	map_object_parameter(ForwardSkewHandle, MonteCarlo, pIn)
	map_object_vector_parameter(Strikes, Strikes, hStrikes)
	map_parameter(endDateIndex, CVector, jDates)
	map_parameter(startDateIndex, CVector, iDates)
	map_optional_parameter(ImpliedVolFlag, long, impliedVolFlag,1L)
	map_optional_object_parameter(Asset, Asset, hAsset)


	if ( iDates.getsize() != jDates.getsize() ){
		throw("size of startDateIndex and EndDateIndex must coincide");
	}

	for ( i = 0 ; i < iDates.getsize(); i++ ){
		if ( iDates[i] >= jDates[i] ){
			throw("startDateIndex must be smaller than endDateIndex");
		}
	}

	pStrikes.resize(hStrikes.size());
	for (i = 0 ; i < hStrikes.size(); i++){
		pStrikes[i] = hStrikes[i]->m_hStrikes;
	}	

	CMatrix		results(iDates.getsize(), pStrikes[0].size());
	int			controlVarFlag = 0;

	if ( hAsset == NULL && pIn->m_nAssets > 1 ){
		throw("optional argument assets must be set in this case");
	}


	for (i = 0; i < iDates.getsize(); i++){

		strikes = getObjectFromVector(pStrikes,i);

		if ( hAsset == NULL ){
			pIn->calculateOptions(results[i], iDates[i],jDates[i], strikes, controlVarFlag, impliedVolFlag);
		}
		else{
			pIn->calculateOptions(results[i], hAsset,iDates[i],jDates[i], strikes, controlVarFlag, impliedVolFlag);																											
		}

	}
	
	// return the result
	CParameterMap pm;
	pm.SetValue(results);
	return pm.GetValue(pResult);
	end_function

}



/****************************************************************
**	 
**	Routine: MLTestHybridModel
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLTestHybridModel(VARIANT ForwardSkewHybridHandle, VARIANT* pResult)
{		

	int	i;	

	// map the inputs
	begin_function			

	map_object_parameter(ForwardSkewHybridHandle, MonteCarlo, hMonteCarlo)
	
	const std::vector<long>& mcDates = 	hMonteCarlo->GetMCDates() ;
	const vector<MlEqAssetHandle>&	hAssets = hMonteCarlo->getAssets();

	int nassets = hAssets.size();

	if ( nassets != 1 ){
		throw("only single asset hybrid supported at present");
	}


	MlEqStrikeHandle stk = new MlEqStrike(0.0);

	vector< MlEqStrikeHandle >  Strikes(1);
	Strikes[0] = stk;

	int ndates = mcDates.size();

	CMatrix res(ndates,5);
	double z,disc;
	CVector result(1);

	MlEqAssetHandle h = hAssets[0];

	long nToday = h->GetCurrentDate();
	double spot = h->GetSpot(nToday);


	for (i = 0; i < ndates; i++)
	{
		res[i][0]	= mcDates[i];
		res[i][1]	= hAssets[0]->GetForward(mcDates[i],false);

		hMonteCarlo->calculateOptions(result, h,0,i, Strikes, 0, 0);

		disc		= hMonteCarlo->GetDiscount(0,i,0);

		res[i][2]	= result[0]*spot;///disc;
		res[i][3]	= disc;
 

		z = 0.0;
		for ( int ipath = 0 ; ipath < hMonteCarlo->m_nPaths; ipath++ ){
			z += hMonteCarlo->GetStochasticDiscount(ipath,i);
		}
		z /= hMonteCarlo->m_nPaths;

		res[i][4]	= z;

	}

	
	// return the result
	CParameterMap pm;
	pm.SetValue(res);
	return pm.GetValue(pResult);
	end_function

}


/****************************************************************
**	 
**	Routine: MLCalculateOptions
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/



HRESULT __stdcall MLTestBarrierBridge(VARIANT ForwardSkewHandle,VARIANT Barrier,VARIANT Rebate, VARIANT* pResult)
{		

	int	ipath;	

	begin_function			

	map_object_parameter(ForwardSkewHandle, MonteCarlo, hForwardSkewMonteCarlo)
	map_parameter(Barrier, double, barrier)
	map_parameter(Rebate, double, rebate)

	int ndates = hForwardSkewMonteCarlo->m_nDates;
	CMatrix pathArray(hForwardSkewMonteCarlo->m_nAssets,ndates+1);

	int idate = 0;
	double L = 0.01;
	double hitProb;

	double payoff = 0.0;
	double prob = 0.0;
	for (ipath = 0; ipath < hForwardSkewMonteCarlo->m_nPaths; ipath++)
	{
		pathArray = hForwardSkewMonteCarlo->GetPathArray( ipath);
		hitProb = hForwardSkewMonteCarlo->GetBrownianBridgeHitProb(ipath, idate, barrier, L);
		payoff += hitProb*pathArray[0][ndates];
		prob += hitProb;
	}
	
	prob /= hForwardSkewMonteCarlo->m_nPaths;
	payoff /= hForwardSkewMonteCarlo->m_nPaths;
	payoff /= prob;
	

	// return the result
	CParameterMap pm;
	pm.SetValue(payoff);
	return pm.GetValue(pResult);

	end_function

}

/****************************************************************
**	 
**	Routine: MLCreateForwardVolatilities
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLCreateForwardVolatilities(VARIANT MarketVolatilityStructure, VARIANT AsymptoticVolatilityStructure, VARIANT DecayRate, VARIANT *pResult)
{
	CComPtr<IVolatilityStructure>		spOut;							// output volatility structure	
	CComBSTR							sHandle;						// output handle
		
	begin_function
	map_object_parameter(MarketVolatilityStructure, VolatilityStructure, hMarket);
	map_object_parameter(AsymptoticVolatilityStructure, VolatilityStructure, hAsymptotic);
	map_parameter(DecayRate, double, fDecayRate);
														
	// create the output volatility structure
	if (spOut.CoCreateInstance(CLSID_VolatilityStructure)) throw "Could not create output volatility structure";			
	map_com_to_analytic(spOut, VolatilityStructure, hOut);		
	hOut->initialize(*hMarket, *hAsymptotic, fDecayRate);
	
	// return a handle to the created volatility structure
	if (_Module.GetSiriusApplication()->InsertObject(spOut, VARIANT_TRUE, &sHandle)) propagate_error;
	return CComVariant(sHandle).Detach(pResult);
	end_function
}



/****************************************************************
**	 
**	Routine: MLGetPV01
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLGetPV01(VARIANT Swap, VARIANT* pResult)
{	
	double								fRate = 0.0;	
	
	begin_function	
	map_com_object_parameter(Swap, Product, spProduct);
					
	CComPtr<IDispatch>					spProductType;
	if (spProduct->GetObject(&spProductType)) propagate_error;
	CComQIPtr<ISwap>					spSwap(spProductType);	
	if (!spSwap) throw "Invalid swap handle returned from product";		
	throw "MLGetPV01 is not supported in this release";
	return CComVariant(fRate).Detach(pResult);	
	end_function
}


/****************************************************************
**	 
**	Routine: MLGetSwapRate
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLGetSwapRate(VARIANT Swap, VARIANT *pResult)
{		
	double								fSwapRate = 0.0;

	begin_function
	map_com_object_parameter(Swap, Product, spProduct);
	
	CComPtr<IDispatch>					spProductType;
	if (spProduct->GetObject(&spProductType)) propagate_error;
	CComQIPtr<ISwap>					spSwap(spProductType);
	if (!spSwap) throw "Invalid swap handle returned from product";		
	throw "MLGetSwapRate is not supported in this release";
	return CComVariant(fSwapRate).Detach(pResult);
	end_function
}


/****************************************************************
**	 
**	Routine: MLInterpolatorGetValue
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLInterpolatorGetValue(VARIANT InterpolatorHandle, VARIANT xValue,VARIANT whichDataOpt, VARIANT* pResult)
{	
	double	f_yValue;
	
	// map the inputs
	begin_function
	map_com_object_parameter(InterpolatorHandle, Interpolator, spInterpolator)
	map_parameter(xValue, double, f_xValue);
	map_optional_parameter(whichDataOpt, long, which,0);

	// do the work
	if (spInterpolator->GetValue(f_xValue, which, &f_yValue)) propagate_error;
	
	// return the value
	return CComVariant(f_yValue).Detach(pResult);
	end_function
}


/****************************************************************
**	 
**	Routine: MLNormalVal
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLNormalVal(VARIANT x, VARIANT* pResult)
{		
	begin_function
	map_parameter(x, double, fx);
	return CComVariant(normal(fx)).Detach(pResult);
	end_function
}	



/****************************************************************
**	 
**	Routine: MLInverseNormalVal
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


	
HRESULT __stdcall MLInverseNormalVal(VARIANT x, VARIANT* pResult)
{		
	begin_function
	map_parameter(x, double, fx);
	return CComVariant(normal_inv(fx)).Detach(pResult);
	end_function
}	
/****************************************************************
**	 
**	Routine: MLSort
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


	
HRESULT __stdcall MLSort(VARIANT x, VARIANT* pResult)
{		
	begin_function
	map_parameter(x, CVector, fx);

	int n = fx.getsize();
	GVector<long> map(n);
	hpsort(n, fx.getPtr()-1,map.getPtr()-1);

	CMatrix res(n,2);
	for ( int i = 0 ; i < n; i++ ){
		res[i][0] = fx[i];
		res[i][1] = map[i];
	}

	CParameterMap pm;
	pm.SetValue(res);
	return pm.GetValue(pResult);

	end_function
}	
	


/****************************************************************
**	 
**	Routine: MLConvertStrike
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


	
HRESULT __stdcall MLConvertStrike(VARIANT OutputStrike,VARIANT InputStrike, VARIANT *pResult)
{		
	begin_function
	
	map_object_parameter(InputStrike, Strikes, hInputStrike)
	map_object_parameter(OutputStrike, Strikes, hOutputStrike)
	
	if ( hInputStrike->m_hStrikes.size() != hOutputStrike->m_hStrikes.size() ){
		throw("number of input strikes must coincide with number of output strikes");
	}

	for ( int i = 0 ; i < hInputStrike->m_hStrikes.size(); i++ ){
		MlEqStrike::convertStrikes(*hOutputStrike->m_hStrikes[i],*hInputStrike->m_hStrikes[i]);
	}

	return CComVariant("OK").Detach(pResult);
	
	end_function
}		

/****************************************************************
**	 
**	Routine: MLHullWhite
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLHullWhite(VARIANT CalibrationDates,VARIANT DateToDouble,VARIANT RateSettingDate,VARIANT RateMaturityDate, VARIANT *pResult)
{	
	begin_function
	map_parameter(CalibrationDates, CVector, calibDates);
	map_object_parameter(DateToDouble, Date, dateToDouble)
	map_parameter(RateSettingDate, long, t);
	map_parameter(RateMaturityDate, long, T);

	MlEqHullAndWhite HW(dateToDouble);
	HW.initialize(calibDates, t, T);

	double res = 1.0;
	return CComVariant(res).Detach(pResult);
	end_function
}	
	
	
/****************************************************************
**	 
**	Routine: ThreeMomentOptionCalculator
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/
	
	
	
HRESULT __stdcall ThreeMomentOptionCalculator(VARIANT Asset, VARIANT AsianWeights, VARIANT AsianDates, VARIANT CorrelationMatrix, VARIANT StrikeHandle,VARIANT DateHandle, VARIANT *pResult)
{		
	begin_function
	map_object_parameter(Asset, Asset, hAsset)
	map_parameter(AsianWeights, CVector, asianWeights)	
	map_parameter(AsianDates, std::vector<DATE>, asianDates)
	map_object_parameter(CorrelationMatrix, CorrelationMatrix, hCorrelationMatrix);			
	map_object_parameter(StrikeHandle, Strikes, hStrike)
	map_object_parameter(DateHandle, Date, dateToDouble)	
		
	MlEqThreeMoment						tM;
	tM.initialize(*hAsset, asianWeights, asianDates, *hCorrelationMatrix, dateToDouble, hStrike->m_hStrikes);
	CVector moments = tM.calculate();
	
	end_function
}		
	
	
/****************************************************************
**	 
**	Routine: ThreeMomentOptionCalculator
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/
	
	

HRESULT __stdcall MLCalcAsianMax(VARIANT Assets, VARIANT PayoffDates, VARIANT CallPut, VARIANT Strikes, VARIANT StrikeFixings,VARIANT SpotFixings,VARIANT Correl,VARIANT nPaths,VARIANT randomNumberFlag, VARIANT *pResult)
{		
	begin_function

	map_object_vector_parameter(Assets, Asset, vAssets);
	map_parameter(CallPut, double, callPut);
	map_parameter(Strikes, CVector, strikes);
	map_parameter(PayoffDates, GVector<long>, payoffDates);
	map_parameter(Correl, CMatrix, correl);
	map_parameter(nPaths, long, npath);
	map_parameter(randomNumberFlag, long, nrandomNumberFlag);
	map_parameter(SpotFixings, CMatrix, fixings);
	map_parameter(StrikeFixings, CVector, strikeFixings);

	

	CAsianMax asianMax;

	asianMax.initialize(callPut,strikes,payoffDates,strikeFixings);

	

	MlEqVolatilityStructureHandle pVol =  vAssets[0]->GetVolatilityStructure();
	MlEqConstDateHandle hToday =  pVol->getDateToDouble();


	CSequentialBlackMC mc(hToday);


	int nassets = vAssets.size();
	int ndates  = payoffDates.getsize();

	mc.initialize(vAssets, asianMax, fixings, correl,npath,nrandomNumberFlag);

//	mc.initialize(CParameterMap::StdVectorToGVector(vAssets), asianMax, fixings, correl,npath);

	CMatrix	results;		
	mc.simulate(results,asianMax);
	
	CParameterMap pm;
	pm.SetValue(results);
	return pm.GetValue(pResult);

	end_function
}		
		

		
/****************************************************************
**	 
**	Routine: CalibrateAssetSkew
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/



HRESULT __stdcall CalibrateAssetSkew(VARIANT AssetHandle, VARIANT Maturity, VARIANT DateHandle, VARIANT Strikes, VARIANT *pResult)
{
	begin_function
	map_object_parameter(AssetHandle, Asset, hAsset)			
	map_parameter(Maturity, long, nMaturity)
	map_object_parameter(DateHandle, Date, hDate)

	vector<vector<MlEqStrikeHandle> >	pStrikes;
	vector< MlEqStrikeHandle >			strikes;

	map_object_vector_parameter(Strikes, Strikes, hStrikes)

	pStrikes.resize(hStrikes.size());
	for (int i = 0 ; i < hStrikes.size(); i++){
		pStrikes[i] = hStrikes[i]->m_hStrikes;
	}	


	CVector vector = hAsset->Calibrate((DATE)nMaturity, *hDate, pStrikes);
	CParameterMap pm;
	pm.SetValue(vector);
	return pm.GetValue(pResult);
	end_function
}



int func (const CVector& x, void* vp, int calc_for_grad, CVector& f)
{

	f[0] = x[2]+x[0]*x[0]+3.0*x[1]-7.0*x[0]*x[1];
	f[1] = x[0]+x[0]*x[1]*x[2]-x[2];


	return 0.0;
}

extern void rootfind_underdetermined_solve(const CVector& initial_x, const CVector& bump_x, const CVector& tolerances, 
								   p_dev_underdet_func fn, void* vp, int max_tries, 
								   int max_restarts, const CMatrix& weights, CVector& found_x);

/****************************************************************
**		 
**	Routine: MlTestUnderdetermined
**	Returns: VARIANT
**	Action : 
**		         
****************************************************************/
	
HRESULT __stdcall MlTestUnderdetermined(VARIANT initialGuess, VARIANT bump, VARIANT Tolerance,VARIANT *pResult)
{	
	begin_function
	map_parameter(initialGuess, CVector, initial_x)			
	map_parameter(bump, CVector, bump_x)			
	map_parameter(Tolerance, CVector, tol)			
	
	int max_tries = 500;
	int max_restarts = 50;
	
	CMatrix weights;
	CVector found_x;
	
	rootfind_underdetermined_solve(initial_x, bump_x, tol, 
								   func, NULL, max_tries, 
								   max_restarts, weights,found_x);
	
	CParameterMap pm;
	pm.SetValue(found_x);
	return pm.GetValue(pResult);
	end_function
}


/****************************************************************
**	 
**	Routine: MLInterpolatorGetValue
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall ML2DimInterpolatorGetValue(VARIANT InterpolatorHandle, VARIANT xValue,VARIANT yValue, VARIANT* pResult)
{	
	begin_function
	
	// map the inputs
	map_object_parameter(InterpolatorHandle, Interpolator, hInterpolator)
	map_parameter(xValue, double, f_xValue);
	map_parameter(yValue, double, f_yValue);

	// do the work
	MlEq2DInterpolator* p = dynamic_cast<MlEq2DInterpolator*>(&*hInterpolator);
	if (!p) throw "You need to pass in a 2D interpolator to use this function";
	double f = p->getValue(f_xValue, f_yValue);

	// return the value
	return CComVariant(f).Detach(pResult);
	end_function
}


 

/****************************************************************
**		 
**	Routine: CapletPrice
**	Returns: VARIANT
**	Action : 
**		         
****************************************************************/
	
HRESULT __stdcall MLCapletPrice(VARIANT HullWhite,VARIANT RateSettingDate,VARIANT RateMaturityDate,VARIANT strike,VARIANT dayCount, VARIANT *pResult)
{	
	begin_function

	map_object_parameter(HullWhite, HullAndWhite, hHullAndWhite);	
	map_parameter(RateSettingDate, long, setdate)
	map_parameter(RateMaturityDate, long, matdate)
	map_parameter(strike, double, xstrike)	
	map_enum_parameter(dayCount, DayCountConventionEnum, dc)

	double res = hHullAndWhite->CapletPrice(setdate, matdate, xstrike, dc);

	CParameterMap pm;
	pm.SetValue(res);
	return pm.GetValue(pResult);
	end_function
}

/****************************************************************
**		 
**	Routine: SwaptionPrice
**	Returns: VARIANT
**	Action : 
**		         
****************************************************************/
	
HRESULT __stdcall MLSwaptionPrice(VARIANT HullWhite,VARIANT SwapFixingDates,VARIANT swaptionStrike,VARIANT dayCount, VARIANT callPut, VARIANT *pResult)
{	
	begin_function
	map_object_parameter(HullWhite, HullAndWhite, hHullAndWhite);	
	map_parameter(SwapFixingDates, GVector<long>, fixDates)
	map_parameter(swaptionStrike, double, xstrike)
	map_parameter(callPut, long, cp)
	map_enum_parameter(dayCount, DayCountConventionEnum, dc)

	double res = hHullAndWhite->SwaptionPrice(fixDates, xstrike, dc, cp);

	return CComVariant(res).Detach(pResult);
	end_function
}

/****************************************************************
**		 
**	Routine: MLDisplayHullWhiteParameters
**	Returns: VARIANT
**	Action : 
**		         
****************************************************************/
	
HRESULT __stdcall MLDisplayHullWhiteParameters(VARIANT HullWhite, VARIANT betaOrGammaFlag, VARIANT *pResult)
{	
	begin_function
	map_object_parameter(HullWhite, HullAndWhite, hHullAndWhite);	
	map_parameter(betaOrGammaFlag, long, flag)			

	CVector res;
	if ( flag == 0 ){
		res = hHullAndWhite->getBetas();
	}else{
		res = hHullAndWhite->getGammas();
	}

	CParameterMap pm;
	pm.SetValue(res);
	return pm.GetValue(pResult);

	end_function
}


/****************************************************************
**		 
**	Routine: MLBlackScholesCapletPrice
**	Returns: VARIANT
**	Action : 
**		         
****************************************************************/
	
HRESULT __stdcall MLBlackScholesCapletPrice(VARIANT ZeroCurve, VARIANT RateSettingDate, VARIANT RateMaturityDate, VARIANT strike, VARIANT vol, VARIANT dayCount, VARIANT callPut, VARIANT *pResult)
{
	begin_function

	map_object_parameter(ZeroCurve, ZeroCurve, hZeroCurve)
	map_parameter(RateSettingDate, long, setdate)
	map_parameter(RateMaturityDate, long, matdate)
	map_parameter(strike, double, xstrike)
	map_parameter(vol, double, xvol)
	map_enum_parameter(dayCount, DayCountConventionEnum, dc)
	map_parameter(callPut, long, cp)			
		
	long valdate = hZeroCurve->GetReferenceDate();
	double res = BlackScholesCapletPricer(valdate, setdate, matdate, xstrike, xvol, dc, *hZeroCurve, cp);

	CParameterMap pm;
	pm.SetValue(res);
	return pm.GetValue(pResult);
	end_function
}

/****************************************************************
**		 
**	Routine: MLBlackScholesSwaptionPricer
**	Returns: VARIANT
**	Action : 
**		         
****************************************************************/
	
HRESULT __stdcall MLBlackScholesSwaptionPricer(VARIANT ZeroCurve, VARIANT SwapFixingDates, VARIANT swaptionStrike, VARIANT vol, VARIANT dayCount, VARIANT callPut, VARIANT *pResult)
{
	begin_function

	map_object_parameter(ZeroCurve, ZeroCurve, hZeroCurve)
	map_parameter(SwapFixingDates, GVector<long>, fixings)			
	map_parameter(swaptionStrike, double, xstrike)			
	map_parameter(vol, double, xvol)			
	map_enum_parameter(dayCount, DayCountConventionEnum, dc)
	map_parameter(callPut, long, cp)			
	
	long valdate = hZeroCurve->GetReferenceDate();	
	double res = BlackScholesSwaptionPricer(valdate, fixings, xstrike, xvol, dc, *hZeroCurve, cp);
	CParameterMap pm;
	pm.SetValue(res);
	return pm.GetValue(pResult);
	end_function
}

/****************************************************************
**	 
**	Routine: CalibrateAssetSkew
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/



/****************************************************************
**	 
**	Routine: MlMathARTrackingPursuit
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/



HRESULT __stdcall MlMathARTrackingPursuit(VARIANT TargetVector, VARIANT ProjectionVectors, VARIANT StoppingTolerance, VARIANT Epsilon,VARIANT MaxIter,VARIANT *pResult)
{
	begin_function
	map_parameter(TargetVector, CVector, targetVec)			
	map_parameter(ProjectionVectors, CMatrix, projectionVecs)			
	map_parameter(Epsilon, double, epsilon)			
	map_parameter(MaxIter, long, maxIter)			
	map_parameter(StoppingTolerance, double, stoppingTol)			

	ARTrackingPursue arTrack;
	arTrack.initialize(targetVec,projectionVecs,epsilon,stoppingTol,maxIter);
	arTrack.project();

	CParameterMap pm;
	pm.SetValue(arTrack.m_weights);
	return pm.GetValue(pResult);
	end_function
}


/****************************************************************
**	 
**	Routine: MLInvertMatrix
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/
	

HRESULT __stdcall MLInvertMatrix(VARIANT InputMatrix, VARIANT *pResult)
{
	begin_function
	map_parameter(InputMatrix, CMatrix, in)

	CMatrix matrixInverse;
	if (!::invertMatrix(matrixInverse, in)) throw "Matrix is not invertible";
	unmap_parameter(matrixInverse, pm)
	return pm.GetValue(pResult);
	end_function
}



/****************************************************************
**	 
**	Routine: MLGetFutureVolatilityFromAsset
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/

HRESULT __stdcall MLGetFutureVolatilityFromAsset(VARIANT Asset,VARIANT StrikeHandle, VARIANT VolStartDate, VARIANT VolEndDate, VARIANT BidOfferOpt, VARIANT ReferenceDateOpt, VARIANT ReturnNaturalVolatilityOpt, VARIANT SkewInvariantStrikeOpt,VARIANT* pResult)
{
	begin_function	

	int									nStrikes;
	vector<MlEqStrikeHandle>			m_pStrikes;
	vector<MlEqStrikeHandle>			m_pSkewInvariantStrike;


	map_object_parameter(Asset, Asset, hAsset);
	map_parameter(VolStartDate, long, nSDate);
	map_parameter(VolEndDate, long, nEDate);

	map_optional_parameter(ReferenceDateOpt, long, nToday, hAsset->GetDateHandle()->GetDate());
	map_optional_parameter(ReturnNaturalVolatilityOpt, bool, bReturnNaturalVolatility, false);
	map_optional_enum_parameter(BidOfferOpt, BidAskMidEnum, nBidOffer, Middle);
	map_object_parameter(StrikeHandle, Strikes, hStrike);

	m_pStrikes = hStrike->m_hStrikes;
	nStrikes = hStrike->m_hStrikes.size();

	map_optional_object_parameter(SkewInvariantStrikeOpt, Strikes, hSkewInvStrike);

	int numberSkewInvStrikes = 0;
	if ( !!hSkewInvStrike ){
		m_pSkewInvariantStrike	= hSkewInvStrike->m_hStrikes;
		numberSkewInvStrikes	= hSkewInvStrike->m_hStrikes.size();

		if ( numberSkewInvStrikes > 1 ){
			throw("onlt one skewInvariantStrike can be entered");
		}
	}


	CVector result(nStrikes);

	MlEqVolatilityStructureHandle vol =  hAsset->GetVolatilityStructure();

	double spot = hAsset->GetSpot(nToday);

	for ( int i = 0 ; i < nStrikes; i++ )
	{
		if (bReturnNaturalVolatility)
		{
			if ( numberSkewInvStrikes == 0 ){
				result[i] = hAsset->GetNaturalVolatility(*(hStrike->m_hStrikes[i]),nSDate,nEDate, nBidOffer);
			}
			else {
				result[i] = hAsset->GetNaturalVolatility(*(hStrike->m_hStrikes[i]),nSDate,nEDate, nBidOffer,hSkewInvStrike->m_hStrikes[0]);
			}

		}
		else 
		{
			// default
			if ( numberSkewInvStrikes == 0 ){
				result[i] = hAsset->GetCompositeVolatility(*(hStrike->m_hStrikes[i]), nSDate,(double)nEDate, nBidOffer);
			}
			else{
				result[i] = hAsset->GetCompositeVolatility(*(hStrike->m_hStrikes[i]), nSDate,(double)nEDate, nBidOffer,hSkewInvStrike->m_hStrikes[0]);
			}
		}
	}

	// return the result
	CParameterMap pm;
	pm.SetValue(result);
	return pm.GetValue(pResult);

	end_function
}


/****************************************************************
**	 
**	Routine: MLGetNaiveFutureVolatilityFromAsset
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/

HRESULT __stdcall MLGetNaiveFutureVolatilityFromAsset(VARIANT Asset,VARIANT StrikeHandleShort,VARIANT StrikeHandleLong, VARIANT VolStartDate, VARIANT VolEndDate, VARIANT BidOfferOpt, VARIANT ReferenceDateOpt, VARIANT ReturnNaturalVolatilityOpt, VARIANT* pResult)
{
	begin_function	

	int									nStrikes;
	vector<MlEqStrikeHandle>			m_pStrikesShort,m_pStrikesLong;

	map_object_parameter(Asset, Asset, hAsset);
	map_parameter(VolStartDate, long, nSDate);
	map_parameter(VolEndDate, long, nEDate);

	map_optional_parameter(ReferenceDateOpt, long, nToday, hAsset->GetDateHandle()->GetDate());
	map_optional_parameter(ReturnNaturalVolatilityOpt, bool, bReturnNaturalVolatility, false);
	map_optional_enum_parameter(BidOfferOpt, BidAskMidEnum, nBidOffer, Middle);		
	map_object_parameter(StrikeHandleShort, Strikes, hStrikeShort);
	map_object_parameter(StrikeHandleLong, Strikes, hStrikeLong);

	m_pStrikesShort = hStrikeShort->m_hStrikes;
	m_pStrikesLong = hStrikeLong->m_hStrikes;

	nStrikes = m_pStrikesShort.size();

	if ( m_pStrikesLong.size() != nStrikes ){
		throw("number of long strikes must coincide with number of short strikes");
	}

	CVector result(nStrikes);

	if ( bReturnNaturalVolatility )
	{
		for ( int i = 0 ; i < nStrikes; i++ ){
			result[i] = hAsset->GetNaturalNaiveFutureVolatility(*(hStrikeShort->m_hStrikes[i]),*(hStrikeLong->m_hStrikes[i]),(double) nSDate,(double) nEDate,nBidOffer) ;
		}
	}
	else
	{
		for ( int i = 0 ; i < nStrikes; i++ ){
			result[i] = hAsset->GetNaiveFutureVolatility(*(hStrikeShort->m_hStrikes[i]),*(hStrikeLong->m_hStrikes[i]),(double) nSDate,(double) nEDate, nBidOffer) ;
		}

	}

	// return the result
	CParameterMap pm;
	pm.SetValue(result);
	return pm.GetValue(pResult);

	end_function
}


/****************************************************************
**	 
**	Routine: MLBasketLookback
**	Returns: VARIANT
**	Action : returns price of option with payoff (max_over_assets [max_over_time [arithm_average_till_time] ] - Strike )+
**  
****************************************************************/


double partie_positive(double a) {if (a<0) return 0;return a;}


HRESULT __stdcall MLExpectationOfMax(VARIANT Dates,VARIANT Forwards,VARIANT Vol,VARIANT MC_Steps, VARIANT* pResult)
{
	begin_function	

	map_parameter(Dates, CVector, vDates);
    map_parameter(Forwards, CVector, forward);
    map_parameter(Vol, CVector, sigma);
	//map_parameter(DiscountFactor, double, vDiscountFactor);
	map_parameter(MC_Steps, long, vMC_Steps);

	//init
	long NT=vDates.getsize();


    CVector s(NT);
	CVector t(NT);

	//find index in array with dates before today
	long TodayDate = MlEqDate::GetCurrentDate();	
	long iSimStart=NT; //start index for simulations, all previous spots are known
	
	
	for(int i=0;i<NT;i++)
		if(vDates[i]>TodayDate) {iSimStart=i;break;}
	//filling up known spots

	for(int i=0;i<iSimStart;i++)// iSimStart if always future date [the first one]
		s[i]=forward[i];
	
	//setting up <double> time from <long> time
	t[0]=0;
	for(int i=1;i<NT;i++)
		t[i]=t[i-1]+(vDates[i]-vDates[i-1])/365.25;


		
	//start caluclation of the expectation
		
	long seed=-123;
	
	double mc_average=0;
		
	for(int mc_i=0;mc_i<vMC_Steps;mc_i++)
	{
	
		double payoff=999999;
		int i; //time counter
		
		for(i=iSimStart;i<NT;i++)  //i is time step
		{
			double dt,sqrt_dt;

			dt=t[i]-t[i-1];
			sqrt_dt=sqrt(dt);

			double n; //standart normal variable
			
			n=sqrt(-2.*log(ran2(&seed)))*cos(2*3.141592653589*ran2(&seed));
			
			s[i]=s[i-1]*forward[i]/forward[i-1]*exp(-.5*sigma[i-1]*sigma[i-1]*dt+sigma[i-1]*sqrt_dt*n);
		
		}
	//find maximal return
	//add if (payoff not 0)
	double a; //current maximum of spot
	a=s[0];
	
	int maxj; //current index of maximum spot
	maxj=0;
	
	for(int j=1;j<NT;j++) { if(s[j]>a) {maxj=j; a=s[j];} }
	
	payoff=s[maxj];
		
	mc_average+=payoff;
	}//end mc
	mc_average/=vMC_Steps;
		
	//put price into Excel

	double price;
	price=mc_average;

	CParameterMap pm;
	pm.SetValue(price);
	return pm.GetValue(pResult);

	end_function
}



/****************************************************************
**	 
**	Routine: MLCalibLocalVol
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLCalibrateLocalVol(VARIANT asset, VARIANT localVolGridIntep,VARIANT calibStrikes, VARIANT calibDates, VARIANT Tolerance, VARIANT nGauss, VARIANT* pResult)
{

	begin_function	



	map_object_vector_parameter(localVolGridIntep, Interpolator, hInterpolator);
	map_object_parameter(asset, Asset, hasset);
	map_parameter(calibDates, GVector < long >, hcalibDates);
	map_parameter(Tolerance, double, hTol);


	map_strikes_vector_parameter(calibStrikes, hStrikes);
	map_parameter(nGauss, long,N );


	int n = hStrikes.size();
	GVector < GVector  < MlEqStrikeHandle > > calibStk;

	calibStk.resize(hStrikes.size());
	for ( int i = 0 ; i < n; i++ )
	{
		calibStk[i].resize(hStrikes[i].size());
		for ( int j = 0 ; j < hStrikes[i].size(); j++ )
		{

			calibStk[i][j] = hStrikes[i][j];
		}
	}

	GVector < MlEqInterpolatorHandle > hInterps(hInterpolator.size());

	for ( int i = 0 ; i < hInterpolator.size(); i++ ){
		hInterps[i] = hInterpolator[i];
	}

	long horizonDate = hcalibDates[hcalibDates.getsize()-1];

	LocalVolGrid lv;

	double hrate=-1,hspacialResolution=-1;//not used any more
	lv.initialise(*hasset,calibStk,hcalibDates,hInterps,horizonDate,hrate,hspacialResolution,N);



	CVector impliedVols;
	CMatrix mostlikelyPaths;

	int idateCalib = 0;
	double accuracy = 0.01;//hTol;
	bool returnMostlikelypath = true;
	
	
	double tolerance = hTol;
	
	
	lv.createLocalVol(impliedVols,tolerance,accuracy);
	
	CParameterMap pm;
	pm.SetValue(impliedVols);
	return pm.GetValue(pResult);
	
	
	end_function
}
	
	
	
/****************************************************************
**	 
**	Routine: MLGetLocalImpliedVol
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/
	
	
HRESULT __stdcall MLGetLocalImpliedVol(VARIANT Strike, VARIANT VolDate,VARIANT localVolGridIntep,VARIANT InterpolatorDates,VARIANT* pResult)
{	


	begin_function	
	
	map_object_vector_parameter(localVolGridIntep, Interpolator, hInterpolator);
	map_parameter(VolDate, long , hDate);
	map_object_parameter(Strike, Strikes, hStrike);
	map_parameter(InterpolatorDates, GVector < long >, hInterpDates);


	GVector  < MlEqStrikeHandle >  cStrikes(hStrike->m_hStrikes.size());

	for ( int i = 0 ; i < cStrikes.getsize(); i++ ){
		cStrikes[i] = hStrike->m_hStrikes[i];
	}


	GVector < MlEqInterpolatorHandle > hInterps(hInterpolator.size());
	for ( int i = 0 ; i < hInterpolator.size(); i++ ){
		hInterps[i] = hInterpolator[i];
	}

	CVector impliedVols;
	CMatrix mostlikelyPaths;

	LocalVolGrid lv;
	bool returnMostlikelypath = false;

	double accuracy = 0.01;

	int index;
	MlEqMaths::Locate(hInterpDates,hDate,index); 
	index++;
	lv.localToImplied(impliedVols,mostlikelyPaths,cStrikes,index,0,accuracy,returnMostlikelypath,hInterps);


	CParameterMap pm;
	pm.SetValue(impliedVols);
	return pm.GetValue(pResult);


	end_function
}


/****************************************************************
**	 
**	Routine: MLCalibLocalVol
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/



HRESULT __stdcall MLLogContractVol(VARIANT asset, VARIANT nMaturity,VARIANT ngaussOpt,VARIANT upperLimitOpt,VARIANT lowerLimitOpt,VARIANT* pResult)
{

	begin_function	

	map_object_parameter(asset, Asset, hasset);
	map_parameter(nMaturity,  long , nmat);
	map_optional_parameter(ngaussOpt,  long , n,50);
	map_optional_parameter(upperLimitOpt,  double , upper,5);
	map_optional_parameter(lowerLimitOpt,  double , lower,-5);

	double vol = impliedLogContractVol(*hasset,nmat,n,upper,lower);


//	vol = expectedVol(hasset,nmat,n,upper,lower);


	CParameterMap pm;
	pm.SetValue(vol);
	return pm.GetValue(pResult);

	end_function

}



/****************************************************************
**	 
**	Routine: MLCalibLocalVol
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/





HRESULT __stdcall MLDeltaHedgedOptionPortfolio(VARIANT asset,VARIANT startDate,VARIANT settlementDates ,VARIANT pastPnl,
											   VARIANT PnlTenor, VARIANT BusinessDayConvention,
											   VARIANT Calendar, VARIANT Portfolio, VARIANT bookingVols,
											   VARIANT bookingRate,VARIANT bookingDiv,
											   VARIANT DayCountConvention, VARIANT nIntegrationPoints,VARIANT useDeltaTerm,
											   VARIANT approximationDateOpt,VARIANT nTimePoints,VARIANT useGaussIntegrator,VARIANT accuracyOpt,VARIANT nstdevOpt,VARIANT* pResult)
{



	begin_function	

	map_object_parameter(asset, Asset, hAsset);
	map_parameter(startDate,  long  , stDate);
	map_parameter(settlementDates, GVector< long > , settleDates);
	map_parameter(pastPnl, CMatrix , pastpnl);
	
	map_parameter(PnlTenor, std::string, szTenor);
	map_enum_parameter(BusinessDayConvention, BusinessDayConventionEnum, bdc);
	
	map_parameter(Calendar, std::string, szCalendar);
	map_parameter(Portfolio, CMatrix, portfolio);
	map_parameter(bookingVols, CVector,vols );
	map_parameter(bookingRate, double, fbookingRate);
	map_parameter(bookingDiv, double , fbookingDiv);
	map_parameter(nIntegrationPoints,  long  , mgauss);
	map_parameter(useDeltaTerm,  long  , fuseDeltaTerm);
	map_optional_parameter(approximationDateOpt,  long  , mapproximationDate,1e10);
	map_parameter(nTimePoints,  long  , tgauss);
	map_parameter(useGaussIntegrator,  bool, useGauss);
	map_optional_parameter(accuracyOpt,  double  , acc,1e-6);
	map_optional_parameter(nstdevOpt,  double  , nstd,4.0);
	map_enum_parameter(DayCountConvention, DayCountConventionEnum, dcc);
	

	DeltaHedgedOptionPortfolio dhP;

	dhP.initialise(	
					hAsset,
					stDate,
					settleDates,
					szTenor,
					bdc,
					szCalendar,
					portfolio,
					vols,
					fbookingRate,
					fbookingDiv  ,
					dcc,
					fuseDeltaTerm,
					useGauss,
					mgauss,
					nstd,
					acc,
					mapproximationDate,
					tgauss

					);


	CMatrix value;

	long approximationDate=0;

	dhP.price(value,hAsset);

	CParameterMap pm;
	pm.SetValue(value);
	return pm.GetValue(pResult);

	end_function

}


/****************************************************************
**	 
**	Routine: MLCalibLocalVol
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/



HRESULT __stdcall MLPatHaagenImpliedVol(VARIANT StrikeHandle,VARIANT maturity,VARIANT forward,VARIANT correl,VARIANT beta,VARIANT volvol, VARIANT shortvol, VARIANT* pResult)
{

	begin_function	


	map_parameter(maturity,  double  , xmat);
	map_parameter(forward,  double  , xfwd);
	map_parameter(correl,  double  , xcorrel);
	map_parameter(beta,  double  , xbeta);
	map_parameter(volvol,  double  , xvolvol);
	map_parameter(shortvol,  double  , xshortvol);
	map_object_parameter(StrikeHandle, Strikes, hStrike);
		
	// do the work


	int nStrikes = hStrike->m_hStrikes.size();


	CVector vols(nStrikes);
	MlEqStrike stk;	

	for (int i = 0; i < nStrikes; i++)
	{
		MlEqStrike::convertStrikes(stk,*(hStrike->m_hStrikes[i]));
		vols[i] =PatHaagenImpliedVol(stk.m_strike,xmat,xfwd,xcorrel,xbeta,xvolvol,xshortvol);
	}



	CParameterMap pm;
	pm.SetValue(vols);
	return pm.GetValue(pResult);

	end_function

}


/****************************************************************
**	 
**	Routine: MLCalibHeston
**	Returns: VARIANT
**	Action : calibration
**           
****************************************************************/


HRESULT __stdcall MLHestonCalibration(VARIANT asset, VARIANT nMaturity, VARIANT calibStrikesOpt, VARIANT* pResult)
{

	begin_function	


	map_object_parameter(asset, Asset, hasset);
	map_parameter(nMaturity,  long , nmat);
	map_optional_object_parameter(calibStrikesOpt, Strikes, hStrike);

	GVector  < MlEqStrikeHandle >  cStrikes;
	if( !!hStrike )
	{
		cStrikes.resize(hStrike->m_hStrikes.size());

		for ( int i = 0 ; i < cStrikes.getsize(); i++ ){
			cStrikes[i] = hStrike->m_hStrikes[i];
		}
	}
	else{
		cStrikes.resize(0);
	}



	RCPtr<heston_fitter> heston = new heston_fitter(hasset, nmat);

	heston->calibrate(cStrikes);

	HestonParameters* params = heston->getParams();

	CVector result(5,0.0);

	result[0]	=	params->v0;
	result[1]	=	params->vas;
	result[2]	=	params->kappa;
	result[3]	=	params->vovol;
	result[4]	=	params->correl;


	CParameterMap pm;
	pm.SetValue(result);
	return pm.GetValue(pResult);

	end_function

}



/****************************************************************
**	 
**	Routine: MLCalibHeston
**	Returns: VARIANT
**	Action : calibration
**           
****************************************************************/


HRESULT __stdcall MLSABRCalibration(VARIANT asset, VARIANT nMaturity, VARIANT calibStrikesOpt, VARIANT* pResult)
{

	begin_function	


 	map_object_parameter(asset, Asset, hasset);
	map_parameter(nMaturity,  long , nmat);
	map_optional_object_parameter(calibStrikesOpt, Strikes, hStrike);

	GVector  < MlEqStrikeHandle >  cStrikes;
	if( !!hStrike )
	{
		cStrikes.resize(hStrike->m_hStrikes.size());
		for ( int i = 0 ; i < cStrikes.getsize(); i++ ){
			cStrikes[i] = hStrike->m_hStrikes[i];
		}
	}
	else{
		cStrikes.resize(0);
	}

	RCPtr<sabr_fitter> sabr = new sabr_fitter(hasset, nmat);

	sabr->calibrate(cStrikes);

	CVector result(1);
	sabr->getParams( result );
	// short vol
	// vol of vol
	// correlation

	CParameterMap pm;
	pm.SetValue(result);
	return pm.GetValue(pResult);

	end_function

}

/****************************************************************
**	 
**	Routine: MLCalibHeston
**	Returns: VARIANT
**	Action : calibration
**           
****************************************************************/


HRESULT __stdcall MLHestonJumpCalibration(VARIANT asset, VARIANT nMaturity, VARIANT calibStrikesOpt, VARIANT* pResult)
{

	begin_function	


	map_object_parameter(asset, Asset, hasset);
	map_parameter(nMaturity,  long , nmat);
	map_optional_object_parameter(calibStrikesOpt, Strikes, hStrike);

	GVector  < MlEqStrikeHandle >  cStrikes;
	if( !!hStrike )
	{
		cStrikes.resize(hStrike->m_hStrikes.size());

		for ( int i = 0 ; i < cStrikes.getsize(); i++ ){
			cStrikes[i] = hStrike->m_hStrikes[i];
		}
	}
	else{
		cStrikes.resize(0);
	}



	RCPtr<heston_jump_fitter> heston = new heston_jump_fitter(hasset, nmat);
	heston->calibrate(cStrikes);

	HestonParameters* params = heston->getParams();

	CVector result(7,0.0);

	result[0]	=	params->v0;
	result[1]	=	params->vas;
	result[2]	=	params->kappa;
	result[3]	=	params->vovol;
	result[4]	=	params->correl;

	result[5]	=	heston->m_lambda;
	result[6]	=	heston->m_jump;


	CParameterMap pm;
	pm.SetValue(result);
	return pm.GetValue(pResult);

	end_function

}




/****************************************************************
**	 
**	Routine: MLCalibHeston
**	Returns: VARIANT
**	Action : check calibration
**           
****************************************************************/


HRESULT __stdcall MLHestonImpliedVol(VARIANT forward, VARIANT maturity, VARIANT StrikeHandle, VARIANT v0, VARIANT vas, VARIANT kappa, VARIANT vovol, VARIANT correl, VARIANT* pResult)
{

	begin_function	


	map_parameter(forward, double, fforward);
	map_object_parameter(StrikeHandle, Strikes, hStrike);
	map_parameter(maturity,  double , fmat);
	map_parameter(v0, double, fv0);
	map_parameter(vas, double, fvas);
	map_parameter(kappa, double, fkappa);
	map_parameter(vovol, double, fvovol);
	map_parameter(correl, double, fcorrel);

	HestonParameters params;

	params.v0 = fv0;
	params.vas = fvas;
	params.kappa = fkappa;
	params.vovol = fvovol;
	params.correl = fcorrel;

	RCPtr<heston_pricer> pricer = new heston_pricer(params);
	
	int nStrikes = hStrike->m_hStrikes.size();
	CVector vols(nStrikes);
	MlEqStrike stk;	

	for (int i = 0; i < nStrikes; i++)
	{
		MlEqStrike::convertStrikes(stk,*(hStrike->m_hStrikes[i]));
		vols[i] = pricer->implied_vol(fforward, fmat, stk.m_strike);
	}


	CParameterMap pm;
	pm.SetValue(vols);
	return pm.GetValue(pResult);

	end_function

}


//////


HRESULT __stdcall MLHestonJumpImpliedVol(VARIANT forward, VARIANT maturity, VARIANT StrikeHandle, VARIANT Params, VARIANT* pResult)
{

	begin_function	


	map_parameter(forward, double, fforward);
	map_object_parameter(StrikeHandle, Strikes, hStrike);
	map_parameter(maturity,  double , fmat);
	map_parameter(Params, CVector, vParams);

	HestonParameters params;

	params.v0		= vParams[0];
	params.vas		= vParams[1];
	params.kappa	= vParams[2];
	params.vovol	= vParams[3];
	params.correl	= vParams[4];

	double lambda	= vParams[5];
	double jump		= vParams[6];
//	jump = log(1.+jump);	// depends on the convention


	RCPtr<heston_jump_pricer> pricer = new heston_jump_pricer(params, lambda, jump);

	int nStrikes = hStrike->m_hStrikes.size();
	CVector vols(nStrikes);
	MlEqStrike stk;	

	for (int i = 0; i < nStrikes; i++)
	{
		MlEqStrike::convertStrikes(stk,*(hStrike->m_hStrikes[i]));
		vols[i] = pricer->implied_vol(fforward, fmat, stk.m_strike);
	}


	CParameterMap pm;
	pm.SetValue(vols);
	return pm.GetValue(pResult);

	end_function

}

/****************************************************************
**	 
**	Routine: MLCalibHestonTermStructure
**	Returns: VARIANT
**	Action : calibration
**           
****************************************************************/


HRESULT __stdcall MLHestonTermCalibration(VARIANT asset, VARIANT nMaturities, VARIANT calibStrikesOpt, VARIANT* pResult)
{

	begin_function	


	map_object_parameter(asset, Asset, hasset);
	map_parameter(nMaturities, GVector < long >, hcalibDates);
	map_optional_object_parameter(calibStrikesOpt, Strikes, hStrike);

	GVector  < MlEqStrikeHandle >  cStrikes;
	if( !!hStrike )
	{
		cStrikes.resize(hStrike->m_hStrikes.size());

		for ( int i = 0 ; i < cStrikes.getsize(); i++ ){
			cStrikes[i] = hStrike->m_hStrikes[i];
		}
	}
	else{
		cStrikes.resize(0);
	}

	int nMat = hcalibDates.getsize();

	RCPtr<heston_term_fitter> heston = new heston_term_fitter(hasset, hcalibDates);
	heston->calibrate( cStrikes );


	CMatrix result(5,nMat);

	for(int t=0; t<nMat; t++)
	{
		HestonParameters* params = heston->getParams(t);

		result[0][t]	=	params->v0;
		result[1][t]	=	params->vas;
		result[2][t]	=	params->kappa;
		result[3][t]	=	params->vovol;
		result[4][t]	=	params->correl;
	}



	CParameterMap pm;
	pm.SetValue(result);
	return pm.GetValue(pResult);

	end_function

}


/****************************************************************
**	 
**	Routine: MLGetDensity
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/



HRESULT __stdcall MLEqPDEBarrier(VARIANT Vol, VARIANT DiscountRate, VARIANT Forward, VARIANT Spot, 
								 VARIANT asOfDateHandle,VARIANT maturityDate,
								 VARIANT numberSpacialPoints,VARIANT numberTimeSteps,VARIANT nstdev,
								 VARIANT UpperBarrier,VARIANT LowerBarrier,VARIANT UpperRebate,VARIANT LowerRebate,
								 VARIANT Callput,VARIANT Strike,

								 VARIANT* pResult)
{		
	
	// map the inputs
	begin_function


	map_parameter(Vol, double, fVol);
	map_parameter(DiscountRate, double, frate);
	map_parameter(Forward, double, fforward);
	map_parameter(Spot, double, fspot);
	map_object_parameter(asOfDateHandle, Date, asof);
	map_parameter(maturityDate, long, matDate);
	map_parameter(numberSpacialPoints, long, nx);
	map_parameter(numberTimeSteps, long, nt);
	map_parameter(nstdev, double, fnstdev);

	map_parameter(UpperBarrier, double, upper);
	map_parameter(LowerBarrier, double, lower);
	map_parameter(UpperRebate, double, upperR);
	map_parameter(LowerRebate, double, lowerR);
	map_parameter(Callput, double, cp);
	map_parameter(Strike, double, strike);


	CVector nothing;
	knockOutBarrierHelperHandle barrier = new knockOutBarrierHelper(upper,lower,
					   upperR,lowerR,
					   fspot,cp,strike,
					   false, nothing);


	pdeBlackScholes pde;

	pde.init(fVol,fforward,frate,fspot,
			 asof,matDate,
			 nx,nt,fnstdev,
			 &(*barrier),lower,upper);


	 CVector result(1);

	double xtime = pde.pde_integrate(result);

	// return the result
	CParameterMap pm;
	pm.SetValue(result);
	return pm.GetValue(pResult);	
	end_function
}


/****************************************************************
**	 
**	Routine: MLGetDensity
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLEqLocalVolPDEBarrier(VARIANT Asset,VARIANT maturityDate,
								 VARIANT numberSpacialPoints,VARIANT numberTimeSteps,VARIANT nstdev,
								 VARIANT UpperBarrier,VARIANT LowerBarrier,VARIANT UpperRebate,VARIANT LowerRebate,
								 VARIANT Callput,VARIANT Strike, VARIANT PayEnd, VARIANT LoanspreadBump,
								 VARIANT* pResult)
{		
	
	// map the inputs
	begin_function


	map_object_parameter(Asset, Asset, hasset);
	map_parameter(maturityDate, long, matDate);
	map_parameter(numberSpacialPoints, long, nx);
	map_parameter(numberTimeSteps, long, nt);
	map_parameter(nstdev, double, fnstdev);

	map_parameter(UpperBarrier, double, upper);
	map_parameter(LowerBarrier, double, lower);
	map_parameter(UpperRebate, double, upperR);
	map_parameter(LowerRebate, double, lowerR);
	map_parameter(Callput, double, cp);
	map_parameter(Strike, double, strike);
	map_parameter(LoanspreadBump, double, dls);
	map_parameter(PayEnd, bool, payEnd);




	MlEqConstDateHandle	dateh = hasset->GetDateHandle();
	int nToday = dateh->GetDate();
	double spot = hasset->GetSpot(nToday);


	double volbump = 0.01;
	DupireLocalVolHandle lv = new DupireLocalVol;
	lv->initialize(hasset, lower, upper, matDate, nt, nx, volbump );


	CVector df(nt,1.);
	if( payEnd )
	{
		GVector<long> date = lv->getDates();

		for(int t=0; t<nt; t++){
			df[t] = hasset->GetPayZeroCurve(true)->GetDiscountFactor(date[t], matDate);
		}
	}

	knockOutBarrierHelperHandle barrier = new knockOutBarrierHelper(upper, lower, upperR, lowerR, spot, cp, strike, payEnd, df);


	
	pdeLocalVol pde;
	pde.init( lv,&(*barrier));

	CVector result(5);
	CVector tmp(1);

	pde.pde_integrate(tmp);
	double price = tmp[0];
	result[0] = price;

	double delta, gamma ;
	barrier->greeks( delta, gamma, 0, 0, 0, 0, &pde );
	result[1] = delta ;
	result[2] = gamma ;
	result[2] *= 0.01 * spot ;	// trader's convention for the gamma

	lv->parallelVegaBump();			// vega
	pde.init( lv,&(*barrier));
	pde.pde_integrate(tmp);
	result[3] = tmp[0] - price ;
	lv->parallelVegaBump();

	lv->bumpLocalDrift( dls );		// loanspread risk
	pde.init( lv,&(*barrier));
	pde.pde_integrate(tmp);
	result[4] = tmp[0] - price ;

	// return the result
	CParameterMap pm;	

	pm.SetValue(result);
	return pm.GetValue(pResult);	
	end_function
}



HRESULT __stdcall MLEqLocalVolPDEAmerican(VARIANT Asset,VARIANT maturityDate,
								 VARIANT numberSpacialPoints,VARIANT numberTimeSteps,VARIANT nstdev,
								 VARIANT Callput,VARIANT Strike, VARIANT LoanspreadBump,
								 VARIANT* pResult)
{		
	
	// map the inputs
	begin_function


	map_object_parameter(Asset, Asset, hasset);
	map_parameter(maturityDate, long, matDate);
	map_parameter(numberSpacialPoints, long, nx);
	map_parameter(numberTimeSteps, long, nt);
	map_parameter(nstdev, double, fnstdev);

	map_parameter(Callput, double, cp);
	map_parameter(Strike, double, strike);
	map_parameter(LoanspreadBump, double, dls);

	MlEqConstDateHandle	dateh = hasset->GetDateHandle();
	int nToday = dateh->GetDate();
	double spot = hasset->GetSpot(nToday);

	RCPtr<americanOptionHelper> american = new americanOptionHelper(spot, cp, strike);

	double volbump = 0.01;
	DupireLocalVolHandle lv = new DupireLocalVol;

	double mat = dateh->GetYearFraction(matDate);
	double fwd = hasset->GetForward(matDate, false);
	double vol = hasset->GetVolatility(MlEqStrike(fwd),matDate);

	double highSpot = spot*exp(fnstdev*vol*sqrt(mat));
	double lowSpot  = spot*exp(- 1.5*fnstdev *vol*sqrt(mat));	// use something asymetric...

	lv->initialize(hasset,lowSpot,highSpot,matDate,nt,nx, volbump);
	
	pdeLocalVol pde;
	pde.init( lv,&(*american));

	CVector result(5);
	CVector tmp(1);

	pde.pde_integrate(tmp);
	double price = tmp[0];
	result[0] = price;

	double delta, gamma ;
	american->greeks( delta, gamma, 0, 0, 0, 0, &pde );
	result[1] = delta ;
	result[2] = gamma ;
	result[2] *= 0.01 * spot ;	// trader's convention for the gamma

	lv->parallelVegaBump();
	pde.init( lv,&(*american));
	pde.pde_integrate(tmp);
	result[3] = tmp[0] - price ;
	lv->parallelVegaBump();

	lv->bumpLocalDrift( dls );
	pde.init( lv,&(*american));
	pde.pde_integrate(tmp);
	result[4] = tmp[0] - price ;

	// return the result
	CParameterMap pm;	

	pm.SetValue(result);
//	pm.Transpose();
	return pm.GetValue(pResult);	
	end_function
}



HRESULT __stdcall MLEqPDEAmerican(VARIANT Asset,VARIANT maturityDate,
								 VARIANT numberSpacialPoints,VARIANT numberTimeSteps,VARIANT nstdev,
								 VARIANT Callput,VARIANT Strike, VARIANT Volatility,
								 VARIANT* pResult)
{
	
	// map the inputs
	begin_function


	map_object_parameter(Asset, Asset, hasset);
	map_parameter(maturityDate, long, matDate);
	map_parameter(numberSpacialPoints, long, nx);
	map_parameter(numberTimeSteps, long, nt);
	map_parameter(nstdev, double, fnstdev);

	map_parameter(Callput, double, cp);
	map_parameter(Strike, double, strike);
	map_parameter(Volatility, double, vol);


	MlEqConstDateHandle	dateh = hasset->GetDateHandle();
	int nToday = dateh->GetDate();
	MlEqDateHandle asof = new MlEqDate(nToday);
	double spot = hasset->GetSpot(nToday);

	RCPtr<americanOptionHelper> american = new americanOptionHelper(spot, cp, strike);

	double volbump = 0.01;
	DupireLocalVolHandle lv = new DupireLocalVol;

	double mat = dateh->GetYearFraction(matDate);
	double fwd = hasset->GetForward(matDate, false);
	double volk = hasset->GetVolatility(MlEqStrike(strike),matDate);

	double upper = spot*exp(fnstdev*vol*sqrt(mat));
	double lower  = spot*exp(- 1.5*fnstdev *vol*sqrt(mat));	// use something asymetric...

	// get average discount / div yield
	double dyld = log( fwd/spot ) / mat;
	double rate = hasset->GetPayZeroCurve(true)->GetDiscountFactor(nToday, matDate);
	rate = -log( rate )/mat;


	pdeBlackScholes pde;

	pde.init(vol,fwd,rate,spot,
			 asof,matDate,
			 nx,nt,fnstdev,
			 &(*american),lower,upper);

	CVector result(4);
	CVector tmp(1);

	pde.pde_integrate(tmp);
	double price = tmp[0];
	result[0] = price;

	double delta, gamma ;
	american->greeks( delta, gamma, 0, 0, 0, 0, &pde );
	result[1] = delta ;
	result[2] = gamma ;
	result[2] *= 0.01 * spot ;	// trader's convention for the gamma

	vol += volbump;
	pde.init(vol,fwd,rate,spot,
			 asof,matDate,
			 nx,nt,fnstdev,
			 &(*american),lower,upper);

	pde.pde_integrate(tmp);
	result[3] = tmp[0] - price ;

	// return the result
	CParameterMap pm;	

	pm.SetValue(result);
//	pm.Transpose();
	return pm.GetValue(pResult);	
	end_function
}



HRESULT __stdcall MLEqLocalVolVanilla(VARIANT Asset,VARIANT maturityDate,
								 VARIANT numberSpacialPoints,VARIANT numberTimeSteps,VARIANT nstdev,
								  VARIANT Callput,VARIANT Strikes, VARIANT Bump, //VARIANT lowBd,VARIANT upBd,
								 VARIANT* pResult)
{		
	
	// map the inputs
	begin_function


	map_object_parameter(Asset, Asset, hasset);
	map_parameter(maturityDate, long, matDate);
	map_parameter(numberSpacialPoints, long, nx);
	map_parameter(numberTimeSteps, long, nt);
	map_parameter(nstdev, double, fnstdev);
	map_parameter(Bump, double, fbump);
//	map_parameter(lowBd, double, lowBound);
//	map_parameter(upBd, double, upBound);

	map_parameter(Callput, CVector, cp);
	map_parameter(Strikes, CVector, strike);

	MlEqConstDateHandle	dateh = hasset->GetDateHandle();
	int nToday = dateh->GetDate();
	double spot = hasset->GetSpot(nToday);


	DupireLocalVolHandle lv = new DupireLocalVol;

	double mat = dateh->GetYearFraction(matDate);
	double fwd = hasset->GetForward(matDate, false);
	double vol = hasset->GetVolatility(MlEqStrike(fwd),matDate);

	double highSpot = spot*exp(fnstdev*vol*sqrt(mat));
	double lowSpot  = spot*exp(- 1.5*fnstdev*vol*sqrt(mat));


//	lv->m_lowBound	= lowBound;
//	lv->m_upBound	= upBound;
	lv->initialize(hasset,lowSpot,highSpot,matDate,nt,nx);

//	RCPtr<FittedLocalVol> fitTest = new FittedLocalVol( *lv );
//	fitTest->fitSurface();

	
	RCPtr<plainVanillaOptionHelper> vanilla = new plainVanillaOptionHelper();

	pdeLocalVol pde;
//	pde.init( &(*fitTest),&(*vanilla));
	pde.init( lv,&(*vanilla));
//	pde.bumpLocalVol( fbump );???

	int nstrike = strike.getsize();
	CMatrix result(nstrike,1);

	for(int k=0; k<nstrike; k++)
	{
		CVector tmp(1);
		
//		vanilla->initialize( cp[k], strike[k], spot );???
		pde.pde_integrate(tmp);
		result[k][0] = tmp[0];

	}



	// return the result
	CParameterMap pm;	

	pm.SetValue(result);
	return pm.GetValue(pResult);	
	end_function
}



/****************************************************************
**	 
**	Routine: MLGetDensity
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/



HRESULT __stdcall MLEqBlackScholesBarrier(VARIANT Vol, VARIANT DiscountRate, VARIANT Forward, VARIANT Spot, 
								 VARIANT asOfDateHandle,VARIANT maturityDate,
								 VARIANT UpperBarrier,VARIANT LowerBarrier,
								 VARIANT Callput,VARIANT Strike,
								 VARIANT* pResult)
{		
	
	// map the inputs
	begin_function


	map_parameter(Vol, double, fVol);
	map_parameter(DiscountRate, double, frate);
	map_parameter(Forward, double, fforward);
	map_parameter(Spot, double, fspot);
	map_object_parameter(asOfDateHandle, Date, asof);
	map_parameter(maturityDate, long, matDate);

	map_parameter(UpperBarrier, double, upper);
	map_parameter(LowerBarrier, double, lower);
	map_parameter(Callput, double, cp);
	map_parameter(Strike, double, strike);


	double mat = asof->GetYearFraction(matDate);


	blackScholesBarrier barrier;

	barrier.initialize(
					 strike,
					 cp,
					 upper,
					 lower,
					 mat,
					 fVol,
					 fforward,
					 fspot,
					 frate
					);


	double res = barrier.price();


	// return the result
	CParameterMap pm;
	pm.SetValue(res);
	return pm.GetValue(pResult);	
	end_function
}


/****************************************************************
**	 
**	Routine: MLCalibLocalVol
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLCreateLocalVolGrid(VARIANT asset,VARIANT lowSpot,VARIANT highSpot,VARIANT maturityDate,VARIANT nx,VARIANT nt, VARIANT* pResult)
{

	begin_function	

	map_object_parameter(asset, Asset, hasset);
	map_parameter(lowSpot, double,sl );
	map_parameter(highSpot, double,sh );
	map_parameter(nx, long,mx );
	map_parameter(nt, long,mt );
	map_parameter(maturityDate, long,matDate );


	DupireLocalVol dv;
	dv.initialize(hasset,sl,sh,matDate,mt,mx);

	dv.createLocalVolGrid();
	
	double res=0.0;

	CParameterMap pm;
	pm.SetValue(res);
	return pm.GetValue(pResult);
	
	
	end_function
}
	
	


/****************************************************************
**	 
**	Routine: MLCalibHeston
**	Returns: VARIANT
**	Action : calibration
**           
****************************************************************/


HRESULT __stdcall MLHestonTermImpliedVol(VARIANT forwardArray, VARIANT maturityArray, VARIANT StrikeHandles, VARIANT Params, VARIANT* pResult)
{

	begin_function	

	map_parameter(forwardArray, CVector, vForward);
	map_strikes_vector_parameter(StrikeHandles, hStrikes);

	map_parameter(maturityArray,  CVector , vMat);
	map_parameter(Params, CMatrix, mParams);

	int nDates = mParams.cols();
	if( mParams.rows() != 6 )	throw "Invalid Heston parameters rows";

	HestonParametersTermStructure params;

	GVector<HestonParameters> hParams(nDates);
	(params.t).resize(nDates);

	for(int t=0; t<nDates; t++)
	{
		hParams[t].v0		= mParams[0][t];
		hParams[t].vas		= mParams[1][t];
		hParams[t].kappa	= mParams[2][t];
		hParams[t].vovol	= mParams[3][t];
		hParams[t].correl	= mParams[4][t];

		params.t[t]			= mParams[5][t];
	}

	params.hParams = hParams;


	RCPtr<heston_term_pricer> pricer = new heston_term_pricer(params);


	nDates	=	vMat.getsize();
	int nStrikes = hStrikes[0].size();

	int n = vForward.getsize();
	if( n != hStrikes.size() )	throw"Number of strike handles and forwards mismatch";

	CMatrix vols(nStrikes, nDates);
	MlEqStrike stk;	

	for(int t=0; t<nDates; t++)
	{
		double forward = vForward[t];
		double mat = vMat[t];
	
		for (int i = 0; i < nStrikes; i++)
		{
			MlEqStrike::convertStrikes(stk,*(hStrikes[t][i]));
			vols[i][t] = pricer->implied_vol(forward, mat, stk.m_strike);
		}
	}

	CParameterMap pm;
	pm.SetValue(vols);
	return pm.GetValue(pResult);

	end_function

}


/****************************************************************
**	 
**	Routine: MLCalibHestonTermStructure
**	Returns: VARIANT
**	Action : calibration
**           
****************************************************************/


HRESULT __stdcall LocalVolTest(VARIANT asset, VARIANT nMaturity, VARIANT StrikesVector, VARIANT nPaths,  VARIANT nSteps, VARIANT* pResult)
{
	begin_function	

	map_object_parameter(asset, Asset, hasset);
	map_parameter(nMaturity, long, matDate);
	map_parameter(StrikesVector, CVector, strikes);
	map_parameter(nPaths, long, npaths);
	map_parameter(nSteps, long, nsteps);


	CVector CallPut(strikes.getsize());
	GVector< MlEqStrikeHandle> xstrikes(strikes.getsize());

	for ( int i = 0 ; i < strikes.getsize(); i++ )
	{
		CallPut[i] = 1.0;
		xstrikes[i] = new MlEqStrike(strikes[i]);
	}

	double spot = hasset->GetSpot(hasset->GetDateHandle()->GetDate());

	plainVanillaOptionHelper plainVanilla;
	plainVanilla.initialize(CallPut,xstrikes,spot);

	GVector<MlEqAnalyticCurveWithTanhWingEdgeHandle> localVolSlice(1);

	CVector coeff(3);

	coeff[0] = 0.2;
	coeff[1] = 0.0;
	coeff[2] = 0.0;

	localVolSlice[0] = new MlEqAnalyticCurveWithTanhWingEdge(-100.0,100.0,0,0,0,0,10,coeff);


	int nd = xstrikes.getsize();
	int nx = 100;
	int nt = 200;

	double vol = hasset->GetVolatility(MlEqStrike(hasset->GetForward(matDate,false)),matDate);
	double T   = hasset->GetDateHandle()->GetYearFraction(matDate);

	double lx = -sqrt(T)*vol*4.0;
	double ux =  sqrt(T)*vol*4.0;
	
	const int keep = 1;

	pdeLocalVolEffectiveHandle pdeDriver = new pdeLocalVolEffective();

	MlEqDateHandle date = new MlEqDate(*hasset->GetDateHandle());

	GVector<long> calibDate(1);

	calibDate[0] = matDate;

	pdeDriver->init(localVolSlice,*hasset,date,matDate,calibDate,&plainVanilla,
					nd, nx, NULL, lx, nt,ux, keep);


	CVector xresult(xstrikes.getsize());
	pdeDriver->pde_integrate(xresult);


/*
	RCPtr<CLocalVolMC> lvmc = new CLocalVolMC( hDate, hasset );

	lvmc->initialize(product, fixings, npaths, rng, nsteps);
	lvmc->CLocalVolMC::simulate( results, product );
*/



	CParameterMap pm;
	pm.SetValue(xresult);
	return pm.GetValue(pResult);

	end_function
}




HRESULT __stdcall LocalVolCalibration00(VARIANT asset, VARIANT* pResult)
{
	begin_function	

	map_object_parameter(asset, Asset, hasset);

//	RCPtr<FittedLocalVol> flv = new FittedLocalVol(hasset);
//	flv->fitSurface();

	CParameterMap pm;
	pm.SetValue(0.0);
	return pm.GetValue(pResult);

	end_function
}






HRESULT __stdcall MLVarianceCall(VARIANT DateHandle, VARIANT dateMaturity, VARIANT StrikesVector, VARIANT Params, VARIANT* pResult)
{

	begin_function	


	map_parameter(StrikesVector, CVector, strikes);
	map_parameter(dateMaturity,  long , dMat);
	map_parameter(Params, CVector, vParams);
	map_object_parameter(DateHandle, Date, hDate);

	
	double mat = hDate->GetYearFraction(dMat);

	VarianceCall ppproduct;
	ppproduct.init(strikes, mat, dMat);

	RCPtr<CStochasticVolMC> hmc;


	if( vParams.getsize() == 2 )
	{
		double v0 = vParams[0];
		double vol = vParams[1];

		hmc = new CSABR_MC(hDate,v0, vol);
	}
	else if( vParams.getsize() == 4 )
	{
		HestonParameters params;

		params.v0		= vParams[0];
		params.vas		= vParams[1];
		params.kappa	= vParams[2];
		params.vovol	= vParams[3];

		hmc = new CHestonMC(hDate, params);
	}
	else if( vParams.getsize() == 6 )
	{
		HestonParameters params;

		params.v0		= vParams[0];
		params.vas		= vParams[1];
		params.kappa	= vParams[2];
		params.vovol	= vParams[3];

		double jump		= vParams[4];
		double lambda	= vParams[5];		

		hmc = new CHestonJumpMC(hDate, params, jump, lambda);
	}
	else 
		throw"Wrong number of parameters";
	

	int npaths = 10000;
	int rng = 5;
	CMatrix results;

	hmc->initialize(ppproduct, dMat, npaths, rng);
	hmc->simulate( results, ppproduct );


	CParameterMap pm;
	pm.SetValue(results);
	return pm.GetValue(pResult);

	end_function

}



/****************************************************************
**	 
**	Routine: MLGetDensity
**	Returns: VARIANT
**	Action : 
**           
****************************************************************/


HRESULT __stdcall MLEqLVCallableNote(VARIANT Asset,VARIANT maturityDate,
								 VARIANT numberSpacialPoints,VARIANT numberTimeSteps,VARIANT nstdev,
								 VARIANT referenceSpot,
								 VARIANT CouponArray, VARIANT CouponDates,
								 VARIANT Callput,VARIANT Strike, VARIANT Gearing, VARIANT FixedCoupon,
								 VARIANT* pResult)
{		
	
	// map the inputs
	begin_function


	map_object_parameter(Asset, Asset, hasset);
	map_parameter(maturityDate, long, matDate);
	map_parameter(numberSpacialPoints, long, nx);
	map_parameter(numberTimeSteps, long, nt);
	map_parameter(nstdev, double, fnstdev);

	map_parameter(CouponDates, GVector<long>, couponDates);
	map_parameter(CouponArray, CVector, coupons);

	map_parameter(Callput, double, cp);
	map_parameter(Strike, double, strike);
	map_parameter(FixedCoupon, double, fixed_coupon);
	map_parameter(Gearing, double, gearing);

	map_parameter(referenceSpot, double, refSpot);




	MlEqConstDateHandle	hDate = hasset->GetDateHandle();
	int nToday = hDate->GetDate();
	double spot = hasset->GetSpot(nToday);


	double volbump = 0.01;
	DupireLocalVolHandle lv = new DupireLocalVol;

	double mat = hDate->GetYearFraction(matDate);
	double fwd = hasset->GetForward(matDate, false);
	double vol = hasset->GetVolatility(MlEqStrike(fwd),matDate);

	double highSpot = spot*exp(fnstdev*vol*sqrt(mat));
	double lowSpot  = spot*exp(- 1.5*fnstdev*vol*sqrt(mat));

//	lv->initialize(hasset, lowSpot, highSpot, matDate, nt, nx, volbump );
	lv->initialize(hasset, lowSpot, highSpot, couponDates, nt, nx, volbump );


	RCPtr<callableNote00Helper> callableNote = new callableNote00Helper();
	GVector<long> pdeDates = lv->getDates();

	callableNote->initialize( spot, refSpot, hDate, coupons, couponDates, pdeDates,
								cp, strike, gearing, fixed_coupon);


	
	pdeLocalVol pde;
	pde.init( lv,&(*callableNote));

	CVector result(4);
	CVector tmp(1);

	pde.pde_integrate(tmp);
	double price = tmp[0];
//	double df = tmp[1];
	result[0] = price;

	double delta, gamma ;
	callableNote->greeks( delta, gamma, 0, 0, 0, 0, &pde );
	result[1] = delta ;
	result[2] = gamma ;
	result[2] *= 0.01 * spot ;

	lv->parallelVegaBump();			// vega
	pde.init( lv,&(*callableNote));
	pde.pde_integrate(tmp);
	result[3] = tmp[0] - price ;
//	lv->parallelVegaBump();


	// return the result
	CParameterMap pm;	

	pm.SetValue(result);
	return pm.GetValue(pResult);	
	end_function
}




HRESULT __stdcall MLEqGenericVanillaPricer(VARIANT Asset,VARIANT maturityDate,
										   VARIANT OptionType, VARIANT Callput,VARIANT Strike, 
										   VARIANT LoanspreadBump,
										   VARIANT UpperBarrier,VARIANT LowerBarrier,VARIANT UpperRebate,VARIANT LowerRebate, VARIANT PayEnd,
										   VARIANT* pResult)
{		
	
	// map the inputs
	begin_function


	map_object_parameter(Asset, Asset, hasset);
	map_parameter(maturityDate, long, matDate);
	map_parameter(OptionType, std::string, sotype);
	map_parameter(Callput, double, cp);
	map_parameter(Strike, double, strike);

	map_optional_parameter(LoanspreadBump, double, dls, 0.0001); // 1bp bump

	map_optional_parameter(UpperBarrier, double, upper, 0);
	map_optional_parameter(LowerBarrier, double, lower,0);
	map_optional_parameter(UpperRebate, double, upperR,0);
	map_optional_parameter(LowerRebate, double, lowerR,0);
	map_optional_parameter(PayEnd, bool, payEnd,false);	


	MlEqConstDateHandle	dateh = hasset->GetDateHandle();
	int nToday = dateh->GetDate();
	double spot = hasset->GetSpot(nToday);
	double fwd = hasset->GetForward(matDate, false);
	double mat = dateh->GetYearFraction(matDate);			

	CVector result(7);
	double volbump = 0.01;

	if( sotype == "European" )	// these greeks don't take the smile in account
	{
		double df = hasset->GetPayZeroCurve(true)->GetDiscountFactor(nToday, matDate);

		double eps = 0.001;
		double vol = hasset->GetVolatility(MlEqStrike(strike),matDate);	

		double price  = Bs(fwd, vol, mat, strike, df, cp);
		double price_ = Bs(fwd*(1+eps), vol, mat, strike, df, cp);
		double _price = Bs(fwd*(1-eps), vol, mat, strike, df, cp);

		result[0] = price;
		result[1] = (price_-_price) / (2.*fwd*eps);
		result[2] = (price_+_price - 2.*price) / (eps*eps*fwd*fwd) ;
		result[3] = Bs(fwd, vol+volbump, mat, strike, df, cp) - price;
		
		result[4] = Bs(fwd*exp(dls*mat), vol, mat, strike, df, cp) - price;
	}
	else 
	{

		double vol = hasset->GetVolatility(MlEqStrike(fwd),matDate);	
		double fnstdev  = (mat < 1.)?6.0 - 2.0 * sqrt(mat) : 4.0;

		double highSpot = fwd*exp(fnstdev*vol*sqrt(mat));
		double lowSpot  = fwd*exp(- 1.5*fnstdev*vol*sqrt(mat));

		double nsteps = (mat < 1.5)?250 - 200 * sqrt(mat/1.5) : 50;
		int nt = int(nsteps * mat);
		int nx = nt ;
		
		DupireLocalVolHandle lv = new DupireLocalVol;
		MlEqPdeHelperHandle payoff;

		if( sotype == "Barrier")
		{
			if( lower < 0.001 * fwd ){
				throw"the barrier is too low";
			}

			double mult = log(upper/lower)/log(highSpot/lowSpot) ;
			nx = int( nx * mult ) ;
	//		nt = nx;
	//		nt = MlEqMaths::Min(nt, matDate - nToday - 10 );			
	
			lv->initialize(hasset, lower, upper, matDate, nt, nx, volbump );
			
			CVector df(nt,1.);
			if( payEnd )
			{
				GVector<long> date = lv->getDates();

				for(int t=0; t<nt; t++){
					df[t] = hasset->GetPayZeroCurve(true)->GetDiscountFactor(date[t], matDate);
				}
			}

			payoff = new knockOutBarrierHelper(	upper, lower, upperR, lowerR, spot, cp, strike, payEnd, df);
		}
		else if( sotype == "American" )
		{		
			upper = highSpot;
			lower = lowSpot;

			lv->initialize(hasset, lower, upper, matDate, nt, nx, volbump );

			payoff = new americanOptionHelper(spot, cp, strike);
		}
	
		pdeLocalVol pde;
		pde.init( lv,&(*payoff));

	
		CVector tmp(1);

		pde.pde_integrate(tmp);
		double price = tmp[0];
		result[0] = price;

		double delta, gamma ;
		payoff->greeks( delta, gamma, 0, 0, 0, 0, &pde );
		result[1] = delta ;
		result[2] = gamma ;

		lv->parallelVegaBump();			// vega
		pde.init( lv,&(*payoff));
		pde.pde_integrate(tmp);
		result[3] = tmp[0] - price ;
		lv->parallelVegaBump();

		lv->bumpLocalDrift( dls );		// loanspread risk
		pde.init( lv,&(*payoff));
		pde.pde_integrate(tmp);
		result[4] = tmp[0] - price ;
	}

	result[5] = fwd;
	result[6] = fwd * exp(dls*mat);	// bumped forward

	result[2] *= 0.01 * spot ;	// trader's convention for the gamma

	// return the result
	CParameterMap pm;	

	pm.SetValue(result);
	pm.Transpose();
	return pm.GetValue(pResult);	
	end_function
}

///



HRESULT __stdcall MLEqTestBarrier(VARIANT Asset,VARIANT maturityDate,
										   VARIANT Callput,VARIANT Strike, 
										   VARIANT LevelType,
										   VARIANT BarrierFeat,
										   VARIANT LoanSpreadBump,
										   VARIANT* pResult)
{		
	
	// map the inputs
	begin_function


	map_object_parameter(Asset, Asset, hasset);
	map_parameter(maturityDate, long, matDate);
	map_parameter(Callput, double, cp);
	map_parameter(Strike, double, strike);
	map_optional_parameter(LoanSpreadBump, double, ls_bump, 1e-4);

	map_parameter(LevelType, std::string, lvtype);
	map_parameter(BarrierFeat, CMatrix, bfeat);

	int nfeat = bfeat.rows();
	int ncol = bfeat.cols();
	if( ncol != 7 && ncol != 6)	{
		throw "Wrong number of columns";
	}


	double spot_multiplier = 1.;
	if( lvtype == "Spot Percentage")
	{
		long nToday = hasset->GetDateHandle()->GetDate();
		spot_multiplier = hasset->GetSpot(nToday);
	}
	else if( lvtype != "Absolute" ){
		throw "Wrong level type - 'Spot Percentage' or 'Absolute' ";
	}

	strike *= spot_multiplier;


	GVector<BarrierFeatures> bFeat(nfeat);

	for(int i=0; i<nfeat; ++i)
	{
		BarrierFeatures features;

		features.startDate	= bfeat[i][0] ;
		features.endDate	= bfeat[i][1] ;

		if(bfeat[i][2] < 0.)	// down barrier
		{
			features.isBarrierDown = true;
			features.isKnockOutDown = (int(bfeat[i][3]) == 0) ;
			features.barrDown = bfeat[i][4] * spot_multiplier;
			features.rebateDown = bfeat[i][5] * spot_multiplier;
			features.payDateDn = 0L;
			if( ncol == 7){
				features.payDateDn = long(bfeat[i][6]);
			}				
		}
		else	// down barrier
		{
			features.isBarrierUp = true;
			features.isKnockOutUp = (int(bfeat[i][3]) == 0) ;
			features.barrUp = bfeat[i][4] * spot_multiplier;
			features.rebateUp = bfeat[i][5] * spot_multiplier;
			features.payDateUp = 0L;
			if( ncol == 7){
				features.payDateUp = long(bfeat[i][6]);
			}	
		}

		bFeat[i] = features;
	}

	
	CVector result(5);
	
	RCPtr<BarrierPricerTest> pricer = new BarrierPricerTest(hasset, bFeat, matDate, strike, cp);
	pricer->m_loanspread_bump = ls_bump;
	pricer->solve( result );

	double fwd = hasset->GetForward(matDate, false);
	double mat = hasset->GetDateHandle()->GetYearFraction(matDate);

	CVector fresult(7);
	for(int i=0; i<5;i++){
		fresult[i] = result[i];
	}
	fresult[5] = fwd;
	fresult[6] = fwd * exp( ls_bump*mat );

// return the result
	CParameterMap pm;	

	pm.SetValue(fresult);
//	pm.Transpose();
	return pm.GetValue(pResult);	
	end_function
}





HRESULT __stdcall MLEqVanillaPricer(VARIANT Asset,VARIANT maturityDate, VARIANT CallPut, VARIANT Strike, VARIANT* pResult)
{		
	
	// map the inputs
	begin_function


	map_object_parameter(Asset, Asset, hasset);
	map_parameter(maturityDate, long, matDate);
//	map_parameter(CallPut, std::string, cptype);
//	double cp = (cptype == "call")? 1.:-1.;
	map_parameter(CallPut, double, cp);
	map_parameter(Strike, double, strike);


	MlEqConstDateHandle	dateh = hasset->GetDateHandle();
	int nToday = dateh->GetDate();
	double spot = hasset->GetSpot(nToday);
	double fwd = hasset->GetForward(matDate, false);
	double mat = dateh->GetYearFraction(matDate);			

	CVector result(4);
	double volbump = 0.01;

	double df = hasset->GetPayZeroCurve(true)->GetDiscountFactor(nToday, matDate);

	double eps = 0.001;
	double vol = hasset->GetVolatility(MlEqStrike(strike),matDate);	

	double price  = Bs(fwd, vol, mat, strike, df, cp);
	double price_ = Bs(fwd*(1+eps), vol, mat, strike, df, cp);
	double _price = Bs(fwd*(1-eps), vol, mat, strike, df, cp);

	result[0] = price;
	result[1] = (price_-_price) / (2.*fwd*eps);
	result[2] = (price_+_price - 2.*price) / (eps*eps*fwd*fwd) ;
	result[3] = Bs(fwd, vol+volbump, mat, strike, df, cp) - price;

	result[2] *= 0.01 * spot ;	// trader's convention for the gamma

	// return the result
	CParameterMap pm;	

	pm.SetValue(result);
	pm.Transpose();
	return pm.GetValue(pResult);	
	end_function
}

///


/****************************************************************
**	 
**	Routine: Asian LocalVol pricer
**	Returns: VARIANT
**	Action : calibration
**           
****************************************************************/


HRESULT __stdcall AsianLocalVolPricer(VARIANT asset, VARIANT AverageDates, VARIANT PastFixings, VARIANT CallPut, VARIANT Strike, VARIANT nPaths,  VARIANT nYearlySteps, VARIANT* pResult)
{
	begin_function	

	map_object_parameter(asset, Asset, hasset);
	map_parameter(AverageDates, GVector<long>, avgDates);
	map_parameter(PastFixings, CVector, past_fixings);

	map_parameter(Strike, double, strike);
	map_parameter(CallPut, double, cp);

	map_parameter(nPaths, long, npaths);
	map_parameter(nYearlySteps, long, nsteps);
	int rng = 1; // sobol...
	
	MlEqConstDateHandle hDate = hasset->GetDateHandle();
	
	CMatrix results;

	AsianPricer product;
	product.initialize( strike, cp, avgDates, past_fixings );

	CMatrix fixings(1, past_fixings.getsize());
	fixings[0] = past_fixings;

	RCPtr<CLocalVolMC> lvmc = new CLocalVolMC( hDate, hasset );
	lvmc->computeGreeks();
	lvmc->initialize(product, fixings, npaths, rng, nsteps);
	lvmc->simulate( results, product );

	long nToday = hasset->GetCurrentDate();
	double spot = hasset->GetSpot( nToday );

	double delta = results[1][0] - results[2][0] ;
	double gamma = results[1][0] + results[2][0] - 2.*results[0][0] ;

	double dSpot = (exp(0.01) - exp(-0.01)) * spot * 0.5 ;
	dSpot = 0.01 * spot ;

	results[1][0] = delta / (2.*dSpot);
	results[2][0] = gamma / (dSpot*dSpot) * 0.01 * spot ;

	results[3][0] -= results[0][0];


	CParameterMap pm;
	pm.SetValue(results);
	pm.Transpose();
	return pm.GetValue(pResult);

	end_function
}


HRESULT __stdcall testBump(VARIANT Asset, VARIANT* pResult)
{		
		// map the inputs
	begin_function

	map_object_parameter(Asset, Asset, hasset);


	MlEqVolatilityStructureHandle hVS = hasset->GetVolatilityStructure();

	double parallelShift = 0.1;
	bool isRelative = true;

	hVS->parallelShiftVol( parallelShift,  isRelative );
//	parallelShift = 0.01;
//	hVS->parallelShiftVol(  parallelShift, !isRelative );


	int shiftSkew = 0;

	if( shiftSkew )
	{
		parallelShift = 0.01;

		MlEqConstDateHandle hDate = hasset->GetDateHandle();
		long nToday = hDate->GetDate();

		int nstrike = 20;
		int nmat = 5;

		std::vector<MlEqStrikeHandle> Strikes;
		CVector Maturities(nmat);
		CMatrix Shift(nmat, nstrike);

		double spot = hasset->GetSpot( nToday );

		for(int k=0; k<nstrike; ++k)
		{
			double Ks = 0.4 + 0.06*k ;
		//	MlEqStrikeHandle hStrike = new MlEqSpotBasedStrike(Ks, spot, hDate);
			MlEqStrikeHandle hStrike = new MlEqForwardBasedStrike(Ks, spot, nToday, hDate);
			Strikes.push_back(hStrike);
		}

	
		for(int t=0; t<nmat; ++t)
		{
			long date = nToday + 500 * t + 100;
			Maturities[t] = date;
			double T = hDate->GetYearFraction( date);

			double fwd = hasset->GetForward( date, false );
			double volatm = hasset->getNaturalATMVol( date );

			for(int k=0; k<nstrike; ++k)
			{
				double Ks = Strikes[k]->m_strike;
				double ns =log( spot/fwd * Ks ) / ( volatm*sqrt(T) );			
				double time_decay = 1.;//exp( -T ) ;

	//			ns = -3.0 + 0.3*k ;

				double bump =  - parallelShift * ns * time_decay;
				Shift[t][k] = bump;
			}
		}
		
		hVS->shiftSkew( Strikes, Maturities, Shift, hasset);
	}

	int reset = 0;

	if( reset ){
		hVS->Reset();
	}



	// return the result
	CParameterMap pm;	

	pm.SetValue(0.0);
	return pm.GetValue(pResult);	
	end_function
}


