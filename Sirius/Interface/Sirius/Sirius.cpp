//	Sirius.cpp : Implementation of DLL Exports.
//
//	Note: Proxy/Stub Information
//		  To build a separate proxy/stub DLL, 
//        run nmake -f Siriusps.mk in the project directory.
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "resource.h"
#include <initguid.h>
#include "dlgsetup.h"
#include "dlginsert.h"
#include "dlgmodulesaveload.h"
#include "Sirius_i.c"
#include "siriusapplication.h"
#include "excelinterface.h"

#include "handle_product.h"
#include "handle_products.h"
#include "handle_zerocurve.h"
#include "handle_zerocurves.h"
#include "handle_position.h"
#include "handle_positions.h"
#include "handle_asset.h"
#include "handle_assets.h"
#include "handle_marketdata.h"
#include "handle_marketdatacollection.h"
#include "handle_montecarlo.h"
#include "handle_montecarlocollection.h"
#include "handle_deal.h"
#include "handle_deals.h"
#include "handle_array.h"
#include "handle_arrays.h"
#include "handle_correlationmatrix.h"
#include "handle_correlationmatrices.h"
#include "handle_dividendschedule.h"
#include "handle_dividendschedules.h"
#include "handle_date.h"
#include "handle_dates.h"
#include "handle_interpolator.h"
#include "handle_interpolators.h"
#include "handle_strikes.h"
#include "handle_strikescollection.h"
#include "handle_volatilitystructure.h"
#include "handle_volatilitystructures.h"
#include "handle_result.h"
#include "handle_results.h"
#include "handle_parameterlist.h"
#include "handle_parameterlists.h"
#include "handle_voldata.h"
#include "handle_voldatacollection.h"
#include "handle_jumpwingvoldata.h"
#include "handle_jumpwingvoldatacollection.h"
#include "handle_dateschedule.h"
#include "handle_dateschedules.h"
#include "handle_spotschedule.h"
#include "handle_spotschedules.h"
#include "handle_svlvoldata.h"
#include "handle_svlvoldatacollection.h"
#include "handle_svlfitvoldata.h"
#include "handle_svlfitvoldatacollection.h"
#include "handle_hullvoldata.h"
#include "handle_hullvoldatacollection.h"
#include "handle_ramvoldata.h"
#include "handle_ramvoldatacollection.h"
#include "handle_ramfitvoldata.h"
#include "handle_ramfitvoldatacollection.h"
#include "handle_arpropvoldata.h"
#include "handle_arpropvoldatacollection.h"
#include "handle_arpropfitvoldata.h"
#include "handle_arpropfitvoldatacollection.h"
#include "objectmanager.h"
#include "product_swap.h"
#include "parameter.h"
#include "parameters.h"
#include "siriusapplication.h"
#include "xmlstreamer.h"
#include "handle_ramvoldata.h"
#include "progress.h"
#include "handle_multipliervoldata.h"
#include "handle_multipliervoldatacollection.h"
#include "handle_hullandwhite.h"
#include "handle_hullandwhitecollection.h"
#include "interface_setup.h"
#include "handle_matrix.h"
#include "handle_matrices.h"
#include "handle_jumpwingfitvoldata.h"
#include "handle_jumpwingfitvoldatacollection.h"
#include "handle_hermitevoldata.h"
#include "handle_hermitevoldatacollection.h"
#include "handle_scenario.h"
#include "handle_scenarios.h"
#include "ResultsCollection.h"
#include "handle_pde.h"
#include "handle_pdecollection.h"

CSiriusComModule						_Module;

BEGIN_OBJECT_MAP(ObjectMap)
	OBJECT_ENTRY(CLSID_ZeroCurve, CZeroCurve)
	OBJECT_ENTRY(CLSID_ZeroCurves, CZeroCurves)
	OBJECT_ENTRY(CLSID_Position, CPosition)
	OBJECT_ENTRY(CLSID_Positions, CPositions)
	OBJECT_ENTRY(CLSID_Asset, CAsset)
	OBJECT_ENTRY(CLSID_Assets, CAssets)
	OBJECT_ENTRY(CLSID_Product, CProduct)
	OBJECT_ENTRY(CLSID_Products, CProducts)
	OBJECT_ENTRY(CLSID_MarketData, CMarketData)
	OBJECT_ENTRY(CLSID_MarketDataCollection, CMarketDataCollection)
	OBJECT_ENTRY(CLSID_MonteCarlo, CMonteCarlo)
	OBJECT_ENTRY(CLSID_MonteCarloCollection, CMonteCarloCollection)
	OBJECT_ENTRY(CLSID_Deal, CDeal)
	OBJECT_ENTRY(CLSID_Deals, CDeals)
	OBJECT_ENTRY(CLSID_Array, CArray)
	OBJECT_ENTRY(CLSID_Arrays, CArrays)
	OBJECT_ENTRY(CLSID_Date, CDate)
	OBJECT_ENTRY(CLSID_Dates, CDates)
	OBJECT_ENTRY(CLSID_Interpolator, CInterpolator)
	OBJECT_ENTRY(CLSID_Interpolators, CInterpolators)
	OBJECT_ENTRY(CLSID_Strikes, CStrikes)
	OBJECT_ENTRY(CLSID_StrikesCollection, CStrikesCollection)
	OBJECT_ENTRY(CLSID_VolatilityStructure, CVolatilityStructure)
	OBJECT_ENTRY(CLSID_VolatilityStructures, CVolatilityStructures)
	OBJECT_ENTRY(CLSID_Result, CResult)
	OBJECT_ENTRY(CLSID_Results, CResults)
	OBJECT_ENTRY(CLSID_ParameterList, CParameterList)
	OBJECT_ENTRY(CLSID_ParameterLists, CParameterLists)
	OBJECT_ENTRY(CLSID_VolData, CVolData)
	OBJECT_ENTRY(CLSID_VolDataCollection, CVolDataCollection)
	OBJECT_ENTRY(CLSID_JumpWingVolData, CJumpWingVolData)
	OBJECT_ENTRY(CLSID_JumpWingVolDataCollection, CJumpWingVolDataCollection)
	OBJECT_ENTRY(CLSID_SVLVolData, CSVLVolData)
	OBJECT_ENTRY(CLSID_SVLVolDataCollection, CSVLVolDataCollection)
	OBJECT_ENTRY(CLSID_SVLFitVolData, CSVLFitVolData)
	OBJECT_ENTRY(CLSID_SVLFitVolDataCollection, CSVLFitVolDataCollection)
	OBJECT_ENTRY(CLSID_HullVolData, CHullVolData)
	OBJECT_ENTRY(CLSID_HullVolDataCollection, CHullVolDataCollection)
	OBJECT_ENTRY(CLSID_Swap, CSwap)
	OBJECT_ENTRY(CLSID_Parameter, CParameter)
	OBJECT_ENTRY(CLSID_Parameters, CParameters)
	OBJECT_ENTRY(CLSID_SiriusApplication, CSiriusApplication)
	OBJECT_ENTRY(CLSID_DateSchedule, CDateSchedule)
	OBJECT_ENTRY(CLSID_DateSchedules, CDateSchedules)
	OBJECT_ENTRY(CLSID_DividendSchedule, CDividendSchedule)
	OBJECT_ENTRY(CLSID_DividendSchedules, CDividendSchedules)
	OBJECT_ENTRY(CLSID_SpotSchedule, CSpotSchedule)
	OBJECT_ENTRY(CLSID_SpotSchedules, CSpotSchedules)
	OBJECT_ENTRY(CLSID_CorrelationMatrix, CCorrelationMatrix)
	OBJECT_ENTRY(CLSID_CorrelationMatrices, CCorrelationMatrices)
	OBJECT_ENTRY(CLSID_RamVolData, CRamVolData)
	OBJECT_ENTRY(CLSID_RamVolDataCollection, CRamVolDataCollection)
	OBJECT_ENTRY(CLSID_RamFitVolData, CRamFitVolData)
	OBJECT_ENTRY(CLSID_RamFitVolDataCollection, CRamFitVolDataCollection)
	OBJECT_ENTRY(CLSID_Progress, CProgress)
	OBJECT_ENTRY(CLSID_ARPropVolData, CARPropVolData)
	OBJECT_ENTRY(CLSID_ARPropVolDataCollection, CARPropVolDataCollection)
	OBJECT_ENTRY(CLSID_ARPropFitVolData, CARPropFitVolData)
	OBJECT_ENTRY(CLSID_ARPropFitVolDataCollection, CARPropFitVolDataCollection)
	OBJECT_ENTRY(CLSID_MultiplierVolData, CMultiplierVolData)
	OBJECT_ENTRY(CLSID_MultiplierVolDataCollection, CMultiplierVolDataCollection)
	OBJECT_ENTRY(CLSID_HullAndWhite, CHullAndWhite)
	OBJECT_ENTRY(CLSID_HullAndWhiteCollection, CHullAndWhiteCollection)
	OBJECT_ENTRY(CLSID_Setup, CSetup)
	OBJECT_ENTRY(CLSID_Matrix, CInterfaceMatrix)
	OBJECT_ENTRY(CLSID_Matrices, CInterfaceMatrices)
	OBJECT_ENTRY(CLSID_JumpWingFitVolData, CJumpWingFitVolData)
	OBJECT_ENTRY(CLSID_JumpWingFitVolDataCollection, CJumpWingFitVolDataCollection)
	OBJECT_ENTRY(CLSID_HermiteVolData, CHermiteVolData)
	OBJECT_ENTRY(CLSID_HermiteVolDataCollection, CHermiteVolDataCollection)
	OBJECT_ENTRY(CLSID_Scenario, CScenario)
	OBJECT_ENTRY(CLSID_Scenarios, CScenarios)
	OBJECT_ENTRY(CLSID_ResultsCollection, CResultsCollection)
	OBJECT_ENTRY(CLSID_Pde, CPde)
	OBJECT_ENTRY(CLSID_PdeCollection, CPdeCollection)
END_OBJECT_MAP()


/////////////////////////////////////////////////////////////////////////////
//	Custom Casts
//
/*static*/ void CObjectManager::Cast(const CParameterMap& pmIn, std::vector<CAdapt<CComQIPtr<IDispatch> > >& asp, std::vector<MlEqVolDataHandle>& apOut)
{
	std::vector<MlEqVolData*> ap;	
	if (pmIn.IsBlank()) return;
	GetObjectVectorStart(pmIn, asp, ap, MlEqVolData)
		GetObjectVectorCast(VolData)
		GetObjectVectorCast(JumpWingVolData)
		GetObjectVectorCast(JumpWingFitVolData)
		GetObjectVectorCast(SVLVolData)
		GetObjectVectorCast(SVLFitVolData)
		GetObjectVectorCast(HullVolData)
		GetObjectVectorCast(RamVolData)
		GetObjectVectorCast(RamFitVolData)
		GetObjectVectorCast(ARPropVolData)
		GetObjectVectorCast(ARPropFitVolData)
		GetObjectVectorCast(MultiplierVolData)
		GetObjectVectorCast(HermiteVolData)
	GetObjectVectorEnd(IDS_INVALID_VOLATILITY_DATA_HANDLE, IID_IVolatilityStructure)
	apOut.resize(ap.size());
	for (long n = 0; n < ap.size(); n++){
		apOut[n] = ap[n];
	}	
}


/////////////////////////////////////////////////////////////////////////////
//	Create
//
//	Creates an object of name szObjectName from the input data, inserts it into
//	the appropriate collection.
//
//	Only ever call this from an MLCreate[...] function since we create an
//  external handle.
//
void Create(const std::string& szObjectName, VARIANT* pResult, int nArgs, /*VARIANT*/...)
{
	CComVariant							vData;
	CComBSTR							sHandle;
	CParameterMap						pm;
	va_list								vl;
	HRESULT								hr;
	CLSID								clsid;
	CComPtr<IDispatch>					spObject;

	if (g_pApplication->GetObjectManager().NameToCLSID(szObjectName, &clsid)) throw "Error creating object '" + szObjectName + "'";	
	va_start(vl, nArgs);
	hr = CParameterMap::VariableArgumentListToArray(&vData, nArgs, vl);
	va_end(vl);
	if (hr) propagate_error;
	if (spObject.CoCreateInstance(clsid)) throw "Error creating object '" + szObjectName + "'";
	if (CComDispatchDriverEx(spObject).PutPropertyByName(L"Value", &vData)) propagate_error;	
	CExcelInterface::Flatten(spObject, pResult);
}


///////////////////////////////////////////////////////////////////////////////
//	DisplayExplorer
//
//	Launches Sirius Explorer.
//
HRESULT __stdcall DisplayExplorer(void)
{
	begin_function
	_Module.DisplayExplorer();
	end_function
}


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
		_Module.Init(ObjectMap, hInstance, &LIBID_Sirius);
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


///////////////////////////////////////////////////////////////////////////////
//	EnumInclude
//
//	This function is used to force various enum definitions to be included in
//  sirius.tlb. It may as well do something slightly better than useless like
//  return the number of parameters passed in!
//
HRESULT __stdcall EnumInclude(PerformanceTypeEnum, PayoffTypeEnum, YesNoEnum, NotionalTypeEnum, MinMaxEnum, SpotOrForwardEnum, long* pn)
{
	*pn = 7L;
	return S_OK;
}


/////////////////////////////////////////////////////////////////////////////
//	InitialiseExcel
//
//	Pass a pointer to the IDispatch interface for the Excel process
//	associated with this DLL. And set up the Excel menu etc.
//
HRESULT __stdcall InitialiseExcel(VARIANT ExcelApplicationInstance)
{
	return CExcelInterface::InitialiseExcel(ExcelApplicationInstance);	
}


///////////////////////////////////////////////////////////////////////////////
//	InsertObject
//
//	handles the Sirius->InsertObject menu option
//
HRESULT __stdcall InsertObject()
{
	CDlgInsert().DoModal(::GetActiveWindow(), CDlgInsert::insert_object);
	return S_OK;
}


///////////////////////////////////////////////////////////////////////////////
//	InsertProduct
//
//	handles the Sirius->InsertProduct menu option
//
HRESULT __stdcall InsertProduct()
{
	CDlgInsert().DoModal(::GetActiveWindow(), CDlgInsert::insert_product);
	return S_OK;
}


/////////////////////////////////////////////////////////////////////////////
//	IsExcelFunctionWizardShowing
//
//	returns TRUE if the excel function wizard is visible
//
HRESULT __stdcall IsExcelFunctionWizardShowing(BOOL* pVal)
{
	*pVal = CParameterMap::IsExcelFunctionWizardVisible();
	return S_OK;
}


//////////////////////////////////////////////////////////////////////////////
//	LoadVBAModule
//
//	impliement the Sirius->LoadVBAModule menu option
//
HRESULT __stdcall LoadVBAModule()
{
	CDlgModuleLoad().DoModal();
	return S_OK;
}


//////////////////////////////////////////////////////////////////////////////
//	Max
//
//	Returns the largest number in a set of values.
//
HRESULT __stdcall Max(VARIANT Number1, VARIANT Number2, double* pVal)
{	
	HRESULT								hr1, hr2;
	CParameterMap						pm1, pm2;
	double								f1, f2;
	
	pm1.SetValue(Number1);
	pm2.SetValue(Number2);
	hr1 = pm1.Max(&f1);
	hr2 = pm2.Max(&f2);

	if (!hr1 && !hr2){
		*pVal = __max(f1, f2);
	} else if (!hr1){
		*pVal = f1;
	} else if (!hr2){
		*pVal = f2;
	} else {
		return CParameterMap::ReturnErrorR(IDS_NO_NUMERIC_DATA_FOUND);
	}
	return S_OK;
}


//////////////////////////////////////////////////////////////////////////////
//	Min
//
//	Returns the smallest number in a set of values.
//
HRESULT __stdcall Min(VARIANT Number1, VARIANT Number2, double* pVal)
{
	HRESULT								hr1, hr2;
	CParameterMap						pm1, pm2;
	double								f1, f2;	

	pm1.SetValue(Number1);
	pm2.SetValue(Number2);
	hr1 = pm1.Min(&f1);
	hr2 = pm2.Min(&f2);

	if (!hr1 && !hr2){
		*pVal = __min(f1, f2);
	} else if (!hr1){
		*pVal = f1;
	} else if (!hr2){
		*pVal = f2;
	} else {
		return CParameterMap::ReturnErrorR(IDS_NO_NUMERIC_DATA_FOUND);
	}
	return S_OK;
}


///////////////////////////////////////////////////////////////////////////////
//	MLAddBusinessDays
//
//	Adds the specified number of business days to given date according to the
//	optional holiday calendar.
//
HRESULT __stdcall MLAddBusinessDays(VARIANT StartDate, VARIANT Days, VARIANT CalendarOpt, VARIANT* pResult)
{
	begin_function

	map_parameter(StartDate, long, nStartDate);
	map_parameter(Days, long, nDays);
	map_optional_parameter(CalendarOpt, estring, szCalendar, "");
	
	MlEqDateHandle hDate = MlEqDate(nStartDate).AddBusinessDays(nDays, szCalendar);
	return CComVariant(hDate->GetDate()).Detach(pResult);
	end_function
}


///////////////////////////////////////////////////////////////////////////////
//	MLCalibrate
//
//	Calibrates an uncalibrated object
//
HRESULT __stdcall MLCalibrate(VARIANT Handle, VARIANT* pResult)
//	Handle - an appropriate object (e.g. MonteCarlo handle or an Asset handle)
{	
	begin_function
	
	MlEqMonteCarloHandle				hMonteCarlo;
	MlEqAssetHandle						hAsset;

	// First try a MonteCarlo handle
	try {
		map_object_parameter(Handle, MonteCarlo, h);
		hMonteCarlo = h;
	} catch (...){
		// Try an Asset handle
		map_object_parameter(Handle, Asset, h);
		hAsset = h;
	}

	if (!!hMonteCarlo){	
		CParameterMap pm;
		CForwardSkewMC* pMC = dynamic_cast<CForwardSkewMC*>(&*hMonteCarlo);
		if (!pMC) throw "Cannot return the calibation results for this type of Monte Carlo";
		CMatrix calibResults = pMC->m_calibResults;
		if (!calibResults.rows()){
			calibResults.resize(1,1);
		}
		pm.SetValue(calibResults);
		return pm.GetValue(pResult);
	} else if (!!hAsset){		
		throw "Direct calibration of an asset is not supported in this release";		
	} else {
		ATLASSERT(false);	// Should never get here
	}
	
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	Temporary place for PDE creation functions
//
HRESULT __stdcall MLCreateBlackScholesPde(VARIANT Vol, VARIANT LoanSpread, VARIANT DiscountRate, VARIANT Spot, VARIANT* pResult)
{
	begin_function
	::Create("Pde", pResult, 5, CComVariant(BlackScholes), Vol, LoanSpread, DiscountRate, Spot);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLCreate[...] - these are the create object stubs
//
HRESULT __stdcall MLCreate(VARIANT ObjectName, VARIANT Arg1, VARIANT Arg2, VARIANT Arg3, VARIANT Arg4, VARIANT Arg5, VARIANT Arg6, VARIANT Arg7, VARIANT Arg8, VARIANT Arg9, VARIANT Arg10, VARIANT* pResult)
{			
	CParameterMap						pmObjectName;
	std::string							szObjectName;
		
	begin_function
	if (pmObjectName.SetValue(ObjectName)) return E_FAIL;	
	if (pmObjectName.GetString(&szObjectName)) CParameterMap::ThrowComErrorRS(IDS_CREATE_OBJECT, szObjectName);
	::Create(szObjectName, pResult, 10, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6, Arg7, Arg8, Arg9, Arg10);
	end_function
}
HRESULT __stdcall MLCreateARPropVolData(VARIANT AtmVolLambda,VARIANT AtmVolS0, VARIANT AtmVolSInfty, VARIANT SkewVolLambda, VARIANT SkewVolS0, VARIANT SkewVolInfty, VARIANT curvature, VARIANT timeCutoff, VARIANT fwds, VARIANT dateToDouble, VARIANT dates,VARIANT* pResult)
{
	begin_function
	::Create("ARPropVolData", pResult, 11, AtmVolLambda, AtmVolS0, AtmVolSInfty, SkewVolLambda, SkewVolS0, SkewVolInfty, curvature, timeCutoff, fwds, dateToDouble, dates);
	end_function
}
HRESULT __stdcall MLCreateARPropFitVolData(VARIANT VolatilityStructure, VARIANT FittingStrikes, VARIANT DateToDouble,VARIANT timeCutoff, VARIANT fittingDates,VARIANT fwds,VARIANT curvature, VARIANT* pResult)
{
	begin_function
	::Create("ARPropFitVolData", pResult, 7, VolatilityStructure, FittingStrikes, DateToDouble, timeCutoff, fittingDates, fwds, curvature);
	end_function
}
HRESULT __stdcall MLCreateAsset(VARIANT ParameterNamesArr, VARIANT ParameterValuesArr, VARIANT* pResult)
{
	begin_function
	::Create("Asset", pResult, 2, ParameterNamesArr, ParameterValuesArr);
	end_function
}
HRESULT __stdcall MLCreateCorrelationMatrix(VARIANT DateOpt, VARIANT DataSource, VARIANT Matrix, VARIANT* pResult)
{
	begin_function
	::Create("CorrelationMatrix", pResult, 3, DateOpt, DataSource, Matrix);
	end_function
}
HRESULT __stdcall MLCreateCubicSplineInterpolator(VARIANT xData, VARIANT yData, VARIANT addTanhWing, VARIANT cL, VARIANT cR, VARIANT yPower, VARIANT* pResult)
{
	begin_function
	::Create("Interpolator", pResult, 7, CComVariant(CubicSpline), xData, yData, addTanhWing, cL, cR, yPower);
	end_function
}
HRESULT __stdcall MLCreateMonotonicSplineInterpolator(VARIANT xData, VARIANT yData, VARIANT addTanhWing, VARIANT cL, VARIANT cR, VARIANT yPower, VARIANT* pResult)
{
	begin_function
	::Create("Interpolator", pResult, 7, CComVariant(MonotonicSpline), xData, yData, addTanhWing, cL, cR, yPower);
	end_function
}
HRESULT __stdcall MLCreateConstantInterpolator(VARIANT xData, VARIANT yData, VARIANT addTanhWingOpt, VARIANT cLOpt, VARIANT cROpt, VARIANT yPowerOpt, VARIANT* pResult)
{
	begin_function
	::Create("Interpolator", pResult, 7, CComVariant(Constant), xData, yData, addTanhWingOpt, cLOpt, cROpt, yPowerOpt);
	end_function
}
HRESULT __stdcall MLCreateFitPolynomialInterpolator(VARIANT xData, VARIANT yData,VARIANT initialGuess,VARIANT xValBounds, VARIANT addTanhWing, VARIANT cL, VARIANT cR, VARIANT yPower,VARIANT finalTolOpt,VARIANT stoppingTolOpt, VARIANT* pResult)
{
	begin_function
	::Create("Interpolator", pResult, 11, CComVariant(FitPolynomial),xData, yData,initialGuess,xValBounds, addTanhWing, cL, cR, yPower,finalTolOpt,stoppingTolOpt); 
	end_function
}
HRESULT __stdcall MLCreateFitPolynomialInterpolatorExplicit(VARIANT coefficients, VARIANT lowerEdge,VARIANT upperEdge,VARIANT addTanhWing,VARIANT cL, VARIANT cR, VARIANT yPower,VARIANT npointsOpt, VARIANT* pResult)
{
	begin_function
	::Create("Interpolator", pResult, 9, CComVariant(FitPolynomialExlicit),coefficients, lowerEdge,upperEdge, addTanhWing, cL, cR, yPower,npointsOpt); 
	end_function
}
HRESULT __stdcall MLCreateTwoDimensionalInterpolator(VARIANT InterpolatorAcrossY, VARIANT InterpolatorAcrossX, VARIANT* pResult)
{
	begin_function
	::Create("Interpolator", pResult, 3, CComVariant(TwoDimensional), InterpolatorAcrossY, InterpolatorAcrossX);
	end_function
}
HRESULT __stdcall MLCreateDate(VARIANT Date, VARIANT DayCountConvention, VARIANT* pResult)
{
	begin_function
	::Create("Date", pResult, 2, Date, DayCountConvention);
	end_function
}
HRESULT __stdcall MLCreateDateSchedule(VARIANT DatesArr, VARIANT ValuesArr, VARIANT* pResult)
{
	begin_function
	::Create("DateSchedule", pResult, 2, DatesArr, ValuesArr);
	end_function
}
HRESULT __stdcall MLCreateDividendSchedule(VARIANT Name, VARIANT DateOpt, VARIANT DataSource, VARIANT DayCountConventionOpt, VARIANT YieldCurveOpt, VARIANT Dates, VARIANT DiscreteDividendsArr, VARIANT ContinuousDividendsArr, VARIANT* pResult)
{
	begin_function
	::Create("DividendSchedule", pResult, 8, Name, DateOpt, DataSource, DayCountConventionOpt, YieldCurveOpt, Dates, DiscreteDividendsArr, ContinuousDividendsArr);
	end_function
}
HRESULT __stdcall MLCreateFixedStrikes(VARIANT Strikes, VARIANT* pResult)
{
	begin_function
	::Create("Strikes", pResult, 2, CComVariant(Fixed), Strikes);
	end_function
}
HRESULT __stdcall MLCreateForwardBasedStrikes(VARIANT Strikes, VARIANT ReferenceForward, VARIANT ReferenceMaturityDate, VARIANT StartDateHandleOpt, VARIANT IsFloatingOpt, VARIANT* pResult)
{
	begin_function
	::Create("Strikes", pResult, 6, CComVariant(ForwardBased), Strikes, ReferenceForward, ReferenceMaturityDate, StartDateHandleOpt, IsFloatingOpt);
	end_function
}
HRESULT __stdcall MLCreateForwardBasedStrikesFromAsset( VARIANT Strikes, VARIANT Asset, VARIANT ReferenceMaturityDate, VARIANT StartDateHandleOpt, VARIANT IsFloatingOpt, VARIANT* pResult)
{
	begin_function
	::Create("Strikes", pResult, 7, CComVariant(CStrikes::ForwardBasedFromAsset), Strikes, Asset, ReferenceMaturityDate, StartDateHandleOpt, IsFloatingOpt);
	end_function
}
HRESULT __stdcall MLCreateForwardSkewMonteCarlo(VARIANT ForwardVolatilities, VARIANT Betas, VARIANT VolVols, VARIANT MeanReversionRates, VARIANT MeanReversionLevels, VARIANT ModelInfo, VARIANT Strikes, VARIANT CalibTimeIndex, VARIANT CalibVolatilities, VARIANT CalibStrikes, VARIANT Spot, VARIANT ControlVariatePricesOpt, VARIANT ExcessVolVolInterpolatorOpt, VARIANT ZeroCurveOpt, VARIANT* pResult)
{
	begin_function
	::Create("MonteCarlo", pResult, 15, CComVariant(ForwardSkew), ForwardVolatilities, Betas, VolVols, MeanReversionRates, MeanReversionLevels, ModelInfo, Strikes, CalibTimeIndex, CalibVolatilities, CalibStrikes, Spot, ControlVariatePricesOpt, ExcessVolVolInterpolatorOpt, ZeroCurveOpt);
	end_function
}
HRESULT __stdcall MLCreateGeneralMonteCarlo(VARIANT Assets, VARIANT Dates, VARIANT ModelParameters, VARIANT HullAndWhiteModels,VARIANT* pResult)
{
	begin_function
	::Create("MonteCarlo", pResult, 5, CComVariant(General), Assets, Dates, ModelParameters,HullAndWhiteModels);
	end_function
}
HRESULT __stdcall MLCreateQuasiMonteCarlo(VARIANT SpotSchedulesArr, VARIANT Dates, VARIANT* pResult)
{
	begin_function
	::Create("MonteCarlo", pResult, 3, CComVariant(Quasi), SpotSchedulesArr, Dates);
	end_function
}
HRESULT __stdcall MLCreateCalibratedForwardSkewMonteCarlo(VARIANT Assets, VARIANT Dates, VARIANT ModelParameters, VARIANT calibTimeIndex, VARIANT calibStrikes, VARIANT* pResult)
{
	begin_function
	::Create("MonteCarlo", pResult, 6, CComVariant(CMonteCarlo::CalibrateForwardSkew), Assets, Dates, ModelParameters,  calibTimeIndex,calibStrikes);
	end_function
}
HRESULT __stdcall MLCreateHybridForwardSkewMonteCarlo(VARIANT Asset, VARIANT Dates, VARIANT ModelParameters, VARIANT HullAndWhiteModels, VARIANT Correlation,VARIANT isMinimalHybrid, VARIANT* pResult)
{
	begin_function
	::Create("MonteCarlo", pResult, 7, CComVariant(HybridForwardSkew), Asset, Dates, ModelParameters, HullAndWhiteModels, Correlation,isMinimalHybrid);
	end_function
}
HRESULT __stdcall MLCreateCalibratedHybridMonteCarlo(VARIANT Asset, VARIANT Dates, VARIANT ModelParameters, VARIANT HullAndWhiteModels, VARIANT Correlation,VARIANT isMinimalHybrid, VARIANT calibTimeIndex, VARIANT calibStrikes,VARIANT* pResult)
{
	begin_function
		::Create("MonteCarlo", pResult, 9, CComVariant(CMonteCarlo::CalibrateHybridModel), Asset, Dates, ModelParameters, HullAndWhiteModels, Correlation,isMinimalHybrid,calibTimeIndex,  calibStrikes);
	end_function
}
HRESULT __stdcall MLCreateLocalVolMonteCarlo(VARIANT Asset, VARIANT Dates, VARIANT ModelParameters,VARIANT* pResult)
{
	begin_function
	::Create("MonteCarlo", pResult, 4, CComVariant(LocalVolMonteCarlo), Asset, Dates, ModelParameters);
	end_function
}
HRESULT __stdcall MLCreateHermiteVolData(VARIANT StartDateHandle, VARIANT DatesArr,  VARIANT HermiteParameters,VARIANT ForwardsArrOpt, VARIANT* pResult)
{
	begin_function
	::Create("HermiteVolData", pResult, 4, StartDateHandle, DatesArr,  HermiteParameters,ForwardsArrOpt);
	end_function
}
HRESULT __stdcall MLCreateHullWhiteModel(VARIANT Dates, VARIANT lambda,VARIANT discountCurve,VARIANT dateToDouble,VARIANT capletVols,VARIANT capletStrikes,VARIANT swaptionVols,VARIANT swaptionStrikes,VARIANT swaptionDates, VARIANT dayCount, VARIANT* pResult)
{
	begin_function
	::Create("HullAndWhite", pResult, 10,  Dates,lambda,discountCurve,dateToDouble,capletVols,capletStrikes,swaptionVols,swaptionStrikes,swaptionDates,dayCount);
	end_function
}
HRESULT __stdcall MLCreateHullVolData(VARIANT Strikes, VARIANT HullParameters, VARIANT* pResult)
{
	begin_function
	::Create("HullVolData", pResult, 2, Strikes, HullParameters);
	end_function
}
HRESULT __stdcall MLCreateRamVolData(VARIANT Strikes, VARIANT BaseVolArr, VARIANT SlopeArr, VARIANT CurvatureArr, VARIANT DownCutoffArr, VARIANT UpCutoffArr, VARIANT* pResult)
{
	begin_function
	::Create("RamVolData", pResult, 6, Strikes, BaseVolArr, SlopeArr, CurvatureArr, DownCutoffArr, UpCutoffArr);
	end_function
}
HRESULT __stdcall MLCreateJumpWingFitVolData(VARIANT MktVols, VARIANT GuessJumpWingVolData,VARIANT FittingDates,VARIANT FitMidPoints,  VARIANT FitSpread,VARIANT UseVega,  VARIANT* pResult)
{
	begin_function
	::Create("JumpWingFitVolData", pResult, 6, MktVols, GuessJumpWingVolData, FittingDates, FitMidPoints, FitSpread, UseVega);
	end_function
}
HRESULT __stdcall MLCreateJumpWingVolData(VARIANT Strikes, VARIANT JumpWingParameters, VARIANT* pResult)
{
	begin_function
	::Create("JumpWingVolData", pResult, 2, Strikes, JumpWingParameters);
	end_function
}
HRESULT __stdcall MLCreateLinearInterpolator(VARIANT xData, VARIANT yData,  VARIANT addTanhWingOpt, VARIANT cLOpt, VARIANT cROpt, VARIANT yPowerOpt,VARIANT* pResult)
{
	begin_function
	::Create("Interpolator", pResult, 7, CComVariant(Linear), xData, yData, addTanhWingOpt, cLOpt, cROpt, yPowerOpt);
	end_function
}
HRESULT __stdcall MLCreateLogBasedStrikes(VARIANT Strikes, VARIANT ForwardValue, VARIANT ReferenceMaturityDate, VARIANT StartDateHandleOpt, VARIANT IsFloatingOpt, VARIANT* pResult)
{
	begin_function
	::Create("Strikes", pResult, 6, CComVariant(LogBased), Strikes, ForwardValue, ReferenceMaturityDate, StartDateHandleOpt, IsFloatingOpt);
	end_function
}
HRESULT __stdcall MLCreateLogBasedStrikesFromAsset(VARIANT Strikes,  VARIANT Asset,VARIANT ReferenceMaturityDate, VARIANT StartDateHandleOpt, VARIANT IsFloatingOpt, VARIANT* pResult)
{
	begin_function
	::Create("Strikes", pResult, 6, CComVariant(CStrikes::LogBasedFromAsset), Strikes, Asset, ReferenceMaturityDate, StartDateHandleOpt, IsFloatingOpt);
	end_function
}
HRESULT __stdcall MLCreateNormalisedStrikes(VARIANT Strikes, VARIANT ReferenceForward, VARIANT ReferenceVolatility, VARIANT ReferenceMaturityDate, VARIANT StartDateHandleOpt, VARIANT IsFloatingOpt, VARIANT* pResult)
{
	begin_function
	::Create("Strikes", pResult, 7, CComVariant(Normalised), Strikes, ReferenceForward, ReferenceVolatility, ReferenceMaturityDate, StartDateHandleOpt, IsFloatingOpt);
	end_function
}
HRESULT __stdcall MLCreateNormalisedStrikesFromAsset(VARIANT Strikes, VARIANT Asset, VARIANT MaturityDate, VARIANT StartDateHandleOpt, VARIANT IsFloatingOpt, VARIANT* pResult)
{
	begin_function
	::Create("Strikes", pResult, 6, CComVariant(CStrikes::NormalisedFromAsset), Strikes, Asset, MaturityDate, StartDateHandleOpt, IsFloatingOpt);
	end_function
}
HRESULT __stdcall MLCreateProduct(VARIANT ProductType, VARIANT DataSourceOpt, VARIANT DateOpt, VARIANT ParameterNamesArr, VARIANT ParameterValuesArr, VARIANT* pResult)
{
	begin_function
	::Create("Product", pResult, 5, ProductType, DataSourceOpt, DateOpt, ParameterNamesArr, ParameterValuesArr);
	end_function
}
HRESULT __stdcall MLCreateRamFitVolData(VARIANT MarketVolatilityStructure, VARIANT FittingStrikes, VARIANT UpperCutoffStrikes,VARIANT LowerCutoffStrikes,VARIANT fittingDates, VARIANT fitCutoffs, VARIANT useVegaWeights,VARIANT fitSpread,VARIANT fitMidVols,VARIANT* pResult)
{
	begin_function
	::Create("RamFitVolData", pResult, 9, MarketVolatilityStructure,FittingStrikes,UpperCutoffStrikes,LowerCutoffStrikes,fittingDates, fitCutoffs,  useVegaWeights,fitSpread,fitMidVols);
	end_function
}
HRESULT __stdcall MLCreateSpotSchedule(VARIANT Name, VARIANT Date, VARIANT DataSource, VARIANT Dates, VARIANT Spots, VARIANT* pResult)
{
	begin_function
	::Create("SpotSchedule", pResult, 5, Name, Date, DataSource, Dates, Spots);
	end_function
}
HRESULT __stdcall MLCreateSpotBasedStrikes(VARIANT Strikes, VARIANT ReferenceSpot, VARIANT StartDateHandleOpt, VARIANT IsFloatingOpt, VARIANT* pResult)
{
	begin_function
	::Create("Strikes", pResult, 5, CComVariant(SpotBased), Strikes, ReferenceSpot, StartDateHandleOpt, IsFloatingOpt);
	end_function
}
HRESULT __stdcall MLCreateSpotBasedStrikesFromAsset(VARIANT Strikes, VARIANT Asset, VARIANT StartDateHandleOpt, VARIANT IsFloatingOpt, VARIANT* pResult)
{
	begin_function
	::Create("Strikes", pResult, 5, CComVariant(CStrikes::SpotBasedFromAsset), Strikes, Asset,StartDateHandleOpt, IsFloatingOpt);
	end_function
}
HRESULT __stdcall MLCreateStochasticBetaVolatilityStructure(VARIANT Name, VARIANT DataSourceAndOrLocation, VARIANT Spot, VARIANT Forwards, VARIANT Discounts, VARIANT Times, VARIANT DateHandle, VARIANT InterpAcrossTime, VARIANT BidVolData, VARIANT MidVolData, VARIANT AskVolData, VARIANT InterpolateType, VARIANT VolatilityDataType, VARIANT AsymptoticVolatilityStructure, VARIANT DecayFactor, VARIANT StochasticBetaSchedules, VARIANT ExcessVolVol2DInterpolatorOpt, VARIANT HermiteElasticityOpt, VARIANT* pResult)
{
	begin_function
	::Create("VolatilityStructure", pResult, 18, Name, DataSourceAndOrLocation, Spot, Forwards, Discounts, Times, DateHandle, InterpAcrossTime, BidVolData, MidVolData, AskVolData, InterpolateType, VolatilityDataType, AsymptoticVolatilityStructure, DecayFactor, StochasticBetaSchedules, ExcessVolVol2DInterpolatorOpt, HermiteElasticityOpt);
	end_function
}
HRESULT __stdcall MLCreateSvlVolData(VARIANT Strikes, VARIANT Parameters, VARIANT SviParametersGivenOpt, VARIANT* pResult)
{
	begin_function
	::Create("SvlVolData", pResult, 3, Strikes, Parameters, SviParametersGivenOpt);
	end_function
}
HRESULT __stdcall MLCreateSvlFitVolData(VARIANT Strikes, VARIANT Bid, VARIANT Mid, VARIANT Ask, VARIANT Times, VARIANT InitialGuess, VARIANT FitIdeals, VARIANT FitMidPoints, VARIANT FitSpread, VARIANT* pResult)
{
	begin_function
	::Create("SVLFitVolData", pResult, 9, Strikes, Bid, Mid, Ask, Times, InitialGuess, FitIdeals, FitMidPoints, FitSpread);
	end_function
}
HRESULT __stdcall MLCreateVolData(VARIANT Strikes, VARIANT Interpolators, VARIANT* pResult)
{
	begin_function
	::Create("VolData", pResult, 2, Strikes, Interpolators);
	end_function
}
HRESULT __stdcall MLCreateMultiplierVolData(VARIANT Strikes, VARIANT Interpolators, VARIANT ReferenceVol, VARIANT* pResult)
{
	begin_function
	::Create("MultiplierVolData", pResult, 3, Strikes, Interpolators, ReferenceVol);
	end_function
}
HRESULT __stdcall MLCreateVolatilityStructure(VARIANT Name, VARIANT DataSourceAndOrLocation, VARIANT Spot, VARIANT Forwards, VARIANT Discounts, VARIANT Times, VARIANT DateHandle, VARIANT InterpAcrossTime, VARIANT BidVolData, VARIANT MidVolData, VARIANT AskVolData, VARIANT InterpolateType, VARIANT VolatilityDataType, VARIANT AsymptoticVolatilityStructureOpt, VARIANT DecayFactorOpt, VARIANT* pResult)
{
	begin_function
	::Create("VolatilityStructure", pResult, 15, Name, DataSourceAndOrLocation, Spot, Forwards, Discounts, Times, DateHandle, InterpAcrossTime, BidVolData, MidVolData, AskVolData, InterpolateType, VolatilityDataType, AsymptoticVolatilityStructureOpt, DecayFactorOpt);
	end_function
}
HRESULT __stdcall MLCreateZeroCurve(VARIANT GdaHandle, VARIANT DataSource, VARIANT* pResult)
{
	begin_function
	::Create("ZeroCurve", pResult, 2, GdaHandle, DataSource);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLCreateBasket
//
//	provides a way of creating an asset basket directly
//
HRESULT __stdcall MLCreateBasket(VARIANT Identifier, VARIANT DataSourceOpt, VARIANT DateOpt, VARIANT AssetsArr, VARIANT CompositeCurrenciesArr, VARIANT WeightsArr, VARIANT PayCurrency, VARIANT* pResult)
{	
	begin_function
	CComPtr<IAsset> spAsset;
	CAsset::CreateBasket(Identifier, DataSourceOpt, DateOpt, AssetsArr, CompositeCurrenciesArr, WeightsArr, PayCurrency, spAsset);
	CExcelInterface::Flatten(spAsset, pResult);	
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLCreateDateScheduleEx
//
//	Returns a handle to an array of dates, or the dates themselves.
//
HRESULT __stdcall MLCreateDateScheduleEx(VARIANT Calendar, VARIANT StartDate, VARIANT EndDate, VARIANT Frequency, VARIANT BusinessDayConvention, VARIANT ReturnHandleOpt, VARIANT* pResult)
{									
	begin_function	
	map_parameter(Calendar, estring, szCalendar);
	map_parameter(StartDate, long, nStartDate);
	map_parameter(EndDate, long, nEndDate);
	map_parameter(Frequency, estring, szFrequency);
	map_enum_parameter(BusinessDayConvention, BusinessDayConventionEnum, bdc);
	map_optional_parameter(ReturnHandleOpt, bool, bReturnHandleOpt, false);
			
	MlEqDateScheduleHandle hDateSchedule(new MlEqDateSchedule);
	hDateSchedule->Create(nStartDate, nEndDate, szFrequency, szCalendar, bdc);
	
	if (bReturnHandleOpt){
		map_analytic_to_com(hDateSchedule, DateSchedule, spDateSchedule);
		CExcelInterface::Flatten(spDateSchedule, pResult);
	} else {		
		std::vector<long> v;
		hDateSchedule->GetDates(v);
		CExcelInterface::Flatten(v, pResult);		
	}
	end_function
}


///////////////////////////////////////////////////////////////////////////////
//	MLDistinct
//
//	returns an array which represents the distinct elements in an input array
//
HRESULT __stdcall MLDistinct(VARIANT InputArray, VARIANT* pResult)
{
	CParameterMap						pmDistinct;
	
	begin_function
	map_parameter(InputArray, CParameterMap, pm);
	pm.GetDistinct(&pmDistinct);
	return pmDistinct.GetValue(pResult);
	end_function
}


///////////////////////////////////////////////////////////////////////////////
//	MLEDate
//
//	Adds a date interval of the form 3d, 1y etc. to an input date and returns
//  the result.
//
HRESULT __stdcall MLEDate(VARIANT StartDate, VARIANT Interval, VARIANT* pResult)
{			
	begin_function	
	map_parameter(StartDate, long, nDate)
	map_parameter(Interval, estring, szInterval)
	szInterval.StripWhiteSpace();	
	MlEqDateHandle hDate = MlEqDate(nDate).AddTenor(szInterval, NoBusinessDayConvention, "", false);
	return CComVariant(hDate->GetDate()).Detach(pResult);
	end_function
}


///////////////////////////////////////////////////////////////////////////////
//	MLEDateInv
//
//	Performs the reverse operation of MLEDate
//
HRESULT __stdcall MLEDateInv(VARIANT StartDate, VARIANT EndDate, VARIANT* pResult)
{
	begin_function
	map_parameter(StartDate, long, nStartDate)
	map_parameter(EndDate, long, nEndDate)
	estring	szTenor = MlEqDate(nStartDate).GetTenor(nEndDate, NoBusinessDayConvention, "", false);
	return szTenor.GetValue(pResult);	
	end_function
}


///////////////////////////////////////////////////////////////////////////////
//	MLEDateInvEx
//
//	Performs the reverse operation of MLEDate on an array of end dates. If at
//	least one of the resulting tenors is of the form [n]d for large n then this
//	function returns back the array of end dates.
//
HRESULT __stdcall MLEDateInvEx(VARIANT ReferenceDate, VARIANT EndDatesArr, VARIANT* pResult)
{
	begin_function
	map_parameter(ReferenceDate, long, nReferenceDate)
	map_parameter(EndDatesArr, std::vector<long>, anEndDates)
	CParameterMap						pm;
	
	std::vector<estring> aszTenors(anEndDates.size());


	for (long n = 0; n < anEndDates.size(); n++){
		aszTenors[n] = MlEqDate(nReferenceDate).GetTenor(anEndDates[n], NoBusinessDayConvention, "", false);
		long nTenor;
		if (aszTenors[n].right(1) == "d" && aszTenors[n].left("d").islong(&nTenor) && nTenor > 6){
			// Return back the array of end dates.
			pm.SetValue(anEndDates);
			pm.GetColumn(0, VT_I4, pResult);
			return S_OK;
		}
	}


	pm.SetValue(aszTenors);
	pm.GetColumn(0, VT_BSTR, pResult);

	
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLEvaluate
//
//	Excel worksheet function to evaluate a product, position or deal
//
HRESULT __stdcall MLEvaluate(VARIANT Handle, VARIANT CalculateOpt, VARIANT ReturnHandleOpt, VARIANT* pResult)
{	
	begin_function
	CComVariant							vResult(VT_DISPATCH);
	CComPtr<IDispatch>					spObject;						// object associated with Handle	
	CParameterMap						pmReturnHandleOpt;
	bool								bReturnHandleOpt = false;


	if (pmReturnHandleOpt.SetValue(ReturnHandleOpt)) CParameterMap::ThrowComErrorRR(IDS_INVALID_PARAMETER, IDS_HEADER_RETURN_HANDLE_OPT);
	if (!pmReturnHandleOpt.IsBlank()){		
		if (pmReturnHandleOpt.GetValue(&bReturnHandleOpt)) CParameterMap::ThrowComErrorRR(IDS_INVALID_PARAMETER, IDS_HEADER_RETURN_HANDLE_OPT);
	}
	// assume the Handle parameter is either a product, position or deal object and test in that order
	if (g_pApplication->GetObjectManager().GetObject(Handle, IID_IProduct, spObject)){
		if (g_pApplication->GetObjectManager().GetObject(Handle, IID_IPosition, spObject)){
			if (g_pApplication->GetObjectManager().GetObject(Handle, IID_IDeal, spObject)){
				CParameterMap::ThrowComErrorRV(IDS_NO_OBJECT_WITH_HANDLE, Handle);
			}
		}
	}
	if (CComDispatchDriverEx(spObject).Invoke1(L"Evaluate", &CalculateOpt, &vResult)) propagate_error;
	if (bReturnHandleOpt){
		CExcelInterface::Flatten(vResult, pResult);
	} else {				
		if (CComQIPtr<IResult>(vResult.pdispVal)->get_Value(pResult)) propagate_error;		
	}
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLEvaluateEx
//
//	Excel worksheet function to evaluate a set of products, positions or deals.
//	Subject to input market data perturbations.
//
//	A result handle is returned for each item.
//
HRESULT __stdcall MLEvaluateEx(VARIANT HandlesArr, VARIANT AssetsArr, VARIANT SpotsArr, VARIANT VolShiftsArr, VARIANT RateShiftsArr, VARIANT CalculateOpt, VARIANT* pResult)
{
	begin_function
	throw "MLEvaluateEx is not supported in this release";
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLGetCollectionIndex
//
//	Returns the 'Name' component of the m_coll variable in the property
//	'CollectionName' of the object associated with 'Collection' corresponding
//	to an input IndexHandle
//
HRESULT __stdcall MLGetCollectionIndex(VARIANT CollectionHandle, VARIANT CollectionName, VARIANT IndexHandle, VARIANT* pResult)
{
	begin_function
	CComPtr<IDispatch>					spObjectContainingCollection;
	CComPtr<IDispatch>					spIndexObject;
	CComVariant							vCollectionObject;
	estring								szIndex;
	
	map_parameter(CollectionName, estring, szCollectionName);
	if (g_pApplication->GetObjectManager().GetObject(CollectionHandle, spObjectContainingCollection)) propagate_error;
	if (g_pApplication->GetObjectManager().GetObject(IndexHandle, spIndexObject)) propagate_error;
	if (g_pApplication->GetObjectManager().GetCollectionIndex(spObjectContainingCollection, spIndexObject, szCollectionName, &szIndex)) propagate_error;
	return szIndex.GetValue(pResult);
	end_function
}


///////////////////////////////////////////////////////////////////////////////
//	MLGetCorrelation
//
//	Returns the correlation number defined between two assets. We make 
//  database calls if necessary,
//
HRESULT __stdcall MLGetCorrelation(VARIANT Asset1, VARIANT Asset2, VARIANT DataSourceOpt, VARIANT DateOpt, VARIANT* pResult)
{
	begin_function
	map_parameter(Asset1, estring, szAsset1);
	map_parameter(Asset2, estring, szAsset2);
	map_optional_parameter(DataSourceOpt, estring, szDataSource, "");
	map_optional_parameter(DateOpt, long, nDate, CComObjectCollectionSerialisableDefaulter::GetDate());
	double f = g_pApplication->GetCorrelation(szAsset1, szAsset2, szDataSource, nDate);
	return CComVariant(f).Detach(pResult);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLGetErrors
//
//	Returns the list of errors currently accumulated in the error buffer.
//
HRESULT __stdcall MLGetErrors(VARIANT* pResult)
{
	begin_function
	if (!_Module.GetSiriusApplication()->get_Errors(pResult)) return S_OK;
	throw CStringResource(IDS_NO_DATA_FOUND);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLGetForward
//
//	Returns the forward for a given spot value.
//
HRESULT __stdcall MLGetForward(VARIANT Spot, VARIANT DomesticZeroCurve, VARIANT Maturity, VARIANT ReferenceDateOpt, VARIANT DividendScheduleOpt, VARIANT ForeignZeroCurveOpt, VARIANT* pResult)
{
	begin_function
	map_parameter(Spot, double, fSpot);
	map_object_parameter(DomesticZeroCurve, ZeroCurve, hZeroCurve);
	map_parameter(Maturity, long, nMaturity);
	map_optional_parameter(ReferenceDateOpt, long, nToday, hZeroCurve->GetReferenceDate());
	map_optional_object_parameter(DividendScheduleOpt, DividendSchedule, hDividendSchedule);	
	map_optional_object_parameter(ForeignZeroCurveOpt, ZeroCurve, hForeignZeroCurve);
	return CComVariant(MlEqDividendSchedule::GetForward(nToday, nMaturity, fSpot, hZeroCurve, hDividendSchedule, hForeignZeroCurve)).Detach(pResult);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLGetForwardFromAsset
//
//	Returns the forward for an asset.
//
HRESULT __stdcall MLGetForwardFromAsset(VARIANT Asset, VARIANT Maturity, VARIANT ReferenceDateOpt, VARIANT ReturnNaturalForwardOpt, VARIANT ReinvestDividendsOpt, VARIANT* pResult)
{
	begin_function	
	map_object_parameter(Asset, Asset, hAsset);
	map_parameter(Maturity, long, nMaturity);
	map_optional_parameter(ReferenceDateOpt, long, nToday, hAsset->GetDateHandle()->GetDate());
	map_optional_parameter(ReturnNaturalForwardOpt, bool, bReturnNaturalForward, false);
	map_optional_parameter(ReinvestDividendsOpt, bool, bReinvestDivs, false);

	if (bReturnNaturalForward){
		return CComVariant(hAsset->GetNaturalForward(nToday, nMaturity, bReinvestDivs)).Detach(pResult);
	} else {
		// default
		return CComVariant(hAsset->GetForward(nToday, nMaturity, bReinvestDivs)).Detach(pResult);
	}
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLGetGDAHandle
//
//	Returns the GDA-style handle, if any, associated with a sirius object
//
HRESULT __stdcall MLGetGDAHandle(VARIANT Handle, VARIANT* pResult)
{
	begin_function
	CComPtr<IDispatch>					spObject;						// object associated with Handle
	std::string							szReturn;
	
	if (!g_pApplication->GetObjectManager().GetObject(Handle, IID_IZeroCurve, spObject)){		
		CZeroCurve* pzc = dynamic_cast<CZeroCurve*>(spObject.p);		
		szReturn = pzc->GetGDAHandle();
	} else {
		estring szHandle(Handle);
		if (szHandle.size()){
			szReturn = "The Sirius handle '" + szHandle + "' does not associate a GDA handle";
		} else {
			szReturn = CStringResource(IDS_NO_OBJECT_WITH_HANDLE);
		}
	}
	return estring::GetValue(szReturn, pResult);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLGetNaturalATMVolatilityFromAsset
//
//	Returns the natural ATM volatility from an asset.
//
HRESULT __stdcall MLGetNaturalATMVolatilityFromAsset(VARIANT Asset, VARIANT Today, VARIANT VolDate, VARIANT SpotOrForward, VARIANT* pResult)
{
	begin_function
	map_object_parameter(Asset, Asset, hAsset);
	map_parameter(VolDate, long, nDate);
	map_parameter(Today, long, nToday);
	map_enum_parameter(SpotOrForward, SpotOrForwardEnum, sf);	
	long nSpotOrFwd;
	switch (sf){
	case Spot:
		nSpotOrFwd = 0;
		break;
	case Forward:
		nSpotOrFwd = 1;
		break;
	default:
		throw "Parameter 'SpotOrForward' could not be set to value '" + CEnumMap::GetString("SpotOrForwardEnum", LIBID_Sirius, sf) + "'";
	}
	double f = hAsset->getNaturalATMVol(nDate, nSpotOrFwd);
	return CComVariant(f).Detach(pResult);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLGetNumBusinessDays
//
//	Returns the number of business days BETWEEN two dates.
//
HRESULT __stdcall MLGetNumBusinessDays(VARIANT StartDate, VARIANT EndDate, VARIANT CalendarOpt, VARIANT* pResult)
{
	begin_function	
	map_parameter(StartDate, long, nStartDate);
	map_parameter(EndDate, long, nEndDate);
	map_optional_parameter(CalendarOpt, estring, szCalendar, "");	
	long n = MlEqDate::GetNumBusinessDays(nStartDate, nEndDate, szCalendar);
	return CComVariant(n).Detach(pResult);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLGetObjectSignature
//
//	Returns the ID of the user and the time when an object was saved to the 
//  database.
//
//	Note that a version of this function that operated on an object handle
//	is conceptually invalid since we cannot assert that an object is
//  unmodified and originates from the database.
//
HRESULT __stdcall MLGetObjectSignature(VARIANT ObjectName, VARIANT Identifier, VARIANT DataSourceOpt, VARIANT DateOpt, VARIANT* pResult)
{
	begin_function	
	map_parameter(ObjectName, estring, szObjectName);
	map_optional_parameter(Identifier, estring, szIdentifier, "");
	map_optional_parameter(DataSourceOpt, estring, szDataSource, CComObjectCollectionSerialisableDefaulter::GetDataSourceStr());
	map_optional_parameter(DateOpt, long, nDate, 0);
	return _Module.GetObjectSignature(szObjectName, szIdentifier, szDataSource, nDate).Detach(pResult);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLGetProperties
//
//	get the available properties associated with a handle
//	
HRESULT __stdcall MLGetProperties(/*[in]*/ VARIANT HandleOrObjectName, /*[out, retval]*/ VARIANT* pResult)
{	
	CComPtr<IDispatch>					spObject;
	IID									iid = IID_NULL;
	HRESULT								hr;
	bool								bIsHandle = true;				// true if a handle is inputted as opposed to an object name
	
	begin_function
	// Get the object correaponding to HandleOrObjectName. We create an object if the user has
	// passed in an object name as opposed to a handle to an existing object.
	if (g_pApplication->GetObjectManager().GetObject(HandleOrObjectName, spObject)){
		bIsHandle = false;
		map_parameter(HandleOrObjectName, estring, szObjectName);
		if (g_pApplication->GetObjectManager().CreateObject(szObjectName, spObject)){
			propagate_error;
		}	
	}
	CParameterMap::GetObjectIID(spObject, &iid);
	if (iid == IID_IProduct && bIsHandle){
		CComQIPtr<IProduct> spProduct(spObject);
		if (!spProduct) throw CStringResource(IDS_UNHANDLED_EXCEPTION);
		if (!(hr = spProduct->GetRequiredParameters(pResult))) return S_OK;
	} else if (iid == IID_IResult && bIsHandle){
		CComQIPtr<IResult> spResult(spObject);
		if (!spResult) throw CStringResource(IDS_UNHANDLED_EXCEPTION);
		if (!(hr = (dynamic_cast<CResult*>(spResult.p))->GetParameters(pResult))) return S_OK;
	} else {
		if (!(hr = g_pApplication->GetProperties(spObject, pResult))) return S_OK;
	}
	ATLASSERT(hr);
	propagate_error;
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLGetSpot
//
//	Returns a spot value according to rules best understood by the
//	function comments.
// 
HRESULT __stdcall MLGetSpot(VARIANT Code, VARIANT CodeTypeOpt, VARIANT* pResult)
{
	begin_function
	map_parameter(Code, estring, szCode);
	map_optional_parameter(CodeTypeOpt, estring, szCodeType, "");

	// Attempt a DDE request from Reuters (if running)
	if (_Module.IsReutersRunning()){
		try {
			double fSpot = _Module.LoadBwlSpot(szCode, szCodeType, "Reuters", "Reuters", "");
			// Just before markets open, the spot is often set to zero.
			// In which case, we want to get the value from Sirius.
			if (fSpot) return CComVariant(fSpot).Detach(pResult);
		} catch (...){
			// do nothing
		}
	}

	// Try obtaining a value from the current native Sirius source
	try {
		CComPtr<ISpotSchedule> spSpotSchedule;	
		CSpotSchedule::Load(szCode, NoDataSource /*use default*/, 0/*use default*/, spSpotSchedule);
		if (spSpotSchedule){
			double fSpot;
			spSpotSchedule->get_LastValue(&fSpot);
			return CComVariant(fSpot).Detach(pResult);
		}
	} catch (...){
		// do nothing
	}

	// Try obtaining a value from Bwl (may have already tried this in the previous step depending on Sirius settings!)
	double fSpot = _Module.LoadBwlSpot(szCode, szCodeType, "", "" /*use default*/, _Module.GetBwlInterpolationRule());
	return CComVariant(fSpot).Detach(pResult);

	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MlGetSpotFromAsset
//
//	Returns a value from the spot schedule on a certain date.
//
HRESULT __stdcall MLGetSpotFromAsset(VARIANT Asset,VARIANT VolDate, VARIANT* pResult)
{
	begin_function	
	map_object_parameter(Asset, Asset, hAsset);
	map_parameter(VolDate, long, nDate);
	double f = hAsset->GetSpot(nDate);
	return CComVariant(f).Detach(pResult);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLGetSystemDate
//
//	Returns the current (GDA) system date
//
HRESULT __stdcall MLGetSystemDate(VARIANT* pResult)
{
	begin_function
	return CComVariant(MlEqDate::GetCurrentDate()).Detach(pResult);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLGetValue
//
//	Abstracts a value from a handle
//
HRESULT __stdcall MLGetValue(VARIANT HandleOrIdentifier, VARIANT PropertyOpt, VARIANT ObjectNameOpt, VARIANT* pResult)
{
	begin_function
	map_parameter(HandleOrIdentifier, estring, szHandle);
	map_optional_parameter(PropertyOpt, estring, szProperty, "");	
	map_optional_parameter(ObjectNameOpt, estring, szObjectName, "");
		
	CComPtr<IDispatch>					spObject;	
	IID									iidObject;	
		
	if (szObjectName.size()){		
		if (g_pApplication->GetObjectManager().GetObject(szHandle, szObjectName, spObject)) propagate_error;
	} else {
		if (g_pApplication->GetObjectManager().GetObject(szHandle.GetValue(), spObject)) propagate_error;
	}	
	if (CParameterMap::GetObjectIID(spObject, &iidObject)) ATLASSERT(false);
	if (iidObject == IID_IResult){
		// Result handle special case.
		CComQIPtr<IResult>				spResult(spObject);
		CComVariant v;
		if (!szProperty.size()) szProperty.LoadString(IDS_HEADER_PRICE);
		if (!szProperty.CompareNoCaseAndSpace("all")){
			if (spResult->get_Value(&v)) propagate_error;			
		} else {
			if (spResult->GetValue(szProperty.GetBSTR(), &v)) propagate_error;
		}
		CExcelInterface::Flatten(v, pResult);		
	} else if (iidObject == IID_IProduct && szProperty.size()){
		// Product special case.
		CComQIPtr<IProduct> spProduct(spObject);
		if (!spProduct) throw CStringResource(IDS_UNHANDLED_EXCEPTION);
		if ((dynamic_cast<CProduct*>(spProduct.p))->GetValue(szProperty, pResult)) propagate_error;
	} else {		
		// General case		
		CComDispatchDriverEx ddObject(spObject);
		
		if (!szProperty.size() || !szProperty.CompareNoCase("value")){
			// The user wants the entire object back.
			CExcelInterface::Flatten(ddObject, pResult);
			return S_OK;
		}
									
		// Try the case where szProperty is the name of an exposed property of spObject.
		CComVariant v;
 		if (!g_pApplication->GetObjectManager().ImplementGetValue(spObject, szProperty, &v)){			
			CExcelInterface::Flatten(v, pResult);
			return S_OK;
		}
		// Try the case where szProperty is a parameter of an MLCreate[...] function corresponding to spObject.			
		if (!ddObject.GetPropertyByName(L"Value", &v) && !g_pApplication->GetOneParameter(ddObject, szProperty, &v)){
			CExcelInterface::Flatten(v, pResult);
			return S_OK;
		}
		// If Handle is one of the schedule types then szProperty could be a date
		if (iidObject == IID_ISpotSchedule || iidObject == IID_IDateSchedule || iidObject == IID_IDividendSchedule){
			CParameterMap pmName(szProperty);
			if (pmName.IsDate()){
				if (ddObject.GetPropertyByName(L"ValueAt", &PropertyOpt, 1, pResult)) propagate_error;
				return S_OK;
			}
		}
		// No property found.
		CParameterMap::ThrowComErrorRS(IDS_PROPERTY_NOT_FOUND, szProperty);
	}
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLGetVolatilityFromAsset
//
//	Returns a volatility value associated with an asset.
//
HRESULT __stdcall MLGetVolatilityFromAsset(VARIANT Asset, VARIANT StrikeHandle, VARIANT VolDate, VARIANT BidOfferOpt, VARIANT ReferenceDateOpt, VARIANT ReturnNaturalVolatilityOpt, VARIANT SkewInvariantStrikeOpt, VARIANT* pResult)
{
	begin_function	
	long								nStrikes;	
	long								nSkewInvStrikes = 0L;
	
	map_object_parameter(Asset, Asset, hAsset);
	map_object_parameter(StrikeHandle, Strikes, hStrike);
	map_parameter(VolDate, long, nDate);
	map_optional_enum_parameter(BidOfferOpt, BidAskMidEnum, nBidOffer, Middle);
	map_optional_parameter(ReferenceDateOpt, long, nToday, hAsset->GetDateHandle()->GetDate());	
	map_optional_parameter(ReturnNaturalVolatilityOpt, bool, bReturnNaturalVolatility, false);
	map_optional_object_parameter(SkewInvariantStrikeOpt, Strikes, hSkewInvStrike);
		
	if (!!hSkewInvStrike){
		nSkewInvStrikes = hSkewInvStrike->m_hStrikes.size();
		if (nSkewInvStrikes > 1){
			throw "The parameter 'SkewInvariantStrikeOpt' cannot contain more than one strike";
		}
	}
	
	nStrikes = hStrike->m_hStrikes.size();

	CVector result(nStrikes);
	for (int n = 0; n < nStrikes; n++){
		if (bReturnNaturalVolatility){
			if (nSkewInvStrikes == 0){
				result[n] = hAsset->GetNaturalVolatility(*hStrike->m_hStrikes[n], nToday, (double)nDate, nBidOffer);
			} else {
				result[n] = hAsset->GetNaturalVolatility(*hStrike->m_hStrikes[n], nToday, (double)nDate, nBidOffer, hSkewInvStrike->m_hStrikes[0]);
			}
		} else {
			// default
			if (nSkewInvStrikes == 0){
				result[n] = hAsset->GetCompositeVolatility(*hStrike->m_hStrikes[n], nToday, (double)nDate, nBidOffer);
			} else{
				result[n] = hAsset->GetCompositeVolatility(*hStrike->m_hStrikes[n], nToday, (double)nDate, nBidOffer, hSkewInvStrike->m_hStrikes[0]);
			}
		}
	}

	// Return the result
	if (nStrikes == 1){		
		return CComVariant(result[0]).Detach(pResult);
	} else {
		unmap_parameter(result, pm);		
		return pm.GetValue(pResult);
	}
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLPutValue
//
//	Inserts a value from a handle
//
HRESULT __stdcall MLPutValue(VARIANT HandleOrIdentifier, VARIANT PropertyOpt, VARIANT ObjectNameOpt, VARIANT Value, VARIANT* pResult)
{
	begin_function
	map_parameter(HandleOrIdentifier, estring, szHandle);
	map_optional_parameter(PropertyOpt, estring, szProperty, "");	
	map_optional_parameter(ObjectNameOpt, estring, szObjectName, "");
		
	CComPtr<IDispatch>					spObject;	
	IID									iidObject;
	CComDispatchDriverEx				ddObject;

	if (szObjectName.size()){		
		if (g_pApplication->GetObjectManager().GetObject(szHandle, szObjectName, spObject)) propagate_error;
	} else {
		if (g_pApplication->GetObjectManager().GetObject(szHandle.GetValue(), spObject)) propagate_error;
	}
	ddObject = spObject;
	if (CParameterMap::GetObjectIID(spObject, &iidObject)) ATLASSERT(false);
	if (iidObject == IID_IResult){
		// Result handle special case.
		CParameterMap::ThrowComErrorR(IDS_PUT_VALUE_RESULT);		
	} else if (iidObject == IID_IProduct && szProperty.size()){
		// Product special case.						
		CComQIPtr<IProduct>				spProduct(spObject);
		if (!spProduct) throw CStringResource(IDS_UNHANDLED_EXCEPTION);
		CProduct*						pProduct = dynamic_cast<CProduct*>(spProduct.p);
		if (!pProduct) throw CStringResource(IDS_UNHANDLED_EXCEPTION);
		if (pProduct->PutValue(szProperty, Value)) propagate_error;
		return estring("OK").GetValue(pResult);
	} else {		
		// General case - all objects other than product and result handles.				
		if (!szProperty.size()) szProperty.assign("Value");
		
		// Try the case where szProperty is the name of an exposed property of spObject
 		if (!g_pApplication->GetObjectManager().ImplementPutValue(spObject, szProperty, Value)){
			return estring("OK").GetValue(pResult);
		}
		
		// If Handle is one of the schedule types then szProperty could be a date
		if (iidObject == IID_ISpotSchedule || iidObject == IID_IDateSchedule || iidObject == IID_IDividendSchedule){						
			map_parameter(Value, double, f);
			map_parameter(PropertyOpt, DATE, date);												
			CComVariant params[] = {f, date};				
			if (ddObject.PutPropertyByName(L"ValueAt", params, 2)) propagate_error;
			return estring("OK").GetValue(pResult);			
		} else {
			// No property found if this point has been reached.
			propagate_error;
		}
	}
	end_function
}


///////////////////////////////////////////////////////////////////////////////
//	MLGetVbaKey
//
//	returns the key for a locked VBA project
//
HRESULT __stdcall MLGetVbaKey(VARIANT FileName, VARIANT* pResult)
{
	CParameterMap		pmFileName;
	std::string			szFileName;
	estring				szKey;

	begin_function
	if (pmFileName.SetValue(FileName) || pmFileName.GetString(&szFileName)) throw CStringResource(IDS_INVALID_FILE);
	CExcelInterface::VbaKey(szFileName, &szKey);
	return szKey.GetValue(pResult);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLGetWorkbookName
//
//	returns the full name (including path) of the current worksheet
//
HRESULT __stdcall MLGetWorkbookName(VARIANT* pResult)
{
	estring								szName;
	
	begin_function
	if (CExcelInterface::GetWorkbookName(&szName)) CParameterMap::ReturnErrorR(IDS_UNHANDLED_EXCEPTION);
	return szName.GetValue(pResult);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLGetXML
//
//	returns the XML associated with an input object
//
HRESULT __stdcall MLGetXML(VARIANT Handle, VARIANT* pResult)
{	
	CComPtr<IDispatch>					spObject;
	xmlstreamer							ssOut;
	std::vector<std::string>			vector;
	CParameterMap						pmOut;
	
	begin_function
	if (g_pApplication->GetObjectManager().GetObject(Handle, spObject)) propagate_error;
	if (CXmlStreamer::GetXML(CComVariant(spObject.p), ssOut)) propagate_error;
	
	// Replace tabs with four spaces since Excel does not display tabs at the start of strings in spreadsheet cells.
	estring								szXml((char*)ssOut);
	szXml.ReplaceStrInStr("\t", "    ");
	szXml.Split("\r\n", &vector);
	if (pmOut.SetValue(vector)) propagate_error;
	return pmOut.GetValue(pResult);	
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLIsBusinessDate
//
//	Returns True if the input date is a business day (i.e. not a weekend or
//  a public holiday.
//
HRESULT __stdcall MLIsBusinessDate(VARIANT QueryDate, VARIANT CalendarOpt, VARIANT* pResult)
{
	begin_function
	map_parameter(QueryDate, long, nDate);
	map_optional_parameter(CalendarOpt, estring, szCalendar, "");

	if (MlEqDate::IsBusinessDay(nDate, szCalendar)){
		return CComVariant("True").Detach(pResult);
	} else {
		return CComVariant("False").Detach(pResult);
	}
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLLoad
//
//	Loads an object from the database. We return a handle to it.
//
HRESULT __stdcall MLLoad(VARIANT ObjectName, VARIANT Identifier, VARIANT DataSourceOpt, VARIANT DateOpt, VARIANT* pResult)
{	
	begin_function	
	CComPtr<IDispatch>					spObject;
	CParameterMap						pm;
	CComBSTR							sHandle;
	
	map_parameter(ObjectName, estring, szObjectName);
	map_optional_parameter(Identifier, estring, szIdentifier, "");
	map_optional_parameter(DataSourceOpt, estring, szDataSource, "");
	map_optional_parameter(DateOpt, long, nDate, 0);
	if (_Module.Load(szObjectName, szIdentifier, szDataSource, nDate, true, spObject)) propagate_error;
	CExcelInterface::Flatten(spObject, pResult);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLPutXML
//
//	Makes an object out of the input Xml and returns a handle to it.
//
HRESULT __stdcall  MLPutXML(VARIANT Xml, VARIANT* pResult)
{		
	begin_function
	CComPtr<IDispatch> spObject;
		
	map_parameter(Xml, std::stringstream, ssXml);
	if (_Module.GetSiriusApplication()->CreateObjectFromXml((CComBSTR)estring(ssXml), &spObject)) propagate_error;
	CExcelInterface::Flatten(spObject, pResult);
	end_function
}


///////////////////////////////////////////////////////////////////////////////
//	MLReinitializeStrikesFromAsset
//
//	Reinitializes strike data with asset
//
HRESULT __stdcall MLReinitializeStrikesFromAsset(VARIANT Strikes, VARIANT Asset, VARIANT* pResult)
{
	std::vector<MlEqStrikeHandle>			strikes;

	begin_function

	map_object_parameter(Strikes, Strikes, hStrike);
	map_object_parameter(Asset, Asset, hAsset);

	strikes = hStrike->m_hStrikes;
	int nstrikes = strikes.size();
	for (long n = 0; n < nstrikes; n++){
		strikes[n]->reinitializeStrike(*hAsset);
	}
	return estring("OK").GetValue(pResult);
	end_function

}


/////////////////////////////////////////////////////////////////////////////
//	MLResize
//
//	Emulates the GdaResize() function
//
HRESULT __stdcall MLResize(VARIANT Value, VARIANT* pResult)
{
	begin_function
	return CExcelInterface::Resize(Value, pResult);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLSave
//
//	saves an object to the database
//
HRESULT __stdcall MLSave(VARIANT Handle, VARIANT* pResult)
{	
	begin_function
	map_parameter(Handle, estring, szHandle);

 	if (g_pApplication->Save((CComBSTR)szHandle)) propagate_error;
	return CComVariant(L"OK").Detach(pResult);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	MLSetToBusinessDay
//
//	Returns the closest business day associated with an input date
//
HRESULT __stdcall MLSetToBusinessDay(VARIANT RefDate, VARIANT BusinessDayConventionOpt, VARIANT CalendarOpt, VARIANT* pResult)
{
	begin_function

	map_parameter(RefDate, long, nRefDate);
	map_optional_enum_parameter(BusinessDayConventionOpt, BusinessDayConventionEnum, bdc, NoChange)
	map_optional_parameter(CalendarOpt, estring, szCalendar, "");
	
	MlEqDateHandle hDate = MlEqDate(nRefDate).SetToBusinessDay(bdc, szCalendar);
	return CComVariant(hDate->GetDate()).Detach(pResult);
	end_function
}


///////////////////////////////////////////////////////////////////////////////
//	MLYearFraction
//
//	Calculates the time in years between two dates subject to the day counting
//  conventions of the first date.
//
HRESULT __stdcall MLYearFraction(VARIANT DateHandle, VARIANT EndDate, VARIANT* pResult)
{			
	begin_function	
	map_object_parameter(DateHandle, Date, hStartDate);
	map_parameter(EndDate, long, nEndDate);		
	double f = hStartDate->GetYearFraction(nEndDate);
	return CComVariant(f).Detach(pResult);
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	Retrieve
//
//	Loads a market data object from the spreadsheet repository.
//
HRESULT __stdcall Retrieve(BSTR ObjectName)
{	
	estring szMessage("Retrieving objects of type '" + estring(ObjectName) + "' is not supported in this release");
	CParameterMap::DisplayError(szMessage, MB_ICONINFORMATION);
	return S_OK;
}


/////////////////////////////////////////////////////////////////////////////
//	ReutersRunning
//
//	returns TRUE if the excel function wizard is visible
//
HRESULT __stdcall ReutersRunning(BOOL* pVal)
{
	*pVal = _Module.IsReutersRunning() ? TRUE : FALSE;
	return S_OK;
}


//////////////////////////////////////////////////////////////////////////////
//	SaveVBAModule
//
//	impliement the Sirius->SaveVBAModule menu option
//
HRESULT __stdcall SaveVBAModule()
{
	CDlgModuleSave().DoModal();
	return S_OK;
}


/////////////////////////////////////////////////////////////////////////////
//	SiriusApplication
//
//	returns a handle to the application object
//
HRESULT __stdcall SiriusApplication(ISiriusApplication** pVal)
{
	return _Module.GetSiriusApplication().CopyTo(pVal);	
}
