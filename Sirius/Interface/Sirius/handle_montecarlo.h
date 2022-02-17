//	handle_montecarlo.h : Declaration of the CMonteCarlo
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __MONTECARLO_H_
#define __MONTECARLO_H_

#include "resource.h"
#include "montecarlo.h"

class ATL_NO_VTABLE CMonteCarlo : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CMonteCarlo, &CLSID_MonteCarlo>,
	public ISupportErrorInfo,
	public IDispatchImpl<IMonteCarlo, &IID_IMonteCarlo, &LIBID_Sirius>
{
public:		
	enum MonteCarloTypeExEnum	// enumerators detailing methods on how to create Monte Carlo handles
								// that are never passed into analytic classes.
								// Make sure that their values do not clash with
								// MonteCarloTypeEnum						
	{
		CalibrateForwardSkew		= -100,	
		CalibrateHybridModel,
	};
	
	
	static const std::string			s_szTimes;
	static const std::string			s_szDiscountFactors;
	static const std::string			s_szForwards;
	static const std::string			s_szVolatilities;
	static const std::string			s_szBetas;
	static const std::string			s_szVolOfVols;
	static const std::string			s_szMeanReversionRates;
	static const std::string			s_szMeanReversionLevels;
	static const std::string			s_szModelInfo;	
	static const std::string			s_sz_cL_Header;
	static const std::string			s_sz_cR_Header;
	static const std::string			s_sz_addTanhWings_Header;
	static const std::string			s_sz_yPower_Header;
	static const std::string			s_sz_seed_Header;
	static const std::string			s_sz_npaths_Header;
	static const std::string			s_sz_calibflag_Header;
	static const std::string			s_sz_numberVolStates_Header;
	static const std::string			s_sz_localVolFlag_Header;
	static const std::string			s_sz_saveVolVolInfo_Header;
	static const std::string			s_sz_numberGridPoints_Header;
	static const std::string			s_sz_saveVolGrid_Header;
	static const std::string			s_sz_randomNumberFlag_Header;
	static const std::string			s_sz_contolVariateFlag_Header;	
	static const std::string			s_sz_globalCalibFlag_Header;
	static const std::string			s_sz_UseHermite_Header;
	
	CMonteCarlo() : m_h(new MlEqMonteCarlo(NULL)){}
	DECLARE_REGISTRY_RESOURCEID(IDR_MONTECARLO)
	DECLARE_NOT_AGGREGATABLE(CMonteCarlo)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CMonteCarlo)
		COM_INTERFACE_ENTRY(IMonteCarlo)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CMonteCarlo, MlEqMonteCarlo);

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
	STDMETHOD(GetForwardValue)(/*[in]*/ long Date, /*[in]*/ long Asset, /*[out, retval]*/ double* pVal);
	STDMETHOD(Average)(/*[in]*/long Path, /*[in]*/  long DateFrom, /*[in]*/ long DateUpTo, /*[in]*/ long Asset,/*[out, retval]*/ double* pVal);
	STDMETHOD(Max)(/*[in]*/long Path, /*[in]*/  long DateFrom, /*[in]*/ long DateUpTo, /*[in]*/ long Asset,/*[out, retval]*/ double* pVal);
	STDMETHOD(BasketMax)(/*[in]*/long Path, /*[in]*/  long DateFrom, /*[in]*/ long DateUpTo, /*[in]*/ IArray* Weights,/*[out, retval]*/ double* pVal);
	STDMETHOD(Basket)(/*[in]*/long Path, /*[in]*/  long Date, /*[in]*/ IArray* Weights,/*[out, retval]*/ double* pVal);
	STDMETHOD(GetSpotValue)(/*[in]*/ long Path, /*[in]*/ long Date, /*[in]*/ long Asset, /*[out, retval]*/ double* pVal);
	STDMETHOD(GetDiscountFactor)(/*[in]*/ long Path, /*[in]*/ long Date,/*[out, retval]*/ double* pVal);
	STDMETHOD(GetZeroBond)(long Path, long Date, long AsOfDate, double* pVal);

	STDMETHOD(GetIndexFromDate)(/*[in]*/long Date,/*[out, retval]*/ long* pVal);
	STDMETHOD(GetDate)(/*[in]*/ long Date, /*[out, retval]*/ long* pVal);
	STDMETHOD(get_Dates)(/*[out, retval]*/ long *pVal);
	STDMETHOD(GetTime)(/*[in]*/ long Date, /*[out, retval]*/ double* pVal);
	STDMETHOD(get_Paths)(/*[out, retval]*/ long *pVal);
	STDMETHOD(get_MonteCarloType)(/*[out, retval]*/ MonteCarloTypeEnum *pVal);
	STDMETHOD(GetVolatility)(long Path, long Date, double strike, long Asset, double* pVal);
	STDMETHOD(get_Assets)(/*[out, retval]*/ long *pVal);
	STDMETHOD(GetImpliedLogContractVol)( long Path, long Date,long Asset,long ngauss,double lowerStdev,double upperStdev, double* pVal);
	STDMETHOD(PutSpotValue)(/*[in]*/ long Path, /*[in]*/ long Date, /*[in]*/ long Asset, /*[in]*/ double newVal);
	
protected:
	HRESULT								GetModelInfoParameter(CParameterMap& pm, const std::string& szHeading, long nRow, double fValue) const;
	HRESULT								put_ValueForwardSkew(std::vector<CParameterMap>& vpm, std::vector<CComVariant>&	vv);
	HRESULT								put_ValueGeneral(std::vector<CParameterMap>& vpm, std::vector<CComVariant>&	vv);
	HRESULT								put_ValueHybridForwardSkew(std::vector<CParameterMap>& vpm, std::vector<CComVariant>& vv);
	HRESULT								put_ValueQuasi(std::vector<CParameterMap>& vpm, std::vector<CComVariant>& vv);
	HRESULT								put_CalibrateMC(std::vector<CParameterMap>& vpm, std::vector<CComVariant>&	vv);
	HRESULT								put_CalibrateHybridMC(std::vector<CParameterMap>& vpm, std::vector<CComVariant>& vv);
	HRESULT								put_ValueLocalVolMC(std::vector<CParameterMap>& vpm, std::vector<CComVariant>& vv);

	HRESULT								SetModelInfoParameter(const CParameterMap& pm, const std::string& szHeading, long nRow, CParameterMap* ppm) const;
	HRESULT								SetModelInfoParameter(const CParameterMap& pmIn, CVector* pv) const;
	
	std::vector<CParameterMap>			m_vpm;
};

#endif