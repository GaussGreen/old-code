//	MlEqAsset.h :			 Asset handle class
//
//	Author :				 David Cuin
/////////////////////////////////////////////////////////////////////////////

#ifndef _MLEQASSET_H_
#define _MLEQASSET_H_

#include "smart.h"
#include "cmatrix.h"
	
class MlEqVolatilityStructure;

class MlEqAsset : public RCObject
{
// forward retrieval functions
public:	
	double									GetForward(long nMaturity, bool bReinvestDivs) const;
	double									GetForward(long nStart, long nMaturity, bool bReinvestDivs) const;
	double									GetNaturalForward(long nMaturity, bool bReinvestDivs) const;
	double									GetNaturalForward(long nStart, long nMaturity, bool bReinvestDivs) const;	
	double									GetQuantoForward(long nMaturity, bool bReinvestDivs) const;
	double									GetQuantoForward(long nStart, long nMaturity, bool bReinvestDivs) const;

	double									GetBumpedForward(long nStart, long nMaturity, double bumpedSpot, bool bReinvestDivs) const;

protected:	
	double									GetNaturalForward(long nStart, long nMaturity, double forwardnStart, bool bReinvestDivs) const;

// other functions


public:

	MlEqAsset();		
	MlEqConstDateHandle						GetDateHandle(void) const;		// Use this in preference to GetInternalDate() since the former takes the handle from the volatility structure.
	long									GetInternalDate(void) const;	// NEVER USE THIS IN QUANT CODE. ALWAYS USE GetDateHandle(). See David Cuin if you intend to violate this.
	void									PutInternalDate(long nDate);

	CVector									Calibrate(DATE maturityDate, MlEqDate& dateToDouble, std::vector<std::vector<MlEqStrikeHandle> >& pVolStrikes);
	bool									IsBasket(void) const {return m_bIsBasket;}
	bool									IsCalibrated(void) const;
	const std::vector<MlEqAssetHandle>&		GetAssets(void) const;	
	virtual double							GetCompositeVolatility(const MlEqStrike& strike,const DATE& date, BidAskMidEnum bidOrAsk = Middle) const;
	virtual double							GetCompositeVolatility(const MlEqStrike& strike,long nStartDate,const DATE& date, BidAskMidEnum bidOrAsk = Middle, MlEqStrikeHandle skewInvariantStrikeHandle = NULL) const;
	virtual double							GetCompositeVolatility(long nStartDate,const DATE& date,double volS, BidAskMidEnum bidOrAsk = Middle) const;
	std::string								GetCountry(void) const;


	DataSourceEnum							GetDataSource(void) const;
	MlEqDividendScheduleHandle				GetDividendSchedule(void) const;	
	const MlEqAssetHandle					GetFxCompositeAsset() const;
	const MlEqAssetHandle					GetFxNaturalAsset() const;
	const MlEqAssetHandle					GetFxPayAsset() const;
	const std::string&						GetName(void) const;		
	double									GetNaturalVolatility(const MlEqStrike& strike, const DATE& date, BidAskMidEnum bidOrAsk = Middle) const;
	double									GetNaturalVolatility(const MlEqStrike& strike, long nStart,const DATE& date, BidAskMidEnum bidOrAsk = Middle, MlEqStrikeHandle skewInvariantStrikeHandle = NULL) const;
	const MlEqCorrelationMatrixHandle		GetCorrelationMatrix(void) const;
	const MlEqZeroCurveHandle				GetCompositeZeroCurve(bool bThrow) const;
	const MlEqZeroCurveHandle				GetPayZeroCurve(bool bThrow) const;
	const MlEqZeroCurveHandle				GetNaturalZeroCurve(bool bThrow) const;
	double									GetNaturalSpot(long nDate) const;
	void									GetNaturalSpots(long nToday, const std::vector<long>& vectorDates, MlEqSpotSchedule& SpotSchedule) const;
	void									GetNaturalSpotsBeforeToday(long nToday, const std::vector<long>& vectorDates, std::vector<double>* pSpots) const;
	std::string								GetRegion(void) const;
	std::string								GetSector(void) const;
	void									GetSpotsBeforeToday(long nToday, const GVector<long>& anDates, CMatrix* pm) const;
	void									GetSpotsBeforeToday(long nToday, const std::vector<long>& anDates, std::vector<double>* pSpots) const;
	void									GetSpotsBeforeToday(long nToday, const std::vector<long>& anDates, CMatrix* pm) const;
	double									GetSpot(long nDate) const;
	void									GetSpots(long nToday, const std::vector<long>& vectorDates, MlEqSpotSchedule& SpotSchedule) const;
	long									GetCurrentDate()const;

	virtual double							GetStockCompositeFxVariance(const MlEqStrike& strike, const DATE& date, BidAskMidEnum bidOrAsk = Middle) const;
	virtual double							GetStockCompositeFxVariance(const MlEqStrike& strike, long nStart, const DATE& date, BidAskMidEnum bidOrAsk = Middle) const;

	virtual double							GetStockNaturalFxVariance(const MlEqStrike& strike, const DATE& date, BidAskMidEnum bidOrAsk = Middle); 
	virtual double							GetStockNaturalFxVariance(const MlEqStrike& strike, long nStart, const DATE& date, BidAskMidEnum bidOrAsk = Middle);
	double									GetVolatility(const MlEqStrike& strike, const DATE& date, BidAskMidEnum bidOrAsk = Middle) const;
	double									GetVolatility(const MlEqStrike& strike, long nStart, const DATE& date, BidAskMidEnum bidOrAsk = Middle, MlEqStrikeHandle skewInvariantStrikeHandle = NULL) const;
	double									GetNaturalNaiveFutureVolatility(const MlEqStrike& shortStrike, const MlEqStrike& longStrike, long nStartDate, const DATE& date, BidAskMidEnum bidOrAsk = Middle) const;
	double									GetNaiveFutureVolatility(const MlEqStrike& shortStrike,const MlEqStrike& longStrike,long nStartDate,const DATE& date, BidAskMidEnum bidOrAsk = Middle) const;
	double									getNaturalATMVol(long nMaturity,int SpotOrForward=1, BidAskMidEnum bidOrAsk = Middle) const;
	double									getNaturalATMVol(long nStart, long nMaturity,int SpotOrForward, BidAskMidEnum bidOrAsk = Middle) const;

	double									GetStockPayFxVariance(const MlEqStrike& strike,const DATE& date, BidAskMidEnum bidOrAsk = Middle) const ;
	double									GetStockPayFxVariance(const MlEqStrike& strike,int nStart,const DATE& date, BidAskMidEnum bidOrAsk = Middle) const;

	MlEqVolatilityStructureHandle           GetVolatilityStructure() const;

	const CVector&							GetWeights(void) const;		

	bool									HasVolatilityStructure(void) const;
		
	void									PutDataSource(DataSourceEnum ds);		
	void									PutCorrelationMatrix(MlEqCorrelationMatrixHandle hCorrelation);	
	void									PutDividendSchedule(MlEqDividendScheduleHandle h)		{m_pdivs = h;}	
	void									PutVolatilityStructure(MlEqVolatilityStructureHandle pVols);	
	void									PutFxCompositeAsset(MlEqAssetHandle h)			{m_fxComposite = h;}
	void									PutFxNaturalAsset(MlEqAssetHandle h)			{m_fxNatural = h;}
	void									PutFxPayAsset(MlEqAssetHandle h)				{m_fxPay = h;}
	void									PutCompositeZeroCurve(MlEqZeroCurveHandle h)	{m_pzcComposite = h;}
	void									PutPayZeroCurve(MlEqZeroCurveHandle h)			{m_pzcPay = h;}
	void									PutName(const std::string& sz)					{m_szName = sz;}
	void									PutNaturalZeroCurve(MlEqZeroCurveHandle h)		{m_pzcNatural = h;}	
	void									PutSpotSchedule(MlEqSpotScheduleHandle h)		{m_pSpots = h;}

public:	
	std::vector<MlEqAssetHandle>			m_assets;
	CVector									m_assetWeights;
	bool									m_bIsBasket;

protected:	
	double									GetCompositeNaturalFxVariance(const MlEqStrike& strike,const DATE& date, BidAskMidEnum bidOrAsk = Middle)const ;
	double									GetCompositeNaturalFxVariance(const MlEqStrike& strike,long nStart,const DATE& date, BidAskMidEnum bidOrAsk = Middle) const;
	double									GetCompositePayFxVariance(const MlEqStrike& strike,const DATE& date, BidAskMidEnum bidOrAsk = Middle) const ;
	double									GetCompositePayFxVariance(const MlEqStrike& strike,long nStart,const DATE& date, BidAskMidEnum bidOrAsk = Middle) const;
	double									GetNaturalPayFxVariance(const MlEqStrike& strike,const DATE& date, BidAskMidEnum bidOrAsk = Middle) const;
	double									GetNaturalPayFxVariance(const MlEqStrike& strike,long nStart, const DATE& date, BidAskMidEnum bidOrAsk = Middle) const;
	
	DataSourceEnum							m_ds;
	long									m_nDate;					// This is kept closely (though not necessarily always) synchronised with the volatility structure date.
	bool									m_bCalibrated;				// True if this represents a calibrated asset.
	MlEqAssetHandle							m_fxNatural;
	MlEqAssetHandle							m_fxPay;
	MlEqAssetHandle							m_fxComposite;
	MlEqCorrelationMatrixHandle				m_hCorrelation;
	MlEqDividendScheduleHandle				m_pdivs;
	MlEqSpotScheduleHandle					m_pSpots;
	MlEqVolatilityStructureHandle			m_pVols;					// NEVER make this public - always use the Get and Set functions to access this pointer
	MlEqZeroCurveHandle						m_pzcComposite;
	MlEqZeroCurveHandle						m_pzcNatural;
	MlEqZeroCurveHandle						m_pzcPay;
	std::string								m_szName;	
};

#endif