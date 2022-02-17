`//	MlEqAsset.cpp :             Implementation of the Asset class
//
//	Author :				    David Cuin
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MlEqAsset.h"
#include "MlEqObjects.h"
#include "threemoments.h"


/////////////////////////////////////////////////////////////////////////////
//	Construction
//
MlEqAsset::MlEqAsset()
{
	m_ds = NoDataSource;
	m_nDate = 0;	
	m_bCalibrated = false;	
	m_pdivs = NULL;
	m_pSpots = NULL;	
	m_pVols = NULL;
	m_pzcComposite = NULL;
	m_pzcNatural = NULL;
	m_pzcPay = NULL;
	m_fxNatural = NULL;
	m_fxPay = NULL;
	m_fxComposite = NULL;
	m_hCorrelation = NULL;
	m_bIsBasket = false;
}


/////////////////////////////////////////////////////////////////////////////
//	GetForward / GetQuantoForward / GetNaturalForward
//
//	These functions implement the retrieval of forwards from the asset:
//	Note that they call call the lowest level forward function
//	MlEqDividendSchedule::GetForward.
//
//	GetForward - these functions call GetQuantoForward (they are provided for convenience)
//
double MlEqAsset::GetForward(long nMaturity, bool bReinvestDivs) const
{
	double f = GetQuantoForward(nMaturity, bReinvestDivs);
	return f;	
}

double MlEqAsset::GetForward(long nStart, long nMaturity, bool bReinvestDivs) const
{
	double f = GetQuantoForward(nStart, nMaturity, bReinvestDivs);
	return f;
}
//
//	GetNaturalForward
//
double MlEqAsset::GetNaturalForward(long nStart, long nMaturity, double forwardnStart, bool bReinvestDivs) const
//  All overloads must call into this one.
{					
	if (!IsBasket()){
		MlEqDividendScheduleHandle	hDivs = NULL;
		MlEqZeroCurveHandle			hzcYield = NULL;
		if (!bReinvestDivs && !!(hDivs = m_pdivs)){
			hzcYield = hDivs->GetForeignCurve();
		}
		double f = MlEqDividendSchedule::GetForward(nStart, nMaturity, forwardnStart, m_pzcNatural, hDivs, hzcYield);
		return f;
	} else {		
		double fTotalForward = 0.0;			
		for (int n = 0; n < m_assets.size(); n++){
			double f = m_assets[n]->GetNaturalForward(nStart, nMaturity, forwardnStart, bReinvestDivs);	// recursion here is deliberate
			fTotalForward += f * m_assetWeights[n];
		}
		return fTotalForward;
	}	
}

double MlEqAsset::GetBumpedForward( long nStart, long nMaturity, double bumpedSpot, bool bReinvestDivs) const
{
	return GetNaturalForward( nStart, nMaturity, bumpedSpot, bReinvestDivs );
}

double MlEqAsset::GetNaturalForward(long nStart, long nMaturity, bool bReinvestDivs) const
{	
	double f = 0.0;

	long nToday = GetDateHandle()->GetDate();

	if (nStart == nToday){
		double f = GetNaturalForward(nMaturity, bReinvestDivs);
		return f;
	} else if (nStart < nToday){
		throw("startdate of forward is lower than current date");
	}
	f = GetNaturalForward(nStart, nMaturity, GetNaturalForward(nToday, nStart, bReinvestDivs), bReinvestDivs);
	return f;
}

double MlEqAsset::GetNaturalForward(long nMaturity, bool bReinvestDivs) const
{	
	double f = 0.0;
	long nToday = GetDateHandle()->GetDate();
	f = GetNaturalForward(nToday, nMaturity, GetNaturalSpot(nToday), bReinvestDivs);
	return f;
}

//
//	GetQuantoForward
//
double MlEqAsset::GetQuantoForward(long nStart, long nMaturity, bool bReinvestDivs) const
//	All GetQuantoForward functions must call into this one.
{	
	double f = 0.0;
	double retval = 0.0;

	if (!IsBasket()){		
		double fNaturalForward = GetNaturalForward(nStart, nMaturity, bReinvestDivs);	// obviously this will never analyse a basket in this case
    		
		if (m_fxPay == m_fxNatural && m_fxPay == m_fxComposite){
			retval =  fNaturalForward;
		} else {

			double fxForward;
			// it doesn't matter if GetNaturalForward has a tendency to analyse basket cases since FX assets are never baskets
			// and we never reinvest dividends in an FX forward
			fxForward = m_fxComposite->GetNaturalForward(nStart, nMaturity, false) /
						m_fxNatural->GetNaturalForward(nStart, nMaturity, false);

			MlEqConstDateHandle d2d = m_pVols->getDateToDouble();
			MlEqStrike atmStrike(fNaturalForward);

			double quantoAdj;

			int SpotOrForward = 0;
			double volfx;

			long nToday = GetDateHandle()->GetDate();

			volfx = m_fxComposite->getNaturalATMVol(nStart,nMaturity,SpotOrForward);

			quantoAdj = GetStockCompositeFxVariance(atmStrike,nStart,nMaturity) 
						- GetStockPayFxVariance(atmStrike,nStart,nMaturity) +
						pow(volfx,2.0)
						- GetCompositePayFxVariance(atmStrike,nStart,nMaturity) 
						- GetCompositeNaturalFxVariance(atmStrike,nStart,nMaturity)
						+ GetNaturalPayFxVariance(atmStrike,nStart,(double)nMaturity);

			double mat = d2d->GetYearFraction(nMaturity) - d2d->GetYearFraction(nStart);
			quantoAdj *= mat;
			
			retval = fNaturalForward * fxForward*exp(quantoAdj);
		}
	} else {
		double fTotalQuantoForward = 0.0;			
		for (int n = 0; n < m_assets.size(); n++){
			double f = m_assets[n]->GetQuantoForward(nStart, nMaturity, bReinvestDivs);	// recursion here is deliberate
			fTotalQuantoForward += f * m_assetWeights[n];
		}
		retval = fTotalQuantoForward;
	}

	return retval;
}

double MlEqAsset::GetQuantoForward(long nMaturity, bool bReinvestDivs) const
{
	long nToday = GetDateHandle()->GetDate();
	double res = GetQuantoForward(nToday, nMaturity, bReinvestDivs);
	return res;
}


/////////////////////////////////////////////////////////////////////////////
//	Other functions from here.
//

DataSourceEnum MlEqAsset::GetDataSource(void) const
{
	return m_ds;
}

void MlEqAsset::PutDataSource(DataSourceEnum ds)
{
	m_ds = ds;
}

const MlEqZeroCurveHandle MlEqAsset::GetCompositeZeroCurve(bool bThrow) const
{
	if (!m_pzcComposite && bThrow) throw "No composite currency zero curve is defined for asset '" + m_szName + "'";	
	return m_pzcComposite;
}
const MlEqZeroCurveHandle MlEqAsset::GetPayZeroCurve(bool bThrow) const
{
	if (!m_pzcPay && bThrow) throw "No pay currency zero curve is defined for asset '" + m_szName + "'";	
	return m_pzcPay;
}
const MlEqZeroCurveHandle MlEqAsset::GetNaturalZeroCurve(bool bThrow) const
{
	if (!m_pzcNatural && bThrow) throw "No natural currency zero curve is defined for asset '" + m_szName + "'";
	return m_pzcNatural;
}

const MlEqAssetHandle MlEqAsset::GetFxNaturalAsset(void) const
{
	return m_fxNatural;
}
const MlEqAssetHandle MlEqAsset::GetFxPayAsset(void) const
{
	return m_fxPay;
}
const MlEqAssetHandle MlEqAsset::GetFxCompositeAsset(void) const
{
	return m_fxComposite;
}

const std::string& MlEqAsset::GetName(void) const
{
	return m_szName;
}

const std::vector<MlEqAssetHandle>& MlEqAsset::GetAssets(void) const
{
	return m_assets;
}

const MlEqCorrelationMatrixHandle MlEqAsset::GetCorrelationMatrix() const
{
	return m_hCorrelation;
}


void MlEqAsset::PutCorrelationMatrix(MlEqCorrelationMatrixHandle hCorrelation)
{
	m_hCorrelation = hCorrelation;
}

MlEqDividendScheduleHandle MlEqAsset::GetDividendSchedule(void) const
{
	return m_pdivs;
}









/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetSpot
**	Returns: double
**	Action : returns  spot in natural currency
**           
****************************************************************/


double MlEqAsset::GetNaturalSpot(long nDate) const
{		
	if (!m_pSpots) throw "no spot schedule defined for asset " + m_szName;
	if (nDate > GetDateHandle()->GetDate()){
		throw "No spot value can be obtained for asset '" + m_szName + "' as the date requested (" + MlEqDate(nDate).GetString() + ") is in the future";
	}
	try {
		double f = m_pSpots->GetValueAt(nDate);
		return f;	
	} catch (std::string){
		throw "The spot schedule for asset '" + m_szName + "' has no value defined on " + MlEqDate(nDate).GetString();
	}
}

/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetSpot
**	Returns: double
**	Action : returns  spot in natural currency
**           
****************************************************************/


double MlEqAsset::GetSpot(long nDate) const
{		
	double spot  = GetNaturalSpot(nDate);
	double Cspot = m_fxComposite->GetNaturalSpot(nDate);
	double Uspot = m_fxNatural->GetNaturalSpot(nDate);

	double res = spot*Cspot/Uspot;
	return res;
}


/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetSpots
**	Returns: double
**	Action : returns  spot in natural currency
**           
****************************************************************/


void MlEqAsset::GetSpots(long nToday, const std::vector<long>& vectorDates, MlEqSpotSchedule& SpotSchedule) const
{
	// ToDo - get this working for baskets. Then we will probably want to cache the schedule
	MlEqSpotSchedule						ssFxComposite, ssFxNatural;	

	GetNaturalSpots(nToday, vectorDates, SpotSchedule);	
	m_fxComposite->GetNaturalSpots(nToday, vectorDates, ssFxComposite);
	m_fxNatural->GetNaturalSpots(nToday, vectorDates, ssFxNatural);
	if (SpotSchedule.size() != ssFxComposite.size() || SpotSchedule.size() != ssFxNatural.size()) throw "Unexpected error in MlEqAsset::GetSpots()";		
	for (std::vector<long>::const_iterator it = vectorDates.begin(); it != vectorDates.end(); it++){
		if (ssFxNatural[*it]) SpotSchedule[*it] *= ssFxComposite[*it] / ssFxNatural[*it];
	}	
}

void MlEqAsset::GetSpotsBeforeToday(long nToday, const std::vector<long>& vectorDates, std::vector<double>* pSpots) const
{
	std::vector<double>					CSpots;
	std::vector<double>					NSpots;
	
	GetNaturalSpotsBeforeToday(nToday, vectorDates, pSpots);
	
	m_fxComposite->GetNaturalSpotsBeforeToday(nToday, vectorDates, &CSpots);
	m_fxNatural->GetNaturalSpotsBeforeToday(nToday, vectorDates, &NSpots);

	for ( int i = 0 ; i < pSpots->size(); i++ ){
		(*pSpots)[i] *= (CSpots)[i]/(NSpots)[i];
	}

}

void MlEqAsset::GetSpotsBeforeToday(long nToday, const std::vector<long>& anDates, CMatrix* pm) const
{
	std::vector<double>					an;
	GetSpotsBeforeToday(nToday, anDates, &an);
	pm->resize(1, an.size(), 0.0);	
	for (long n = 0; n < an.size(); n++){
		(*pm)[0][(int)n] = an[n];
	}
}

void MlEqAsset::GetSpotsBeforeToday(long nToday, const GVector<long>& anDates, CMatrix* pm) const
{
	std::vector<long>					an;
	an.resize(anDates.getsize());
	for (long n = 0; n < an.size(); n++){
		an[n] = anDates[n];
	}
	GetSpotsBeforeToday(nToday, an, pm);
}

void MlEqAsset::GetNaturalSpots(long nToday, const std::vector<long>& vectorDates, MlEqSpotSchedule& SpotSchedule) const
{
	// ToDo - get this working for baskets
	if (!m_pSpots) throw "no spot schedule defined for asset " + m_szName;
	SpotSchedule.clear();
	for (std::vector<long>::const_iterator it = vectorDates.begin(); it != vectorDates.end(); it++){
		double f = 0.0;
		try {
			f = m_pSpots->GetValueAt(*it);
		} catch (std::string& sz) {		
			// only a problem if *it is before Today
			if (*it < nToday) throw sz;
		}
		SpotSchedule[*it] = f;		
	}
}

void MlEqAsset::GetNaturalSpotsBeforeToday(long nToday, const std::vector<long>& vectorDates, std::vector<double>* pSpots) const
{
	// ToDo - get this working for baskets
	// ToDo - use MlEqSpotSchedule::GetValuesBeforeDate
	if (!m_pSpots) throw "no spot schedule defined for asset '" + m_szName + "'";
	pSpots->clear();
	for (std::vector<long>::const_iterator it = vectorDates.begin(); it != vectorDates.end(); it++){		
		if (*it >= nToday) break;
		try {
			if (m_szName == "USD.USD"){
				pSpots->push_back(1.0);
			} else {
				pSpots->push_back(m_pSpots->GetValueAt(*it));
			}
		} catch (const std::string&){
			throw "No value is defined at " + MlEqDate(*it).GetString() + " for asset '" + m_szName + "'";
		}
	}
}


/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetVolatilityStructure
**	Returns: MlEqVolatilityStructure
**	Action : returns  volstructure
**           
****************************************************************/

MlEqVolatilityStructureHandle MlEqAsset::GetVolatilityStructure() const
{
	if (!m_pVols) throw "No volatility structure defined for asset '" + m_szName + "'";

	return m_pVols;
}


/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetInternalDate / GetDateHandle
**	Action : returns  currentDate
**           
****************************************************************/

long MlEqAsset::GetInternalDate(void) const
{
	return m_nDate;	// Never use GetDateHandle here.
}

MlEqConstDateHandle MlEqAsset::GetDateHandle(void) const
{
	if (m_bIsBasket){
		// Basket case - use the volatility structure associated with the first basket constituent.
		return m_assets[0]->GetDateHandle();
	} else if (!m_pVols){
		throw "Cannot get the date for asset '" + m_szName + "' because it does not contain a volatility structure.";
	} else {
		MlEqConstDateHandle hDate = m_pVols->getDateToDouble();
		if (m_nDate){
			// Check the date member of the asset is consistent with h.
			if (m_nDate != hDate->GetDate()) throw "Cannot get the date for asset '" + m_szName + "' because its date is " + MlEqDate(m_nDate).GetString() + " but its volatility structure date is " + hDate->GetString();
		}			
		return hDate;
	}
}

const CVector& MlEqAsset::GetWeights(void) const
{
	return m_assetWeights;
}

bool MlEqAsset::HasVolatilityStructure(void) const
{
	return !!m_pVols ? true : false;
}

bool MlEqAsset::IsCalibrated(void) const
{
	return m_bCalibrated;
}

void MlEqAsset::PutVolatilityStructure(MlEqVolatilityStructureHandle pVols)
{
	m_pVols = pVols;	
}


/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetNaturalVolatility
**	Returns: double
**	Action : vol in natural currency
**           
****************************************************************/


double	MlEqAsset::GetNaturalVolatility(const MlEqStrike& strike,  const DATE& date, BidAskMidEnum bidOrAsk) const 
{
	MlEqVolatilityStructureHandle  pVolstruct  = GetVolatilityStructure();
	double volS = pVolstruct->getVol(strike,date,bidOrAsk);
	return volS;
}

/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetStockCompositeFxVariance
**	Returns: double
**	Action : vol in composite currency
**           
****************************************************************/

double	MlEqAsset::GetStockCompositeFxVariance(const MlEqStrike& strike,long nStart,const DATE& date,BidAskMidEnum bidOrAsk) const
{
	long nToday = GetDateHandle()->GetDate();

	if ( nStart < nToday ){
		throw("startDate is smaller than valuation date");
	}

	MlEqVolatilityStructureHandle  pVolstruct  = GetVolatilityStructure();

	double shortFwd = GetNaturalForward(nStart, false);
	double volS = m_pVols->getFutureVol(strike,nStart,(double) date,shortFwd);
	MlEqVolatilityStructureHandle vol =   m_fxComposite->GetVolatilityStructure();
	int SpotOrForward = 1;
	double volC = m_fxComposite->getNaturalATMVol(nStart,date,SpotOrForward,bidOrAsk) ;

//  start calculating correlations

	double cSC = m_hCorrelation->GetCorrelation(m_szName, m_fxComposite->GetName());

	double quantoCompVar;
	quantoCompVar = volS*volC*cSC;
	return quantoCompVar;

}


/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetStockCompositeFxVariance
**	Returns: double
**	Action : vol in composite currency
**           
****************************************************************/

double	MlEqAsset::GetStockCompositeFxVariance(const MlEqStrike& strike,const DATE& date,BidAskMidEnum bidOrAsk) const
{
	long nToday = GetDateHandle()->GetDate();
	double	res = GetStockCompositeFxVariance(strike,nToday,date,bidOrAsk);
	return res;

}

/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetCompositeNaturalFxVariance
**	Returns: double
**	Action : vol in composite currency
**           
****************************************************************/

double	MlEqAsset::GetCompositeNaturalFxVariance(const MlEqStrike& strike,long nStart,const DATE& date,BidAskMidEnum bidOrAsk) const
{

	long nToday = GetDateHandle()->GetDate();

	if ( nStart < nToday ){
		throw("startDate is smaller than valuation date");
	}

	MlEqVolatilityStructureHandle vol =   m_fxNatural->GetVolatilityStructure();

	double shortFwd = GetNaturalForward(nStart, false);
	double volS = vol->getFutureVol(strike,nStart,(double) date,shortFwd) ;

	vol =   m_fxComposite->GetVolatilityStructure();
	int SpotOrForward = 1;
	double volC = m_fxComposite->getNaturalATMVol(nStart,date,SpotOrForward,bidOrAsk) ;

// start calculating correlations

	double cSC = m_hCorrelation->GetCorrelation(m_fxNatural->GetName(), m_fxComposite->GetName());	

	double quantoCompVar;
	quantoCompVar = volS*volC*cSC;
	return quantoCompVar;


}

/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetCompositeNaturalFxVariance
**	Returns: double
**	Action : vol in composite currency
**           
****************************************************************/

double	MlEqAsset::GetCompositeNaturalFxVariance(const MlEqStrike& strike,const DATE& date,BidAskMidEnum bidOrAsk) const
{

	long nToday = GetDateHandle()->GetDate();
	double res = GetCompositeNaturalFxVariance(strike,date,bidOrAsk);
	return res;
}

/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetCompositePayFxVariance
**	Returns: double
**	Action : vol in composite currency
**           
****************************************************************/

double	MlEqAsset::GetCompositePayFxVariance(const MlEqStrike& strike,long nStart,const DATE& date,BidAskMidEnum bidOrAsk) const
{

	long nToday = GetDateHandle()->GetDate();

	if ( nStart < nToday ){
		throw("startDate is smaller than valuation date");
	}

	MlEqVolatilityStructureHandle vol =   m_fxPay->GetVolatilityStructure();

	double shortFwd = GetNaturalForward(nStart, false);
	double volS = vol->getFutureVol(strike,nStart,(double) date,shortFwd) ;

	vol =   m_fxComposite->GetVolatilityStructure();
	int SpotOrForward = 0;
	double volC = m_fxComposite->getNaturalATMVol(nStart,date,SpotOrForward,bidOrAsk) ;

// start calculating correlations

	double cSC = m_hCorrelation->GetCorrelation(m_fxPay->GetName(), m_fxComposite->GetName());	

	double quantoCompVar;
	quantoCompVar = volS*volC*cSC;
	return quantoCompVar;


}

/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetCompositePayFxVariance
**	Returns: double
**	Action : vol in composite currency
**           
****************************************************************/

double	MlEqAsset::GetCompositePayFxVariance(const MlEqStrike& strike,const DATE& date,BidAskMidEnum bidOrAsk) const
{

	long nToday = GetDateHandle()->GetDate();
	double res = GetCompositePayFxVariance(strike,date,bidOrAsk);
	return res;
}


/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetNatualPayFxVariance
**	Returns: double
**	Action : vol in composite currency
**           
****************************************************************/


double	MlEqAsset::GetNaturalPayFxVariance(const MlEqStrike& strike,long nStart, const DATE& date,BidAskMidEnum bidOrAsk) const
{

	long nToday = GetDateHandle()->GetDate();

	if ( nStart < nToday ){
		throw("startDate is smaller than valuation date");
	}

	MlEqVolatilityStructureHandle vol =   m_fxPay->GetVolatilityStructure();

	double shortFwd = GetNaturalForward(nStart, false);

	double volS = vol->getFutureVol(strike,nStart,(double) date,shortFwd) ;

	vol =   m_fxNatural->GetVolatilityStructure();

	int SpotOrForward = 0;
	double volC = m_fxNatural->getNaturalATMVol(nStart,date,SpotOrForward,bidOrAsk) ;

// start calculating correlations

	double cSC = m_hCorrelation->GetCorrelation(m_fxPay->GetName(), m_fxNatural->GetName());	

	double quantoCompVar;
	quantoCompVar = volS*volC*cSC;
	return quantoCompVar;

}

/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetNatualPayFxVariance
**	Returns: double
**	Action : vol in composite currency
**           
****************************************************************/


double	MlEqAsset::GetNaturalPayFxVariance(const MlEqStrike& strike,const DATE& date,BidAskMidEnum bidOrAsk) const
{

	long nToday = GetDateHandle()->GetDate();
	double res = GetNaturalPayFxVariance(strike,nToday,date,bidOrAsk);
	return res;
}

/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetStockPayFxVariance
**	Returns: double
**	Action : vol in composite currency
**           
****************************************************************/

double	MlEqAsset::GetStockPayFxVariance(const MlEqStrike& strike,int nStart,const DATE& date,BidAskMidEnum bidOrAsk) const
{
	long nToday = GetDateHandle()->GetDate();

	if ( nStart < nToday ){
		throw("startDate is smaller than valuation date");
	}

	MlEqVolatilityStructureHandle  vol  = GetVolatilityStructure();

	double shortFwd = GetNaturalForward(nStart, false);

	double volS = vol->getFutureVol(strike,nStart,(double) date,shortFwd) ;

	vol =   m_fxPay->GetVolatilityStructure();
	
	int SpotOrForward = 0;
	double volC = m_fxPay->getNaturalATMVol(nStart,date,SpotOrForward,bidOrAsk) ;

// start calculating correlations

	double cSC = m_hCorrelation->GetCorrelation(m_szName, m_fxPay->GetName());	

	double quantoCompVar;
	quantoCompVar = volS*volC*cSC;
	return quantoCompVar;

}

/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetStockPayFxVariance
**	Returns: double
**	Action : vol in composite currency
**           
****************************************************************/

double	MlEqAsset::GetStockPayFxVariance(const MlEqStrike& strike,const DATE& date,BidAskMidEnum bidOrAsk) const
{
	long nToday = GetDateHandle()->GetDate();
	double res = GetStockPayFxVariance(strike,date,bidOrAsk);
	return res;
}
/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetStockNaturalFxVariance
**	Returns: double
**	Action : vol in composite currency
**           
****************************************************************/

double	MlEqAsset::GetStockNaturalFxVariance(const MlEqStrike& strike,long nStart,const DATE& date,BidAskMidEnum bidOrAsk) 
{
	long nToday = GetDateHandle()->GetDate();

	if ( nStart < nToday ){
		throw("startDate is smaller than valuation date");
	}

	MlEqVolatilityStructureHandle  vol  = GetVolatilityStructure();


	double shortFwd = GetNaturalForward(nStart, false);
	double volS = vol->getFutureVol(strike,nStart,(double) date,shortFwd) ;

	vol =   m_fxNatural->GetVolatilityStructure();
	int SpotOrForward = 0;
	double volU = m_fxNatural->getNaturalATMVol(nStart,date,SpotOrForward,bidOrAsk) ;

// start calculating correlations

	double cSU = m_hCorrelation->GetCorrelation(m_szName, m_fxNatural->GetName());	

	double quantoCompVar;
	quantoCompVar = volS*volU*cSU;
	return quantoCompVar;

}

/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetStockNaturalFxVariance
**	Returns: double
**	Action : vol in composite currency
**           
****************************************************************/

double	MlEqAsset::GetStockNaturalFxVariance(const MlEqStrike& strike,const DATE& date,BidAskMidEnum bidOrAsk) 
{

	long nToday = GetDateHandle()->GetDate();
	double res = GetStockNaturalFxVariance(strike,nToday,date,bidOrAsk) ;
	return res;

}
	
double	MlEqAsset::GetCompositeVolatility(long nStartDate,const DATE& date,double volS ,BidAskMidEnum bidOrAsk) const
{
//  inputs should be composite strikes!

    if ( m_pzcComposite == m_pzcNatural && m_pzcComposite == m_pzcPay){
		return volS;
	}

	MlEqStrike atmStrike(1.0);// this is atm spot

	long nToday = GetDateHandle()->GetDate();


	int SpotOrForward = 1;
	double volC = m_fxComposite->getNaturalATMVol(nStartDate,date,SpotOrForward,bidOrAsk);
	double volP = m_fxPay->getNaturalATMVol(nStartDate,date,SpotOrForward,bidOrAsk);
	double volU = m_fxNatural->getNaturalATMVol(nStartDate,date,SpotOrForward,bidOrAsk);

//  start calculating correlations

	double cSC = m_hCorrelation->GetCorrelation(m_szName, m_fxComposite->GetName());	
	double cSP = m_hCorrelation->GetCorrelation(m_szName, m_fxPay->GetName());	
	double cSU = m_hCorrelation->GetCorrelation(m_szName, m_fxNatural->GetName());	

	double cCP = m_hCorrelation->GetCorrelation(m_fxComposite->GetName(),m_fxPay->GetName());	
	double cUC = m_hCorrelation->GetCorrelation(m_fxNatural->GetName(),m_fxComposite->GetName());	
	double cUP = m_hCorrelation->GetCorrelation(m_fxNatural->GetName(),m_fxPay->GetName());	

//  do something here if correlations are timedependent

	double quantoCompVar;

	quantoCompVar = volS*volS - volS*volU*cSU - volS*volU*cSU +
					volU*volU + volS*volC*cSC - volU*volC*cUC +
					volS*volC*cSC - volU*volC*cUC + volC*volC;

	if ( quantoCompVar < 0.0 ){
		throw("negative composite variance encountered in getNaturalVolatility");
	}

	quantoCompVar = sqrt(quantoCompVar);
	return quantoCompVar;


}

/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetStockNaturalFxVariance
**	Returns: double
**	Action : vol in composite currency
**           
****************************************************************/

double	MlEqAsset::GetNaiveFutureVolatility(const MlEqStrike& shortStrike,const MlEqStrike& longStrike,long nStartDate,const DATE& date,BidAskMidEnum bidOrAsk) const
{
// this function should be checked sos

	double	volS = GetNaturalNaiveFutureVolatility(shortStrike,longStrike,nStartDate,date,bidOrAsk);
	double	vol = GetCompositeVolatility(nStartDate,date,volS ,bidOrAsk);
	return vol;
}

/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetStockNaturalFxVariance
**	Returns: double
**	Action : vol in composite currency
**           
****************************************************************/

double	MlEqAsset::GetNaturalNaiveFutureVolatility(const MlEqStrike& shortStrike,const MlEqStrike& longStrike,long nStartDate,const DATE& date,BidAskMidEnum bidOrAsk) const
{
// this function should be checked sos

	MlEqVolatilityStructureHandle pVol =  GetVolatilityStructure();

	double volS = pVol->getNaiveFutureVol(shortStrike,longStrike,nStartDate,date,bidOrAsk) ;
	return volS;
}



/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetCompositeVolatility
**	Returns: double
**	Action : vol in composite currency
**           
****************************************************************/

double	MlEqAsset::GetCompositeVolatility(const MlEqStrike& strike,long nStartDate,const DATE& date,BidAskMidEnum bidOrAsk /*= MlEqVolatilityStructure::MidVol*/,MlEqStrikeHandle skewInvariantStrikeHandle) const
{

//  inputs should be composite strikes !

	long nToday = GetDateHandle()->GetDate();
	MlEqStrike fStrike;
	MlEqStrike::convertStrikes(fStrike,strike);
	double spotFx = m_fxComposite->GetSpot(nToday)/m_fxNatural->GetSpot(nToday);
	fStrike.m_strike /= spotFx;

	double refSpot = GetNaturalSpot(nToday);

	MlEqVolatilityStructureHandle pVol =  GetVolatilityStructure();

	MlEqConstDateHandle		d		= GetDateHandle();
	DayCountConventionEnum  conv	= d->GetDayCountConvention() ;

	MlEqDateHandle startDateHandle	= new MlEqDate(nStartDate,conv);

	double fwdShort = GetNaturalForward(nStartDate, false); 

	MlEqStrikeHandle strikeInv;
	strikeInv = MlEqVolatilityStructure::setUpDefaultStrikeInvariant(skewInvariantStrikeHandle,*this,startDateHandle,date);

	double volS =  pVol->getFutureVol(*this,fStrike,*startDateHandle,date,strikeInv,bidOrAsk);
	double res	= GetCompositeVolatility(nStartDate,date,volS, bidOrAsk);

	return res;
}



/****************************************************************
**	Class  : MlEqAsset 
**	Routine: getATMVol
**	Returns: double
**	Action : returns atm vol
**           
****************************************************************/

double	MlEqAsset::getNaturalATMVol( long nMaturity,int SpotOrForward,BidAskMidEnum bidOrAsk) const
{
	double strike;
	long nToday = GetDateHandle()->GetDate();
	if ( SpotOrForward != 1 ){
		strike = GetNaturalSpot(nToday);
	} else{
		strike = GetNaturalForward(nMaturity, false);
	}

	MlEqStrike hstrike(strike);
	double  vol = GetNaturalVolatility(hstrike,(double)nMaturity,bidOrAsk);
	return vol;
}

/****************************************************************
**	Class  : MlEqAsset 
**	Routine: getATMVol
**	Returns: double
**	Action : returns atm vol
**           
****************************************************************/

long	MlEqAsset::GetCurrentDate()const
{

	long nToday = GetDateHandle()->GetDate();
	return nToday;
}

/****************************************************************
**	Class  : MlEqAsset 
**	Routine: getATMVol
**	Returns: double
**	Action : returns atm vol
**           
****************************************************************/

double	MlEqAsset::getNaturalATMVol(long nStart, long nMaturity,int SpotOrForward,BidAskMidEnum bidOrAsk) const
{	
	double strike;
	long nToday = GetDateHandle()->GetDate();

	if ( nStart == nToday ){
		double res = getNaturalATMVol(nMaturity,SpotOrForward,bidOrAsk);
		return res;
	}
	else if ( nStart < nToday ){
		throw("StartDate is smaller than current date in getNaturalVol");
	}

	if ( SpotOrForward != 1 ){
		strike = GetNaturalForward(nStart, false);
//		this is because the forward to nStart is our best guess for the spot
	}
	else{
		strike = GetNaturalForward(nMaturity, false);
	}
	
	MlEqVolatilityStructureHandle  pVol  = GetVolatilityStructure();
	MlEqStrike hstrike(strike);

	double shortFwd = GetNaturalForward(nStart, nMaturity, false);
	double  vol = pVol->getFutureVol(hstrike,(double)nStart,(double)nMaturity,shortFwd,bidOrAsk) ;
	return  vol;
}


/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetCompositeVolatility
**	Returns: double
**	Action : vol in composite currency
**           
****************************************************************/

double	MlEqAsset::GetCompositeVolatility(const MlEqStrike& strike,const DATE& date,BidAskMidEnum bidOrAsk) const 
{

//  inputs should be composite strikes!

	long nToday = GetDateHandle()->GetDate();
	double res = GetCompositeVolatility(strike,nToday,date,bidOrAsk);
	return res;
}


/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetVolatility
**	Returns: double
**	Action : vol 
**           
****************************************************************/

double	MlEqAsset::GetVolatility(const MlEqStrike& strike,const DATE& date,BidAskMidEnum bidOrAsk) const
{

	double volS = GetNaturalVolatility(strike,date,bidOrAsk) ;

	if ( !IsCalibrated() ){
		return volS ;
	}

	return GetCompositeVolatility(strike,date,bidOrAsk);

}

/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetVolatility
**	Returns: double
**	Action : vol 
**           
****************************************************************/

double	MlEqAsset::GetNaturalVolatility(const MlEqStrike& strike, long nStart,const DATE& date, BidAskMidEnum bidOrAsk ,MlEqStrikeHandle skewInvariantStrikeHandle ) const
{
	MlEqVolatilityStructureHandle  pVolstruct  = GetVolatilityStructure();

	

	MlEqDate StartDate(nStart, GetDateHandle()->GetDayCountConvention());


	double volS = pVolstruct->getFutureVol(*this,strike,StartDate,(long)date, skewInvariantStrikeHandle,bidOrAsk) ;	
	

	return volS;
}

/****************************************************************
**	Class  : MlEqAsset 
**	Routine: GetVolatility
**	Returns: double
**	Action : vol 
**           
****************************************************************/

double	MlEqAsset::GetVolatility(const MlEqStrike& strike,long nStart, const DATE& date,BidAskMidEnum bidOrAsk,MlEqStrikeHandle skewInvariantStrikeHandle ) const
{
	double volS = GetNaturalVolatility(strike,nStart,date,bidOrAsk) ;

	if ( !IsCalibrated() ){
		return volS ;
	}

	return GetCompositeVolatility(strike,nStart,date,bidOrAsk,skewInvariantStrikeHandle );
}



void MlEqAsset::PutInternalDate(long nDate)
{
	m_nDate = nDate;
}


/****************************************************************
**	Class  : MlEqAsset 
**	Routine: Calibrate
**	Returns: void
**	Action : calibrates asset 
**           
****************************************************************/


CVector MlEqAsset::Calibrate(DATE maturityDate, MlEqDate& dateToDouble, vector < vector < MlEqStrikeHandle > >& pVolStrikes)
{

//	vector < vector < MlEqStrikeHandle > >& pVolStrikes : [iasset][istrike]
//	iasset = 0 corresponds to the basket strikes

	const std::vector<MlEqAssetHandle>		assets = GetAssets();
	const CVector	assetWeights = 	GetWeights();		

	return CalculateEuropeanThreeMomentBasket(assets,assetWeights,maturityDate,pVolStrikes);

	
/*	int nassets = m_assetWeights.getsize();

	if ( pVolStrikes.size() != nassets+1) {
		throw("dimension of strikevector must be number of assets +1");
	}

	for ( int i = 0 ; i < nassets+1; i++ )
	{
		if (  pVolStrikes[0].size() != pVolStrikes[i].size() ){
			throw("same number of strikes must be entered per asset");
		}
	}

	int nstrikes = pVolStrikes[0].size() ;

	vector < DATE > asianDates;
	asianDates.resize(1);

	CVector asianWeights(1);
	asianWeights[0]		= 1.0;
	asianDates[0]		= maturityDate;

	CVector res;
	CVector impliedVols(nstrikes);

	double mat = dateToDouble.GetYearFraction(maturityDate);


	double fwd = GetQuantoForward(dateToDouble.GetDate(),maturityDate, false);

	vector < MlEqStrikeHandle >  pStrikes(nassets);

	vector < MlEqStrikeHandle > basketStrike(1);

	for ( int istrike = 0 ; istrike < pVolStrikes[0].size(); istrike++ )
	{

		MlEqThreeMoment tm;

		for ( int i = 0 ; i < nassets; i++ ){
			pStrikes[i] = pVolStrikes[i+1][istrike];
		}



		tm.initialize(
					  *this,
					  asianWeights,
					  asianDates,
					  *hCorrelation,
					  dateToDouble,
					  pStrikes);



		MlEqStrike xstrike;

		basketStrike[0] = pVolStrikes[0][istrike];
		res =  tm.calculatePrice(basketStrike);
		
		MlEqStrike::convertStrikes(xstrike,*basketStrike[0]);

		impliedVols[istrike] =  MlEqBSImpliedVol(res[0],fwd,mat,xstrike.m_strike,1.0,1);

	}


	m_bCalibrated = true;

	return impliedVols;
*/
}

// do this over many slices and create volcurve

// create volcurve

