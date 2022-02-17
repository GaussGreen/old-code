//	NapoleonPricer.cpp : Implementation of CNapoleonPricer
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MlEqObjects.h"
#include "NapoleonPricer.h"
#include "MonteCarlo.h"
#include "MlEqPde.h"

#undef max
#undef min











void MlEqNapoleonPricer::initialize(MlEqAssetHandle				hUnderlying,
									const GVector<long>&		anPayoffDates,
									const GVector<bool>&		abStartsPeriod,			// [date]
									const GVector<long>&		anCallPuts,				// [date]
									const CVector&				afStrikes,				// [date]
									const CVector&				afLocalFloors,			// [date]
									const CVector&				afLocalCaps,			// [date]
									const CVector&				afWeights,				// [date]
									const CVector&				afCoupons,				// [period]
									const CVector&				afGearings,				// [period]									
									const CVector&				afPeriodFloors,			// [period]
									const CVector&				afPeriodCaps,			// [period]					
									GVector<long>				anPayDatesArr,			// [period, optional]
									bool						bDelayPaymentsToEnd)
{
	m_payoffDates = anPayoffDates;
	m_abStartsPeriod = abStartsPeriod;
	m_anCallPuts = anCallPuts;
	m_afStrikes = afStrikes;
	m_afLocalFloors = afLocalFloors;
	m_afLocalCaps = afLocalCaps;
	m_afWeights = afWeights;
	m_afCoupons = afCoupons;
	m_afGearings = afGearings;	
	m_afPeriodFloors = afPeriodFloors;
	m_afPeriodCaps = afPeriodCaps;	

	// All the arrays indexed [date] should have one element less than the number of dates (although we relax this case for a truncated abStartsPeriod).		
	if (m_abStartsPeriod.getsize() > m_payoffDates.getsize() - 1) throw "The number of period starts cannot be greater than the number of payoff dates";
	if (m_anCallPuts.getsize() != m_payoffDates.getsize() - 1) throw "The number of call/puts is inconsistent with the number of payoff dates";
	if (m_afStrikes.getsize() != m_payoffDates.getsize() - 1) throw "The number of strikes is inconsistent with the number of payoff dates";		
	if (m_afLocalFloors.getsize() != m_payoffDates.getsize() - 1) throw "The number of local floors is inconsistent with the number of payoff dates";
	if (m_afLocalCaps.getsize() != m_payoffDates.getsize() - 1) throw "The number of local caps is inconsistent with the number of payoff dates";
	if (m_afWeights.getsize() != m_payoffDates.getsize() - 1) throw "The number of weights is inconsistent with the number of payoff dates";
		
	// Set m_nPeriods
	m_nPeriods = 0;
	for (int n = 0; n < m_abStartsPeriod.getsize(); n++){
		if (m_abStartsPeriod[n]) m_nPeriods++;
	}
	if (!m_abStartsPeriod.getsize() || !m_nPeriods) throw "You have not declared any Napoleon periods";
	if (!m_abStartsPeriod[0]) throw "The first element in the periods array must be 1";

		
	// set m_afPayDiscountFactors (this and only this uses hUnderlying, anPayDatesArr and bDelayPaymentsToEnd)	
	if (anPayDatesArr.getsize() && anPayDatesArr.getsize() != m_nPeriods) throw "The number of pay dates must match the number of periods";
	if (bDelayPaymentsToEnd){
		// Set all the discount factors to the same value.
		// (The payment date defaults to the final payoff date.)
		long nPayDate = anPayDatesArr.getsize() ? anPayDatesArr[anPayDatesArr.getsize() - 1] : m_payoffDates[m_payoffDates.getsize() - 1];
		double f = hUnderlying->GetPayZeroCurve(true)->GetDiscountFactor(hUnderlying->GetDateHandle()->GetDate(), nPayDate);
		m_afPayDiscountFactors.resize(m_nPeriods, f);
	} else {		
		if (!anPayDatesArr.getsize()){
			// Assume that we receive payment for the previous Napoleon period at the start of the next period,
			// with the exception of the final period which we receive at the last date.
			int nPeriod = 0;
			anPayDatesArr.resize(m_nPeriods);
			for (int n = 1/*first element always starts a period but we want the start of the next period*/; n < m_abStartsPeriod.getsize(); n++){
				if (m_abStartsPeriod[n]) anPayDatesArr[nPeriod++] = m_payoffDates[n];
			}
			anPayDatesArr[m_nPeriods - 1] = m_payoffDates[m_payoffDates.getsize() - 1];	// get the last date
		}						
		// now set the discount factors from anPayDatesArr
		m_afPayDiscountFactors.resize(m_nPeriods);
		for (int nPeriod = 0; nPeriod < m_nPeriods; nPeriod++){
			double f = hUnderlying->GetPayZeroCurve(true)->GetDiscountFactor(hUnderlying->GetDateHandle()->GetDate(), anPayDatesArr[nPeriod]);
			m_afPayDiscountFactors[nPeriod] = f;
		}
	}

	// All the arrays indexed [period] should have the same number of elements as there are true elements in m_abStartsPeriod	
	if (m_afCoupons.getsize() != m_nPeriods) throw "The number of coupons must match the number of periods";
	if (m_afGearings.getsize() != m_nPeriods) throw "The number of gearings must match the number of periods";	
	if (m_afPeriodFloors.getsize() != m_nPeriods) throw "The number of period floors must match the number of periods";
	if (m_afPeriodCaps.getsize() != m_nPeriods) throw "The number of period caps must match the number of periods";
	if (m_afPayDiscountFactors.getsize() != m_nPeriods) throw "Internal error: The pay discount factors have not been set up correctly";
}

void MlEqNapoleonPricer::payout(CMatrix& values, CMatrix& pathArray, CVector& discountFactors, MlEqMonteCarlo& mc)
//	pathArray -	[0][Date]
//	discountFactors - [0][Date]
{			
	int									nDate = 1;						// date index (indexed s.t. the return is spot[nDate] / spot[nDate - 1])
	double								fPayoff = 0.0;
	long								nPeriod = 0;

	while (nDate < m_payoffDates.getsize())
	{
		if (!m_abStartsPeriod[nDate - 1]) throw "Unhandled exception - starts period array is inconsistent with the rest of the payoff data";
		
		// initialise to start a new period
		double		fMinReturn = std::numeric_limits<double>::max();
		nPeriod++;

		// process one period - get fMinReturn for this period
		do {		
			double fReturn = m_afWeights[nDate - 1] * std::min(m_afLocalCaps[nDate - 1], std::max(m_afLocalFloors[nDate - 1], m_anCallPuts[nDate - 1] * (pathArray[0][nDate] / pathArray[0][nDate - 1] - m_afStrikes[nDate - 1])));
			fMinReturn = std::min(fMinReturn, fReturn);
			nDate++;
		}
		while (nDate < m_payoffDates.getsize() /*finished*/ && (nDate - 1 >= m_abStartsPeriod.getsize() || !m_abStartsPeriod[nDate - 1]));
		
		// get the contribution of this to the payoff...
		double		fR = fMinReturn * m_afGearings[nPeriod - 1] + m_afCoupons[nPeriod - 1];
		double		fPayoff_period = std::min(m_afPeriodCaps[nPeriod - 1], std::max(m_afPeriodFloors[nPeriod - 1], fR));				
		fPayoff_period *= m_afPayDiscountFactors[nPeriod - 1];
		
		// ...and add it
		fPayoff += fPayoff_period;
	}
	
	values[0][0] = fPayoff;
}

void MlEqNapoleonPricer::setUp(CMatrix& value, MlEqMonteCarlo& mc)
{
	value.resize(1, 2);
	m_results.resize(1, 2);
	product::setUp(value,mc);
};











void DeltaHedgedOptionPortfolio::initialise(	
												MlEqAssetHandle hAsset,
												long startDate,
												GVector<long>& settlementDates,
												std::string	&szTenor,
												BusinessDayConventionEnum	bdc,
												std::string					szCalendar,
												CMatrix	&					portfolio,
												CVector	&					bookingVols,
												double						bookingRate,
												double						bookingDiv,
												DayCountConventionEnum		dcc,
												long						includeDeltaTerm,
												bool						gaussIntegrator,
												int							ngauss,
												double						nstdev,
												double						eps,
												long						approximationDate,
												int							ntimeGauss
											)
{
		m_szTenor			= szTenor;
		m_bdc				= bdc;
		m_szCalendar		= szCalendar;
		m_portfolio			= portfolio;
		m_bookingVols		= bookingVols;
		m_bookingRate		= bookingRate;
		m_bookingDiv		= bookingDiv ;
		m_settlementDates	= settlementDates;
		m_dcc				= dcc;
		m_startDate			= startDate;
		m_hAsset			= hAsset;
		m_includeDeltaTerm	= includeDeltaTerm;

		int n = m_bookingVols.getsize();
		if ( n != m_portfolio.rows()){
			throw("incorrect dimesnioning of input encountered in setting up deltahedged portfolio");
		};



	MlEqConstDateHandle	currentDateHandle	= hAsset->GetDateHandle();	
	long currentDate = currentDateHandle->GetDate();

	if ( m_startDate >= m_settlementDates[m_settlementDates.getsize()-1] ){
		throw("what is going on horsey: the startdate is greater or equal to the last settlement date");
	}


//  create pnl recording dates

	MlEqDateHandle hFirstPnlDate;
	hFirstPnlDate = new MlEqDate(m_startDate,m_dcc);
	
	long lastDate = m_settlementDates[m_settlementDates.getsize() - 1];
	vector<double> vRecordingDates;	
	for (MlEqDateHandle hDate = hFirstPnlDate, hNextDate; hDate->GetDate() < lastDate; hDate = hNextDate)
	{	
		hNextDate = hDate->AddTenor(m_szTenor, m_bdc);
		if ( hNextDate->GetDate() > lastDate ){
			break;
		}

		if ( hNextDate->GetDate() == hDate->GetDate() ){
			throw("got stuck in the generation of PnL schedule: try different daycount convention horse");
		}

		vRecordingDates.push_back(hNextDate->GetDate());					
	}	
	

	int i = 0;
	m_recordingDates.resize(vRecordingDates.size());
	for (MlEqDateHandle hDate = hFirstPnlDate, hNextDate; hDate->GetDate() < lastDate; hDate = hNextDate)
	{	
		hNextDate = hDate->AddTenor(m_szTenor, m_bdc);
		if ( hNextDate->GetDate() > lastDate ){
			break;
		}
		m_recordingDates[i] = hNextDate;
		i++;
	}	

	if ( m_recordingDates.getsize() < 1 ){
		throw("what's going on horsey: U must have at list 1 PnL date");
	}

	int ndate = m_recordingDates.getsize();
	m_nDates = ndate;
	m_nMaturity = -9999;//m_settlementDates[m_settlementDates.getsize()-1];

	m_mapToRecordingDate.resize(m_settlementDates.getsize());

	int k = 0 ;
	m_mapToSettlement.resize(ndate);// this maps recoring date to settlement date
	m_mapToSettlement[0] = 0;
	for ( i = 1 ; i < ndate; i++ )
	{
		if ( m_recordingDates[i]->GetDate() > m_settlementDates[k] )
		{
			m_mapToRecordingDate[k] = i;
			k++;
		}
		m_mapToSettlement[i]	= k;
	}

	m_mapToRecordingDate[m_settlementDates.getsize()-1] = ndate;

	m_gaussIntegrator	= gaussIntegrator;
	m_ngauss			= ngauss;
	m_nstdev			=	nstdev;
	m_eps				= eps;
	m_approximationDate	= approximationDate;
	m_ntimeGauss		= ntimeGauss;

	m_sq = 1.0/sqrt(2.0*3.141592654);	
};







class deltaPortContainer
{
public:

	DeltaHedgedOptionPortfolio* vp;
	double						mat;
	long						nMaturity;
	MlEqAssetHandle				hAsset;
	int							indexStart;
	double						voleps;
	double						spoteps;
	bool						calcDeltaTerm;
	double						mktFwd;
	double						bookingDiv;
	double						bookingRate;
	int							idate;			   
	double						refSpot;
	bool						calcAddOn;
};








void integratePortfolio(double x, double z[]  ,void *vp)
{
	
	deltaPortContainer* dc  = (deltaPortContainer*) vp;

	DeltaHedgedOptionPortfolio* op	= dc->vp;
	double mat						= dc->mat;
	int nMaturity					= dc->nMaturity;
	MlEqAssetHandle pAsset			= dc->hAsset;
	int indexStart					= dc->indexStart;
	CVector& vols					= op->m_bookingVols;
	double   voleps					= dc->voleps;
	double   spoteps				= dc->spoteps;
	double	 mktFwd					= dc->mktFwd;
	bool     calcDeltaTerm			= dc->calcDeltaTerm;
	double bookingDiv				= dc->bookingDiv;
	double bookingRate				= dc->bookingRate;
	double spot						= x;
	double fwd						= spot*exp((bookingRate-bookingDiv)*mat);
	double df						= exp(bookingRate*mat);


	double value = 0.0;
	for ( int i = 0 ; i < op->m_portfolio.rows(); i++ ){
		value += op->m_portfolio[i][0]*
				 Bs(fwd,vols[i],mat,op->m_portfolio[i][1],df,op->m_portfolio[i][2]);
	}
	

	MlEqStrike strike(x);
	double vol = pAsset->GetVolatility(strike,nMaturity);

	MlEqStrike strike_up(x*(1.0+voleps));
	double vol_up = pAsset->GetVolatility(strike_up,nMaturity);

	MlEqStrike strike_do(x*(1.0-voleps));
	double vol_do = pAsset->GetVolatility(strike_do,nMaturity);


	double density = 


	( Bs(mktFwd,vol_up,mat,x*(1.0+voleps),1.0,1) +
	  Bs(mktFwd,vol_do,mat,x*(1.0-voleps),1.0,1) -
	  2.0*Bs(mktFwd,vol,mat,x,1.0,1)
	  );

	density /=  pow(x*voleps,2.0);

	z[1] = density*value;

	if ( !calcDeltaTerm ){
		z[2] = 0.0;
		return;
	}

	double deltaTerm = 0.0;
	for ( int i = 0 ; i < op->m_portfolio.rows(); i++ )
	{
		deltaTerm += op->m_portfolio[i][0]*
					 Bs(fwd*(1.0+spoteps),vols[i],mat,op->m_portfolio[i][1],df,op->m_portfolio[i][2]);
	}

	deltaTerm	=	(deltaTerm-value)/(x*spoteps);
	deltaTerm	*=	x*density;

	z[2] = deltaTerm;
}







double DeltaHedgedOptionPortfolio::price(CMatrix& values,MlEqAssetHandle hAsset)
{

	double value = 0.0;

	MlEqConstDateHandle	currentDateHandle	= hAsset->GetDateHandle();	
	long currentDate = currentDateHandle->GetDate();

	if ( m_startDate >= m_settlementDates[m_settlementDates.getsize()-1] ){
		throw("what is going on horsey: the startdate is greater or equal to the last settlement date");
	}

//  generate pnl recording dates

	int approximationDateIndex,i;
	int ndate = m_recordingDates.getsize();

	for ( i = 1 ; i < ndate; i++ )
	{
		if ( m_recordingDates[i]->GetDate() > m_approximationDate ){
			approximationDateIndex = i;
			break;
		}
	}
	approximationDateIndex = i;

	long aDate;
	long lDate = m_recordingDates[m_recordingDates.getsize()-1]->GetDate(); 

	int nlimit;
	if ( approximationDateIndex == ndate ){
		aDate	= lDate;
		nlimit	= ndate-1;
	}
	else{
		aDate	= m_recordingDates[approximationDateIndex]->GetDate();
		nlimit	= approximationDateIndex-1;
	}

//  start pricing via summation

	int idateStart = 0;
	int idateEnd   = nlimit;
	bool calcAddOn = true;

	if ( int n = m_ntimeGauss*m_settlementDates.getsize() >= ndate  )
	{
//		you can as well do the entire summation explicitly without continuous approximation in this case
		idateEnd = ndate-1;
		approximationDateIndex = ndate;
	}

	CMatrix vals;
	value = calcExpectedPnl(vals,calcAddOn,m_gaussIntegrator,idateStart,idateEnd,hAsset,m_nstdev,m_eps);

	if ( approximationDateIndex == ndate ){
		values.resize(vals.rows()+1,vals.cols());

		values[0][0] = value;
		for ( int i = 0 ; i < vals.rows(); i++ ){
			for ( int j = 0 ; j < vals.cols() ; j++ ){
				values[i+1][j] = vals[i][j];
			}
		}
		return value;
	}

//  do the last one also

	vector < vector < double > > vvals;
	for ( int i = 0 ; i < vals.rows(); i++ ) {
		vector <double> vRow; 
		for ( int j = 0 ; j < vals.cols(); j++ ){
			vRow.push_back( vals[i][j] ) ; 
		}

		vvals.push_back(vRow); 
	}
	


	CVector	 lowerL(m_settlementDates.getsize());
	CVector	 upperL(m_settlementDates.getsize());


	vector <double> vRow(2); 

	double res;
	lowerL[0] = idateEnd+1;
	for ( int isettlement = 0 ; isettlement < m_settlementDates.getsize(); isettlement++ )
	{
		i = m_mapToRecordingDate[isettlement]-1;
		if ( i <= idateEnd ){
			continue;
		}

		res = calcExpectedPnl(vals,true,m_gaussIntegrator,i,i,hAsset,m_nstdev,m_eps);	
		value += res;

		vRow[0] = m_recordingDates[i]->GetDate(); 
		vRow[1] = res; 
		vvals.push_back(vRow); 


		upperL[isettlement] = i;

		if ( isettlement )
		{
			i = m_mapToRecordingDate[isettlement-1];
			res		= calcExpectedPnl(vals,true,m_gaussIntegrator,i,i,hAsset,m_nstdev,m_eps);	
			value	+= res;

			vRow[0] = m_recordingDates[i]->GetDate(); 
			vRow[1] = res; 
			vvals.push_back(vRow); 

			lowerL[isettlement] = i;
		}
	}



//  continue with continuous integral approximation to the summation:

	int ngauss		= m_ntimeGauss;
	CVector gaussPoint(ngauss),gaussWeight(ngauss);

	for (int isettlement = 0 ; isettlement < m_settlementDates.getsize(); isettlement++ )
	{
		int lowerDate	=	m_recordingDates[lowerL[isettlement]+1]->GetDate();
		int upperDate	=	m_recordingDates[upperL[isettlement]-1]->GetDate();
		double Mat		=	currentDateHandle->GetYearFraction(upperDate); 

		double mat		=	currentDateHandle->GetYearFraction(lowerDate);
		double lower	=	0.0;
		double upper	=	Mat-mat;

		if ( fabs(Mat-mat)<1e-3 ){
	//		there is just one point missing

			res = calcExpectedPnl(vals,true,m_gaussIntegrator,lowerL[isettlement],lowerL[isettlement],hAsset,m_nstdev,m_eps);
			value += res;


			for ( int i = 0 ; i < vals.rows(); i++ ) 
			{
				vRow[0] = vals[i][0]; 
				vRow[1] = vals[i][1]; 
				vvals.push_back(vRow);
			}


			continue;
		}

		MlEqMaths::dGauleg(lower,upper,gaussPoint,gaussWeight,ngauss,true);			
		
		int idate = 0 ;
		double t;
		GVector<long> mapToIdate(ngauss);
		long nDate;
		for (int igauss = 0 ; igauss < ngauss; igauss++ )
		{	
			t = gaussPoint[igauss];		
			nDate = lowerDate + (upperDate-lowerDate)*t/(Mat-mat);
		
			while ( m_recordingDates[idate]->GetDate() < nDate ) 
			{
				idate++;
				if ( idate >= ndate ){
					throw("what's going on horse?: trouble in setting up time integration");
				}
			}

			if ( idate == ndate-1 ){
				mapToIdate[igauss] = idate-1;
			}
			else{
				mapToIdate[igauss] = idate;		
			}
		}

		double deltaT;
		double integralValue = 0;
		for ( int igauss = 0 ; igauss < ngauss; igauss++ )
		{	
			idate = mapToIdate[igauss];
			if ( idate < ndate-1 ){
				deltaT	= m_recordingDates[idate]->GetYearFraction(m_recordingDates[idate+1]->GetDate());
			}
			else{
				deltaT	= m_recordingDates[idate-1]->GetYearFraction(m_recordingDates[idate]->GetDate());
			}

			res		= gaussWeight[igauss]*calcExpectedPnl(values,true,m_gaussIntegrator,idate,idate,hAsset,m_nstdev,m_eps)/deltaT;
			integralValue += res;

			vRow[0] = m_recordingDates[idate]->GetDate(); 
			vRow[1] = res; 
			vvals.push_back(vRow);

		}	

		value	+= integralValue;


	}

	values.resize(vvals.size()+1,2);
	values[0][0] = value;

	for ( int i = 0 ; i < vvals.size(); i++ ){
		values[i+1][0] = vvals[i][0];
		values[i+1][1] = vvals[i][1];
	}

	return value;
}		





double DeltaHedgedOptionPortfolio::calcExpectedPnl(CMatrix& values,bool calcAddOn,bool gaussIntegrator,int idateStart,int idateEnd,MlEqAssetHandle hAsset,double nstdev=4,double eps=1e-6)
{
	double result;
	if ( gaussIntegrator ){
		result = calcExpectedPnlGauss(values,calcAddOn,idateStart,idateEnd,hAsset,nstdev);
	}
	else{
		result = calcExpectedPnlRungeKutta(calcAddOn,idateStart,idateEnd,hAsset,nstdev,eps);
	}

	return result;

}




void integrateExpectedPnlRungeKutta(double x, double z[]  ,void *vp)
{

	deltaPortContainer* dC			= (deltaPortContainer*)vp;

	DeltaHedgedOptionPortfolio* dP = dC->vp;

	MlEqAssetHandle hAsset = dP->m_hAsset;
	double refSpot		   = 1.0;
	
	int idate			   = dC->idate;

	double spot	= x;
	double res = dP->calcPnlB(dC->calcAddOn,idate,hAsset,spot,refSpot);
	z[1] = res;
}




double DeltaHedgedOptionPortfolio::calcExpectedPnlRungeKutta(bool calcAddOn,int idateStart,int idateEnd,MlEqAssetHandle hAsset,double nstdev,double eps=1e-6)
{

	MlEqConstDateHandle	currentDateHandle	= hAsset->GetDateHandle();	

	deltaPortContainer dc;
	dc.vp			= this;
	dc.nMaturity	= m_nMaturity;
	dc.voleps		= 0.01;
	dc.spoteps		= 0.01;

	double hInit=0.1,hMin=1e-2,bump=0.0008,scale = 10;

	CVector integrationLimits(3);
	double value = 0.0;
	double integral;
	for ( int idate = idateStart; idate <= idateEnd; idate++ )
	{

		dc.idate				= idate;

		long aDate				= m_recordingDates[idate]->GetDate();
		double mat				= currentDateHandle->GetYearFraction(aDate);
		double fwd				= hAsset->GetForward(aDate, false);

		MlEqStrike	stk(fwd);
		double volref			= hAsset->GetVolatility(stk,aDate);

		if ( m_useControlVariate ){
			volref = MlEqMaths::Max(volref,m_bookingVols[0]);
		}

		integrationLimits[0]	= fwd*exp(-nstdev*volref*mat);
		integrationLimits[1]	= fwd;
		integrationLimits[2]	= fwd*exp(nstdev*volref*mat);
		

		double addOnVal = calcPnlB(calcAddOn,idate,hAsset,0.0,0.0);

		if ( calcAddOn ){
			dc.calcAddOn = false;
		}

		Runge_Kutta_Integrate_new(integral,integrateExpectedPnlRungeKutta , integrationLimits, hMin, eps, bump, &dc, hInit, scale);

		value += integral+addOnVal;
	}
	return value;
}







double DeltaHedgedOptionPortfolio::calcExpectedPnlGauss(CMatrix& values,bool calcAddOn,int idateStart,int idateEnd,MlEqAssetHandle hAsset,double nstdev=4)
{

	MlEqConstDateHandle	currentDateHandle	= hAsset->GetDateHandle();	

	values.resize(idateEnd-idateStart+1,2);

	double fwd,mat,volref,lower,upper;
	long cDate;

	double val;
	double value = 0.0;
	for ( int idate = idateStart; idate <= idateEnd; idate++ )
	{
		m_gaussPoint.resize(m_ngauss);
		m_gaussWeight.resize(m_ngauss);
		
		cDate	= m_recordingDates[idate]->GetDate();
		fwd		= hAsset->GetForward(cDate, false);
		mat		= currentDateHandle->GetYearFraction(cDate); 
		volref	= hAsset->GetVolatility(MlEqStrike(fwd),cDate);

		values[idate-idateStart][0] = cDate;

		if ( m_useControlVariate ){
			volref = MlEqMaths::Max(volref,m_bookingVols[0]);
		}
		
		lower	= 1e-10;
		upper   = 1.0;
		
		MlEqMaths::dGauleg(lower,upper,m_gaussPoint,m_gaussWeight,m_ngauss,true);			
		
		for ( int i = 0 ; i < m_gaussWeight.getsize(); i++ ){
			m_gaussWeight[i] *= fwd;
		}
		
		val		=	calcExpectedPnlGauss(calcAddOn,idate,hAsset,fwd);
		
		lower   = 1.0;
		upper	= exp( nstdev*volref*sqrt(mat));
		MlEqMaths::dGauleg(lower,upper,m_gaussPoint,m_gaussWeight,m_ngauss,true);			
		
		for ( int i = 0 ; i < m_gaussWeight.getsize(); i++ ){
			m_gaussWeight[i] *= fwd;
		}
		
		val		+=	calcExpectedPnlGauss(false,idate,hAsset,fwd);
		values[idate-idateStart][1] = val;

		value	+=	val;

	}	
	
	return value;
}	




double DeltaHedgedOptionPortfolio::calcExpectedPnlGauss(bool calcAddOn,int idate,MlEqAssetHandle hAsset,double refSpot)
{

	double result = calcPnlC(calcAddOn,idate,hAsset,refSpot);
	return result;
}






double DeltaHedgedOptionPortfolio::calcPortfolioB(double& deltaTerm,bool calcAddOn,int idate,MlEqAssetHandle hAsset,double spot,double refSpot)
{

	MlEqConstDateHandle		hDate   = hAsset->GetDateHandle();
	MlEqDate	hiDate(m_recordingDates[idate]->GetDate(),hDate->GetDayCountConvention());

	double   mat					= hiDate.GetYearFraction(m_nMaturity);
	long currentDate				= hiDate.GetDate();

	double	 mktFwd					= hAsset->GetForward(currentDate, false);

	double payoffFwdFac				= exp((m_bookingRate-m_bookingDiv)*mat);
	double payoffFwd				= spot*payoffFwdFac;
	double payoffDf					= exp(-m_bookingRate*mat);

//	double	 TestMktFwd					= hAsset->GetForward(m_nMaturity, false)/hAsset->GetForward(m_recordingDates[idate]->GetDate(),false);

	double asOfmat					= hDate->GetYearFraction(currentDate);

	MlEqStrike strike(spot);
	double vol = hAsset->GetVolatility(strike,currentDate);

	double bsdensity,d1;

	double volsqrt;
	double result = 0.0;
	double value;

	deltaTerm = 0.0;

	if ( calcAddOn )		// this is the first term before the integral...
	{

		double addonVal;
		for ( int i = 0 ; i < m_portfolio.rows(); i++ )
		{
			addonVal	=	Bs(mktFwd*payoffFwdFac,m_bookingVols[i],mat,m_portfolio[i][1],payoffDf,m_portfolio[i][2]);
			result		+=	m_portfolio[i][0]*addonVal;

			if ( fabs( mat) < 1e-3 )	
			{
				spot = m_portfolio[i][1];
				vol = hAsset->GetVolatility(MlEqStrike(spot),currentDate);

				if ( spot < mktFwd ){
					result += Bs(mktFwd,vol,asOfmat,spot,1.0,-1);

					if( m_portfolio[i][2] == 1 ) // call case
						deltaTerm += payoffFwdFac*mktFwd;
				}
				else{
					result += Bs(mktFwd,vol,asOfmat,spot,1.0,1);

					if( m_portfolio[i][2] == -1 ) // put case
						deltaTerm += - payoffFwdFac*mktFwd;
				}
			}
			else
			{
				double volsqrt = sqrt(mat)*m_bookingVols[i];
				double d1 = (log(mktFwd*payoffFwdFac/m_portfolio[i][1])+0.5*volsqrt*volsqrt)/volsqrt;
			
				if( m_portfolio[i][2] == 1 ) // call case
					deltaTerm += normal(d1)*payoffFwdFac*mktFwd;
				else	// put case
					deltaTerm += (normal(d1) - 1.)*payoffFwdFac*mktFwd;
			}
		}

		return result;
	}

	if ( fabs( mat < 1e-3 )){
		return 0.0;
	}

	double deltaDensity,put,call;
	for ( int i = 0 ; i < m_portfolio.rows(); i++ )
	{

		volsqrt = m_bookingVols[i]*sqrt(mat);
		d1 = (log(payoffFwd/m_portfolio[i][1])+0.5*volsqrt*volsqrt)/volsqrt;

		bsdensity = m_sq*payoffFwdFac*exp(-0.5*d1*d1)/(spot*volsqrt);

		deltaDensity= -bsdensity/spot*(1.0-d1/volsqrt)+2.0*bsdensity;

		if ( spot < mktFwd )
		{
			put			= Bs(mktFwd,vol,asOfmat,spot,1.0,-1);
			value		= bsdensity*put;
			deltaTerm	+= deltaDensity*put;
		}
		else
		{
			call		= Bs(mktFwd,vol,asOfmat,spot,1.0,1);
			value		= bsdensity*call;
			deltaTerm	+= deltaDensity*call;
		}

		value *= m_portfolio[i][0];
		result += value;
	}


//	result *= refSpot;
//	deltaTerm* = refSpot;

	return result;
}


double DeltaHedgedOptionPortfolio::calcPnlB(bool calcAddOn,int idate,MlEqAssetHandle hAsset,double spot,double refSpot)
{


	double result = 0.0;
	double deltaT,df;

	const MlEqZeroCurveHandle	zc = hAsset->GetPayZeroCurve(true) ;

	double payoffFwdFactor,mktFwdFactor;

	bool j=false;
	if ( idate < m_nDates-1 )
	{
		if ( m_mapToSettlement[idate+1] == m_mapToSettlement[idate] )
		{
			j		= true;
			deltaT	= m_recordingDates[idate]->GetYearFraction(m_recordingDates[idate+1]->GetDate());
			df		= exp(m_bookingRate*deltaT);

			payoffFwdFactor		= exp((m_bookingRate-m_bookingDiv)*deltaT);
			mktFwdFactor		= hAsset->GetForward(m_recordingDates[idate+1]->GetDate(),false)/hAsset->GetForward(m_recordingDates[idate]->GetDate(),false);

/*df = zc->GetDiscountFactor(m_settlementDates[m_mapToSettlement[idate+1]])/
		zc->GetDiscountFactor(m_settlementDates[m_mapToSettlement[idate]]);

df = 1.0/df;
*/
		}
	}
	else{
		df = 1.0;// is this corrdct?sos
	}

	double deltaTerm;
	double Val = calcPortfolioB(deltaTerm,calcAddOn,idate,hAsset,spot,refSpot);

	if ( idate < m_nDates-1 ){
		deltaTerm = 0.0;
	}
	else{
		deltaTerm *= (mktFwdFactor-payoffFwdFactor);
	}



//	if ( idate == 0 )

	if ( idate == m_nDates-1 || m_mapToSettlement[idate+1] > m_mapToSettlement[idate] )
	{
		// this is the last term
		//  do nothing
		int i=1;
	}
	else if ( idate == 0 || m_mapToSettlement[idate] > m_mapToSettlement[idate-1] )
	{
// this is the first term in the series
		Val *= -df;
	}
	else 
	{// this is the term in the middle
			Val *= (1.0-df);
	}

	result = Val;
		
//  this gets paid at the next setting Date

//	const MlEqZeroCurveHandle	zc = hAsset->GetPayZeroCurve(true) ;
	double disc = zc->GetDiscountFactor(m_settlementDates[m_mapToSettlement[idate]]);
	result *= disc;

	return result;
}



double DeltaHedgedOptionPortfolio::calcPnlC(bool calcAddOn,int idate,MlEqAssetHandle hAsset,double refSpot)
{

	double spot,result,deltaT,df,payoffFwdFactor=0.0,mktFwdFactor=0.0,deltaTerm,addOnVal,bsdensity;
	double d1,volsqrt,value,deltaDensity,put,call,addOnDeltaTerm,vol,payoffFwd;

	const MlEqZeroCurveHandle	zc = hAsset->GetPayZeroCurve(true) ;

	MlEqConstDateHandle		hDate   = hAsset->GetDateHandle();
	MlEqDate				hiDate(m_recordingDates[idate]->GetDate(),hDate->GetDayCountConvention());
	

	long currentDate				= hiDate.GetDate();
	double	 mktFwd					= hAsset->GetForward(currentDate, false);

	double asOfmat					= hDate->GetYearFraction(currentDate);

	if ( idate < m_nDates-1 )
	{
		if ( m_mapToSettlement[idate+1] == m_mapToSettlement[idate] )
		{
			deltaT	= m_recordingDates[idate]->GetYearFraction(m_recordingDates[idate+1]->GetDate());
			df		= exp(m_bookingRate*deltaT);

			payoffFwdFactor		= exp((m_bookingRate-m_bookingDiv)*deltaT);
			mktFwdFactor		= hAsset->GetForward(m_recordingDates[idate+1]->GetDate(),false)/hAsset->GetForward(m_recordingDates[idate]->GetDate(),false);

		}
	}
	else{
			df = 1.0;
	}

// calculate constant terms first
	
	
	addOnVal = 0.0;
	addOnDeltaTerm = 0.0;

	if ( calcAddOn )
	{
		for ( int i = 0 ; i < m_portfolio.rows(); i++ )
		{	
			double   mat					= hiDate.GetYearFraction(m_portfolio[i][3]);
			double payoffFwdFac				= exp((m_bookingRate-m_bookingDiv)*mat);
			double payoffDf					= exp(-m_bookingRate*mat);

			addOnVal		+=	m_portfolio[i][0]*Bs(mktFwd*payoffFwdFac,m_bookingVols[i],mat,m_portfolio[i][1],payoffDf,m_portfolio[i][2]);
			if ( fabs( mat) < 1e-3 )	
			{
				double bookingStrike = m_portfolio[i][1];

				vol = hAsset->GetVolatility(MlEqStrike(bookingStrike),currentDate);
				payoffFwd	= bookingStrike*payoffFwdFac;

				if ( bookingStrike < mktFwd )
				{
						double vanilla = Bs(mktFwd,vol,asOfmat,bookingStrike,1.0,-1);
						addOnVal += vanilla ;

						addOnDeltaTerm += 0.5 * vanilla ;

						if( m_portfolio[i][2] == 1 )
						{ // call case
							addOnDeltaTerm += mktFwd;
						}
				}
				else
				{
						double vanilla = Bs(mktFwd,vol,asOfmat,bookingStrike,1.0,1);
						addOnVal += vanilla ;

						addOnDeltaTerm += 0.5*vanilla ;

						if( m_portfolio[i][2] == -1 )
						{ // put case
							addOnDeltaTerm += - mktFwd;
						}
				}
			}
			else
			{
				double volsqrt = sqrt(mat)*m_bookingVols[i];
				double d1 = (log(mktFwd*payoffFwdFac/m_portfolio[i][1])+0.5*volsqrt*volsqrt)/volsqrt;
				
				if( m_portfolio[i][2] == 1 ){ // call case
					addOnDeltaTerm += exp(-m_bookingRate*mat)*normal(d1)*payoffFwdFac*mktFwd;
				}
				else{	// put case
					addOnDeltaTerm += exp(-m_bookingRate*mat)*(normal(d1) - 1.0)*payoffFwdFac*mktFwd;
				}
			}
		}	
	}
	


//  start integration


	double deltaTermFactor;
	if ( idate == m_nDates-1 ){
		deltaTermFactor = 0.0;
	}
	else{
		deltaTermFactor = -m_includeDeltaTerm*(mktFwdFactor-payoffFwdFactor);

	}

	result = 0.0;
	for ( int igauss = 0 ; igauss < m_ngauss; igauss++ )
	{	
		
		spot = m_gaussPoint[igauss]*refSpot;

		MlEqStrike strike(spot);
		double vol = hAsset->GetVolatility(strike,currentDate);

		for ( int i = 0 ; i < m_portfolio.rows(); i++ )
		{

			double   mat					= hiDate.GetYearFraction(m_portfolio[i][3]);

			if ( fabs( mat < 1e-3 )){
				continue;
			}


			double payoffFwd			= spot*exp((m_bookingRate-m_bookingDiv)*mat);
			double payoffFwdFac				= exp((m_bookingRate-m_bookingDiv)*mat);
			double payoffDf					= exp(-m_bookingRate*mat);

			volsqrt = m_bookingVols[i]*sqrt(mat);
			d1 = (log(payoffFwd/m_portfolio[i][1])+0.5*volsqrt*volsqrt)/volsqrt;
			bsdensity = exp(-m_bookingRate*mat)*m_sq*payoffFwdFac*exp(-0.5*d1*d1)/(spot*volsqrt);
			deltaDensity =	bsdensity*(1.0-d1/volsqrt);

			if ( spot < mktFwd )
			{
				put			=	Bs(mktFwd,vol,asOfmat,spot,1.0,-1);
				value		=	bsdensity*put;
				deltaTerm	=	deltaDensity*put;
			}
			else
			{
				call		=	Bs(mktFwd,vol,asOfmat,spot,1.0,1);
				value		=	bsdensity*call;
				deltaTerm	=	deltaDensity*call;
			}

			value		*= m_portfolio[i][0];
			deltaTerm	*= m_portfolio[i][0];

			result += (value+deltaTerm*deltaTermFactor)*m_gaussWeight[igauss];
		}

	}

	result += addOnVal+addOnDeltaTerm*deltaTermFactor;

	if ( idate == m_nDates-1 || m_mapToSettlement[idate+1] > m_mapToSettlement[idate] )
	{
	//		this is the last term
	//		do nothing
			int i=1;
	}
	else if ( idate == 0 || m_mapToSettlement[idate] > m_mapToSettlement[idate-1] )
	{
	//		this is the first term in the series
			result *= -df;
	}
	else 
	{
	//		this is the term in the middle
			result *= (1.0-df);
	}


//  this gets paid at the next setting Date

//	const MlEqZeroCurveHandle	zc = hAsset->GetPayZeroCurve(true) ;
	double disc = zc->GetDiscountFactor(m_settlementDates[m_mapToSettlement[idate]]);
	result *= disc;

	return result;
}




double DeltaHedgedOptionPortfolio::calcPnl(int idate,MlEqAssetHandle hAsset,double spot,double refSpot)
{
	double result = 0.0;
	double nextVal=0;
	double deltaT,df;

	bool j=false;
	if ( idate < m_nDates-1 )
	{
		if ( m_mapToSettlement[idate+1] == m_mapToSettlement[idate] )
		{
			j		= true;
			deltaT	= m_recordingDates[idate]->GetYearFraction(m_recordingDates[idate+1]->GetDate());
			df		= exp(m_bookingRate*deltaT);
		}
	}
	else{
		df = 1.0;
	}

	double Val = calcPortfolio(idate,hAsset,spot);
	Val *= refSpot;
	if ( j ){
		if ( m_mapToSettlement[idate+1] == m_mapToSettlement[idate] ){
			nextVal = calcPortfolio(idate+1,hAsset,spot);
			nextVal *= refSpot;
		}
	}
	else{
		nextVal = 0.0;
	}

	result = (nextVal - df*Val);
		
//  this gets paid at the next setting Date

	const MlEqZeroCurveHandle	zc = hAsset->GetPayZeroCurve(true) ;
	double disc = zc->GetDiscountFactor(m_settlementDates[m_mapToSettlement[idate]]);
	result *= disc;

	return result;
}


double DeltaHedgedOptionPortfolio::calcPortfolio(int idate,MlEqAssetHandle hAsset,double spot)
{

	double   voleps					= 0.01;
	double   spoteps				= 0.01;

	MlEqConstDateHandle		hDate   = hAsset->GetDateHandle();
//	long	 valDate				= hAsset->GetDateHandle()->GetDate();

	double   mat					= m_recordingDates[idate]->GetYearFraction(m_nMaturity);
	double	 mktFwd					= hAsset->GetForward(m_recordingDates[idate]->GetDate(), false);
	double	 mktSpot				= hAsset->GetSpot(hAsset->GetDateHandle()->GetDate());

	double payoffFwd				= spot*exp((m_bookingRate-m_bookingDiv)*mat);
	double payoffDf					= exp(-m_bookingRate*mat);

	double asOfmat					= hDate->GetYearFraction(m_recordingDates[idate]->GetDate());
	double shortPayoffFwd			= mktSpot*exp((m_bookingRate-m_bookingDiv)*asOfmat);

	double value = 0.0;
	for ( int i = 0 ; i < m_portfolio.rows(); i++ )
	{
		value += m_portfolio[i][0]*
				 Bs(payoffFwd,m_bookingVols[i],mat,m_portfolio[i][1],payoffDf,m_portfolio[i][2]);
	}

	
	double bsdensity=0.0;

	long currentDate = m_recordingDates[idate]->GetDate();
	double x = spot;
	MlEqStrike strike(x);
	double vol = hAsset->GetVolatility(strike,currentDate);

	MlEqStrike strike_up(x*(1.0+voleps));
	double vol_up = hAsset->GetVolatility(strike_up,currentDate);

	MlEqStrike strike_do(x*(1.0-voleps));
	double vol_do = hAsset->GetVolatility(strike_do,currentDate);


	double density = 


	( Bs(mktFwd,vol_up,asOfmat,x*(1.0+voleps),1.0,1) +
	  Bs(mktFwd,vol_do,asOfmat,x*(1.0-voleps),1.0,1) -
	  2.0*Bs(mktFwd,vol,asOfmat,x,1.0,1)
	  );

	density /=  pow(x*voleps,2.0);

	double bsVol = m_bookingVols[0];

	if ( m_useControlVariate )
	{
		 bsdensity = 

		( Bs(shortPayoffFwd,bsVol,asOfmat,x*(1.0+voleps),1.0,1) +
		  Bs(shortPayoffFwd,bsVol,asOfmat,x*(1.0-voleps),1.0,1) -
		  2.0*Bs(shortPayoffFwd,bsVol,asOfmat,x,1.0,1)
		  );

		bsdensity /=  pow(x*voleps,2.0);
	}

	value = (density-bsdensity)*value;

	return value;

}


