

#include "stdafx.h"
#include "StaticReplication.h"


CStaticReplication::CStaticReplication( MlEqAssetHandle hUdly,
										DATE dateMaturity,
										int nPoints):
m_hUnderlying(hUdly)
{	
	m_nToday				= m_hUnderlying->GetDateHandle()->GetDate();
	m_dateMaturity			= dateMaturity;
	MlEqDateHandle	hDate	= new MlEqDate(*m_hUnderlying->GetDateHandle());
	m_mat					= hDate->GetYearFraction(m_dateMaturity) ;
	m_forward				= m_hUnderlying->GetQuantoForward(m_nToday, m_dateMaturity, false);	
	m_spot					= m_hUnderlying->GetSpot(m_nToday);		
	m_hVolatilityStructure	= m_hUnderlying->GetVolatilityStructure();


	initializeIntegrator(nPoints);

	m_splitLimit = m_forward;

	m_integrationLimits = CVector(3, 0.0);
	m_integrationLimits[0] = 0.0;
	m_integrationLimits[1] = m_forward;
	m_integrationLimits[2] = 2.* m_forward;

	m_infiniteIntegrate = true;	
}

CStaticReplication::CStaticReplication( MlEqAssetHandle hUdly,
										DATE dateMaturity,
										double splitLimit,
										int nPoints):
m_hUnderlying(hUdly)
{	
	m_nToday				= m_hUnderlying->GetDateHandle()->GetDate();
	m_dateMaturity			= dateMaturity;
	MlEqDateHandle	hDate	= new MlEqDate(*m_hUnderlying->GetDateHandle());
	m_mat					= hDate->GetYearFraction(m_dateMaturity) ;
	m_forward				= m_hUnderlying->GetQuantoForward(m_nToday, m_dateMaturity, false);	
	m_hVolatilityStructure	= m_hUnderlying->GetVolatilityStructure();


	initializeIntegrator(nPoints);

	m_splitLimit = splitLimit;

	m_integrationLimits = CVector(3, 0.0);
	m_integrationLimits[0] = 0.0;
	m_integrationLimits[1] = m_splitLimit;
	m_integrationLimits[2] = 2.* m_forward;

	m_infiniteIntegrate = true;	
}

void CStaticReplication::initializeIntegrator( int nPoints )
{		
	m_nPoints		= nPoints;
	m_gaussWeights	= CVector(m_nPoints, 0.0);
	m_gaussPoints	= CVector(m_nPoints, 0.0);

	MlEqMaths::dGauleg(0.0, 1.0, m_gaussPoints, m_gaussWeights, m_nPoints, false);	// true ?
}

double CStaticReplication::vanillaPrice( double strike, double cp)
{
	double vol  = m_hVolatilityStructure->getFutureVol(MlEqStrike(strike), m_nToday, m_dateMaturity, m_spot);
	vol = m_hUnderlying->GetCompositeVolatility(long(m_nToday), double(m_dateMaturity), vol);

	return ::Bs(m_forward, vol, m_mat, strike, 1.0, cp);	// non discounted vanilla price...
}

double CStaticReplication::integrate( double down, double up)
{
	double dx = up - down ;
	double cp = (up <= m_splitLimit)? -1.0:1.0 ;

	double res = 0.0;
	for ( int i = 0 ; i < m_nPoints; i++ )
	{
		double strike = down + dx * m_gaussPoints[i];
		double eval = vanillaPrice(strike, cp);
		eval *= ddf( strike ); 
	
		res += m_gaussWeights[i] * eval ;		
	}
	res *= dx;
	return res ;
}


double CStaticReplication::integrate( double finite )
{
	double res = 0.0;
	for ( int i = 0 ; i < m_nPoints; i++ )
	{
		double u = m_gaussPoints[i];
		double strike = finite + 1./(u*u) -1.  ;// change of variables if upper integration limit is infinite: u : [0,1]
		double du = 2./(u*u*u) ;
		double eval = vanillaPrice(strike, 1.0);
		eval *= ddf( strike ); 
	
		res += m_gaussWeights[i] * eval * du ;		
	}
	return res ;
}


double CStaticReplication::fComputePrice()	// should eventually be discounted
{
// start with the "principal value"
	double value = f( m_splitLimit ); 

// second term
	value += df( m_splitLimit ) * (m_forward - m_splitLimit);

// add the integral parts...
	int nslice = m_integrationLimits.getsize();

	for(int i=0; i<nslice-1; i++){
		value += integrate( m_integrationLimits[i], m_integrationLimits[i+1] );
	}
	
	if( m_infiniteIntegrate )
		value += integrate( m_integrationLimits[nslice-1] );

	return value;
}


double testPayoff::f(double x)
{ 
	return -log(x/m_splitLimit);
//	return 1/x;
}


double testPayoff::df(double x)
{
	return 1/x ;
//	return 2/(x*x);
}

double testPayoff::ddf(double x)
{
	return 1/(x*x);
//	return 2/(x*x*x);
}

  