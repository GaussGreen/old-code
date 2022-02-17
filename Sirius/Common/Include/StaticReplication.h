
#pragma once
#include "mleqobjects.h"


class CStaticReplication
{
public:
	CStaticReplication( MlEqAssetHandle hUdly,
						DATE dateMaturity,
						int nPoints = 30);

	CStaticReplication( MlEqAssetHandle hUdly,
						DATE dateMaturity,
						double splitLimit,
						int nPoints = 30);


	virtual ~CStaticReplication(){};

	virtual double fComputePrice();

protected:							// f and df should be at least continuous...to check
	virtual double f(double x) = 0;	
	virtual double df(double x) {return 0.0;}
	virtual double ddf(double x){return 0.0;}

protected:
	MlEqAssetHandle		m_hUnderlying;
	long				m_dateMaturity;

	double				m_splitLimit;
	CVector				m_integrationLimits;
	void				initializeIntegrator( int nPoints );
	double				integrate( double down, double up);
	double				integrate( double finite );	// integrate from finite to infinity...
	bool				m_infiniteIntegrate;

	CVector				m_gaussWeights;
	CVector				m_gaussPoints;
	int					m_nPoints;

	virtual double		vanillaPrice(double strike, double cp);
	double				m_mat;	
	MlEqVolatilityStructureHandle	m_hVolatilityStructure;
	double				m_forward;
	double				m_spot;
	long				m_nToday;
};


class testPayoff	:	public CStaticReplication
{
public:
	testPayoff(MlEqAssetHandle hUdly, DATE dateMaturity, int nPoints)
		:CStaticReplication(hUdly, dateMaturity, nPoints){};

	~testPayoff(){};

	double f(double x);
	double df(double x);
	double ddf(double x);
};

