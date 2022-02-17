
#pragma once


#include "resource.h"       // main symbols
#include "mleqobjects.h"



class MlEqMonteCarlo;
class CMatrix;

class IguanaPricer : public product
{
public:
	double m_Strike;
	double m_Barrier;

	void initialize(double Strike,
					double Barrier,
					GVector< long >& payoffDates,
					CVector& spotFixings);

	void payout(CMatrix& values,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& mc);
	void setUp(CMatrix& value,MlEqMonteCarlo& mc);


protected:
	double m_pastMin;
	double m_pastMax;
	int m_start;
	int m_scheduleSize;
};

