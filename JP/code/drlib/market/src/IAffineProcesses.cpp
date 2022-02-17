#include "edginc/config.hpp"
#include "edginc/IAffineProcesses.hpp"
#include "edginc/Algorithm.hpp"

DRLIB_BEGIN_NAMESPACE
 

static void IAffineProcessload(CClassSP& clazz) {
	REGISTER_INTERFACE(IAffineProcess, clazz);
	EXTENDS(IObject);
}

CClassConstSP const IAffineProcess::TYPE = 
CClass::registerInterfaceLoadMethod("IAffineProcess", 
									typeid(IAffineProcess), 
									IAffineProcessload);


void ricatti(IAffineProcess *ricattiEq,
			 int index, 
			 double bigT, 
			 double t,
			 double &alpha, 
			 double &beta,
			 double h)
{
	if (t>=bigT) return; 

	double k1a,k2a,k3a,k4a, k1b,k2b,k3b,k4b,hh;
	double tn=bigT;
	while (tn>t+1e-10)
	{
		hh = min(tn-t,h);
		ricattiEq->ricattiVF(index, tn,	     beta,              k1a, k1b);
		ricattiEq->ricattiVF(index, tn-0.5*hh,  beta - 0.5*hh*k1b, k2a, k2b);
		ricattiEq->ricattiVF(index, tn-0.5*hh,  beta - 0.5*hh*k2b, k3a, k3b);
		ricattiEq->ricattiVF(index, tn-hh,	     beta - hh*k3b,	    k4a, k4b);
		alpha -= 0.1666666666666*hh*(k1a+2*k2a+2*k3a+k4a);
		beta  -= 0.1666666666666*hh*(k1b+2*k2b+2*k3b+k4b);
		tn -= hh;
	}
}

double minusLogSurvProbability(IAffineProcess *ricattiEq,
						 double bigT, 
						 double t, 
						 double lambdat)
{
	DoubleArraySP pieceT = ricattiEq->getPiecewiseTimes();
	QLIB_VERIFY(t<=(*pieceT)[0], "t in the past");
	double h = ricattiEq->getRungeKuttaDiscretization();
	if (t>=bigT) return 0.0;

	double alpha = 0.0, beta = 0.0;

	int index = lower_bound(pieceT->begin(),pieceT->end(),bigT) - pieceT->begin();  // pieceT[index-1] < bigT <= pieceT[index]
	if (index==pieceT->size()) --index;

	double t1, t2 = bigT;
	while ( (index>0) && ( (*pieceT)[index-1]>t ) )
	{
		t1 = (*pieceT)[index-1];
		ricatti(ricattiEq,index, t2, t1, alpha, beta,h);
		t2 = t1;
		--index;
	}
	ricatti(ricattiEq,index, t2, t, alpha, beta,h);
	return -alpha-lambdat*beta;
}


static void IJumpDistributionload(CClassSP& clazz) {
	REGISTER_INTERFACE(IJumpDistribution, clazz);
	EXTENDS(IObject);
}

CClassConstSP const IJumpDistribution::TYPE = 
CClass::registerInterfaceLoadMethod("IJumpDistribution", 
									typeid(IJumpDistribution), 
									IJumpDistributionload);



DRLIB_END_NAMESPACE
