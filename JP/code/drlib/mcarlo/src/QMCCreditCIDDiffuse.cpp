//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : QMCCreditCIDDiffuse.cpp
//
//   Description : Extension CIR credit path generation
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/QMCCreditCIDDiffuse.hpp"

DRLIB_BEGIN_NAMESPACE


/** constructor */
QMCCreditCIDDiffuse::QMCCreditCIDDiffuse(IQMCDiffusibleInterestRateSP srmRatesDiffuse, bool isFullMC) 
: SRMCreditCIRDiffuse(srmRatesDiffuse, isFullMC) {}

/** destructor */
QMCCreditCIDDiffuse::~QMCCreditCIDDiffuse(){}

void QMCCreditCIDDiffuse::setQMCCreditCIDDiffuse(
    int                    _randomIndex,
    const DateTime&        _today,
	const vector<double>&  _prob,
	double				   _sigma,
	double				   _meanReversionSpeed,
	double				   _initialValue,
	const DateTimeArray&   _datesForTheta,
	const vector<double>&  _theta)
{
    static const string method("QMCCreditCIDDiffuse::setQMCCreditCIDDiffuse");

    randomIndex=_randomIndex;
    crFxCorr=0.0;
    today=_today;
    probStart=-1; // FIXME probStart : see QMCRatesDiffuse how to calc it
    prob=_prob;

	setModelParametersExceptTheta(_sigma,0,_meanReversionSpeed);
	vector<double> SRMtheta(_theta);  // SRM CIR USED A DIFFERENT DEFINITION OF THETA THAN CID !!!!!!!!!
	for (size_t i=0; i<_theta.size(); i++) SRMtheta[i] *= _meanReversionSpeed;
	setTheta(_datesForTheta, _initialValue, SRMtheta);
}


void QMCCreditCIDDiffuse::setQMCCreditCIDDiffuse(
    int                    _randomIndex,
    const DateTime&        _today,
	const vector<double>&  _prob,
	double				   _sigma,
	double				   _meanReversionSpeed,
	const DateTimeArray&   _datesSurvProb,
	const vector<double>&  _survProba)
{
    static const string method("QMCCreditCIDDiffuse::setQMCCreditCIDDiffuse");

    randomIndex=_randomIndex;
    crFxCorr=0.0;
    today=_today;
    probStart=-1; // FIXME probStart : see QMCRatesDiffuse how to calc it
    prob=_prob;

	setModelParametersExceptTheta(_sigma,0,_meanReversionSpeed);
	calibrateTheta(_datesSurvProb, _survProba);
}



DRLIB_END_NAMESPACE

