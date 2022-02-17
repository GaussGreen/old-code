//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : QMCCreditCIDDiffuse.hpp
//
//   Description : CIR credit path generation
//
//
//----------------------------------------------------------------------------

#ifndef QMCCreditCIDDiffuse_HPP
#define QMCCreditCIDDiffuse_HPP

#include "edginc/SRMCreditCIRDiffuse.hpp"

DRLIB_BEGIN_NAMESPACE

/** class for CIR credit path generation modified for CID */
class QMCCreditCIDDiffuse: public SRMCreditCIRDiffuse {
public:

    /** constructor */
    QMCCreditCIDDiffuse(IQMCDiffusibleInterestRateSP srmRatesDiffuse, bool isFullMC);

    /** destructor */
    virtual ~QMCCreditCIDDiffuse();

    /** initialization */
    void setQMCCreditCIDDiffuse(
            int                    _randomIndex,
            const DateTime&        _today,
			const vector<double>&  _prob,
	        double				   _sigma,
	        double				   _meanReversionSpeed,
	        double				   _initialValue,
	        const DateTimeArray&   _datesForTheta,
	        const vector<double>&  _theta);
		

    void setQMCCreditCIDDiffuse(
		    int                    _randomIndex,
		    const DateTime&        _today,
			const vector<double>&  _prob,
			double				   _sigma,
			double				   _meanReversionSpeed,
			const DateTimeArray&   _datesSurvProb,
			const vector<double>&  _survProba);

};

/** declare smartPtr versions */
DECLARE(QMCCreditCIDDiffuse);

DRLIB_END_NAMESPACE

#endif // QMCCreditCIDDiffuse_HPP
