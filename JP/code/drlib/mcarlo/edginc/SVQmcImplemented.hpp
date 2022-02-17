//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : SVQmcImplemented.hpp
//
//   Description : QMC implementation of State variables for accessing information
//                 in various IQMCDiffusibleAsset objects
//
//----------------------------------------------------------------------------

#ifndef SVQmcImplemented_HPP
#define SVQmcImplemented_HPP

#include "edginc/IQMCStateVariableBase.hpp"
#include "edginc/IQMCDiffusibleAsset.hpp"

#include "edginc/SVAllInterfaces.hpp"

#include "edginc/QMCHelperBoundedDiffusion.hpp"
#include <cstdio> // FIXME remove debug printings
#include <cassert>

DRLIB_BEGIN_NAMESPACE


/*******************************************************************
 * Asset specific state variables with some function names aliased *
 *******************************************************************/

// note: there is a bunch of other inherited functions in every class here -- like
// advance();  -- for stateless pricing
// reset();    -- for stateless pricing
// path();     -- accessing SVPath structure, etc
//
//  they are not to be overridden or specified in a regular case. The classes
//  below add one or two additional asset-type-specific methods of information
//  retrieval and also define a correct function for calculating a value
//
// element(int i);  -- accessing proper element


/** Interest-rate specific function names realized in these state variables */

class SVQmcDiscFactor : public IQMCStateVariableSpot, public SVDiscFactor
{
public:
    virtual ~SVQmcDiscFactor() {}

    SVQmcDiscFactor(
        IQMCDiffusibleInterestRate*  _pAsset,
        const DateTimeArray&        _measurementDates);

    /** asset specific name of one of the generic functions */
    virtual double getSpotDF() const {
    	return getValue();
    }

    /** asset specific name of one of the generic functions */
    virtual double firstDF() const {
    		return (*this)[0];
    }

    virtual double getDF(int i) const {
        return (*this)[i];
    }

    virtual double element(int dtIndex) const
    {
        return pAsset->getDiscFactor(measurementDatesIdx[dtIndex]);
    }

    virtual void prepare(bool mm);

protected:

    SVQmcDiscFactor();
    IQMCDiffusibleInterestRate* pAsset;

};

//DECLARE(SVQmcDiscFactor);


class SVQmcDiscFactorMM : public SVQmcDiscFactor
{
public:
    virtual ~SVQmcDiscFactorMM() {}

    SVQmcDiscFactorMM(
        IQMCDiffusibleInterestRate*  _pAsset,
        const DateTimeArray&        _measurementDates) : SVQmcDiscFactor(_pAsset, _measurementDates), nAcc(0) {}

    // these are the only functions overridden in the class compared to its immediate parent
    virtual double element(int dtIndex) const
    {
        return pAsset->getDiscFactor(measurementDatesIdx[dtIndex]) * correction.at(dtIndex);
    }
    virtual void prepare(bool mm); // mm is redundant here -- if this class was initialized it is already MM
    void resetMMCorrection();
    void accumulateMMCorrection();
    void setMMCorrection();

private:

    SVQmcDiscFactorMM();

    // for Moment Matching
    int nAcc;
    vector<double> correction;  // initialized in prepare()
    vector<double> accumulation;
    vector<double> curvedf;

};

//DECLARE(SVQmcDiscFactorMM);



class SVQmcExpectedDiscFactor : public IQMCStateVariableExpected, public SVExpectedDiscFactor
{
public:
    virtual ~SVQmcExpectedDiscFactor() {}

    SVQmcExpectedDiscFactor(
        IQMCDiffusibleInterestRate*    _pAsset,
        const DateTime&                _measurementDate,
        const DateTimeArray&           _futureDates,
        bool                           _doLog,
        YieldCurveConstSP			   _yc);

    /** asset specific name of one of the generic functions */
    inline double getFwdDF( int idx ) const {
        return element(idx);
    }

    /** asset specific name of one of the generic functions */
    /*inline*/ double firstDF() const {
        return element(0);
    }

    virtual double element(int cfIndex) const
    {
        return  doLog ?  pAsset->getLnExpectedDiscFactor(ycIdx,
                        								measurementDateIdx,
                                                        futureDatesIdx[cfIndex]):
                        pAsset->getExpectedDiscFactor(ycIdx,
                        						      measurementDateIdx,
                                                      futureDatesIdx[cfIndex]);
    }

    virtual void prepare(bool mm);

protected:
    SVQmcExpectedDiscFactor();
    IQMCDiffusibleInterestRate* pAsset;
    size_t ycIdx; // choosing a proper index YC for EDF

};

//DECLARE(SVQmcExpectedDiscFactor);

// MM extension
class SVQmcExpectedDiscFactorMM : public SVQmcExpectedDiscFactor
{
public:
    virtual ~SVQmcExpectedDiscFactorMM() {}

    SVQmcExpectedDiscFactorMM(
        IQMCDiffusibleInterestRate*    _pAsset,
        const DateTime&                _measurementDate,
        const DateTimeArray&           _futureDates,
        bool                           _doLog,
        YieldCurveConstSP			   _yc) ;


    virtual double element(int cfIndex) const
    {
        return  doLog ?  pAsset->getLnExpectedDiscFactor(ycIdx,
                        								measurementDateIdx,
                                                        futureDatesIdx[cfIndex]) + correction.at(cfIndex):
                        pAsset->getExpectedDiscFactor(ycIdx,
                        						      measurementDateIdx,
                                                      futureDatesIdx[cfIndex]) * correction.at(cfIndex);
    }

    virtual void prepare(bool mm);

    void resetMMCorrection();
    void accumulateMMCorrection();
    void setMMCorrection();

private:
    SVQmcExpectedDiscFactorMM();

    // for Moment Matching
    vector<double> correction;  // initialized in prepare()
    double sumDF;
    SpotIdx futureDateFXIdx;
    SpotIdx futureDateDomDFIdx;
    vector<double> curveLnEDF;
    vector<double> accumulation;

};

//DECLARE(SVQmcExpectedDiscFactorMM);




//////////////////////////////////////////////////////////////////////////




/** Credit-spread specific function names realized in these state variables */

class SVQmcSurvivalDiscFactor : public IQMCStateVariableSpot, public SVSurvivalDiscFactor
{
public:
    virtual ~SVQmcSurvivalDiscFactor() {}

    SVQmcSurvivalDiscFactor(
        IQMCDiffusibleCreditSpreadBase*  _pAsset,
        const DateTimeArray&        _measurementDates);

    /** asset specific name of one of the generic functions */
    virtual double firstSDF() const { return element(0); }
    virtual double getSDF(int idx) const { return element(idx); }

    virtual double element(int dtIndex) const
    {
        return pAsset->getSurvivalDiscFactor(measurementDatesIdx[dtIndex]);
    }

    virtual double getRecoveryRate(int dtIndex) const 
    { 
        return pAsset->getRecoveryRate(measurementDatesIdx[dtIndex]); 
    }

    virtual void prepare(bool mm);

protected:

    SVQmcSurvivalDiscFactor();
    IQMCDiffusibleCreditSpreadBase* pAsset;
};


class SVQmcSurvivalDiscFactorMM : public SVQmcSurvivalDiscFactor
{
public:
    virtual ~SVQmcSurvivalDiscFactorMM() {}

    SVQmcSurvivalDiscFactorMM(
        IQMCDiffusibleCreditSpreadBase*  _pAsset,
        const DateTimeArray&             _measurementDates);

    virtual double element(int dtIndex) const
    {
        return pAsset->getSurvivalDiscFactor(measurementDatesIdx[dtIndex]) * correction.at(dtIndex);;
    }

    virtual void prepare(bool mm);

    void resetMMCorrection();
    void accumulateMMCorrection();
    void setMMCorrection();

private:

    SVQmcSurvivalDiscFactorMM();

      // for Moment Matching
    vector<SpotIdx> measurementDFDatesIdx;// not same as measurementDatesIdx, initialized in prepare()
    vector<double> correction;
    vector<double> accumulation;
    vector<double> accumulationDF;
    vector<double> curveSDF;
};



class SVQmcExpSurvDiscFactor : public IQMCStateVariableExpected, public SVExpSurvDiscFactor
{
public:
    virtual ~SVQmcExpSurvDiscFactor() {}

    SVQmcExpSurvDiscFactor(
        IQMCDiffusibleCreditSpreadBase*     _pAsset,
        const DateTime&                _measurementDate,
        const DateTimeArray&           _futureDates,
        bool                           _doLog);

    /** asset specific name of one of the generic functions */
    virtual double firstExpSDF() const { return element(0); }
    virtual double getExpSDF(int idx) const { return element(idx); }
    virtual double getExpRecoveryRate(int idx) const;

    virtual double element(int cfIndex) const
    {
            return doLog ?  pAsset->getLnExpectedSurvivalDiscFactor(
                                                        measurementDateIdx,
                                                        futureDatesIdx[cfIndex]):
                            pAsset->getExpectedSurvivalDiscFactor(
                                                        measurementDateIdx,
                                                        futureDatesIdx[cfIndex]) ;
    }

    virtual void prepare(bool mm);

protected:

    SVQmcExpSurvDiscFactor();
    IQMCDiffusibleCreditSpreadBase* pAsset;

};



class SVQmcExpSurvDiscFactorMM : public SVQmcExpSurvDiscFactor
{
public:
    virtual ~SVQmcExpSurvDiscFactorMM() {}

    SVQmcExpSurvDiscFactorMM(
        IQMCDiffusibleCreditSpreadBase*     _pAsset,
        const DateTime&                _measurementDate,
        const DateTimeArray&           _futureDates,
        bool                           _doLog);

    virtual double element(int cfIndex) const
    {
            return doLog ?  pAsset->getLnExpectedSurvivalDiscFactor(
                                                        measurementDateIdx,
                                                        futureDatesIdx[cfIndex]) + correction.at(cfIndex) :
                            pAsset->getExpectedSurvivalDiscFactor(
                                                        measurementDateIdx,
                                                        futureDatesIdx[cfIndex]) * correction.at(cfIndex) ;
    }

    virtual void prepare(bool mm);

    void resetMMCorrection();
    void accumulateMMCorrection();
    void setMMCorrection();

private:

    SVQmcExpSurvDiscFactorMM();

    vector<double> correction;  // initialized in prepare()

    SpotIdx measureDateSDFIdx;
    SpotIdx measureDateDFIdx;
    FwdIdx measureDateEDFIdx;
    vector<FwdIdx> futureDatesDFIdx;

    double accumulationDF;
    double accumulationDFSP;
    vector<double> accumulation;
    vector<double> accumulationDFEDF;
    vector<double> curveESDF;

};



class SVQmcAggregatedSurvDiscFactor : public IQMCStateVariableBase, public SVAggregatedSurvDiscFactor
{

    public:
        SVQmcAggregatedSurvDiscFactor(
                IQMCDiffusibleCreditSpreadBase*     _pAsset,
                DateTimeArraySP   _sdfDates,            // set of dates for discount factors
                SpotIdxArraySP    _sdfIdxSP,
                DateTimeArraySP   _esdfRequestedDates,  //  union of {t_i} for all expected factors
                FwdIdxArraySP     _esdfReqIdxSP,
                DateTimeArraySP   _esdfForwardDates,    // union of all {T_j} for all expected factors
                FwdIdxArraySP     _esdfForIdxSP,
                const DateTime&   maxDiffDate,
                const DateTime&   maxCurveMat,
                bool _doLog);

        virtual void prepare(bool mm);

        IQMCDiffusibleCreditSpreadBase*          getAsset(void) const {return pAsset;}
    private:
        DateTimeArraySP   sdfDates;            // set of dates for discount factors
        SpotIdxArraySP    sdfIdxSP;

        DateTimeArraySP   esdfRequestedDates;  //  union of {t_i} for all expected factors
        FwdIdxArraySP     esdfReqIdxSP;

        DateTimeArraySP   esdfForwardDates;    // union of all {T_j} for all expected factors
        FwdIdxArraySP     esdfForIdxSP;

        SVQmcAggregatedSurvDiscFactor();
        IQMCDiffusibleCreditSpreadBase* pAsset;
};

class SVQmcDateOfDefault : public IQMCStateVariableSpot, public SVDateOfDefault
{
public:
    virtual ~SVQmcDateOfDefault() {}

    SVQmcDateOfDefault(
        IQMCDiffusibleDefaultableCreditSpread*  _pAsset);

    /** asset specific name of one of the generic functions */
    virtual DateTime getDateOfDefault() const { return pAsset->getDateOfDefault(); }

    virtual double   getRecoveryRateAtDefault() const { return pAsset->getRecoveryRateAtDefault(); }

    virtual double element(int dtIndex) const
    {
        throw ModelException("SVDateOfDefault::element",
            "incorrect way to access this information. Use getDateOfDefault() instead");
        return 0.0;
    }

    virtual void prepare(bool mm);

private:

    SVQmcDateOfDefault();
    IQMCDiffusibleDefaultableCreditSpread* pAsset;

//  Moment Matching not supported
};








//////////////////////////////////////////////////////////////////////////

/** FX specific function names realized in these state variables */
class SVSpotFX : public IQMCStateVariableSpot
{
public:
    virtual ~SVSpotFX() {}

    SVSpotFX(
        IQMCDiffusibleFX*            _pAsset,
        const DateTimeArray&        _measurementDates);


    /** asset specific name of one of the generic functions */
    inline double getSpotFX() const { return element(0); }

    /** asset specific name of one of the generic functions */
    inline double getSpotPrice() const { return element(0);} // legacy

    virtual double element(int dtIndex) const
    {
        return pAsset->getSpotFX(measurementDatesIdx[dtIndex]);
    }

    virtual void prepare(bool mm);

protected:
    SVSpotFX();
    IQMCDiffusibleFX* pAsset;

};

DECLARE(SVSpotFX);

/** moment-matched FX specific function names realized in these state variables */
class SVSpotFXMM : public SVSpotFX
{
public:
    virtual ~SVSpotFXMM() {}

    SVSpotFXMM(
        IQMCDiffusibleFX*            _pAsset,
        const DateTimeArray&        _measurementDates);


    virtual double element(int dtIndex) const
    {
        return pAsset->getSpotFX(measurementDatesIdx[dtIndex]) * correction.at(dtIndex);
    }

    virtual void prepare(bool mm);

    // Moment matching condition
    // Expectation { FX(t) DF_dom(t) } = FX(0) DF_for(t) = \bar{FX(t)} \bar{DF_for(t)}
    // where \bar{.} mean for forward value of.
    void resetMMCorrection();
    void accumulateMMCorrection();
    void setMMCorrection();

private:
    SVSpotFXMM();

    // for Moment Matching
    vector<SpotIdx> measurementDFDatesIdx;// not same as measurementDatesIdx, initialized in prepare()
    vector<double> correction;            // initialized in prepare()
    vector<double> accumulation;
    vector<double> expFX;
    vector<double> accumulationDomDF;

};


class SVQmcExpectedFX : public IQMCStateVariableExpected, public SVExpectedFX
{
public:
    virtual ~SVQmcExpectedFX() {}

    SVQmcExpectedFX(
        IQMCDiffusibleFX*               _pAsset,
        const DateTime&                _measurementDate,
        const DateTimeArray&           _futureDates,
        bool                           _doLog);


    /** asset specific name of one of the generic functions */
    inline double getFwdFX(int idx) const { return element(idx); }

    /** asset specific name of one of the generic functions */
    inline double getFwd(int idx) const { return element(idx); } // legacy name

    virtual double element(int cfIndex) const
    {
        return doLog ?  pAsset->getLnExpectedFX(measurementDateIdx, futureDatesIdx[cfIndex]) :
                        pAsset->getExpectedFX(measurementDateIdx, futureDatesIdx[cfIndex]) ;
    }

    virtual void prepare(bool mm);

protected:

    SVQmcExpectedFX();
    IQMCDiffusibleFX* pAsset;
};


class SVQmcExpectedFXMM : public SVQmcExpectedFX
{
public:
    virtual ~SVQmcExpectedFXMM() {}

    SVQmcExpectedFXMM(
        IQMCDiffusibleFX*              _pAsset,
        const DateTime&                _measurementDate,
        const DateTimeArray&           _futureDates,
        bool                           _doLog);


    virtual double element(int cfIndex) const
    {
        return doLog ?  pAsset->getLnExpectedFX(measurementDateIdx, futureDatesIdx[cfIndex])
                                                + correction.at(cfIndex):
                        pAsset->getExpectedFX(measurementDateIdx, futureDatesIdx[cfIndex])
                                                * correction.at(cfIndex);
    }

    virtual void prepare(bool mm);

    void resetMMCorrection();
    void accumulateMMCorrection();
    void setMMCorrection();

private:

    SVQmcExpectedFXMM();

    SpotIdx        dfDateIdx;
    vector<FwdIdx> edfDomDateIdx;
    FwdIdx         edfDomMeasureDateIdx;
    vector<double> correction;
    vector<double> accumulationDF;
    vector<double> accumulation;
    vector<double> originalExpectedFX;
};



//////////////////////////////////////////////////////////////////////////

/** EQ specific function names realized in these state variables */
class SVSpotEQ : public IQMCStateVariableSpot
{
public:
    virtual ~SVSpotEQ() {}

    SVSpotEQ(
        IQMCDiffusibleEQ*            _pAsset,
        const DateTimeArray&        _measurementDates);


    /** asset specific name of one of the generic functions */
    inline double getSpotPrice() const { return element(0); } // legacy

    virtual double element(int dtIndex) const
    {
        return pAsset->getSpotPrice(measurementDatesIdx[dtIndex]);
    }

    virtual void prepare(bool mm);

protected:

    SVSpotEQ();
    IQMCDiffusibleEQ* pAsset;
};

DECLARE(SVSpotEQ);

/** EQ specific function names realized in these state variables */
class SVSpotEQMM : public SVSpotEQ
{
public:
    virtual ~SVSpotEQMM() {}

    SVSpotEQMM(
        IQMCDiffusibleEQ*            _pAsset,
        const DateTimeArray&        _measurementDates);


    virtual double element(int dtIndex) const
    {
        return pAsset->getSpotPrice(measurementDatesIdx[dtIndex]) * correction.at(dtIndex);
    }

    virtual void prepare(bool mm);

    void resetMMCorrection();
    void accumulateMMCorrection();
    void setMMCorrection();

private:

    SVSpotEQMM();

    vector<SpotIdx> measurementDFIdx;
    vector<SpotIdx> measurementFXIdx;
    vector<double> correction;
    vector<double> accumulation;
    vector<double> accumulationDF;
    vector<double> origFwdPrice;
    IQMCDiffusibleInterestRateSP domPtr;
    IQMCDiffusibleFX* fxPtr;
};



class SVQmcExpectedEQ : public IQMCStateVariableExpected, public SVExpectedEQ
{
public:
    virtual ~SVQmcExpectedEQ() {}

    SVQmcExpectedEQ(
        IQMCDiffusibleEQ*               _pAsset,
        const DateTime&                _measurementDate,
        const DateTimeArray&           _futureDates,
        bool                           _doLog);

    /** asset specific name of one of the generic functions */
    inline double getFwdPrice(int idx) const { return element(idx); }

    /** asset specific name of one of the generic functions */
    inline double getFwd(int idx) const { return element(idx); } // legacy name

    virtual double element(int cfIndex) const
    {
        return doLog ?  pAsset->getLnExpectedPrice(measurementDateIdx,
                                                   futureDatesIdx[cfIndex]):
                        pAsset->getExpectedPrice(measurementDateIdx,
                                                 futureDatesIdx[cfIndex]);
    }

    virtual void prepare(bool mm);

protected:
    SVQmcExpectedEQ();
    IQMCDiffusibleEQ* pAsset;

};


class SVQmcExpectedEQMM : public SVQmcExpectedEQ
{
public:
    virtual ~SVQmcExpectedEQMM() {}

    SVQmcExpectedEQMM(
        IQMCDiffusibleEQ*               _pAsset,
        const DateTime&                _measurementDate,
        const DateTimeArray&           _futureDates,
        bool                           _doLog);

    virtual double element(int cfIndex) const
    {
        return doLog ?  pAsset->getLnExpectedPrice(measurementDateIdx,
                                                   futureDatesIdx[cfIndex]) + correction.at(cfIndex):
                        pAsset->getExpectedPrice(measurementDateIdx,
                                                 futureDatesIdx[cfIndex]) * correction.at(cfIndex);
    }

    virtual void prepare(bool mm);

    void resetMMCorrection();
    void accumulateMMCorrection();
    void setMMCorrection();

private:
    SVQmcExpectedEQMM();

    SpotIdx dfDomSpotIdx;
    SpotIdx fxSpotIdx;
    FwdIdx  edfForSpotIdx;
    vector<FwdIdx> edfForFwdIdx;

    vector<double> correction;

    vector<double> origExpectedFwdPrice;
    vector<double> accumulation;
    vector<double> accumulationDF;
    IQMCDiffusibleInterestRateSP domPtr; // domestic ccy ptr
    IQMCDiffusibleInterestRateSP undPtr; // underlying ccy ptr, whether it is domestic or foreign
    IQMCDiffusibleFX* fxPtr;    // or NULL if priced in domestic ccy
};


/////////// Energy specific state variables ////////////////////////////////////
// Generic expected energy SV, can be used for future and spread:
class SVQmcExpEnergyFuture : public IQMCStateVariableExpected, public SVExpEnergyFuture
{
public:
    virtual ~SVQmcExpEnergyFuture() {}

    SVQmcExpEnergyFuture(
        IQMCDiffusibleEnergy * _pAsset,
        const DateTime & _measurementDate,
        const DateTimeArray & _futureDates,
        bool _doLog ) ;

    const DateTimeArray & getFutureDates() const { return futureDates; }

    /** asset specific name of one of the generic functions */
    inline double getExpFuturePrice( int idx ) const {
		return (*this)[idx];	// return element(idx);
    }

    virtual double element(int cfIndex) const
    {
        assert(prepared);
		if (doLog)
			return pAsset->getLnExpectedPrice(measurementDateIdx, futureDatesIdx[cfIndex]);
		else
			return pAsset->getExpectedPrice(measurementDateIdx, futureDatesIdx[cfIndex]);
    }

    virtual void prepare(bool mm);

protected:

    SVQmcExpEnergyFuture();
    IQMCDiffusibleEnergy * pAsset;

    bool prepared; // for debug purposes, to detect access to SV before prepare() method was called
};

/////////// Energy specific state variables ////////////////////////////////////
// Generic expected energy SV, can be used for future and spread:
class SVQmcExpEnergyFutureMM : public SVQmcExpEnergyFuture
{
public:
    virtual ~SVQmcExpEnergyFutureMM() {}

    SVQmcExpEnergyFutureMM(
        IQMCDiffusibleEnergy * _pAsset,
        const DateTime & _measurementDate,
        const DateTimeArray & _futureDates,
        bool _doLog ) ;

    virtual double element(int cfIndex) const
    {
        assert(prepared);
		if (doLog)
			return pAsset->getLnExpectedPrice(measurementDateIdx, futureDatesIdx[cfIndex]) + correction.at(cfIndex);
		else
			return pAsset->getExpectedPrice(measurementDateIdx, futureDatesIdx[cfIndex]) * correction.at(cfIndex);
    }

    virtual void prepare(bool mm);

    void resetMMCorrection();
    void accumulateMMCorrection();
    void setMMCorrection();

private:

    SVQmcExpEnergyFutureMM();

    int nAcc;
    vector<double> correction;
    vector<double> accumulation;
    vector<double> origLnFuturePrice;


protected:
    bool prepared; // for debug purposes, to detect access to SV before prepare() method was called
};

/////////////////////////////////////////////////////////////////////////

class SVQmcExpectedBasisFwdSpread : public IQMCStateVariableExpected, public SVExpectedBasisFwdSpread
{
public:
    virtual ~SVQmcExpectedBasisFwdSpread() {}

    // pAsset will compute the TMat from resetDate and the basis index curve's
    // ref tenor (3ML, etc.).  The basis index curve is independently collected
    // by the pathgenerator and used in setting up pAsset so it is not
    // needed here.

    SVQmcExpectedBasisFwdSpread(
        IQMCDiffusibleBasisIndex* _pAsset,
        const DateTime& _measurementDate,
        const DateTimeArray& _resetDate);

    // asset specific name of one of the generic functions
    virtual double getFwdSpread( int idx ) const {
        return element(idx);
    }

    // asset specific name of one of the generic functions
    virtual double firstSpread() const {
        return element(0);
    }

    virtual double element(int fwdSpreadIndex) const
    {
        assert(prepared);
        double result = pAsset->getExpectedBasisFwdSpread(measurementDateIdx, futureDatesIdx[fwdSpreadIndex]);
        return result;
    }

    virtual void prepare(bool mm);

private:

    SVQmcExpectedBasisFwdSpread();
    IQMCDiffusibleBasisIndex* pAsset;

    // choosing a proper basis curve for expected basis coupons
protected:
    bool prepared; // for debug purposes, to detect access to SV before prepare() method was called
};

class MCPathConfigSRMGen;

// this is a simple implementation of time-independent weight
class SVQmcPathWeightTimeIndependent : public IQMCStateVariableSpot, public SVPathWeight
{
public:
    virtual ~SVQmcPathWeightTimeIndependent() {}

    SVQmcPathWeightTimeIndependent(MCPathConfigSRMGen*  _pPathConfig) : 
        IQMCStateVariableSpot(DateTimeArray()), 
        pPathConfig(_pPathConfig) {}

    virtual void prepare(bool /*mm*/) {}

    virtual double getWeight(int i) const
    { return element(i); }

    virtual double element(int /*idx*/ ) const;

private:
    SVQmcPathWeightTimeIndependent();
    MCPathConfigSRMGen* pPathConfig;

};




/*
class SVExpectedBasisCoupon : public IQMCStateVariableExpected
{
public:
virtual ~SVExpectedBasisCoupon() {}

SVExpectedBasisCoupon(
SVGenIRSwap::IStateVarSP _swapYieldSV, // To do: only defined in cpp!!
SVExpectedBasisFwdSpreadSP _fwdSpreadSV );


// asset specific name of one of the generic functions
inline double getExpBasisCoupon( int idx ) const {
return element(idx);
}

// asset specific name of one of the generic functions
inline double firstDF() const {
return element(0);
}

virtual double element(int cfIndex) const;

virtual void prepare();

private:

SVExpectedBasisCoupon();
SVGenIRSwap::IStateVarSP swapYieldSV;
SVExpectedBasisFwdSpreadSP fwdSpreadSV;

protected:
bool prepared; // for debug purposes, to detect access to SV before prepare() method was called
};

DECLARE(SVExpectedBasisCoupon);
*/

DRLIB_END_NAMESPACE
#endif // EDR_IQMCSTATEVARIABLE_HPP

