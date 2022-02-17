//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CDSParSpreadsLegalBasis.cpp
//
//   Description : Legal Basis adjustment for CDSParSpreads
//   
//   Author      : Antoine Gregoire 
//
//   Date        : Jan 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_CDSPARSPREADSLEGALBASIS_CPP
#include "edginc/RiskProperty.hpp"
#include "edginc/CDSParSpreadsLegalBasis.hpp"


DRLIB_BEGIN_NAMESPACE

/** Constructor */
CDSParSpreadsLegalBasis::CDSParSpreadsLegalBasis() : CDSParSpreadsAdjustment(TYPE) {}

CDSParSpreadsLegalBasis:: ~CDSParSpreadsLegalBasis() {}

/** Used for reflection */
IObject* CDSParSpreadsLegalBasis::defaultConstructor() {
     return new CDSParSpreadsLegalBasis();
}

/** Invoked when Class is 'loaded', used for reflection */
void CDSParSpreadsLegalBasis::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(CDSParSpreadsLegalBasis, clazz);
    SUPERCLASS(CDSParSpreadsAdjustment);
    IMPLEMENTS(TweakableWith<LegalBasisAdditiveParallelTweak>);
    IMPLEMENTS(ITweakableWithRespectTo<LegalBasisAbsolutePointwise>);
    IMPLEMENTS(TweakableWith<LegalBasisMultiplierParallelTweak>);
    IMPLEMENTS(ITweakableWithRespectTo<LegalBasisAdditiveRelativeParallelTweak>);
    IMPLEMENTS(ITweakableWithRespectTo<LegalBasisRelativePointwise>);
    IMPLEMENTS(TweakableWith<LegalBasisAdditiveRecoveryTweak>);
    IMPLEMENTS(TweakableWith<LegalBasisMultiplierRecoveryTweak>);
    IMPLEMENTS(QuasiContractualBaseCorrelation::IShift);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(name,
        "Name of that market object");
    FIELD(expiries,
        "Expiries array. Must match those in the par spread curve");
    FIELD(parSpreadsMultiplierCoefficients,
        "Multiplier coefficients to apply to par spread curve");
    FIELD(parSpreadsAdditiveCoefficients, 
        "Additive coefficients to apply to par spread curve");
    FIELD(recoveryMultiplierCoefficient,
        "Multiplier coefficient to apply to recovery");
    FIELD(recoveryAdditiveCoefficient,
        "Additive coefficient to apply to recovery");
}

/** Returns the name of this object. This is the name with which
 * it is stored in the market data cache and is the name with
 * which results (eg tweaks) should be reported against */
string CDSParSpreadsLegalBasis::getName() const {
    return name;
}


void CDSParSpreadsLegalBasis::validatePop2Object() {
    static const string method = "CDSParSpreadsLegalBasis::validatePop2Object";

    // checks lengths of arrays are the same
    if (parSpreadsMultiplierCoefficients.size() != expiries->size()
        || parSpreadsAdditiveCoefficients.size() != expiries->size())
    {
        throw ModelException(method,
            "Wrong size : "
            "Expiries (size=" +
            Format::toString(expiries->size()) +
            "), multiplier coefficients (size=" +
            Format::toString(parSpreadsMultiplierCoefficients.size()) +
            ") and additive coefficients (size=" +
            Format::toString(parSpreadsAdditiveCoefficients.size()) +
            ") have different sizes.");    
    }
}


bool CDSParSpreadsLegalBasis::expiriesValid(const ICDSBootstrappableWrapper& cdsToAdjust) const {
    return Expiry::equals(
                      cdsToAdjust->getValueDate(),
                      cdsToAdjust->getParSpreadsExpiries().get(),
                      expiries.get());
}

/** Returns the expiries themselves */
ExpiryArrayConstSP CDSParSpreadsLegalBasis::getAdjustmentExpiries() const
{
    return expiries;
}

double CDSParSpreadsLegalBasis::getAdjustedRecovery (const ICDSBootstrappableWrapper& cdsToAdjust) const {
    double recovery = recoveryMultiplierCoefficient * cdsToAdjust->getRecovery() +
                        recoveryAdditiveCoefficient;

    // Adjust recovery so that it never goes below 0 or above 1 (100%)
    if (recovery < 0) {
        recovery = 0;
    }
    else if (recovery > 1) {
        recovery = 1;
    }
    return recovery;
}



DoubleArrayConstSP CDSParSpreadsLegalBasis::getAdjustedParSpreads (const ICDSBootstrappableWrapper& cdsToAdjust) const {

    DoubleArrayConstSP unadjustedSpreads = cdsToAdjust->getParSpreads();
    DoubleArraySP adjustedParSpreads = DoubleArraySP(
        new DoubleArray(unadjustedSpreads->size()));
    
    for (int i = 0; i < adjustedParSpreads->size(); ++i) {
        (*adjustedParSpreads)[i] = parSpreadsAdditiveCoefficients[i] +
                                   (parSpreadsMultiplierCoefficients[i] * 
                                    (*unadjustedSpreads)[i]);
    }
    return adjustedParSpreads;
}

/** Returns a single adjusted rate, corresponding to expiry */
double CDSParSpreadsLegalBasis::getAdjustedParSpread(const ICDSBootstrappable* cdsToAdjust,
                                                     ExpiryConstSP expiry) const
{
    static const string method = "CDSParSpreadsLegalBasis::getAdjustedParSpread";
    try
    {
        //find the expiry - exception thrown in search method if not found
        int idx = expiry->search(cdsToAdjust->getParSpreadsExpiries().get());

        //get and adjusted spread
        DoubleArrayConstSP unadjustedSpreads = cdsToAdjust->getParSpreads();
        double spd = parSpreadsAdditiveCoefficients[idx] +
                     (parSpreadsMultiplierCoefficients[idx] *
                     (*unadjustedSpreads)[idx]);
        return spd;
    }
    catch(exception& e)
    {
       throw ModelException(e, method);
    }
}

/** Returns the adjusted maximum tweak sizes */
const DoubleArrayConstSP CDSParSpreadsLegalBasis::adjustMaxTweakSizes(const ICDSBootstrappableWrapper& cdsToAdjust,
                                                                      DoubleArrayConstSP calculatedMaxTwkSizes,
                                                                      IExpiryRiskPropertyConstSP withRespectTo) const
{
    DoubleArray adjTwks(calculatedMaxTwkSizes->size());

    //adjustments will depend upon what is being tweaked
    if (RiskProperty<ParSpreadPointwise>::TYPE->isInstance(withRespectTo))
    {
        //the spread
        for (int i=0; i<calculatedMaxTwkSizes->size(); i++)
        {
            adjTwks[i] = (*calculatedMaxTwkSizes)[i]
                         / 
                         parSpreadsMultiplierCoefficients[i];
        }
    }
    else if (RiskProperty<LegalBasisAbsolutePointwise>::TYPE->isInstance(withRespectTo))
    {
        //the additive component
        for (int i=0; i<calculatedMaxTwkSizes->size(); i++)
        {
            adjTwks[i] = (*calculatedMaxTwkSizes)[i];
        }
    }
    else if (RiskProperty<LegalBasisRelativePointwise>::TYPE->isInstance(withRespectTo))
    {
        //the multiplicative component
        DoubleArrayConstSP parSpreads = cdsToAdjust->getParSpreads();
        for (int i=0; i<calculatedMaxTwkSizes->size(); i++)
        {
            adjTwks[i] = (*calculatedMaxTwkSizes)[i] / (*parSpreads)[i];
        }
    }
    else
    {
        throw ModelException("CDSParSpreadsLegalBasis::adjustMaxTweakSizes",
                             "unknown tweak associated with maxTweakSizes.");
    }

    DoubleArrayConstSP maxTweaks(new DoubleArray(adjTwks));
    return DoubleArrayConstSP(maxTweaks);
}

/*****************************************
 * Additive curve coefficients - Parallel
 *****************************************/

/** Returns name identifying Legal Basis for a Parallel Tweak to the additive coefficients */
string CDSParSpreadsLegalBasis::sensName(LegalBasisAdditiveParallelTweak* shift) const{
    return name;
}


/** Shifts the additive coefficients using given shift */
bool CDSParSpreadsLegalBasis::sensShift(LegalBasisAdditiveParallelTweak* shift){
    static const string method = "CDSParSpreadsLegalBasis::sensShift";
	
    try {
        double shiftSize = shift->getShiftSize();

        if (!Maths::isZero(shiftSize)) {
            for (int i = 0; i < parSpreadsAdditiveCoefficients.size(); ++i) {
                parSpreadsAdditiveCoefficients[i] += shiftSize;
            }
        }
    } 
    catch (exception& e) {
       throw ModelException(e, method,
							"LegalBasisAdditiveParallelTweak failed for " + getName());
   }
   return false; // No sub-components need tweaking
}


/******************************************
 * Additive curve coefficients - Pointwise
 ******************************************/

/** Returns name identifying Legal Basis for a Pointwise Tweak to the additive coefficients */
string CDSParSpreadsLegalBasis::sensName(const LegalBasisAbsolutePointwise*) const {
    return name;
}

/** Shifts the additive coefficients using given TweakSizes */
TweakOutcome CDSParSpreadsLegalBasis::sensShift(const PropertyTweak<LegalBasisAbsolutePointwise>& tweak) {
    try {
        if (!Maths::isZero(tweak.coefficient)) {
            int i = tweak.qualifier->expiry->search(expiries.get());
            parSpreadsAdditiveCoefficients[i] += tweak.coefficient;
        }

        return TweakOutcome(tweak.coefficient,
                            false); // none of our components has this sensitivity
    }
    catch (exception& e) {
        throw ModelException(e, "CDSParSpreadsLegalBasis::sensShift()",
                             "LegalBasisAdditivePointwiseTweak failed for " + getName());
    }
}

ExpiryWindowArrayConstSP CDSParSpreadsLegalBasis::sensQualifiers(const LegalBasisAbsolutePointwise*) const {
    return ExpiryWindow::series(expiries);
}


/*******************************************
 * Multiplier curve coefficients - Parallel
 *******************************************/

/** Returns name identifying Legal Basis for a Parallel Tweak to the additive coefficients */
string CDSParSpreadsLegalBasis::sensName(LegalBasisMultiplierParallelTweak*) const {
    return name;
}


/** Shifts the additive coefficients using given shift */
bool CDSParSpreadsLegalBasis::sensShift(LegalBasisMultiplierParallelTweak* shift) {
    static const string method = "CDSParSpreadsLegalBasis::sensShift";

    try {
        double shiftSize = shift->getShiftSize();

        if (!Maths::isZero(shiftSize)) {
            for (int i = 0; i < parSpreadsMultiplierCoefficients.size(); ++i) {
                parSpreadsMultiplierCoefficients[i] += shiftSize;
            }
        }
    }
    catch (exception& e) {
       throw ModelException(e, method,
							"LegalBasisMultiplierParallelTweak failed for " + getName());
   }
   return false; // No sub-components need tweaking
}



/********************************************
 * Additive curve coefficients - Parallel and relative
 ********************************************/

//// Returns name identifying Legal Basis for a parallel Tweak to the additive coefficients 
string CDSParSpreadsLegalBasis::sensName(const LegalBasisAdditiveRelativeParallelTweak*) const {
    return name;
}

//// Shifts the additive coefficients using given TweakSizes 
TweakOutcome CDSParSpreadsLegalBasis::sensShift(const PropertyTweak<LegalBasisAdditiveRelativeParallelTweak>& tweak) {
    try {
        if (!Maths::isZero(tweak.coefficient)) {
            for (int i = 0; i < parSpreadsAdditiveCoefficients.size(); ++i) {
                parSpreadsAdditiveCoefficients[i] *= (1 + tweak.coefficient);
            }
        }

        return TweakOutcome(tweak.coefficient, false); 
        // none of our components has this sensitivity
    }
    catch (exception& e) {
       throw ModelException(e, "CDSParSpreadsLegalBasis::sensShift()",
                           "LegalBasisAdditiveRelativeParallel tweak failed for " + getName());
    }
}



/********************************************
 * Multiplier curve coefficients - Pointwise
 ********************************************/

/** Returns name identifying Legal Basis for a Pointwise Tweak to the additive coefficients */
string CDSParSpreadsLegalBasis::sensName(const LegalBasisRelativePointwise*) const {
    return name;
}

/** Shifts the additive coefficients using given TweakSizes */
TweakOutcome CDSParSpreadsLegalBasis::sensShift(const PropertyTweak<LegalBasisRelativePointwise>& tweak) {
    try {
        if (!Maths::isZero(tweak.coefficient)) {
            int i = tweak.qualifier->expiry->search(expiries.get());
            parSpreadsMultiplierCoefficients[i] += tweak.coefficient;
        }

        return TweakOutcome(tweak.coefficient,
                            false); // none of our components has this sensitivity
    }
    catch (exception& e) {
       throw ModelException(e, "CDSParSpreadsLegalBasis::sensShift()",
                           "LegalBasisRelativePointwise tweak failed for " + getName());
    }
}

ExpiryWindowArrayConstSP CDSParSpreadsLegalBasis::sensQualifiers(const LegalBasisRelativePointwise*) const {
    return ExpiryWindow::series(expiries);
}


/********************************
 * Additive recovery coefficient
 ********************************/

/** Returns name identifying Legal Basis for a Parallel Tweak to the additive recovery coefficient */
string CDSParSpreadsLegalBasis::sensName(LegalBasisAdditiveRecoveryTweak* shift) const{
    return name;
}


/** Shifts the additive recovery coefficient using given shift */
bool CDSParSpreadsLegalBasis::sensShift(LegalBasisAdditiveRecoveryTweak* shift){
    static const string method = "CDSParSpreadsLegalBasis::sensShift";

    try {
        double shiftSize = shift->getShiftSize();

        if (!Maths::isZero(shiftSize)) {
            recoveryAdditiveCoefficient += shiftSize;
        }
    }
    catch (exception& e) {
        throw ModelException(e, method,
							"LegalBasisAdditiveRecoveryTweak failed for " + getName());
   }
   return false; // No sub-components need tweaking
}


/**********************************
 * Multiplier recovery coefficient
 **********************************/

/** Returns name identifying Legal Basis for a Parallel Tweak to the additive recovery coefficient */
string CDSParSpreadsLegalBasis::sensName(LegalBasisMultiplierRecoveryTweak* shift) const{
    return name;
}


/** Shifts the additive recovery coefficient using given shift */
bool CDSParSpreadsLegalBasis::sensShift(LegalBasisMultiplierRecoveryTweak* shift){
    static const string method = "CDSParSpreadsLegalBasis::sensShift";

    try {
        double shiftSize = shift->getShiftSize();

        if (!Maths::isZero(shiftSize)) {
            recoveryMultiplierCoefficient += shiftSize;
        }
    }
    catch (exception& e) {
       throw ModelException(e, method,
							"LegalBasisMultiplierRecoveryTweak failed for " + getName());
   }
   return false; // No sub-components need tweaking
}


/////////////


// Remove the legal basis adjustment, ie, set the additive coefficients to
// zero and the multiplicative coefficients to 1.
bool CDSParSpreadsLegalBasis::sensShift(QuasiContractualBaseCorrelation* shift) {
    static const string method = "CDSParSpreadsLegalBasis::sensShift(QCBC)";

    try {
        for (int i=0; i < parSpreadsMultiplierCoefficients.size(); ++i) {
            parSpreadsMultiplierCoefficients[i] = 1.0;
            parSpreadsAdditiveCoefficients[i]   = 0.0;
        }
    }
    catch (exception& e) {
       throw ModelException(e, method,
                            "QuasiContractualBaseCorrelation failed for " + 
                            getName());
   }
   return false; // No sub-components need tweaking
}



CClassConstSP const CDSParSpreadsLegalBasis::TYPE = CClass::registerClassLoadMethod(
    "CDSParSpreadsLegalBasis", typeid(CDSParSpreadsLegalBasis), CDSParSpreadsLegalBasis::load);

DEFINE_TEMPLATE_TYPE(CDSParSpreadsLegalBasisWrapper);

DRLIB_END_NAMESPACE
