/**
 * @file FieldSensitivityDefinition.hpp
 */

#ifndef QLIB_FieldSensitivityDefinition_H
#define QLIB_FieldSensitivityDefinition_H

#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(OutputName)
FORWARD_DECLARE(IResultsFunction)
FORWARD_DECLARE(IScalarDerivative)
FORWARD_DECLARE(RiskQuantityFactorySensitivity)
FORWARD_DECLARE(IFieldRiskPropertyDefinition)
FORWARD_DECLARE(FieldSensitivityDefinition)

// should just be static methods on FieldSensitivityDefinition
// but that crashes VC6

template <class T>
struct RQFFactory_1packet {
    static RiskQuantityFactorySensitivitySP newOne(
            const double* shiftSize,
            const string* packetName1, const string* packetName2) {
        return RiskQuantityFactorySensitivitySP(
            new T(shiftSize ? *shiftSize : T::DEFAULT_SHIFT,
                  packetName1 ? *packetName1 : string(T::NAME)));
    }
};

template <class T>
struct RQFFactory_2packet {
    static RiskQuantityFactorySensitivitySP newOne(
            const double* shiftSize,
            const string* packetName1, const string* packetName2) {
        return RiskQuantityFactorySensitivitySP(
            new T(shiftSize ? *shiftSize : T::DEFAULT_SHIFT,
                  packetName1 ? *packetName1 : string(T::NAME),
                  packetName2 ? *packetName2 : string(T::NAME2)));
    }
};

/**
 * Client-friendly vehicle for dynamically defining a RiskPropertySensitivity
 *
 * This class is specific to the "flexible scenarios and greeks" facility: it's
 * a flat, text-based way of defining at least the most common types of greek.
 *
 * See FlexibleSensitivity for an overview of the framework, and
 * IRiskQuantityFactory for background on 'declarative' greeks.
 *
 *
 * <H3>Usage</H3>
 *
 * Internally, FieldSensitivityDefinition is a factory for
 * RiskQuantityFactorySensitivity's: see sensitivity().  It's used by
 * FlexibleSensitivity to designate the greek to be computed.
 *
 *
 * <H3>Interface</H3>
 *
 * [This section is reproduced from "Flexible scenarios and greeks" in the QLib
 * doc db:
 *
 * -     Notes://PPUSMC006/8525703B0051D832/55B689E36191AD7285256DFD00470F6B/272C882E688B88B3852570CB00403228
 *
 * ].
 *
 * A Sensitivity represents a greek, i.e. a derivative of (typically)
 * instrument value with respect to something.  For example, QLib's built-in
 * Delta reports a two-sided estimate of the derivative of value with respect
 * to spot, for each object carrying the property "spot".
 *
 * To define a new Sensitivity through the flexible interface, you need to name
 * the quantity to be differentiated (derivand below) and the kind of
 * derivative (derivative1 below), designate the RiskProperty with respect to
 * which the derivative is to be taken (property below), and provide a name
 * under which it's to appear in the results (packetName1 below).  You can also
 * specify a constant factor to change the units in which the greek is reported
 * (unit below: 0.0001 means "per basis point").
 *
 * The defaultShiftSize governs the magnitude of the property tweaks used to
 * estimate the derivative, unless overridden when the sensitivity is requested
 * for a particular pricing run (see FlexibleSensitivity below).
 *
 * <H4>Scalar or pointwise resolution</H4>
 *
 * If the property with respect to which the greek is calculated is
 * term-structured (see <I>Vector and matrix properties</I> in the
 * IFieldRiskPropertyDefinition docs), you can choose to have it reported either
 * with respect to each expiry separately, or with respect to a parallel tweak
 * to all entries at once (resolution below = POINTWISE or SCALAR).
 *
 * For matrix-valued properties, POINTWISE resolution is still with respect to
 * expiries, i.e. the greek is reported per expiry, using parallel tweaks to
 * each row (or column) of the matrix in turn.  There is no mechanism for
 * reporting greeks with respect to individual entries in the matrix.
 *
 * <H4>Reporting both two-sided and second derivatives</H4>
 *
 * Mostly in QLib, we report two-sided first derivatives and second derivatives
 * together, since they're obtaining using the same repricings.  In the
 * flexible interface you can get the same effect by specifying another
 * derivative (derivative2 below) and a name under which it's to appear in the
 * results (packetName2 below).  The operator and unit will be the same as for
 * derivative1.
 */

class RISKMGR_DLL FieldSensitivityDefinition: public CObject {

    FieldSensitivityDefinition(const FieldSensitivityDefinition& rhs);
    FieldSensitivityDefinition& operator=(
        const FieldSensitivityDefinition& rhs);

    FieldSensitivityDefinition();
    static IObject* emptyShell();
    static void load(CClassSP& clazz);
    void validatePop2Object();

public:

    static CClassConstSP const TYPE;

private:

    /**
     * A name for the sensitivity (not the same as packetName1)
     *
     * See registerBuiltin().
     */

    string name; // $required(Unique name for the sensitivity)

    /**
     * Whether to use a hard-coded implementation of the sensitivity,
     * or follow this "flexible" definition
     *
     * FLEXIBLE, INTERNAL or DEFAULT.  See registerBuiltin().
     */

    string implementation; // $required(FLEXIBLE, INTERNAL, DEFAULT)

    /**
     * What the greek is a derivative of
     *
     * PRICE or an OutputRequest name.
     */

    string derivand; // $required(What the greek is a derivative of, e.g. PRICE)

    /**
     * What the greek is a derivative with respect to
     *
     * IFieldRiskPropertyDefinition is a factory for RiskPropertyDefinition's.
     */

    // $required(What the greek is a derivative with respect to)
    IFieldRiskPropertyDefinitionConstSP property;

    /**
     * Optional second property, making the greek a cross-gamma
     */

    // $optional(Second property with respect to which to take deriv, for cross-gammas)
    IFieldRiskPropertyDefinitionConstSP property2;

    /**
     * Whether to compute the greek for all expiries
     *
     * SCALAR, POINTWISE, ELEMENTWISE, POINTWISExPOINTWISE, or SCALARxSCALAR.
     */

    // $required(SCALAR, POINTWISE, ELEMENTWISE, POINTWISExPOINTWISE, or SCALARxSCALAR)
    string resolution;

    /**
     * Units in which to report greek
     */

    double unit; // $required(Units in which to report greek, e.g. 0.0001 for "per bp")

    /**
     * What derivative the greek is
     *
     * "ONE_SIDED", "TWO_SIDED" or "SECOND".
     */

    string derivative1; // $required(Deriv to compute: ONE_SIDED, TWO_SIDED or SECOND)

    /**
     * The name under which the greek gets reported
     *
     * A packet name in the Results.
     */

    string packetName1; // $required(Packet name under which to report the greek)

    /**
     * Optional other derivative to compute
     *
     * "ONE_SIDED", "TWO_SIDED", "SECOND" or absent (in which case value is "").
     */

    string derivative2; // $optional(Deriv to compute in addition to 'derivative1')

    /**
     * If derivative2 present, name under which it gets reported
     *
     * A packet name in the Results.
     */

    string packetName2; // $optional(Packet name under which to report derivative2)

    /**
     * Shift size to use in estimating derivative1
     */

    double defaultShiftSize;  // $required(Shift size for estimating derivative1)

    /**
     * Shift size to use in estimating derivative2
     */

    double defaultShiftSize2; // $optional(Shift size for estimating derivative2)

    IScalarDerivativeConstSP _derivative1;  // $transient
    IScalarDerivativeConstSP _derivative2;  // $transient
    IResultsFunctionConstSP _derivand;      // $transient
    bool _useOverride, _requireOverride;    // $transient

public:

    ~FieldSensitivityDefinition();

    /**
     * A RiskQuantityFactorySensitivity constructed according to this
     * definition
     *
     * Called by FlexibleSensitivity::validatePop2Object().
     *
     * Returns either a ScalarRiskPropertySensitivity (if resolution == SCALAR)
     * or an ExpiryRiskPropertySensitivity (if resolution == POINTWISE).
     *
     * If given, @a shiftSize, @a packetName and @a packetName2 override
     * respectively defaultShiftSize, packetName1 and packetName2.
     */

    RiskQuantityFactorySensitivitySP sensitivity(
        OutputNameArrayConstSP toTweak,
        const double* shiftSize,
        const string* packetName,
        const string* packetName2) const;

    /**
     * Register a hard-coded sensitivity implementation to override a flexible
     * one
     *
     * If a FieldSensitivityDefinition comes in with name = the @a name you give
     * here and implementation == "DEFAULT" or "INTERNAL", then its
     * sensitivity() factory method will ignore the rest of the definition and
     * call instead your @a builtinConstructor.
     *
     * If implementation == "FLEXIBLE", however, then the flexible version will
     * be returned even if you've registered a builtin version.
     */

    static void registerBuiltin(
        const string& name,
        RiskQuantityFactorySensitivitySP (*builtinConstructor)(
            const double* shiftSize,
            const string* packetName,
            const string* packetName2));
};

DRLIB_END_NAMESPACE

#endif
