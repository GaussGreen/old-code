/**
 * @file IFieldRiskPropertyDefinition.hpp
 */

#ifndef QLIB_IFieldRiskPropertyDefinition_H
#define QLIB_IFieldRiskPropertyDefinition_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Range.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/IRiskProperty.hpp"
#include "edginc/IFieldTweak.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(FieldPath)
FORWARD_DECLARE(IFieldRiskPropertyDefinition)

/**
 * Client-friendly vehicle for dynamically defining an IRiskProperty
 *
 * This class is specific to the "flexible scenarios and greeks" facility: it's
 * a flat, text-based way of defining a FieldRiskProperty, i.e. a named field
 * of a named class which can be tweaked in some way in order to construct a
 * scenario or estimate a greek.
 *
 * See FlexibleSensitivity for an overview of the framework, and
 * IRiskQuantityFactory for background on 'declarative' greeks.
 *
 *
 * <H3>Usage</H3>
 *
 * Internally, IFieldRiskPropertyDefinition is a factory for
 * FieldRiskProperty's: see scalarProperty(), parallelProperty(),
 * pointwiseProperty().  It is used by ScalarFlexiblePerturbation and
 * VectorFlexiblePerturbation classes to designate a property to be adjusted to
 * construct a scenario, and by FieldSensitivityDefinition (and hence
 * FlexibleSensitivity) to designate a property to be tweaked to compute a
 * greek.
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
 * A FieldRiskPropertyDefinition describes a field, its term structure (if
 * any), and how it is to be tweaked.
 *
 * [An "IRiskProperty" in QLib represents a property which a market object can
 * carry, and with respect to which risk can be estimated or scenarios
 * constructed.  For instance, QLib's built-in Delta estimates first-order risk
 * with respect to RiskProperty<Spot>.]
 *
 * You need to name the type of object which carries the property (type below),
 * the field in which the property's value is stored (field below), and the
 * kind of operation to which it's subjected (additive or multiplicative ---
 * operation below).
 *
 * Whenever you mention a FieldRiskPropertyDefinition in defining a scenario or
 * greek, it's understood that results will be reported for all market objects
 * of the given type, tagged with each object's market name.
 *
 * <H4>Type</H4>
 *
 * Technically type names a QLib class.  The class need not be externally
 * visible and does not have to have any correlate in Pyramid, although QR will
 * generally be cautious about approving RiskProperty's defined against classes
 * which are not already part of the QLib public interface.
 *
 * Most commonly, type will be a descendant of QLib's internal MarketObject
 * base class for named objects provided by the client in the market data.
 * Other types are supported, but only for tweaks to all objects of the type at
 * once (since if they don't have a predefined name it's impossible to
 * distinguish between them).
 *
 * <H4>Scalar, vector and matrix properties</H4>
 *
 * Some properties are scalar-valued, i.e. are one-dimensional per market name.
 * The field in which they're stored may be of type double, or of type
 * DoubleArray or CDoubleMatrix, in which case a tweak to the property will
 * move all its elements together.
 *
 * Other properties are array-valued and are indexed by "qualifier" as well as
 * market name.  Most commonly, "pointwise" (aka vector) properties are
 * qualified by ExpiryArray.  For such a property, the field you name must be
 * of type DoubleArray, and you need to provide the name of another field, of
 * type ExpiryArray, which holds the expiries for which the property is defined
 * on any given object (expiryQualifier below).  Each entry in the ExpiryArray
 * provides an Expiry label for the corresponding entry in the DoubleArray.
 *
 * A few "matrix" properties are qualified jointly by expiry and something else
 * (typically strike).  In these cases, the field you name must be of type
 * CDoubleMatrix, and you need to provide an expiryQualifier as for pointwise
 * fields.  You also need to specify which dimension corresponds to expiries
 * (expiryDimension below: by QLib internal convention 0 is called "column" and
 * 1 "row").
 *
 * Because the interface at present only supports tweaking of whole
 * rows/columns of matrix-valued fields, as opposed to individual entries,
 * there is no need to provide a field for labelling the other dimension; the
 * subQualifier field is therefore unused and is included only for possible
 * future expansion.
 *
 * <H4>Field paths</H4>
 *
 * In specfying field, expiryQualifier and subQualifier you can provide either
 * a single string, naming a field directly attached to objects of type type,
 * or an array of strings which give a "path" to the field.  For instance if
 * type is "Foo" and field is the array [ "a", "b" ], then you've defined a
 * property stored in the field b of the object referred to by field a of every
 * object of type Foo.
 *
 * <H4>Floor and cap</H4>
 *
 * Many fields must be positive, or between -1 and 1, ?  The optional floor and
 * cap (below) allow you to specify a range outside which QLib will not attempt
 * to shift the field's value.  They default to -8 and +8 respectively if
 * omitted.
 *
 * floorInclusive and capInclusive say whether the floor and cap endpoints are
 * themselves included in the valid range.
 *
 * clip specifies what should happen if a tweak takes the field's value (or in
 * the case of an array value, one of its entries) outside its valid range.  If
 * clip is "false", the tweak simply fails.  Otherwise, the value (or entry) is
 * set equal to the boundary it tried to pass, if that's included in the valid
 * range, or to halfway between the original value and the boundary if it's
 * not.
 *
 * <H4>Additive/multiplicative tweaks</H4>
 *
 * Some fields are most naturally adjusted using additive tweaks of a given
 * magnitude; sometimes multiplicative (aka relative) tweaks make more sense.
 * For example, QLib's built-in Delta sensitivity interprets shiftSize as a
 * proportion of current spot.  In the flexible interface, you can choose
 * either (operator_ below).  AUTO means "whatever is most appropriate given
 * floor and cap"? if neither are given, that's an additive tweak; if one is
 * given, so that the field is bounded on one side, it's a multiplicative tweak
 * relative to the boundary; if both are given, so that the field is bounded to
 * an interval, it's a "sigmoidal" tweak which behaves like a multiplicative
 * tweak at each end, with a smooth changeover in the middle.
 *
 * <H4>Internal or flexible implementation</H4>
 *
 * We envisage that we may sometimes find that tweaks defined through the
 * flexible interface turn out to be inadequate under some unforseen
 * circumstances.  It's proposed that these cases should be handled by
 * including, in a subsequent QLib release, some code which intercepts the
 * flexibly-defined tweak and transparently substitutes a hard-coded one.  All
 * the flexible greeks defined in terms of that tweak will automatically pick
 * up the enhanced implementation.
 *
 * The implementation flag is usually set to DEFAULT, meaning use the
 * hard-wired version if available, else the flexible one.  To ensure use of
 * the former (and fail if it's not there) you can use INTERNAL.  To override
 * the internal implementation (if e.g. it turns out to be buggy, or for
 * comparison) you can use FLEXIBLE.
 */

class RISKMGR_DLL IFieldRiskPropertyDefinition:
        public virtual IObject { // $public

    IFieldRiskPropertyDefinition(const IFieldRiskPropertyDefinition& rhs);
    IFieldRiskPropertyDefinition& operator=(
        const IFieldRiskPropertyDefinition& rhs);

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

    IFieldRiskPropertyDefinition();
    ~IFieldRiskPropertyDefinition();

    /**
     * A ScalarFieldRiskProperty constructed from this definition
     *
     * Used ScalarFlexiblePerturbation and by FieldSensitivityDefinition when
     * its 'resolution' is "SCALAR".
     */

    virtual IScalarRiskPropertyConstSP scalarProperty(
        IObjectConstSP argument) const = 0;

    /**
     * A ParallelFieldRiskProperty constructed from this definition
     *
     * Used by ExpiryFlexiblePerturbation.  Will fail if no expiryQualifier is
     * provided and no getGetExpiriesMethod is registered for the field through
     * Calibrator::IAdjustable::registerBootstrappableField.
     */

    virtual IScalarRiskPropertyConstSP parallelProperty(
        IArrayConstSP qualifiers,
        IArrayConstSP arguments) const = 0;

    /**
     * A PointwiseFieldRiskProperty constructed from this definition
     *
     * Used by FieldSensitivityDefinition when its 'resolution' is POINTWISE.
     * Will fail if no expiryQualifier is provided and no getGetExpiriesMethod
     * is registered for the field through
     * Calibrator::IAdjustable::registerBootstrappableField.
     */

    virtual IExpiryRiskPropertyConstSP pointwiseProperty(IObjectConstSP argument)
        const = 0;

    /**
     * An ExpiryAndStrikeRiskProperty constructed from this definition
     *
     * Used by FieldSensitivityDefinition when its 'resolution' is
     * EXPIRYANDSTRIKE.
     */

    virtual IExpiryAndStrikeRiskPropertyConstSP expiryAndStrikewiseProperty(
        IObjectConstSP argument) const = 0;

    /**
     * An ElementsFieldRiskProperty constructed from this definition
     *
     * Used by FieldSensitivityDefinition when its 'resolution' is ELEMENTS.
     * Will fail if the field is not a DoubleArray or DoubleMatrix.
     */

    virtual IScalarRiskPropertyConstSP elementsProperty(
        IntArrayConstSP indices,
        IArrayConstSP arguments) const = 0;

    /**
     * An ElementwiseFieldRiskProperty constructed from this definition
     *
     * Used by FieldSensitivityDefinition when its 'resolution' is ELEMENTWISE.
     * Will fail if the field is not a DoubleArray or DoubleMatrix.
     */

    virtual IIntRiskPropertyConstSP elementwiseProperty(IObjectConstSP argument)
        const = 0;

    /**
     * See IPerturbation::applyBeforeGetMarket() --- whether scenarios specified
     * using this property definition should be applied to the model/instrument
     * before the "getMarket" phase or after
     */

    virtual bool applyBeforeGetMarket() const = 0;

    /**
     * Register a hard-coded property implementation to override a flexible one
     * (scalar version)
     *
     * If a FieldRiskPropertyDefinition comes in with name == the @a name you
     * give here and implementation == "DEFAULT" or "INTERNAL", then its
     * scalarProperty() factory method will ignore the rest of the definition
     * and return instead your @a builtin IRiskProperty<Void>.
     *
     * If implementation == "FLEXIBLE", however, then the flexible version will
     * be returned even if you've registered a builtin version.
     */

    static void registerScalarBuiltin(
        const string& name,
        IScalarRiskPropertyConstSP (*builtin)(IObjectConstSP argument));

    /**
     * Register a hard-coded property implementation to override a flexible one
     * (parallel version)
     *
     * If a FieldRiskPropertyDefinition comes in with name == the @a name you
     * give here and implementation == "DEFAULT" or "INTERNAL", then its
     * parallelProperty() factory method will ignore the rest of the definition
     * and call instead your @a builtinConstructor.
     *
     * If implementation == "FLEXIBLE", however, then the flexible version will
     * be returned even if you've registered a builtin version.
     */

    static void registerParallelBuiltin(
        const string& name,
        IScalarRiskPropertyConstSP (*builtin)(IArrayConstSP indices,
                                              IArrayConstSP arguments));

    /**
     * Register a hard-coded property implementation to override a flexible one
     * (term-structured version)
     *
     * If a FieldRiskPropertyDefinition comes in with name == the @a name you
     * give here and implementation == "DEFAULT" or "INTERNAL", then its
     * pointwiseProperty() factory method will ignore the rest of the
     * definition and return instead your @a builtin
     * IRiskProperty<ExpiryWindow>.
     *
     * If implementation == "FLEXIBLE", however, then the flexible version will
     * be returned even if you've registered a builtin version.
     */

    static void registerPointwiseBuiltin(
        const string& name,
        IExpiryRiskPropertyConstSP (*builtin)(IObjectConstSP argument));

    /**
     * Register a hard-coded property implementation to override a flexible one
     * (term-structured version)
     *
     * If a FieldRiskPropertyDefinition comes in with name == the @a name you
     * give here and implementation == "DEFAULT" or "INTERNAL", then its
     * pointwiseProperty() factory method will ignore the rest of the
     * definition and return instead your @a builtin
     * IRiskProperty<ExpiryAndStrike>.
     *
     * If implementation == "FLEXIBLE", however, then the flexible version will
     * be returned even if you've registered a builtin version.
     */

    static void registerExpiryAndStrikewiseBuiltin(
        const string& name,
        IExpiryAndStrikeRiskPropertyConstSP (*builtin)(IObjectConstSP argument));

    /**
     * Register a hard-coded property implementation to override a flexible one
     * (version for changes to multiple array elements at once)
     *
     * If a FieldRiskPropertyDefinition comes in with name == the @a name you
     * give here and implementation == "DEFAULT" or "INTERNAL", then its
     * elementsProperty() factory method will ignore the rest of the definition
     * and return instead your @a builtin IScalarRiskProperty.
     *
     * If implementation == "FLEXIBLE", however, then the flexible version will
     * be returned even if you've registered a builtin version.
     */

    static void registerElementsBuiltin(
        const string& name,
        IScalarRiskPropertyConstSP (*builtin)(IntArrayConstSP indices,
                                              IArrayConstSP arguments));
    /**
     * Register a hard-coded property implementation to override a flexible one
     * (version for changes to each array element in turn)
     *
     * If a FieldRiskPropertyDefinition comes in with name == the @a name you
     * give here and implementation == "DEFAULT" or "INTERNAL", then its
     * elementwiseProperty() factory method will ignore the rest of the definition
     * and return instead your @a builtin IRiskProperty<BoxedInt>.
     *
     * If implementation == "FLEXIBLE", however, then the flexible version will
     * be returned even if you've registered a builtin version.
     */

    static void registerElementwiseBuiltin(
        const string& name,
        IIntRiskPropertyConstSP (*builtin)(IObjectConstSP argument));

};

DRLIB_END_NAMESPACE

#endif
