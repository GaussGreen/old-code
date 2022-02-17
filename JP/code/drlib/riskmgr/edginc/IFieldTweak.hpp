/**
 * @file IFieldTweak.hpp
 */

#ifndef QLIB_IFieldTweak_H
#define QLIB_IFieldTweak_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(FieldPath)
FORWARD_DECLARE(TweakFunction)
FORWARD_DECLARE(IFieldTweak)
class Range;

FORWARD_DECLARE(RiskAxis)
FORWARD_DECLARE(OutputName)
FORWARD_DECLARE(IAbstractRiskProperty)

/**
 * An operator for "tweaking" an object (for calculating greeks or scenarios)
 * by adjusting the value of a named field
 *
 * This is the core of the "flexible scenarios" system (see
 * FlexibleSensitivity, ScalarFlexiblePerturbation,
 * PerEntryFlexiblePerturbation) --- the method tweak() is used by
 * FieldTweakHypothesis to make the necessary edits to the objects it wants to
 * zap.
 *
 * The methods you should use are
 *
 *    -  bulk() for a scalar operator (or parallel array operator)
 *    -  elementwise() for a selective array operator
 *
 *    -  IOperator::numeric() for tweaking int or double (array) fields
 *    -  IOperator::setter() for changing the value of non-numeric fields
 *
 *    -  IIndices::direct() for selecting array elements literally by index
 *    -  IIndices::lookup() for selecting array indices indirectly by "key"
 *       (e.g. Expiry)
 * 
 * It's quite tightly designed but looks complicated because it supports the
 * following features:
 *
 *    -  Array fields can be tweaked all-at-once or per-entry.
 *
 *    -  Entries in array fields may be objects of which it's actually a
 *       (sub-)field which gets tweaked.
 *
 *    -  Array fields can be addressed either directly by index,
 *
 *    -  or via "keys" --- typically Expiry's --- listed in some other array,
 *
 *    -  or in a different sub-field of the same array,
 *
 *    -  or in a field registered with the Calibrator.
 *
 *    -  The tweak can be a numeric operation (+, *, etc.),
 *
 *    -  or a straight replacement of a non-numeric value,
 *
 *    -  or replacement with a value constructed from literal XML with optional
 *       interpolation of a user-supplied argument (!)
 */

class RISKMGR_DLL IFieldTweak: public virtual IObject {
public:

    static CClassConstSP const TYPE;

    IFieldTweak(); ~IFieldTweak();

    /**
     * Tweak some field of an object, returning the absolute amount by which it
     * was changed
     */

    virtual double apply(IObjectSP object) const = 0;

    /**
     * A tweak which edits an object the same way we do but @a coeff times more
     */

    virtual IFieldTweakConstSP scaled(double coeff) const = 0;

    /**
     * Whether this->scaled(0) is always a "no-op"
     */

    virtual bool zeroIsNoop() const = 0;

    /**
     * The keys (typically ExpiryWindow's) specifying the elements we edit
     * in the (array-valued) field which we will tweak; or NULL if we are
     * bulk() tweak
     */

    virtual IArrayConstSP qualifiers() const = 0;

    FORWARD_DECLARE(IIndices)

    /**
     * Different ways of specifying the elements to be tweaked in an
     * array-valued field
     */

    class RISKMGR_DLL IIndices: public virtual IObject {
    public:

        static CClassConstSP const TYPE;

        IIndices(); ~IIndices();

        /**
         * The indices to be tweaked of a field on a given object
         */

        virtual IntArrayConstSP indices(IObjectConstSP object) const = 0;

        /**
         * The keys specifying the elements to be tweaked
         */

        virtual IArrayConstSP qualifiers() const = 0;

        /**
         * Indices specified directly
         */

        static IIndicesConstSP direct(IntArrayConstSP indices);

        /**
         * Indices (typically ExpiryWindow's) specified indirectly
         *
         * The indices that get tweaked on a given object x are
         *
         * -     { i : ((x.keysField)[i]).keyField is in keys }
         */

        static IIndicesConstSP lookup(FieldPathConstSP keysField,
                                      FieldPathConstSP keyField,
                                      IArrayConstSP keys);

        /**
         * Indices (typically ExpiryWindow's) specified indirectly, possibly
         * via an 'expiries' field registered with the Calibrator
         *
         * If @a keysField is given: as lookup()
         *
         * If not: @a field must have been registered as having a
         * <tt>getExpiriesMethod</tt> via
         * Calibrator::IAdjustable::registerBootstrappableField(), and the
         * lookup keys are obtained from that method
         */

        static IIndicesConstSP lookupOrCalibratorExpiries(
            FieldPathConstSP field,
            FieldPathConstSP keysField,
            FieldPathConstSP keyField,
            IArrayConstSP indices);
    };

    FORWARD_DECLARE(IOperator)

    /**
     * Different kinds of tweak
     *
     * This is the function used to compute a new value for a field from its
     * current value and a user-supplied argument.
     */

    class RISKMGR_DLL IOperator: public virtual IObject {
    public:

        static CClassConstSP const TYPE;

        /**
         * The new value for a @a field, given its current @a value and a
         * user-supplied @a argument
         */

        virtual IObjectConstSP operator ()(IObjectConstSP argment,
                                           IObjectConstSP value,
                                           CFieldConstSP field) const = 0;

        /**
         * Whether supplying 0 as operator()()'s @a argument makes it a no-op
         */

        virtual bool zeroIsNoop() const = 0;

        IOperator(); ~IOperator();

        /**
         * A numeric operation (+, *, etc.) --- just an adaptor for a
         * TweakFunction
         */

        static IOperatorConstSP numeric(TweakFunctionConstSP function,
                                        const Range& range,
                                        bool useCalibratorRange,
                                        bool clip);

        /**
         * A polymorphic (non-numeric) operator which just
         * sets the field's value to the user-supplied 'argument'
         *
         * For int- or double-valued fields you're better of using numeric(),
         * with TweakFunction::setter(), since then you'll get range checking.
         */

        static IOperatorConstSP setter();

        /**
         * A polymorphic (non-numeric) operator which sets the field's value to
         * an object constructed from a literal XML string, optionally
         * interpolated with a user-supplied argument
         *
         * E.g. if you give @a xml as
         *
         * <pre>
         *   <dependenceMaker TYPE='DependenceMaker@'></dependenceMaker>
         * </pre>
         *
         * (note the @ sign) and 'argument' is given as <tt>Gauss</tt> then the
         * new value will be a freshly-built DependenceMakerGauss.  See
         * <tt>flexiblegreeksinp/CorrSwapBasis.xml</tt>.
         */

        static IOperatorConstSP xmlSetter(const string& xml);
    };

    /**
     * An operator for tweaking a named scalar field of an object (or tweaking
     * all elements of an array field in parallel)
     *
     * @param field          FieldPath naming the field to be tweaked
     *
     * @param subField       Optionally, if @a field is array-valued: the
     *                       field within each array element to be tweaked
     *
     * @param argument       The new value is obtained by passing the old
     *                       value along with @a argument to @a op's
     *                       IFieldTweak::IOperator::operator()() method
     *
     * @param op             IFieldTweak::IOperator which gives the field's
     *                       new value
     */

    static IFieldTweakConstSP bulk(FieldPathConstSP field,
                                   FieldPathConstSP subField,
                                   IObjectConstSP argument,
                                   IOperatorConstSP op);

    /**
     * An operator for individually tweaking the elements of a named array
     * field of an object
     *
     * @param field          FieldPath naming the field to be tweaked
     *
     * @param subField       Optionally, if @a field is array-valued: the
     *                       field within each array element to be tweaked
     *
     * @param argument       The new value is obtained by passing the old
     *                       value along with @a argument to @a op's
     *                       IFieldTweak::IOperator::operator()() method
     *
     * @param op             IFieldTweak::IOperator which gives the field's
     *                       new value
     */

    static IFieldTweakConstSP elementwise(FieldPathConstSP field,
                                          FieldPathConstSP subField,
                                          IIndicesConstSP indices,
                                          int indexDimension,
                                          IArrayConstSP arguments,
                                          IOperatorConstSP op);

    /**
     * How a IRiskAxis defined via this 'tweak' is represented externally
     * (we used to have a simpler "flexible greeks" scheme which was used
     * for RiskMapping and hence leaked out, and I'm not yet game to get
     * the Pyramid guys to change their datamodel)
     */

    RiskAxisSP legacyRepresentation(IAbstractRiskPropertyConstSP property,
                                    OutputNameConstSP name,
                                    bool absoluteDistance) const;
};

DRLIB_END_NAMESPACE

#endif
