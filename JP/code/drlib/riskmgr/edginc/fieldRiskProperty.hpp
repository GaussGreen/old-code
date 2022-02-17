/**
 * @file fieldRiskProperty.hpp
 */

#ifndef QLIB_fieldRiskProperty_H
#define QLIB_fieldRiskProperty_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IRiskProperty.hpp"
#include "edginc/IFieldTweak.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(FieldPath)

/**
 * Constructors for IRiskProperty's corresponding to a named field of some
 * named class
 *
 * These functions allow you to construct IRiskProperty's --- <i>i.e.</i> ways
 * of "tweaking" the world for greeks or scenarios --- dynamically from class
 * and field names.  They are used to implement the "flexible greeks and
 * scenarios" facility (see FlexibleSensitivity, ScalarFlexiblePerturbation,
 * PerEntryFlexiblePerturbation for more info).
 */

struct RISKMGR_DLL fieldRiskProperty {

/**
 * An IRiskProperty under which instances of a specified class have a named
 * field tweaked ("in parallel", if it's an array-valued field)
 *
 * You can use this to create a "scalar" greek
 * (PerNameRiskPropertySensitivity), an AllNamesRiskPropertySensitivity,
 * etc. without having to program any data classes to be
 * ITweakableWithRespectTo some hardwired property.  It's also how scalar
 * FlexibleSensitivity's and ScalarFlexiblePerturbation's are implemented (via
 * IFieldRiskPropertyDefinition::scalarProperty()).
 *
 * @param clazz             The type of object on which the property is defined
 * @param tweak             An IFieldTweak saying how instances of @a clazz get
 *                          tweaked
 * @param absoluteDistance  When calculating greeks, 'true' means divide the
 *                          value difference by the nominal shift size, while
 *                          'false' means divide by the actual change in the
 *                          property (e.g. the former may be in relative terms
 *                          while the latter is always absolute). See
 *                          also ITweakableWithRespectTo<TAG>::sensShift().
 */

static IScalarRiskPropertyConstSP scalar(CClassConstSP clazz,
                                         IFieldTweakConstSP tweak,
                                         bool absoluteDistance);

/**
 * An IRiskProperty under which instances of a specified class have a named
 * array-valued field tweaked once for each entry
 *
 * You can use this to create a "elementwise" greek
 * (PerNameRiskPropertySensitivity) without having to program any data classes
 * to be ITweakableWithRespectTo some hardwired property.  It's also how
 * elementwise FlexibleSensitivity's and PerEntryFlexiblePerturbation's are
 * implemented (via IFieldRiskPropertyDefinition::elementwiseProperty()).
 *
 * @param clazz             The type of object on which the property is defined
 * @param field             The array field (on @a clazz) holding the values to
 *                          be tweaked
 * @param subField          Optionally, the field within each value that
 *                          gets tweaked (if they're objects)
 * @param indexDimension    For matrix-valued fields, 0 means tweak whole
 *                          columns at a time, 1 means tweak rows
 * @param op                The tweaking operator to be applied
 * @param arg               Argument (typically shift size or new value) to
 *                          be passed to @a op along with old value
 * @param absoluteDistance  When calculating greeks, 'true' means divide the
 *                          value difference by the nominal shift size, while
 *                          'false' means divide by the actual change in the
 *                          property (e.g. the former may be in relative terms
 *                          while the latter is always absolute). See
 *                          also ITweakableWithRespectTo<TAG>::sensShift().
 */

static IIntRiskPropertyConstSP elementwise(CClassConstSP clazz,
                                           FieldPathConstSP field,
                                           FieldPathConstSP subField,
                                           int indexDimension,
                                           IFieldTweak::IOperatorConstSP op,
                                           IObjectConstSP arg,
                                           bool absoluteDistance);

/**
 * An IRiskProperty under which instances of a specified class have a named
 * term-structured field tweaked once for each relevant tenor
 *
 * You can use this to create a "pointwise" greek
 * (PerNameRiskPropertySensitivity) without having to program any data classes
 * to be ITweakableWithRespectTo some hardwired property.  It's also how
 * pointwise FlexibleSensitivity's and PerEntryFlexiblePerturbation's are
 * implemented (via IFieldRiskPropertyDefinition::pointwiseProperty()).
 *
 * @param clazz             The type of object on which the property is defined
 * @param field             The array field (on @a clazz) holding the values to
 *                          be tweaked
 * @param subField          Optionally, the field within each value that
 *                          gets tweaked (if they're objects)
 * @param keyListField      The array field holding a corresponding list of
 *                          Expiry's --- or of some object type of which ...
 * @param keyField          ... this names the field holding the actual Expiry
 * @param indexDimension    For matrix-valued fields, 0 means tweak whole
 *                          columns at a time, 1 means tweak rows
 * @param op                The tweaking operator to be applied
 * @param arg               Argument (typically shift size or new value) to
 *                          be passed to @a op along with old value
 * @param absoluteDistance  When calculating greeks, 'true' means divide the
 *                          value difference by the nominal shift size, while
 *                          'false' means divide by the actual change in the
 *                          property (e.g. the former may be in relative terms
 *                          while the latter is always absolute). See
 *                          also ITweakableWithRespectTo<TAG>::sensShift().
 */

static IExpiryRiskPropertyConstSP pointwise(CClassConstSP clazz,
                                            FieldPathConstSP field,
                                            FieldPathConstSP subField,
                                            FieldPathConstSP keyListField,
                                            FieldPathConstSP keyField,
                                            int indexDimension,
                                            IFieldTweak::IOperatorConstSP op,
                                            IObjectConstSP arg,
                                            bool absoluteDistance);

/**
 * An IRiskProperty under which instances of a specified class have a named
 * doubly-term-structured field tweaked once for each relevant combination of
 * tenors
 *
 * You can use this to create an ExpiryPair-wise greek
 * (PerNameRiskPropertySensitivity) without having to program any data classes
 * to be ITweakableWithRespectTo some hardwired property.
 *
 * @param clazz             The type of object on which the property is defined
 * @param field             The array field (on @a clazz) holding the values to
 *                          be tweaked
 * @param subField          Optionally, the field within each value that
 *                          gets tweaked (if they're objects)
 * @param keyListField      The array field holding a corresponding list of
 *                          ExpiryPair's --- or of some object type of which ...
 * @param keyField          ... this names the field holding the actual ExpiryPair
 * @param indexDimension    For matrix-valued fields, 0 means tweak whole
 *                          columns at a time, 1 means tweak rows
 * @param op                The tweaking operator to be applied
 * @param arg               Argument (typically shift size or new value) to
 *                          be passed to @a op along with old value
 * @param absoluteDistance  When calculating greeks, 'true' means divide the
 *                          value difference by the nominal shift size, while
 *                          'false' means divide by the actual change in the
 *                          property (e.g. the former may be in relative terms
 *                          while the latter is always absolute). See
 *                          also ITweakableWithRespectTo<TAG>::sensShift().
 */

static IExpiryPairRiskPropertyConstSP expiryPairwise(
    CClassConstSP clazz,
    FieldPathConstSP field,
    FieldPathConstSP subField,
    FieldPathConstSP keyListField,
    FieldPathConstSP keyField,
    int indexDimension,
    IFieldTweak::IOperatorConstSP op,
    IObjectConstSP arg,
    bool absoluteDistance);

static IExpiryAndStrikeRiskPropertyConstSP expiryAndStrikewise(
    CClassConstSP clazz,
    FieldPathConstSP field,
    FieldPathConstSP subField,
    FieldPathConstSP keyListField,
    FieldPathConstSP keyField,
    int indexDimension,
    IFieldTweak::IOperatorConstSP op,
    IObjectConstSP arg,
    bool absoluteDistance);

};

DRLIB_END_NAMESPACE

#endif
