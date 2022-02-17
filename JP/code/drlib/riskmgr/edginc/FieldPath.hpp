/**
 * @file FieldPath.hpp
 */

#ifndef QLIB_FieldPath_H
#define QLIB_FieldPath_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(FieldPath)

/**
 * A chain of field names which can be used to get and set values of fields
 * attached to a given object either directly, or indirectly through a
 * series of intermediate objects
 *
 * For instance, if p is FieldPath([ "a", "b", "c" ]) then
 *
 *    -  p.get(x) is the value of the [QLib-registered] field "c" on the
 *       object which is the value of the field "b" on the object which is
 *       the value of the field "a" on x.
 *
 *    -  p.set(x, v) sets the value of that "c" field to v
 *
 * Fields typed as MarketObjectWrapper are transparently resolved to the
 * referenced MarketObject's, for both reading and writing (and for
 * intermediate values in the chain).
 *
 * FieldPath is a deliberately defined in terms of untyped field names,
 * mapped to actual CField's on an object-by-object basis, rather than
 * more "statically" in terms of particular classes and fields specified
 * up front.  The latter is cleaner but means you can't follow links
 * typed as abstract interfaces even if you know in practice what the
 * object's concrete type is going to be at run-time.
 *
 * This class is used in the "flexible scenarios and greeks" framework: see
 * FieldRiskProperty / FieldRiskAxis / FieldTweakHypothesis, and (for
 * background info) FlexibleSensitivity.
 */

class RISKMGR_DLL FieldPath: public CObject {

    FieldPath();
    static IObject* defaultOne();
    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

private:

    CStringArrayConstSP _path;

public:

    /**
     * Constructor
     *
     * @param path     Chain of field names to followed from any given "anchor"
     *                 object.  May be empty, in which case the chain has
     *                 length zero and the result of calling get() on an object
     *                 is the object itself [but an exception is thrown if you
     *                 call set() or objectChain() with an empty path].
     */

    FieldPath(CStringArrayConstSP path);

    /**
     * Constructor returning a smartPtr.
     */

    static FieldPathSP SP(CStringArrayConstSP path);

    /**
     * Constructor for a length-1 chain, returning a smartPtr.
     */

    static FieldPathSP SP(const string& name);

    /**
     * Constructor for a length-0 chain, returning a smartPtr.
     */

    static FieldPathSP SP();

    ~FieldPath();

    /**
     * The names of the fields in the chain
     */

    CStringArrayConstSP path() const;

    /**
     * The value of the field at the end of the chain, starting from a given
     * anchor object
     */

    IObjectConstSP get(IObjectConstSP root, CClassConstSP type = 0) const;

    /**
     * The value of the field at the end of the chain, starting from a given
     * anchor object (non-const version)
     */

    IObjectSP get(IObjectSP root, CClassConstSP type = 0) const;

    /**
     * Return value of objectChain()
     */

    struct ObjectChain {
        ObjectArraySP objects;
        CFieldArray fields;

        void fieldsUpdated() const;

        ObjectChain();
        ~ObjectChain();
    };

    /**
     * The chain of fields and values followed by our series of field names
     * starting from a given object
     *
     * In the returned ObjectChain,
     *
     *   -  objects[0] is @a root
     *   -  fields[0] is the field on @a root's class or one of its ancestor
     *      classes called path()[0]
     *
     *   -  objects[1] is the value of fields[0] on objects[0]
     *   -  fields[1] is the field field on objects[1]'s class or one
     *      of its ancestors called path()[1]
     *        
     * and so on.  Note that 'objects' is one element longer than 'fields',
     * and objects.back() is the leaf field value.
     *
     * For convenience, the returned chain is guaranteed to contain at least
     * one field (and hence two objects, @a root and its field value).
     * An exception is thrown is path() is of zero length.
     */

    ObjectChain objectChain(IObjectSP root) const;

    /**
     * Set the value of the field at the end of the chain, starting from a
     * given anchor object
     */

    void set(IObjectSP root, IObjectSP value) const;

    /**
     * String representation of the FieldPath
     */

    string toString(int highlight = -1) const;

    /**
     * Message suitable for use in error messages arising from dealing with
     * the chain w.r.t. a given anchor object
     */

    string contextMessage(IObjectConstSP obj, int highlight = -1) const;
};

DRLIB_END_NAMESPACE

#endif
