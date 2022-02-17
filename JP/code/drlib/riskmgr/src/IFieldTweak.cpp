/**
 * @file IFieldTweak.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/Range.hpp"
#include "edginc/Format.hpp"
#include "edginc/XMLReader.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/ExpiryWindow.hpp"
#include "edginc/FieldPath.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/TweakFunction.hpp"
#include "edginc/RiskAxis.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/IFieldTweak.hpp"
#include "edginc/fieldRiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

// 
// **********
//  IIndices
// **********
// 

IFieldTweak::IIndices::IIndices() {}
IFieldTweak::IIndices::~IIndices() {}

// 
// ========
//  direct
// ========
// 

struct DirectIndices: CObject,
                      virtual IFieldTweak::IIndices {

    static CClassConstSP const TYPE;
    static IObject* emptyShell(); static void load(CClassSP&);

    IntArrayConstSP _indices; // $required

    DirectIndices(IntArrayConstSP indices = IntArrayConstSP()):
        CObject(TYPE),
        _indices(indices)
    {
        ASSERT(!!indices);
    }

    IntArrayConstSP indices(IObjectConstSP) const {
        return _indices;
    }

    IArrayConstSP qualifiers() const {
        return _indices;
    }
};

IFieldTweak::IIndicesConstSP IFieldTweak::IIndices::direct(IntArrayConstSP indices) {
    return IIndicesConstSP(new DirectIndices(indices));
}

// 
// ========
//  lookup
// ========
// 

struct LookupIndices: CObject,
                      virtual IFieldTweak::IIndices {

    static CClassConstSP const TYPE;
    static IObject* emptyShell(); static void load(CClassSP&);

    FieldPathConstSP keyListField;  // $required
    FieldPathConstSP keyField;      // $required
    IArrayConstSP indexKeys;        // $required

    bool matches(IObjectConstSP keyA, IObjectConstSP keyB) const {
        return !keyA ? !keyB :
               (keyA->equalTo(keyB.get()) ||
                    ExpiryWindow::TYPE->isInstance(keyA) &&
                    ExpiryWindowConstSP::dynamicCast(keyA)->expiry->equalTo(
                                                                keyB.get()));
    }

    IntArrayConstSP indices(IObjectConstSP root) const {
        try {
            IntArraySP indices(new IntArray(indexKeys->getLength()));

            IArrayConstSP keyList;

            IObjectConstSP v = keyListField->get(root);
            try {
                keyList = IArrayConstSP::dynamicCast(v);
            }
            catch (exception& e) {
                throw ModelException(e, "Field is not an array");
            }

            for (int k = 0; k < indexKeys->getLength(); ++k) {
                IObjectConstSP key = indexKeys->get(k);
                try {
                    int i;
                    for (i = 0; i < keyList->getLength(); ++i) {
                        IObjectConstSP keyI = keyField->get(keyList->get(i));
                        if (matches(key, keyI)) {
                            (*indices)[k] = i;
                            break;
                        }
                    }

                    if (i == keyList->getLength()) {
                        throw ModelException("Not found");
                    }
                }
                catch (exception& e) {
                    throw ModelException(e,
                        "Index " + (!key ? string("NULL") : key->toString()) +
                        " (#" + Format::toString(k) + ")");
                }
            }

            return indices;
        }
        catch (exception& e) {
            throw ModelException(e, __FUNCTION__,
                                 "Looking up indices in " +
                                 keyListField->contextMessage(root, -1));
        }
     }

    IArrayConstSP qualifiers() const {
        return indexKeys;
    }

    LookupIndices(FieldPathConstSP keyListField,
                  FieldPathConstSP keyField,
                  IArrayConstSP indexKeys):
        CObject(TYPE),
        keyListField(keyListField),
        keyField(!keyField ? FieldPath::SP() : keyField),
        indexKeys(indexKeys)
    {
        ASSERT(!!keyListField && !keyListField->path()->empty());
        ASSERT(!!indexKeys);
    }

    LookupIndices(): CObject(TYPE) {}
};

IFieldTweak::IIndicesConstSP IFieldTweak::IIndices::lookup(
        FieldPathConstSP keyListField,
        FieldPathConstSP keyField,
        IArrayConstSP indexKeys) {
    return IIndicesConstSP(new LookupIndices(keyListField, keyField,
                                             indexKeys));
}

// 
// ============================
//  lookupOrCalibratorExpiries
// ============================
// 

struct CalibratorExpiriesIndices: CObject,
                                  virtual IFieldTweak::IIndices {

    static CClassConstSP const TYPE;
    static IObject* emptyShell(); static void load(CClassSP&);

    FieldPathConstSP field;                    // $required
    ExpiryWindowArrayConstSP indexKeys;        // $required

    IntArrayConstSP indices(IObjectConstSP root) const {
        try {
            ExpiryArrayConstSP expiries;

            try {
                FieldPath::ObjectChain chain =
                    field->objectChain(IObjectSP::constCast(root));
                expiries = Calibrator::IAdjustable::
                    getGetExpiriesMethod(chain.fields.back())(
                        (*chain.objects)[chain.fields.size() - 1].get());
            }
            catch (exception& e) {
                throw ModelException(e,
                    "Looking up Calibrator-registered expiries",
                    "Did you forget to specify an explicit expiryQualifier "
                    "in your FieldRiskPropertyDefinition?");
            }

            if (!expiries) {
                throw ModelException("Expiries are NULL");
            }

            IntArraySP indices(new IntArray(indexKeys->getLength()));

            for (int k = 0; k < indexKeys->getLength(); ++k) {
                (*indices)[k] = (*indexKeys)[k]->expiry->search(expiries.get());
            }

            return indices;
        }
        catch (exception& e) {
            throw ModelException(e, __FUNCTION__,
                                 "Looking up Calibrator-registered expiries "
                                 "for " + field->contextMessage(root, -1));
        }
     }

    IArrayConstSP qualifiers() const {
        return indexKeys;
    }

    CalibratorExpiriesIndices(FieldPathConstSP field,
                              ExpiryWindowArrayConstSP indexKeys):
        CObject(TYPE),
        field(field),
        indexKeys(indexKeys)
    {
        ASSERT(!!field && !field->path()->empty());
        ASSERT(!!indexKeys);
    }

    CalibratorExpiriesIndices(): CObject(TYPE) {}
};

IFieldTweak::IIndicesConstSP IFieldTweak::IIndices::lookupOrCalibratorExpiries(
        FieldPathConstSP field,
        FieldPathConstSP keysField,
        FieldPathConstSP keyField,
        IArrayConstSP indexKeys) {
    if (!!keysField) {
        return IIndicesConstSP(new LookupIndices(keysField, keyField, indexKeys));
    }
    else {
        ExpiryWindowArraySP exps(new ExpiryWindowArray(indexKeys->getLength()));
        for (int i = 0; i < exps->size(); ++i) {
            IObjectSP k = IObjectSP::constCast(indexKeys->get(i));
            try {
                if (!k) throw ModelException("Value is NULL");
                exps->set(i, Expiry::TYPE->isInstance(k) ?
                                 ExpiryWindow::SP(ExpirySP(),
                                                  ExpirySP::dynamicCast(k),
                                                  ExpirySP()) :
                                 ExpiryWindowSP::dynamicCast(k));
            }
            catch (exception& e) {
                throw ModelException(e, __FUNCTION__,
                    "You have not specified a 'qualifier' field, "
                    "so I am trying to fall back on any Expiry array "
                    "fields registered with the Calibrator for " +
                    field->toString() + ", but the supplied qual #" +
                    Format::toString(i) + " is not an Expiry but rather " +
                    (!k ? string("NULL") : k->getClass()->getName() +
                                           " \"" + k->toString() + "\""));
            }
        }

        return IIndicesConstSP(new CalibratorExpiriesIndices(field, exps));
    }
}

// 
// ***********
//  IOperator
// ***********
// 

IFieldTweak::IOperator::IOperator() {}
IFieldTweak::IOperator::~IOperator() {}

// 
// --------------------------------
//  convertedTo, convertedToDouble
// --------------------------------
// 

static const double notDouble = -1.234e56;

static double convertedToDouble(IObjectConstSP x, const string& desc) {
    const CDouble* d;
    const CInt* i;

    if (!!(d = dynamic_cast<const CDouble*>(x.get()))) {
        return d->doubleValue();
    }
    else if (!!(i = dynamic_cast<const CInt*>(x.get()))) {
        return i->intValue();
    }
    else if (desc.empty()) {
        return notDouble;
    }
    else {
        throw ModelException(
            __FUNCTION__,
            "Can't convert " + desc + " " +
              (!x ? string("NULL") :
                    x->getClass()->getName() + " \"" + x->toString() + "\"") +
              " to double");
    }
}

static int convertedToInt(IObjectConstSP x, const string& desc) {
    const CDouble* d;
    const CInt* i;

    if (!!(d = dynamic_cast<const CDouble*>(x.get()))) {
        return int(d->doubleValue() + (d->doubleValue() > 0 ? .5 : -.5));
    }
    else if (!!(i = dynamic_cast<const CInt*>(x.get()))) {
        return i->intValue();
    }
    else {
        throw ModelException(
            __FUNCTION__,
            "Can't convert " + desc + " " +
              (!x ? string("NULL") :
                    x->getClass()->getName() + " \"" + x->toString() + "\"") +
              " to int");
    }
}

static IObjectConstSP convertedTo(CClassConstSP type,
                             IObjectConstSP x,
                             const string& desc) {
    if (!x || type->isInstance(x)) {
        return x;
    }
    else if (type == CInt::TYPE) {
        return IObjectConstSP(CInt::create(convertedToInt(x, desc)));
    }
    else if (type == CDouble::TYPE) {
        return CDouble::SP(convertedToDouble(x, desc));
    }
    else {
        throw ModelException(
            "Don't know how to convert " + desc + " " +
             x->getClass()->getName() + " \"" + x->toString() + "\" "
            "to " + type->getName());
    }
}

// 
// ---------------------
//  numericOrEqDistance
// ---------------------
// 

static double numericOrEqDistance(IObjectConstSP a, IObjectConstSP b) {
    if (!a || !b) {
        return (!a && !b) ? 0. : 1.;
    }

    double da = convertedToDouble(a, "");
    if (da != notDouble) {
        double db = convertedToDouble(b, "");
        if (db != notDouble) {
            return db - da;
        }
    }

    return a->equalTo(b.get()) ? 0. : 1.;
}

// 
// =========
//  numeric
// =========
// 

struct NumericOperator: CObject,
                        virtual IFieldTweak::IOperator {

    static CClassConstSP const TYPE;
    static IObject* emptyShell(); static void load(CClassSP&);

    TweakFunctionConstSP function;   // $required

    const bool hasFloor;             // $required
    const double floor;              // $optional
    const bool floorInclusive;       // $optional

    const bool hasCap;               // $required
    const double cap;                // $optional
    const bool capInclusive;         // $optional

    const bool useCalibratorRange;   // $required

    const bool clip;                 // $optional

    mutable CFieldConstSP _field;    // $unregistered
    mutable Range _range;            // $unregistered

    Range range(CFieldConstSP field) const {
        if (field != _field) {
            _range = Range(hasFloor ? floor : -HUGE_VAL, floorInclusive,
                           hasCap ? cap : HUGE_VAL, capInclusive);

            if (!!field && useCalibratorRange) {
                try {
                    _range = Calibrator::IAdjustable::getRange(field);
                }
                catch (...) {}
            }
        }

        return _range;
    }

    IObjectConstSP operator ()(IObjectConstSP argument,
                               IObjectConstSP x,
                               CFieldConstSP field) const {
        return CDouble::SP(
            (*function)(convertedToDouble(argument, "argument"),
                        convertedToDouble(x, "current value"),
                        range(field),
                        clip));
    }

    bool zeroIsNoop() const {
        return function->zeroIsNoop();
    }

    NumericOperator():
        CObject(TYPE),
        hasFloor(false),
        floor(0),
        floorInclusive(true),
        hasCap(false),
        cap(0),
        capInclusive(true),
        useCalibratorRange(true),
        clip(false),
        _field((CFieldConstSP)0xDEADBEEF),
        _range(InfiniteRange())
    {}

    NumericOperator(TweakFunctionConstSP function,
                    const Range& range,
                    bool useCalibratorRange,
                    bool clip):
        CObject(TYPE),
        function(function),
        hasFloor(!range.getLower().isInfinite()),
        floor(range.getLower().isInfinite() ? 0 : range.getLower().getValue()),
        floorInclusive(range.getLower().isClosedBracket()),
        hasCap(!range.getUpper().isInfinite()),
        cap(range.getUpper().isInfinite() ? 0 : range.getUpper().getValue()),
        capInclusive(range.getUpper().isClosedBracket()),
        useCalibratorRange(useCalibratorRange),
        clip(clip),
        _field((CFieldConstSP)0xDEADBEEF),
        _range(InfiniteRange())
    {}

    string toString() const {
        return function->toString();
    }
};

IFieldTweak::IOperatorConstSP IFieldTweak::IOperator::numeric(
        TweakFunctionConstSP function,
        const Range& range,
        bool useCalibratorRange,
        bool clip) {
    return IOperatorConstSP(new NumericOperator(function, range,
                                                useCalibratorRange, clip));
}

// 
// ========
//  setter
// ========
// 

struct SetterOperator: CObject,
                       virtual IFieldTweak::IOperator {

    static CClassConstSP const TYPE;
    static IObject* emptyShell(); static void load(CClassSP&);

    IObjectConstSP operator ()(IObjectConstSP argument, IObjectConstSP x,
                               CFieldConstSP field) const {
        return argument;
    }

    bool zeroIsNoop() const {
        return true;
    }

    SetterOperator():
        CObject(TYPE)
    {}

    string toString() const {
        return "=";
    }
};

IFieldTweak::IOperatorConstSP IFieldTweak::IOperator::setter() {
    return IOperatorConstSP(new SetterOperator());
}

// 
// ===========
//  xmlSetter
// ===========
// 

struct XMLSetterOperator: CObject,
                          virtual IFieldTweak::IOperator {

    static CClassConstSP const TYPE;
    static IObject* emptyShell(); static void load(CClassSP&);

    string xml; // $required

    IObjectConstSP operator ()(IObjectConstSP argument, IObjectConstSP x,
                               CFieldConstSP field) const {
        try {
            string argstring;

            if (!!argument && CString::TYPE->isInstance(argument)) {
                argstring = CStringConstSP::dynamicCast(argument)->
                                stringValue();
            }
            else {
                argstring = IObject::stringOf(argument);
            }

            string axml;

            for (size_t a = 0; a < xml.size();) {
                int b = a;
                a = xml.find('@', a);
                if (a == string::npos) a = xml.size();
                axml.append(xml, b, a - b);
                if (a != xml.size()) {
                    axml.append(argstring);
                    ++a;
                }
            }

            XMLReader reader(axml, false);
            return reader.read(reader.root());
        }
        catch (exception& e) {
            throw ModelException(e, __FUNCTION__,
                "Generating object from xml to replace existing value");
        }
    }

    bool zeroIsNoop() const {
        return true;
    }

    XMLSetterOperator(const string& xml):
        CObject(TYPE),
        xml(xml)
    {}

    string toString() const {
        return "=";
    }
};

IFieldTweak::IOperatorConstSP IFieldTweak::IOperator::xmlSetter(
        const string& xml) {
    return IOperatorConstSP(new XMLSetterOperator(xml));
}

// 
// *************
//  IFieldTweak
// *************
// 

IFieldTweak::IFieldTweak()
{}

IFieldTweak::~IFieldTweak() {}

// 
// ============
//  FieldTweak
// ============
// 

static string s(IObjectConstSP it) {
    return !it ? string("NULL") :
            CString::TYPE->isInstance(it) ? CStringConstSP::dynamicCast(it)->stringValue() :
            CDouble::TYPE->isInstance(it) ? Format::toString(CDoubleConstSP::dynamicCast(it)->doubleValue()) :
            CInt::TYPE->isInstance(it) ? Format::toString(CIntConstSP::dynamicCast(it)->intValue()) : it->toString();
}

struct FieldTweak: CObject,
                   virtual IFieldTweak {

    static CClassConstSP const TYPE;
    static void load(CClassSP&);

    FieldPathConstSP field;     // $required
    FieldPathConstSP subField;  // $required
    IOperatorConstSP op;        // $required
    int indexDimension;         // $required
    
    virtual IObjectConstSP argument(int i) const = 0;
    virtual int numArguments() const = 0;

    double applyArray(const FieldPath::ObjectChain& chain,
                      IntArrayConstSP inds) const;

    virtual double _apply(const FieldPath::ObjectChain& chain) const = 0;

    double apply(IObjectSP root) const {
        try {
            FieldPath::ObjectChain chain = field->objectChain(root);
            double distance = _apply(chain);
            chain.fieldsUpdated();
            return distance;
        }
        catch (exception& e) {
            throw ModelException(e, __FUNCTION__,
                                 "Tweaking " + field->contextMessage(root));
        }
    }

    bool zeroIsNoop() const {
        return op->zeroIsNoop();
    }

    FieldTweak(CClassConstSP type,
               FieldPathConstSP field,
               FieldPathConstSP subField,
               IOperatorConstSP op,
               int indexDimension):
        CObject(type),
        field(field),
        subField(!subField ? FieldPath::SP() : subField),
        op(op),
        indexDimension(indexDimension)
    {
        ASSERT(!!field && !field->path()->empty());
    }

    FieldTweak(CClassConstSP type):
        CObject(type),
        indexDimension(0)
    {}

    string toString() const {
        IArrayConstSP quals = qualifiers();
        int nargs = numArguments();

        ostringstream it;

        it << field->toString();

        if (!!quals && quals->getLength() != 0) {
            it << "[";
            for (int a = 0; a < quals->getLength(); ++a) {
                if (a) it << ", ";
                it << s(quals->get(a));
            }
            it << "]";
        }

        if (!subField->path()->empty()) {
            if (!quals || quals->getLength() == 0) it << "[*]";
            it << subField->toString();
        }

        it << " (" << op->toString() << " ";

        if (nargs == 1) {
            it << s(argument(0));
        }
        else {
            it << "[";
            for (int a = 0; a < nargs; ++a) {
                if (a) it << ", ";
                it << s(argument(a));
            }
            it << "]";
        }

        it << ")";

        return it.str();
    }
};

static IntArrayConstSP iota(int n) {
    IntArraySP it(new IntArray(n));
    for (int i = 0; i < n; ++i) (*it)[i] = i;
    return it;
}

static const double notArray = -1.234e56;

static IObjectSP maybeCloned(IObjectConstSP y) {
    return !y ? IObjectSP() :
            y->getRefCount() == 1 ? IObjectSP::constCast(y) :
           IObjectSP(y->clone());
}

double FieldTweak::applyArray(const FieldPath::ObjectChain& chain,
                              IntArrayConstSP reqInds) const {
    double dist = 0;
    int distTerms = 0;

    IObjectSP v = chain.objects->back();
    if (!v) {
        throw ModelException("Field's value is NULL");
    }

    CDoubleMatrix* m;
    IArray* a;
    if (!!(m = dynamic_cast<CDoubleMatrix*>(v.get()))) {

        // 
        // --------------
        //  DoubleMatrix
        // --------------
        // 

        if (!subField->path()->empty()) {
            throw ModelException(
                "Doesn't make sense to tweak field " +
                subField->toString() + " of the elements of a DoubleMatrix");
        }

        int d = indexDimension;
        int ne = d == 0 ? m->numCols() : m->numRows();
        int ns = d == 0 ? m->numRows() : m->numCols();

        IntArrayConstSP inds = !reqInds ? iota(ne) : reqInds;

        for (int j = 0; j < inds->size(); ++j) {
            int t = (*inds)[j];
            try {
                if (!(0 <= t && t < ne)) {
                    throw ModelException(
                        string("Matrix only has ") +
                        Format::toString(ne) +
                        (d == 0 ? "columns" : "rows"));
                }

                for (int k = ns - 1; k >= 0; --k) {
                    double& y = d == 0 ? (*m)[t][k] : (*m)[k][t];
                    double x = y;
                    y = convertedToDouble(
                            (*op)(argument(j), CDouble::SP(y),
                                  chain.fields.back()),
                            "result");

                    dist += y - x;
                    ++distTerms;
                }
            }
            catch (exception& e) {
                throw ModelException(e,
                string("Tweaking ") +
                (d == 0 ? "column " : "row ") + 
                Format::toString(t) + " of a matrix");
            }
        }
    }
    else if (!!(a = dynamic_cast<IArray*>(v.get()))) {

        IntArrayConstSP inds = !reqInds ? iota(a->getLength()) : reqInds;

        // 
        // --------------------------
        //  Array (inc. DoubleArray)
        // --------------------------
        // 

        for (int j = 0; j < inds->size(); ++j) {
            int i = (*inds)[j];
            try {
                if (!(0 <= i && i < a->getLength())) {
                    throw ModelException(
                        "Index is out of range: array is of length " + 
                        Format::toString(a->getLength()));
                }

                if (subField->path()->empty()) {
                    // Array element
                    
                    IObjectConstSP x = a->get(i);
                    IObjectConstSP y = convertedTo(
                        chain.fields.back()->getType()->getComponentType(),
                        (*op)(argument(j), x, chain.fields.back()),
                        "result");

                    a->set(i, maybeCloned(y));

                    double d = numericOrEqDistance(x, y);
                    dist += d;
                    ++distTerms;
                }
                else {
                    // Field of array element

                    try {
                        FieldPath::ObjectChain ch =
                            subField->objectChain(a->get(i));

                        IObjectConstSP y = convertedTo(ch.fields.back()->getType(),
                                                  (*op)(argument(j),
                                                        ch.objects->back(),
                                                        ch.fields.back()),
                                                  "result");

                        ch.fields.back()->set((*ch.objects)[ch.fields.size()-1],
                                              maybeCloned(y));
                        ch.fieldsUpdated();

                        dist += numericOrEqDistance(ch.objects->back(), y);
                        ++distTerms;
                    }
                    catch (exception& e) {
                        throw ModelException(e,
                         "Tweaking " +
                         subField->contextMessage(a->get(i)));
                    }
                }
            }
            catch (exception& e) {
                throw ModelException(e, "Tweaking element #" +
                                           Format::toString(i));
            }
        }
    }
    else {
        return notArray;
    }

    return dist == 0 ? 0. : dist / distTerms;
}

// 
// ======
//  bulk
// ======
//

struct BulkFieldTweak: FieldTweak {

    static CClassConstSP const TYPE;
    static IObject* emptyShell(); static void load(CClassSP&);

    IObjectConstSP _argument;   // $required

    IObjectConstSP argument(int i) const {
        return _argument;
    }

    int numArguments() const {
        return 1;
    }

    double _apply(const FieldPath::ObjectChain& chain) const {
        double distance = applyArray(chain, IntArrayConstSP());

        if (distance == notArray) {
            IObjectConstSP newVal = convertedTo(chain.fields.back()->getType(),
                                           (*op)(argument(0),
                                                 chain.objects->back(),
                                                 chain.fields.back()),
                                           "result");

            chain.fields.back()->set((*chain.objects)[chain.fields.size() - 1],
                                     maybeCloned(newVal));

            distance = numericOrEqDistance(chain.objects->back(), newVal);
        }

        return distance;
    }

    IFieldTweakConstSP scaled(double coeff) const {
         try {
             return coeff == 1. ?
                 IFieldTweakConstSP::attachToRef(this) :
                 IFieldTweakConstSP(new BulkFieldTweak(
                     field, subField,
                     CDouble::SP(convertedToDouble(_argument,
                                                   "argument") * coeff),
                     op));
         }
         catch (exception& e) {
             throw ModelException(e, __FUNCTION__,
                                  "Scaling argument");
         }
    }

    IArrayConstSP qualifiers() const {
        return IArrayConstSP();
    }

    BulkFieldTweak(FieldPathConstSP field,
                   FieldPathConstSP subField,
                   IObjectConstSP argument,
                   IOperatorConstSP op):
        FieldTweak(TYPE, field, subField, op, 0),
        _argument(argument)
    {}

    BulkFieldTweak(): FieldTweak(TYPE) {}
};

IFieldTweakConstSP IFieldTweak::bulk(FieldPathConstSP field,
                                     FieldPathConstSP subField,
                                     IObjectConstSP argument,
                                     IOperatorConstSP op) {
    return IFieldTweakConstSP(new BulkFieldTweak(
        field, subField, argument, op));
}

// 
// =============
//  elementwise
// =============
// 

struct ElementwiseFieldTweak: FieldTweak {

    static CClassConstSP const TYPE;
    static IObject* emptyShell(); static void load(CClassSP&);

    IIndicesConstSP indices;    // $required
    IArrayConstSP arguments;    // $required

    IObjectConstSP argument(int i) const {
        if (i < 0 || arguments->getLength() <= i) {
            throw ModelException("Number of arguments specified for tweak "
                                 "doesn't match number of indices (" +
                                 Format::toString(i) + " isn't in [0, " +
                                 Format::toString(arguments->getLength()) + "))");
        }

        return arguments->get(i);
    }

    int numArguments() const {
        return arguments->getLength();
    }

    double _apply(const FieldPath::ObjectChain& chain) const {
        double distance = applyArray(
           chain, indices->indices((*chain.objects)[0]));

        if (distance == notArray) {
            const IObject* v = chain.objects->back().get();
            throw ModelException(
                "Field's value is not an array but rather " +
                (!v ? string("NULL") : "\"" + v->toString() + "\""));
        }

        return distance;
    }

    IFieldTweakConstSP scaled(double coeff) const {
        if (coeff == 1.) return IFieldTweakConstSP::attachToRef(this);

        DoubleArraySP scaledArgs(new DoubleArray(arguments->getLength()));
        for (int i = 0; i < scaledArgs->size(); ++i) {
            try {
                (*scaledArgs)[i] = convertedToDouble(arguments->get(i),
                                                     "argument") * coeff;
            }
            catch (exception& e) {
                throw ModelException(e, __FUNCTION__,
                                     "Scaling argument #" + Format::toString(i));
            }
        }

        return IFieldTweakConstSP(new ElementwiseFieldTweak(
            field, subField, indices, indexDimension, scaledArgs, op));
    }

    IArrayConstSP qualifiers() const {
        return indices->qualifiers();
    }

    ElementwiseFieldTweak(FieldPathConstSP field,
                          FieldPathConstSP subField,
                          IIndicesConstSP indices,
                          int indexDimension,
                          IArrayConstSP arguments,
                          IOperatorConstSP op):
        FieldTweak(TYPE, field, subField, op, indexDimension),
        indices(indices),
        arguments(arguments)
    {
        ASSERT(!!indices);
        ASSERT(!!arguments);
    }

    ElementwiseFieldTweak(): FieldTweak(TYPE) {}
};

IFieldTweakConstSP IFieldTweak::elementwise(FieldPathConstSP field,
                                            FieldPathConstSP subField,
                                            IIndicesConstSP indices,
                                            int indexDimension,
                                            IArrayConstSP arguments,
                                            IOperatorConstSP op) {
    return IFieldTweakConstSP(new ElementwiseFieldTweak(
        field, subField, indices, indexDimension, arguments, op));
}






static void IFieldTweak_load(CClassSP& clazz) {
  REGISTER_INTERFACE(IFieldTweak, clazz);
  EXTENDS(IObject);
}

CClassConstSP const IFieldTweak::TYPE = CClass::registerInterfaceLoadMethod(
  "IFieldTweak", typeid(IFieldTweak), IFieldTweak_load);

void IFieldTweak_IOperator_load(CClassSP& clazz) {
  REGISTER_INTERFACE(IFieldTweak::IOperator, clazz);
  EXTENDS(IObject);
}

typedef IFieldTweak::IOperator IFieldTweak_IOperator;

CClassConstSP const IFieldTweak_IOperator::TYPE = CClass::registerInterfaceLoadMethod(
  "IFieldTweak::IOperator", typeid(IFieldTweak_IOperator), IFieldTweak_IOperator_load);

static void IFieldTweak_IIndices_load(CClassSP& clazz) {
  REGISTER_INTERFACE(IFieldTweak::IIndices, clazz);
  EXTENDS(IObject);
}

typedef IFieldTweak::IIndices IFieldTweak_IIndices;

CClassConstSP const IFieldTweak_IIndices::TYPE = CClass::registerInterfaceLoadMethod(
  "IFieldTweak::IIndices", typeid(IFieldTweak_IIndices), IFieldTweak_IIndices_load);



IObject* SetterOperator::emptyShell() {
  return new SetterOperator();
}

void SetterOperator::load(CClassSP& clazz) {
  REGISTER(SetterOperator, clazz);
  SUPERCLASS(CObject);
  IMPLEMENTS(IFieldTweak::IOperator);
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const SetterOperator::TYPE = CClass::registerClassLoadMethod(
  "SetterOperator", typeid(SetterOperator), SetterOperator::load);

IObject* XMLSetterOperator::emptyShell() {
  return new XMLSetterOperator("");
}

void XMLSetterOperator::load(CClassSP& clazz) {
  REGISTER(XMLSetterOperator, clazz);
  SUPERCLASS(CObject);
  IMPLEMENTS(IFieldTweak::IOperator);
  EMPTY_SHELL_METHOD(emptyShell);
  FIELD(xml, "xml");
}

CClassConstSP const XMLSetterOperator::TYPE = CClass::registerClassLoadMethod(
  "XMLSetterOperator", typeid(XMLSetterOperator), XMLSetterOperator::load);

void FieldTweak::load(CClassSP& clazz) {
  REGISTER(FieldTweak, clazz);
  SUPERCLASS(CObject);
  IMPLEMENTS(IFieldTweak);
  FIELD(field, "field");
  FIELD(subField, "subField");
  FIELD(op, "op");
  FIELD(indexDimension, "indexDimension");
}

CClassConstSP const FieldTweak::TYPE = CClass::registerClassLoadMethod(
  "FieldTweak", typeid(FieldTweak), FieldTweak::load);

IObject* NumericOperator::emptyShell() {
  return new NumericOperator();
}

void NumericOperator::load(CClassSP& clazz) {
  REGISTER(NumericOperator, clazz);
  SUPERCLASS(CObject);
  IMPLEMENTS(IFieldTweak_IOperator);
  FIELD(function, "function");
  FIELD(hasFloor, "hasFloor");
  FIELD(floor, "floor");
  FIELD_MAKE_OPTIONAL(floor);
  FIELD(floorInclusive, "floorInclusive");
  FIELD_MAKE_OPTIONAL(floorInclusive);
  FIELD(hasCap, "hasCap");
  FIELD(cap, "cap");
  FIELD_MAKE_OPTIONAL(cap);
  FIELD(capInclusive, "capInclusive");
  FIELD_MAKE_OPTIONAL(capInclusive);
  FIELD(useCalibratorRange, "useCalibratorRange");
  FIELD(clip, "clip");
  FIELD_MAKE_OPTIONAL(clip);
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const NumericOperator::TYPE = CClass::registerClassLoadMethod(
  "NumericOperator", typeid(NumericOperator), NumericOperator::load);

IObject* DirectIndices::emptyShell() {
  return new DirectIndices();
}

void DirectIndices::load(CClassSP& clazz) {
  REGISTER(DirectIndices, clazz);
  SUPERCLASS(CObject);
  IMPLEMENTS(IFieldTweak::IIndices);
  FIELD(_indices, "_indices");
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const DirectIndices::TYPE = CClass::registerClassLoadMethod(
  "DirectIndices", typeid(DirectIndices), DirectIndices::load);

IObject* LookupIndices::emptyShell() {
  return new LookupIndices();
}

void LookupIndices::load(CClassSP& clazz) {
  REGISTER(LookupIndices, clazz);
  SUPERCLASS(CObject);
  IMPLEMENTS(IFieldTweak::IIndices);
  FIELD(keyListField, "keyListField");
  FIELD(keyField, "keyField");
  FIELD(indexKeys, "indexKeys");
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const LookupIndices::TYPE = CClass::registerClassLoadMethod(
  "LookupIndices", typeid(LookupIndices), LookupIndices::load);

IObject* BulkFieldTweak::emptyShell() {
  return new BulkFieldTweak();
}

void BulkFieldTweak::load(CClassSP& clazz) {
  REGISTER(BulkFieldTweak, clazz);
  SUPERCLASS(FieldTweak);
  FIELD(_argument, "_argument");
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const BulkFieldTweak::TYPE = CClass::registerClassLoadMethod(
  "BulkFieldTweak", typeid(BulkFieldTweak), BulkFieldTweak::load);

IObject* ElementwiseFieldTweak::emptyShell() {
  return new ElementwiseFieldTweak();
}

void ElementwiseFieldTweak::load(CClassSP& clazz) {
  REGISTER(ElementwiseFieldTweak, clazz);
  SUPERCLASS(FieldTweak);
  FIELD(indices, "indices");
  FIELD(arguments, "arguments");
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const ElementwiseFieldTweak::TYPE = CClass::registerClassLoadMethod(
  "ElementwiseFieldTweak", typeid(ElementwiseFieldTweak), ElementwiseFieldTweak::load);

IObject* CalibratorExpiriesIndices::emptyShell() {
  return new CalibratorExpiriesIndices();
}

void CalibratorExpiriesIndices::load(CClassSP& clazz) {
  REGISTER(CalibratorExpiriesIndices, clazz);
  SUPERCLASS(CObject);
  IMPLEMENTS(IFieldTweak_IIndices);
  FIELD(field, "field");
  FIELD(indexKeys, "indexKeys");
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const CalibratorExpiriesIndices::TYPE = CClass::registerClassLoadMethod(
  "CalibratorExpiriesIndices", typeid(CalibratorExpiriesIndices), CalibratorExpiriesIndices::load);



// 
// *************
//  Legacy junk
// *************
//
// This a workaround to preserve the existing interface for storage of
// FieldRiskAxis's in Pyramid --- unfortunately they got out into the wild in a
// previous form as part of RiskMappingMatrix.
// 

struct LegacyFieldRiskAxis: CObject,
                            virtual RiskAxis {
    static void load(CClassSP& clazz);
    static CClassConstSP const TYPE;

    // These give us a flat (string) representation of 'property' for RiskAxis
    // storage -- a nasty hack; they're filled in by
    // Scalar/ExpiryFieldRiskAxis:validatePop2Object

    string className;                    // $required

    CStringArrayConstSP fieldPath;       // $required
    string tweakType;                    // $required
    OutputNameConstSP _marketDataName;   // $optional

    bool hasFloor;                       // $required
    double floor;                        // $required
    bool floorInclusive;                 // $required

    bool hasCap;                         // $required
    double cap;                          // $required
    bool capInclusive;                   // $required;

    bool clip;                           // $required

    bool useCalibratorRange;             // $required

    bool absoluteDistance;               // $required

    IFieldTweak::IOperatorConstSP operat0r() const {
        return IFieldTweak::IOperator::numeric(
                   TweakFunction::fromRepr(tweakType),
                   Range(hasFloor ? floor : -HUGE_VAL, floorInclusive,
                         hasCap ? cap : HUGE_VAL, capInclusive),
                   useCalibratorRange,
                   clip);
    }

    LegacyFieldRiskAxis(CClassConstSP type):
        CObject(type)
    {}

    ~LegacyFieldRiskAxis() {}
};

DECLARE(LegacyFieldRiskAxis)

struct ScalarFieldRiskAxis: LegacyFieldRiskAxis {  // $public

    static CClassConstSP const TYPE;
    static IObject* emptyShell();

    static void load(CClassSP& clazz);

    IScalarRiskPropertyConstSP prop;     // $transient

    void validatePop2Object() {
        try {
            prop = fieldRiskProperty::scalar(
                CClass::forName(className),
                IFieldTweak::bulk(
                    FieldPath::SP(fieldPath), FieldPathSP(),
                    CDouble::SP(1.), operat0r()),
                absoluteDistance);
        }
        catch (exception& e) {
            throw ModelException(e, __FUNCTION__);
        }
    }

    IRiskAxisConstSP thawed() const {
        return prop->axisFor(_marketDataName);
    }

    ScalarFieldRiskAxis():
        LegacyFieldRiskAxis(TYPE)
    {}

    ~ScalarFieldRiskAxis() {}
};

DECLARE(ScalarFieldRiskAxis)

struct ExpiryFieldRiskAxis: LegacyFieldRiskAxis {  // $public

    static CClassConstSP const TYPE;
    static IObject* emptyShell();
    static void load(CClassSP& clazz);

    void validate();

    CStringArrayConstSP expiryPath;                // $optional
    ExpiryConstSP previous;                        // $optional
    ExpiryArraySP expiries;                        // $required
    ExpiryConstSP next;                            // $optional
    int expiryDimension;                           // $required
    CDoubleArraySP weights;                        // $optional

    IScalarRiskPropertyConstSP scalarProp;         // $transient
    IExpiryRiskPropertyConstSP pointwiseProp;      // $transient

    void validatePop2Object() {
        try {
            if (!!weights && weights->size() != expiries->size()) {
                throw ModelException(
                    "'weights' are a different length than 'expiries'");
            }

            // horrible hack to "simplify" storage in pyramid

            if (!weights) {
                pointwiseProp = fieldRiskProperty::pointwise(
                   CClass::forName(className),
                   FieldPath::SP(fieldPath), FieldPathSP(),
                   !expiryPath ? FieldPathSP() :
                                 FieldPath::SP(expiryPath),
                   FieldPathSP(),
                   expiryDimension,
                   operat0r(),
                   CDouble::SP(1.),
                   absoluteDistance);
            }
            else {
                scalarProp = fieldRiskProperty::scalar(
                    CClass::forName(className),
                    IFieldTweak::elementwise(
                        FieldPath::SP(fieldPath), FieldPathSP(),
                        IFieldTweak::IIndices::lookupOrCalibratorExpiries(
                            FieldPath::SP(fieldPath),
                            !expiryPath ? FieldPathSP() :
                                          FieldPath::SP(expiryPath),
                            FieldPathSP(),
                            expiries),
                        expiryDimension,
                        weights,
                        operat0r()),
                    absoluteDistance);
            }
        }
        catch (exception& e) {
            throw ModelException(
                e, "ExpiryFieldRiskAxis::validatePop2Object()");
        }
    }

    IRiskAxisConstSP thawed() const {
        if (!pointwiseProp) {
            ASSERT(!!scalarProp);
            return scalarProp->axisFor(_marketDataName);
        }
        else {
            ASSERT(!!expiries && expiries->size() == 1);
            return pointwiseProp->axisFor(
                _marketDataName,
                ExpiryWindow::SP(previous, (*expiries)[0], next));
        }
    }

    ExpiryFieldRiskAxis():
        LegacyFieldRiskAxis(TYPE)
    {}

    ~ExpiryFieldRiskAxis() {}
};

DECLARE(ExpiryFieldRiskAxis)

RiskAxisSP IFieldTweak::legacyRepresentation(
        IAbstractRiskPropertyConstSP property,
        OutputNameConstSP name,
        bool absoluteDistance) const {
    try {
        const FieldTweak* t = dynamic_cast<const FieldTweak*>(this);
        if (!t) throw ModelException("Don't know how to translate this class");

        const NumericOperator* n = dynamic_cast<const NumericOperator*>(
                                       t->op.get());
        if (!n) throw ModelException("Can only translate NumericOperator, not " +
                                     t->op->getClass()->getName());

        LegacyFieldRiskAxisSP a;

        const ElementwiseFieldTweak* e;
        if (!!(e = dynamic_cast<const ElementwiseFieldTweak*>(this))) {
            ExpiryFieldRiskAxisSP ea(new ExpiryFieldRiskAxis());
            a = ea;

            IArrayConstSP indexKeys;

            {
                const LookupIndices* l;
                const CalibratorExpiriesIndices* c;

                if (!!(l = dynamic_cast<const LookupIndices*>(e->indices.get()))) {
                    if (!l->keyField->path()->empty()) {
                        throw ModelException("Cannot translate sub-fields");
                    }

                    ea->expiryPath = !l->keyListField ? CStringArrayConstSP() :
                                                        l->keyListField->path();

                    indexKeys = l->indexKeys;
                }
                else if (!!(c = dynamic_cast<const CalibratorExpiriesIndices*>(
                                    e->indices.get()))) {
                    indexKeys = c->indexKeys;
                }
                else {
                    throw ModelException("Can only translate LookupIndices "
                                         "or CalibratorExpiriesIndices, not " +
                                         e->indices->getClass()->getName());
                }
            }

            ea->expiryDimension = t->indexDimension;
            if (IExpiryRiskProperty::TYPE->isInstance(property)) {
                if (e->arguments->getLength() != 1) {
                    throw ModelException("Cannot translate pointwise axes for "
                                         "multiple expiries");
                }

                if (convertedToDouble(e->arguments->get(0), "argument") != 1.) {
                    throw ModelException("Cannot translate pointwise axes with "
                                         "non-unit arguments");
                }

                ASSERT(indexKeys->getLength() == 1); // like arguments
                    
                ExpiryWindowConstSP w;
                try {
                    w = ExpiryWindowConstSP::dynamicCast(indexKeys->get(0));
                }
                catch (exception& e) {
                    throw ModelException(e, "Cannot translate qualifier");
                }

                ea->expiries = ExpiryArray::SP(1, ExpirySP::constCast(w->expiry));
                ea->previous = w->previous;
                ea->next = w->next;
            }
            else {
                ea->weights = DoubleArray::SP(e->arguments->getLength());
                for (int w = 0; w < ea->weights->size(); ++w) {
                    try {
                        ea->weights->set(w, IObjectSP::constCast(e->arguments->get(w)));
                    }
                    catch (exception& q) {
                        throw ModelException(q,
                            "Cannot translate argument #" +
                            Format::toString(w) + " \"" +
                            IObject::stringOf(e->arguments->get(w)) + "\"");
                    }
                }

                ea->expiries = ExpiryArray::SP(indexKeys->getLength());
                for (int x = 0; x < ea->expiries->size(); ++x) {
                    try {
                        ea->expiries->set(x, IObjectSP::constCast(indexKeys->get(x)));
                    }
                    catch (exception& e) {
                        throw ModelException(e,
                            "Cannot translate qualifier #" +
                            Format::toString(x) + " \"" +
                            IObject::stringOf(indexKeys->get(x)) + "\"");
                    }
                }
            }
        }
        else if (dynamic_cast<const BulkFieldTweak*>(this)) {
            if (convertedToDouble(t->argument(0), "argument") != 1.) {
                throw ModelException("Cannot translate pointwise axes with "
                                     "non-unit arguments");
            }

            a.reset(new ScalarFieldRiskAxis());
        }
        else {
            throw ModelException("Don't know how to translate this class");
        }

        a->className = property->subjectInterface()->getName();
        a->fieldPath = t->field->path();
        if (!t->subField->path()->empty()) {
            throw ModelException("Cannot translate sub-fields");
        }

        a->tweakType = n->function->repr();
        a->hasFloor = n->hasFloor;
        a->floor = n->floor;
        a->floorInclusive = n->floorInclusive;
        a->hasCap = n->hasCap;
        a->cap = n->cap;
        a->capInclusive = n->capInclusive;
        a->useCalibratorRange = n->useCalibratorRange;
        a->clip = n->clip;
        a->absoluteDistance = absoluteDistance;
        a->_marketDataName = name;

        return a;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__,
            "Can't make a Pyramid representation for the " +
            getClass()->getName() + " '" + toString() + "'");
    }
}




void LegacyFieldRiskAxis::load(CClassSP& clazz) {
  REGISTER(LegacyFieldRiskAxis, clazz);
  SUPERCLASS(CObject);
  IMPLEMENTS(RiskAxis);
  FIELD(className, "className");
  FIELD(fieldPath, "fieldPath");
  FIELD(tweakType, "tweakType");
  FIELD(_marketDataName, "_marketDataName");
  FIELD_MAKE_OPTIONAL(_marketDataName);
  FIELD(hasFloor, "hasFloor");
  FIELD(floor, "floor");
  FIELD(floorInclusive, "floorInclusive");
  FIELD(hasCap, "hasCap");
  FIELD(cap, "cap");
  FIELD(capInclusive, "capInclusive");
  FIELD(clip, "clip");
  FIELD(useCalibratorRange, "useCalibratorRange");
  FIELD(absoluteDistance, "absoluteDistance");
}

CClassConstSP const LegacyFieldRiskAxis::TYPE = CClass::registerClassLoadMethod(
  "LegacyFieldRiskAxis", typeid(LegacyFieldRiskAxis), LegacyFieldRiskAxis::load);

IObject* ScalarFieldRiskAxis::emptyShell() {
  return new ScalarFieldRiskAxis();
}

void ScalarFieldRiskAxis::load(CClassSP& clazz) {
  REGISTER(ScalarFieldRiskAxis, clazz);
  clazz->setPublic();
  SUPERCLASS(LegacyFieldRiskAxis);
  FIELD(prop, "prop");
  FIELD_MAKE_TRANSIENT(prop);
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const ScalarFieldRiskAxis::TYPE = CClass::registerClassLoadMethod(
  "ScalarFieldRiskAxis", typeid(ScalarFieldRiskAxis), ScalarFieldRiskAxis::load);

IObject* ExpiryFieldRiskAxis::emptyShell() {
  return new ExpiryFieldRiskAxis();
}

void ExpiryFieldRiskAxis::load(CClassSP& clazz) {
  REGISTER(ExpiryFieldRiskAxis, clazz);
  clazz->setPublic();
  SUPERCLASS(LegacyFieldRiskAxis);
  FIELD(expiryPath, "expiryPath");
  FIELD_MAKE_OPTIONAL(expiryPath);
  FIELD(previous, "previous");
  FIELD_MAKE_OPTIONAL(previous);
  FIELD(expiries, "expiries");
  FIELD(next, "next");
  FIELD_MAKE_OPTIONAL(next);
  FIELD(expiryDimension, "expiryDimension");
  FIELD(weights, "weights");
  FIELD_MAKE_OPTIONAL(weights);
  FIELD(scalarProp, "scalarProp");
  FIELD_MAKE_TRANSIENT(scalarProp);
  FIELD(pointwiseProp, "pointwiseProp");
  FIELD_MAKE_TRANSIENT(pointwiseProp);
  EMPTY_SHELL_METHOD(emptyShell);
}

CClassConstSP const ExpiryFieldRiskAxis::TYPE = CClass::registerClassLoadMethod(
  "ExpiryFieldRiskAxis", typeid(ExpiryFieldRiskAxis), ExpiryFieldRiskAxis::load);


DRLIB_END_NAMESPACE
