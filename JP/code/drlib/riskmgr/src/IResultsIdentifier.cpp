/**
 * @file IResultsIdentifier.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/Void.hpp"
#include "edginc/BoxedInt.hpp"
#include "edginc/ExpiryWindow.hpp"
#include "edginc/Results.hpp"
#include "edginc/ExpiryResult.hpp"
#include "edginc/ExpiryPairResult.hpp"
#include "edginc/MatrixResult.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/IResultsIdentifier.hpp"

DRLIB_BEGIN_NAMESPACE

IResultsIdentifier::IResultsIdentifier() {}

IResultsIdentifier::~IResultsIdentifier() {}

CClassConstSP const IResultsIdentifier::TYPE = CClass::registerInterfaceLoadMethod(
    "IResultsIdentifier", typeid(IResultsIdentifier), 0);

/**
 * Base class for IResultsIdentifier implementations.
 *
 * Exists simply to keep IResultsIdentifier pure virtual.
 */

class ResultsIdentifier: public CObject,
                         public virtual IResultsIdentifier {

    static void load(CClassSP& clazz) {
        REGISTER(ResultsIdentifier, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IResultsIdentifier);
        FIELD(_packet, "_packet");
        FIELD(_entry, "_entry");
    }

public:

    static CClassConstSP const TYPE;

    ResultsIdentifier(CClassConstSP type, string packet, string entry):
        CObject(type), _packet(packet), _entry(new OutputName(entry))
    {}

    ResultsIdentifier(CClassConstSP type, string packet, OutputNameConstSP entry):
        CObject(type), _packet(packet), _entry(entry)
    {}

    ~ResultsIdentifier() {}

protected:

    string _packet;
    OutputNameConstSP _entry;

public:

    // FIXME it would be nice not to have these but we need them for
    // e.g. filtering BasketDelta results in Delta.cpp

    string packet() const { return _packet; }
    OutputNameConstSP entry() const { return _entry; }

    void storeUntweakable(CResultsSP results, UntweakableConstSP oops) const {
        results->storeGreek(IObjectSP((IObject*)oops.get()), _packet,
                            !_entry ? OutputNameSP(new OutputName(_packet))
                                    : _entry);
    }

    void storeNotApplicableToName(CResultsSP results) const {
        results->storeNotApplicable(_packet,
                                    OutputNameSP((OutputName*)_entry.get()),
                                    IObjectSP(new NotApplicable()));
    }

    void storeNotApplicableToInstrument(CResultsSP results) const {
        results->storeNotApplicable(_packet);
    }

    string toString() const {
        return !_entry ? _packet : _packet + " " + _entry->toString();
    }

    void store(ResultsSP results, IObjectSP value) const {
        results->storeGreek(value, _packet,  entry());
    }
};

CClassConstSP const ResultsIdentifier::TYPE = CClass::registerClassLoadMethod(
    "ResultsIdentifier", typeid(ResultsIdentifier), load);

// 
// ============================================
//  IResultsIdentifier::SP(string, OutputName)
// ============================================
// 

/**
 * Names a scalar entry in a Results dictionary
 *
 * This is a pretty thin wrapper round Results::storeScalarGreek().
 */

class ScalarResultsIdentifier: public ResultsIdentifier {

    static IObject* emptyShell() {
        return new ScalarResultsIdentifier("", "");
    }

    static void load(CClassSP& clazz) {
        REGISTER(ScalarResultsIdentifier, clazz);
        SUPERCLASS(ResultsIdentifier);
        EMPTY_SHELL_METHOD(emptyShell);
    }

public:

    static CClassConstSP const TYPE;

    ScalarResultsIdentifier(string packet, string entry):
        ResultsIdentifier(TYPE, packet, entry)
    {}

    ScalarResultsIdentifier(string packet, OutputNameConstSP entry):
        ResultsIdentifier(TYPE, packet, entry)
    {}

    ~ScalarResultsIdentifier() {}

    bool exists(ResultsConstSP results) const {
        return !_entry ? results->packetExists(_packet) :
                         results->exists(_packet, _entry);
    }

    void store(ResultsSP results, double value) const {
        results->storeScalarGreek(value, _packet, entry());
    }
};

CClassConstSP const ScalarResultsIdentifier::TYPE = CClass::registerClassLoadMethod(
    "ScalarResultsIdentifier", typeid(ScalarResultsIdentifier), load);

IResultsIdentifierSP IResultsIdentifier::SP(string packet,
                                            const string& entry) {
    return IResultsIdentifierSP(new ScalarResultsIdentifier(
                                    packet, OutputName::SP(entry)));
}

IResultsIdentifierSP IResultsIdentifier::SP(string packet,
                                            OutputNameConstSP entry,
                                            VoidConstSP, VoidConstSP) {
    return IResultsIdentifierSP(new ScalarResultsIdentifier(packet, entry));
}

// 
// ********************************
//  Qualified IResultsIdentifier's
// ********************************
// 

/**
 * Names an entry in some container in a Results dictionary
 *
 * This class handles the logic of storing numbers in a Results dictionary
 * which are qualified not only by packet and name, but also by something else.
 * It can be used uniformly with ScalarResultsIdentifier.
 */

template <class QUALIFIER, class CONTAINER>
struct QualifiedResultsIdentifier: ResultsIdentifier {

    static void load(CClassSP& clazz) {
        REGISTER(QualifiedResultsIdentifier, clazz);
        SUPERCLASS(ResultsIdentifier);
        FIELD(qualifier, "qualifier");
    }

    typedef QUALIFIER Qualifier; 
    DECLARE(Qualifier)

    typedef CONTAINER Container;
    DECLARE(Container)

    static CClassConstSP const TYPE;

    QualifierConstSP qualifier;

    virtual bool existsInContainer(const Container& container) const = 0;
    virtual void storeInContainer(Container& container,
                                  double value) const = 0;

    void validatePop2Object() {
        if (!qualifier) {
            throw ModelException(
                "QualifiedResultsIdentifier::validate()",
                "qualifier must be present");
        }
    }

    QualifiedResultsIdentifier(CClassConstSP type,
                               string packet, OutputNameConstSP entry,
                               QualifierConstSP qualifier):
        ResultsIdentifier(type, packet, entry),
        qualifier(qualifier)
    {
        validatePop2Object();
    }

    ~QualifiedResultsIdentifier() {}

    /**
     * Whether a number has been stored against this name in @a results.
     *
     * True if there's an entry under the right packet and name, it's
     * an ExpiryResultArray, and it contains an entry for the right
     * expiry.
     */

    bool exists(ResultsConstSP results) const {
        if (!_entry) return results->packetExists(_packet);

        if (!results->exists(_packet, _entry)) return false;

        // Actually it could still be a NotApplicable for the whole packet

        if (results->exists(_packet, Results::emptyName) &&
            !Container::TYPE->isInstance(
                results->retrieveGreek(_packet, Results::emptyName))) {
            return true;
        }

        const Container* them = dynamic_cast<const Container*>(
            results->retrieveGreek(_packet, _entry).get());

        return !them /* probably an Untweakable */ || existsInContainer(*them);
    }

    /**
     * Store a @a value against this name in @a results.
     *
     * Stores an ExpiryResultArray under our given packet and name, if it's not
     * there already.  Sets (if necessary inserting) the entry for our expiry
     * to @a value.
     */

    void store(ResultsSP results, double value) const {
        OutputNameConstSP e = !_entry ? OutputNameConstSP(new OutputName(_packet))
                                      : _entry;

        ContainerSP us;

        if (!results->exists(_packet, e)) {
            us.reset(new Container());
        }
        else {
            const IObject* it = results->retrieveGreek(_packet, e).get();

            // Actually it could still be a NotApplicable for the whole packet

            if (results->exists(_packet, Results::emptyName) &&
                !Container::TYPE->isInstance(
                    results->retrieveGreek(_packet, Results::emptyName))) {
                return;
            }

            const Container* them =
                    dynamic_cast<const Container*>(it);

            if (!them) {
                if (dynamic_cast<const Untweakable*>(it)) {
                    return;
                }
                else {
                    throw ModelException("QualifiedResultsIdentifier::store",
                        "Expected either " + Container::TYPE->getName() + 
                        " or Untweakable in results, but found a " +
                        it->getClass()->getName());
                }
            }

            // Cloning the container to avoid upsetting anyone in case of
            // aliasing was a nice defensive touch but is actually unnecessary
            // since in practice we create it ourselves (new Container()
            // above), and expensive now we can handle MatrixResult's

            us.reset(IObject::copy(them));
            // us.reset(const_cast<Container*>(them));
        }

        storeInContainer(*us, value);
        results->storeGreek(us, _packet, e);
    }

    string toString() const {
        return ResultsIdentifier::toString() + "[" +
               (!qualifier ? "all" : qualifier->toString()) + "]";
    }
};

template <>
CClassConstSP const QualifiedResultsIdentifier<Expiry, ExpiryResultArray>::TYPE = CClass::registerClassLoadMethod(
    "QualifiedResultsIdentifier<Expiry, ExpiryResultArray>", typeid(QualifiedResultsIdentifier<Expiry, ExpiryResultArray>), load);

template <>
CClassConstSP const QualifiedResultsIdentifier<ExpiryPair, ExpiryPairResultArray>::TYPE = CClass::registerClassLoadMethod(
    "QualifiedResultsIdentifier<ExpiryPair, ExpiryPairResultArray>", typeid(QualifiedResultsIdentifier<ExpiryPair, ExpiryPairResultArray>), load);

template <>
CClassConstSP const QualifiedResultsIdentifier<ExpiryAndStrike, MatrixResult>::TYPE = CClass::registerClassLoadMethod(
    "QualifiedResultsIdentifier<ExpiryAndStrike, MatrixResult>", typeid(QualifiedResultsIdentifier<ExpiryAndStrike, MatrixResult>), load);

template <>
CClassConstSP const QualifiedResultsIdentifier<BoxedInt, DoubleArray>::TYPE = CClass::registerClassLoadMethod(
    "QualifiedResultsIdentifier<BoxedInt, DoubleArray>", typeid(QualifiedResultsIdentifier<BoxedInt, DoubleArray>), load);

// 
// ====================================================
//  IResultsIdentifier::SP(string, OutputName, Expiry)
// ====================================================
// 

/**
 * Names an entry in an ExpiryResultsArray in a Results dictionary
 *
 * This class handles the logic of storing numbers in a Results dictionary
 * which are qualified not only by packet and name, but also by Expiry.  It can
 * be used uniformly with ScalarResultsIdentifier.
 */

class ExpiryResultsIdentifier:
    public QualifiedResultsIdentifier<Expiry, ExpiryResultArray> {

    static IObject* emptyShell() {
        return new ExpiryResultsIdentifier("", OutputNameConstSP(),
                                           ExpiryConstSP());
    }

    static void load(CClassSP& clazz) {
        REGISTER(ExpiryResultsIdentifier, clazz);
        SUPERCLASS(QualifiedResultsIdentifier<Expiry _COMMA_ ExpiryResultArray>);
        EMPTY_SHELL_METHOD(emptyShell);
    }

public:

    static CClassConstSP const TYPE;

    ExpiryResultsIdentifier(string packet, OutputNameConstSP entry,
                            ExpiryConstSP expiry):
        QualifiedResultsIdentifier<Expiry, ExpiryResultArray>(
            TYPE, packet, entry, expiry)
    {}

    ~ExpiryResultsIdentifier() {}

    bool existsInContainer(const ExpiryResultArray& them) const {
        for (int e = 0; e < them.size(); ++e) {
            if (them[e].getExpiry()->equals(qualifier.get())) return true;
        }

        return false;
    }

    void storeInContainer(ExpiryResultArray& us, double value) const {
        int i;
        for (i = 0; i < us.size(); ++i) {
            if (us[i].getExpiry()->equals(qualifier.get())) break;
        }

        if (i == us.size()) {
            us.push_back(ExpiryResult(qualifier, value));
        }
        else {
            us[i] = ExpiryResult(qualifier, value);
        }
    }
};

CClassConstSP const ExpiryResultsIdentifier::TYPE = CClass::registerClassLoadMethod(
    "ExpiryResultsIdentifier", typeid(ExpiryResultsIdentifier), load);

IResultsIdentifierSP IResultsIdentifier::SP(
        string packet, OutputNameConstSP entry, ExpiryConstSP expiry) {
    return IResultsIdentifierSP(
        new ExpiryResultsIdentifier(packet, entry, expiry));
}

IResultsIdentifierSP IResultsIdentifier::SP(
        string packet, OutputNameConstSP entry, ExpiryWindowConstSP expiry) {
    return IResultsIdentifierSP(
        new ExpiryResultsIdentifier(packet, entry,
                                    !expiry ? ExpirySP() : expiry->expiry));
}

// 
// ========================================================
//  IResultsIdentifier::SP(string, OutputName, ExpiryPair)
// ========================================================
// 

/**
 * Names an entry in an ExpiryPairResultsArray in a Results dictionary
 *
 * This class handles the logic of storing numbers in a Results dictionary
 * which are qualified not only by packet and name, but also by ExpiryPair.  It
 * can be used uniformly with ScalarResultsIdentifier.
 */

class ExpiryPairResultsIdentifier:
        public QualifiedResultsIdentifier<ExpiryPair, ExpiryPairResultArray> {

    static IObject* emptyShell() {
        return new ExpiryPairResultsIdentifier("", OutputNameConstSP(),
                                               ExpiryPairConstSP());
    }

    static void load(CClassSP& clazz) {
        REGISTER(ExpiryPairResultsIdentifier, clazz);
        SUPERCLASS(QualifiedResultsIdentifier<ExpiryPair _COMMA_ ExpiryPairResultArray>);
        EMPTY_SHELL_METHOD(emptyShell);
    }

public:

    static CClassConstSP const TYPE;

    ExpiryPairResultsIdentifier(string packet, OutputNameConstSP entry,
                                ExpiryPairConstSP expiry):
        QualifiedResultsIdentifier<ExpiryPair, ExpiryPairResultArray>(
            TYPE, packet, entry, expiry)
    {}

    ~ExpiryPairResultsIdentifier() {}

    bool existsInContainer(const ExpiryPairResultArray& them) const {
        for (int e = 0; e < them.size(); ++e) {
            if (them[e].getExpiryPair()->equals(qualifier.get())) return true;
        }

        return false;
    }

    void storeInContainer(ExpiryPairResultArray& us, double value) const {
        int i;
        for (i = 0; i < us.size(); ++i) {
            if (us[i].getExpiryPair()->equals(qualifier.get())) break;
        }

        if (i == us.size()) {
            us.push_back(ExpiryPairResult(qualifier, value));
        }
        else {
            us[i] = ExpiryPairResult(qualifier, value);
        }
    }
};

CClassConstSP const ExpiryPairResultsIdentifier::TYPE = CClass::registerClassLoadMethod(
    "ExpiryPairResultsIdentifier", typeid(ExpiryPairResultsIdentifier), load);

IResultsIdentifierSP IResultsIdentifier::SP(
        string packet, OutputNameConstSP entry, ExpiryPairConstSP expiry) {
    return IResultsIdentifierSP(
        new ExpiryPairResultsIdentifier(packet, entry, expiry));
}

IResultsIdentifierSP IResultsIdentifier::SP(
        string packet, OutputNameConstSP entry,
        ExpiryWindowConstSP expiry0, ExpiryWindowConstSP expiry1) {
    return IResultsIdentifierSP(
        new ExpiryPairResultsIdentifier(packet, entry,
                                        ExpiryPair::SP(expiry0->expiry,
                                                       expiry1->expiry)));
}

// 
// =============================================================
//  IResultsIdentifier::SP(string, OutputName, ExpiryAndStrike)
// =============================================================
// 

/**
 * Names an entry in a MatrixResult in a Results dictionary
 *
 * This class handles the logic of storing numbers in a Results dictionary
 * which are qualified not only by packet and name, but also by
 * ExpiryAndStrike.  It can be used uniformly with ScalarResultsIdentifier.
 */

class ExpiryAndStrikeResultsIdentifier:
        public QualifiedResultsIdentifier<ExpiryAndStrike, MatrixResult> {

    static IObject* emptyShell() {
        return new ExpiryAndStrikeResultsIdentifier("", OutputNameConstSP(),
                                                    ExpiryAndStrikeConstSP());
    }

    static void load(CClassSP& clazz) {
        REGISTER(ExpiryAndStrikeResultsIdentifier, clazz);
        SUPERCLASS(QualifiedResultsIdentifier<ExpiryAndStrike _COMMA_ MatrixResult>);
        EMPTY_SHELL_METHOD(emptyShell);
    }

public:

    static CClassConstSP const TYPE;

    ExpiryAndStrikeResultsIdentifier(string packet, OutputNameConstSP entry,
                                     ExpiryAndStrikeConstSP expiryAndStrike):
        QualifiedResultsIdentifier<ExpiryAndStrike, MatrixResult>(
            TYPE, packet, entry, expiryAndStrike)
    {}

    ~ExpiryAndStrikeResultsIdentifier() {}

    bool existsInContainer(const MatrixResult& them) const {
        // This is "imperfect": if the expiry/strike exist in the matrix, but
        // the actual value is zero, then that may be because it hasn't yet
        // been set (in which case we should return false) or because it has
        // been set to zero (ditto true).  Hence we may end up "wrongly"
        // overwriting values previously set to zero by some other greek.
        // Arguably this shouldn't matter since like-labelled values should be
        // the same whoever computed them, and in practice this is true apart
        // from a few buggy cases such as when the tree gets confused about
        // Delta & BasketDelta.

        return them.getOr0(*qualifier) != 0;
    }

    void storeInContainer(MatrixResult& us, double value) const {
        us.set(*qualifier, value);
    }
};

CClassConstSP const ExpiryAndStrikeResultsIdentifier::TYPE = CClass::registerClassLoadMethod(
    "ExpiryAndStrikeResultsIdentifier", typeid(ExpiryAndStrikeResultsIdentifier), load);

IResultsIdentifierSP IResultsIdentifier::SP(
        string packet, OutputNameConstSP entry,
        ExpiryAndStrikeConstSP expiryAndStrike) {
    return IResultsIdentifierSP(
        new ExpiryAndStrikeResultsIdentifier(packet, entry, expiryAndStrike));
}

// 
// =================================================
//  IResultsIdentifier::SP(string, OutputName, int)
// =================================================
// 

/**
 * Names an entry in a DoubleArray in a Results dictionary
 *
 * This class handles the logic of storing numbers in a Results dictionary
 * which are qualified not only by packet and name, but also by int.  It can
 * be used uniformly with ScalarResultsIdentifier.
 */

class ElementResultsIdentifier:
    public QualifiedResultsIdentifier<BoxedInt, DoubleArray> {

    static IObject* emptyShell() {
        return new ElementResultsIdentifier("", OutputNameConstSP(),
                                           BoxedIntConstSP());
    }

    static void load(CClassSP& clazz) {
        REGISTER(ElementResultsIdentifier, clazz);
        SUPERCLASS(QualifiedResultsIdentifier<BoxedInt _COMMA_ DoubleArray>);
        EMPTY_SHELL_METHOD(emptyShell);
    }

public:

    static CClassConstSP const TYPE;

    void validatePop2Object() {
        QualifiedResultsIdentifier<BoxedInt, DoubleArray>::validatePop2Object();

        if (qualifier->intValue() < 0) {
            throw ModelException(
                "ElementResultsIdentifier::validatePop2Object()",
                "Index qualifier must not be negative but is " +
                Format::toString(qualifier->intValue()));
        }
    }

    ElementResultsIdentifier(string packet, OutputNameConstSP entry,
                             BoxedIntConstSP index):
        QualifiedResultsIdentifier<BoxedInt, DoubleArray>(
            TYPE, packet, entry, index)
    {
        validatePop2Object();
    }

    ~ElementResultsIdentifier() {}

    bool existsInContainer(const DoubleArray& them) const {
        return qualifier->intValue() < them.size();
    }

    void storeInContainer(DoubleArray& us, double value) const {
        int index = qualifier->intValue();
        ASSERT(index >= 0);

        if (index >= us.size()) us.resize(index + 1, 0.);
        us[index] = value;
    }
};

CClassConstSP const ElementResultsIdentifier::TYPE = CClass::registerClassLoadMethod(
    "ElementResultsIdentifier", typeid(ElementResultsIdentifier), load);

IResultsIdentifierSP IResultsIdentifier::SP(
        string packet, OutputNameConstSP entry, int i) {
    return SP(packet, entry, BoxedInt::SP(i));
}

IResultsIdentifierSP IResultsIdentifier::SP(
        string packet, OutputNameConstSP entry, BoxedIntConstSP i) {
    return IResultsIdentifierSP(new ElementResultsIdentifier(packet, entry, i));
}

DRLIB_END_NAMESPACE
