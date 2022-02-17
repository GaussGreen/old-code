/**
 * @file IResultsFunction.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/Results.hpp"
#include "edginc/NotApplicableException.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/TRACE.hpp"

DRLIB_BEGIN_NAMESPACE

CClassConstSP const IResultsFunction::TYPE = CClass::registerInterfaceLoadMethod(
    "IResultsFunction", typeid(IResultsFunction), 0);

/**
 * IResultsFunction::price()
 */

//@{

class PriceResultsFunction: public CObject,
                            public virtual IResultsFunction {

    static void load(CClassSP& clazz) {
        REGISTER(PriceResultsFunction, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IResultsFunction);
        EMPTY_SHELL_METHOD(DefaultConstructor<PriceResultsFunction>::iObject);
    }

public:
    static CClassConstSP const TYPE;

    PriceResultsFunction(): CObject(TYPE) {}

    double operator ()(const Results *results) const {
        TRACE_METHOD;
        TRACE_SHOW(results->retrievePrice());

        return results->retrievePrice();
    }

    OutputRequestArrayConstSP outputRequests() const {
        return OutputRequestArrayConstSP();
    }

    string toString() const { return "price"; }
};

CClassConstSP const PriceResultsFunction::TYPE = CClass::registerClassLoadMethod(
    "PriceResultsFunction", typeid(PriceResultsFunction), load);

IResultsFunctionConstSP IResultsFunction::price() {
    static IResultsFunctionSP it(new PriceResultsFunction());
    return it;
}

//@}

/**
 * IResultsFunction::outputRequest()
 */

//@{

class OutputRequestResultsFunction: public CObject,
                                    public virtual IResultsFunction {

    static void load(CClassSP& clazz) {
        REGISTER(OutputRequestResultsFunction, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IResultsFunction);
        EMPTY_SHELL_METHOD(DefaultConstructor<OutputRequestResultsFunction>::iObject);
        FIELD(request, "request");
    }

public:

    static CClassConstSP const TYPE;

    OutputRequestConstSP request;

    OutputRequestResultsFunction(OutputRequestConstSP request =
                                     OutputRequestConstSP()):
        CObject(TYPE),
        request(request)
    {}

    double operator ()(const Results *results) const {
        TRACE_METHOD;

        const string& requestPacketName = request->getPacketName();
        const string& requestName = request->getRequestName();
        OutputNameSP requestOutName(new OutputName(requestName));

        IObjectConstSP result;

        if (results->exists(request.get())){
            result = results->retrieveGreek(requestPacketName,
                                            requestOutName);
        }

        if (!result) {
            throw ModelException("IResultsFunction::operator ()()",
                                 requestName + " output not requested");
        }
        else if (NotApplicable::TYPE->isInstance(result)) {
            throw NotApplicableException();
        }
        else if (!CDouble::TYPE->isInstance(result)) {
            throw ModelException("IResultsFunction::operator ()()",
                                 requestName + " calculation failed");
        }
        else {
            double it = CDoubleConstSP::dynamicCast(result)->doubleValue();
            TRACE("results->retrieveGreek(\"" << requestPacketName <<
                  "\", \"" << requestName << "\") = " << it);
            return it;
        }
    }

    OutputRequestArrayConstSP outputRequests() const {
        return OutputRequestArray::SP(1, OutputRequestSP::constCast(request));
    }

    string toString() const { return request->getRequestName(); }
};

CClassConstSP const OutputRequestResultsFunction::TYPE = CClass::registerClassLoadMethod(
    "OutputRequestResultsFunction", typeid(OutputRequestResultsFunction), load);

IResultsFunctionConstSP IResultsFunction::outputRequest(
        OutputRequestConstSP request) {
    return IResultsFunctionConstSP(new OutputRequestResultsFunction(request));
}

IResultsFunctionConstSP IResultsFunction::outputRequest(
        const string& request) {
    return IResultsFunctionConstSP(new OutputRequestResultsFunction(
               OutputRequestSP(new OutputRequest(request))));
}

//@}

/**
 * IResultsFunction::zero()
 */

//@{

class ZeroResultsFunction: public CObject,
                           public virtual IResultsFunction {

    static void load(CClassSP& clazz) {
        REGISTER(ZeroResultsFunction, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IResultsFunction);
        EMPTY_SHELL_METHOD(DefaultConstructor<PriceResultsFunction>::iObject);
    }

public:
    static CClassConstSP const TYPE;

    ZeroResultsFunction(): CObject(TYPE) {}

    double operator ()(const Results *results) const {
        return 0;
    }

    OutputRequestArrayConstSP outputRequests() const {
        return OutputRequestArrayConstSP();
    }

    string toString() const { return "zero"; }
};

CClassConstSP const ZeroResultsFunction::TYPE = CClass::registerClassLoadMethod(
    "ZeroResultsFunction", typeid(ZeroResultsFunction), load);

IResultsFunctionConstSP IResultsFunction::zero() {
    static IResultsFunctionSP it(new ZeroResultsFunction());
    return it;
}

//@}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
