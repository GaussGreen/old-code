/**
 * @file RiskMappingMatrix.cpp
 */

#include "edginc/config.hpp"
#include "edginc/TRACE.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/BoxedInt.hpp"
#include "edginc/Maths.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/ExpiryPair.hpp"
#include "edginc/ExpiryAndStrike.hpp"
#include "edginc/ExpiryWindow.hpp"
#include "edginc/OutputName.hpp"
#include "edginc/IQualifiedRiskAxis.hpp"
#include "edginc/CompoundHypothesis.hpp"
#include "edginc/AbstractPropertyTweakHypothesis.hpp"
#include "edginc/RiskAxis.hpp"
#include "edginc/RiskMappingMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

RiskMappingMatrix::RiskMappingMatrix():
    MarketObject(TYPE)
{}

static RiskAxisArrayConstSP frozen(IRiskAxisArrayConstSP axes,
                                   const string& inOrOut) {
    RiskAxisArraySP them(new RiskAxisArray(axes->size()));

    for (int i = 0; i < axes->size(); ++i) {
        if (!(*axes)[i]) {
            throw ModelException("Null pointer in '" + inOrOut + "Basis'");
        }

        (*them)[i] = RiskAxisSP::constCast((*axes)[i]->frozen());

        if (!(*them)[i]) {
            throw ModelException("Null pointer returned from frozen() in "
                                 "'" + inOrOut + "Basis'");
        }
    }

    return them;
}

RiskMappingMatrix::RiskMappingMatrix(
        const string &name,
        IRiskAxisArrayConstSP inputBasis,
        CDoubleArrayConstSP inputDistanceScaling,
        IRiskAxisArrayConstSP outputBasis,
        const DoubleMatrix& coefficients,
        CDoubleArrayConstSP errorBars):
    MarketObject(TYPE),
    name(name),
    inputDistanceScaling(new DoubleArray(*inputDistanceScaling)),
    coefficients(new DoubleMatrix(coefficients)),
    errorBars(new DoubleArray(*errorBars)),
    _inputBasis(inputBasis),
    _outputBasis(outputBasis)
{
    try {
        this->inputBasis = frozen(inputBasis, "input");
        this->outputBasis = frozen(outputBasis, "output");
        validate();
    }
    catch (exception& e) {
        throw ModelException(e, "RiskMappingMatrix::RiskMappingMatrix()");
    }
}

RiskMappingMatrix::RiskMappingMatrix(
        const string &name,
        IRiskAxisArrayConstSP inputBasis,
        CDoubleArrayConstSP inputDistanceScaling,
        IRiskAxisArrayConstSP outputBasis,
        const DoubleMatrix& coefficients,
		const DoubleMatrix& coefficients2,
        CDoubleArrayConstSP errorBars):
    MarketObject(TYPE),
    name(name),
    inputDistanceScaling(new DoubleArray(*inputDistanceScaling)),
    coefficients(new DoubleMatrix(coefficients)),
	coefficients2(new DoubleMatrix(coefficients2)),
    errorBars(new DoubleArray(*errorBars)),
    _inputBasis(inputBasis),
    _outputBasis(outputBasis)
{
    try {
        this->inputBasis = frozen(inputBasis, "input");
        this->outputBasis = frozen(outputBasis, "output");
        validate();
    }
    catch (exception& e) {
        throw ModelException(e, "RiskMappingMatrix::RiskMappingMatrix()");
    }
}

RiskMappingMatrixSP RiskMappingMatrix::SP(
        const string &name,
        IRiskAxisArrayConstSP inputBasis,
        CDoubleArrayConstSP inputDistanceScaling,
        IRiskAxisArrayConstSP outputBasis,
        const DoubleMatrix& coefficients,
        CDoubleArrayConstSP errorBars) {
    return RiskMappingMatrixSP(new RiskMappingMatrix(
        name, inputBasis, inputDistanceScaling,
        outputBasis, coefficients, errorBars));
}

RiskMappingMatrixSP RiskMappingMatrix::SP(
        const string &name,
        IRiskAxisArrayConstSP inputBasis,
        CDoubleArrayConstSP inputDistanceScaling,
        IRiskAxisArrayConstSP outputBasis,
        const DoubleMatrix& coefficients,
		const DoubleMatrix& coefficients2,
        CDoubleArrayConstSP errorBars) {
    return RiskMappingMatrixSP(new RiskMappingMatrix(
        name, inputBasis, inputDistanceScaling,
        outputBasis, coefficients, coefficients2, errorBars));
}

RiskMappingMatrix::~RiskMappingMatrix() {}

void RiskMappingMatrix::validate() {
    try {
        _inputProperties.reset(
            new IAbstractRiskPropertyArray(_inputBasis->size()));

        for (int i = 0; i < _inputBasis->size(); ++i) {
            try {
                if (!(*_inputBasis)[i]) {
                    throw ModelException("Entry is null");
                }

                if (!(*_inputBasis)[i]->marketDataName()) {
                    throw ModelException(
                        "Axis " + (*_inputBasis)[i]->getClass()->getName() +
                        " won't say what market name it's for");
                }

                (*_inputProperties)[i] = IAbstractRiskPropertySP::constCast(
                    (*_inputBasis)[i]->abstractProperty());

                if (!(*_inputProperties)[i]) {
                    throw ModelException(
                        "Axis " + (*_inputBasis)[i]->getClass()->getName() +
                        " won't say what property it parameterises");
                }
            }
            catch (exception& e) {
                throw ModelException(
                    e, "inputBasis entry #" + Format::toString(i));
            }
        }

        for (int o = 0; o < _outputBasis->size(); ++o) {
            if (!(*_outputBasis)[o]) {
                throw ModelException("outputBasis entry #" +
                                     Format::toString(o) + " is null");
            }
        }

        if (coefficients->numRows() != outputBasis->size()) {
            throw ModelException("'outputBasis' must have length = num rows "
                                 "in 'coefficients' matrix");
        }

        if (coefficients->numCols() != inputBasis->size()) {
            throw ModelException("'inputBasis' must have length = num cols "
                                 "in 'coefficients' matrix");
        }

        if (coefficients->numCols() != inputDistanceScaling->size()) {
            throw ModelException("'inputDistanceScaling' must have length = "
                                 "num cols in 'coefficients' matrix");
        }

        if (coefficients->numCols() != errorBars->size()) {
            throw ModelException("'errorBars' must have length = "
                                 "num cols in 'coefficients' matrix");
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

static IRiskAxisArrayConstSP thawed(RiskAxisArrayConstSP axes,
                                    const string& inOrOut) {
    IRiskAxisArraySP them(new IRiskAxisArray(axes->size()));

    for (int i = 0; i < axes->size(); ++i) {
        if (!(*axes)[i]) {
            throw ModelException("Null pointer in '" + inOrOut + "Basis'");
        }

        (*them)[i] = IRiskAxisSP::constCast((*axes)[i]->thawed());

        if (!(*them)[i]) {
            throw ModelException("Null pointer returned from thawed() in "
                                 "'" + inOrOut + "Basis'");
        }
    }

    return them;
}

void RiskMappingMatrix::validatePop2Object() {
    try {
        _inputBasis = thawed(inputBasis, "input");
        _outputBasis = thawed(outputBasis, "output");
        validate();
    }
    catch (exception& e) {
        throw ModelException(e, "RiskMappingMatrix::validatePop2Object()");
    }
}

string RiskMappingMatrix::getName() const {
    return name;
}

IHypothesisConstSP RiskMappingMatrix::mapped(IHypothesisConstSP hypothesis)
        const {
    TRACE_METHOD;
    TRACE_SHOW(this->getName());
    TRACE_SHOW(*hypothesis);

    AtomicHypothesisArraySP outs(new AtomicHypothesisArray());

    DoubleArray v(coefficients->numRows(), 0.);
    DoubleArray dists(hypothesis->numAtomics(), 0.);

    for (int a = 0; a < hypothesis->numAtomics(); ++a) {
        AtomicHypothesisConstSP atomic = hypothesis->atomic(a);
        TRACE_BLOCK("Trying to map " << *atomic);

        outs->push_back(AtomicHypothesisSP::constCast(atomic));

        IRiskAxisConstSP axis = atomic->axis();
        if (!!axis) {
            for (int j = 0; j < _inputBasis->size(); ++j) {
                if ((*_inputBasis)[j]->equalTo(axis.get())) {
                    TRACE_BLOCK("Column #" << j << " " << *(*_inputBasis)[j] <<
                                " matches");
                    double k = atomic->axisCoefficient();

                    for (int i = 0; i < _outputBasis->size(); ++i) {
                        TRACE("Want coeff[" << i << ", " << j << "] = " <<
                              (*coefficients)[j][i] << " * " << k <<
                              " of row #" << i << " "  << *(*_outputBasis)[i]);

						if (!!coefficients2 && (atomic->getApproxOrder()==2)){
							v[i] += (*coefficients)[j][i] * k + 0.5*(*coefficients2)[j][i] * k * k;
						} else {
							v[i] += (*coefficients)[j][i] * k;
						}
                    }

                    TRACE("Equivalent 'distance' denominator is " <<
                          (*inputDistanceScaling)[j] << " * " << k);

                    dists[a] += (*inputDistanceScaling)[j] * k;
                }
            }
        }
    }

    for (int i = 0; i < v.size(); ++i) {
        if (!Maths::isZero(v[i])) {
            IHypothesisConstSP hyp = (*_outputBasis)[i]->hypothesis(v[i]);
            for (int a = 0; a < hyp->numAtomics(); ++a) {
                outs->push_back(AtomicHypothesisSP::constCast(hyp->atomic(a)));
            }
        }
    }

    return outs->size() == hypothesis->numAtomics() ?
        hypothesis :
        CompoundHypothesis::SP(outs,
                               IHypothesis::IDistanceMetric::constant(
                                   (*hypothesis->distanceMetric())(dists)))->
            grouped();
}

void RiskMappingMatrix::appendMappedNames(
        IAbstractRiskPropertyConstSP property,
        OutputNameArraySP names) const {

    TRACE_METHOD;
    TRACE_BLOCK("Finding names for which RiskMappingMatrix " << getName() <<
                " can map " << *property);

    for (int j = 0; j < _inputBasis->size(); ++j) {
        if (property->equalTo((*_inputProperties)[j].get())) {
            OutputNameConstSP name = (*_inputBasis)[j]->marketDataName();
            if (!!name) {
                // duplicates are OK
                TRACE(*name << " from column #" << j << " " <<
                      *(*_inputBasis)[j]);
                names->push_back(OutputNameSP::constCast(name));
            }
        }
    }
}

template <class Q>
bool RiskMappingMatrix::appendMappedQualifiers(
        smartConstPtr<IRiskProperty<Q> > property,
        OutputNameConstSP name,
        smartPtr<array<smartPtr<Q>, Q> > qualifiers) const {

    TRACE_METHOD;
    TRACE_BLOCK("Finding qualifiers which RiskMappingMatrix " << getName() <<
                " can map for " << *property << " on name '" <<
                name->toString() << "'");

    bool any = false;

    for (int j = 0; j < _inputBasis->size(); ++j) {
        const IQualifiedRiskAxis<Q>* axis =
            dynamic_cast<IQualifiedRiskAxis<Q> *>((*_inputBasis)[j].get());
        if (axis &&
                property->equalTo((*_inputProperties)[j].get()) &&
                name->equals(axis->marketDataName().get())) {
            smartConstPtr<Q> w = axis->qualifier();
            if (!!w) {
                TRACE(*w << " from column #" << j << " " << *(*_inputBasis)[j]);
                qualifiers->push_back(smartPtr<Q>::constCast(w));
                any = true;
            }
        }
    }

    return any;
}

template bool RiskMappingMatrix::appendMappedQualifiers(
    smartConstPtr<IRiskProperty<ExpiryWindow> > property,
    OutputNameConstSP name,
    smartPtr<array<smartPtr<ExpiryWindow>, ExpiryWindow> > qualifiers) const;

template bool RiskMappingMatrix::appendMappedQualifiers(
    smartConstPtr<IRiskProperty<ExpiryPair> > property,
    OutputNameConstSP name,
    smartPtr<array<smartPtr<ExpiryPair>, ExpiryPair> > qualifiers) const;

template bool RiskMappingMatrix::appendMappedQualifiers(
    smartConstPtr<IRiskProperty<ExpiryAndStrike> > property,
    OutputNameConstSP name,
    smartPtr<array<smartPtr<ExpiryAndStrike>, ExpiryAndStrike> > qualifiers) const;

template bool RiskMappingMatrix::appendMappedQualifiers(
    smartConstPtr<IRiskProperty<BoxedInt> > property,
    OutputNameConstSP name,
    smartPtr<array<smartPtr<BoxedInt>, BoxedInt> > qualifiers) const;

IObject* RiskMappingMatrix::defaultRiskMappingMatrix() {
    return new RiskMappingMatrix();
}

void RiskMappingMatrix::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(RiskMappingMatrix, clazz);
    SUPERCLASS(MarketObject);
    EMPTY_SHELL_METHOD(defaultRiskMappingMatrix);
    FIELD(name, "name");
    FIELD(inputBasis, "inputBasis");
    FIELD(inputDistanceScaling, "inputDistanceScaling");
    FIELD(outputBasis, "outputBasis");
    FIELD(coefficients, "coefficients");
	FIELD(coefficients2, "coefficients2");
	FIELD_MAKE_OPTIONAL(coefficients2);
    FIELD(errorBars, "errorBars");

    FIELD(_inputBasis, "_inputBasis");
    FIELD_MAKE_TRANSIENT(_inputBasis);
    FIELD(_outputBasis, "_outputBasis");
    FIELD_MAKE_TRANSIENT(_outputBasis);
    FIELD(_inputProperties, "_inputProperties");
    FIELD_MAKE_TRANSIENT(_inputProperties);
}

CClassConstSP const RiskMappingMatrix::TYPE = CClass::registerClassLoadMethod(
    "RiskMappingMatrix", typeid(RiskMappingMatrix), load);

DEFINE_TEMPLATE_TYPE(RiskMappingMatrixArray);

DRLIB_END_NAMESPACE
