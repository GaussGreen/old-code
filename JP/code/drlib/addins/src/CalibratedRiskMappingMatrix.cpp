/**
 * @file CalibratedRiskMappingMatrix.cpp
 *
 * See class docs on CalibratedRiskMappingMatrix for information.
 */

#include "edginc/config.hpp"
#include "edginc/TRACE.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/ExpiryResult.hpp"
#include "edginc/ExpiryWindow.hpp"
#include "edginc/ExpiryPairResult.hpp"
#include "edginc/ExpiryPair.hpp"
#include "edginc/Model.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/RiskQuantity.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskQuantity.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/AbstractPropertyTweakHypothesis.hpp"
#include "edginc/HypotheticalQuantity.hpp"
#include "edginc/FieldPath.hpp"
#include "edginc/TweakFunction.hpp"
#include "edginc/fieldRiskProperty.hpp"
#include "edginc/VanillaGrid.hpp"
#include "edginc/VanillaGridInstrumentCollection.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"
#include "edginc/RiskMappingMatrix.hpp"
#include "edginc/RiskQuantityEvaluator.hpp"
#include "edginc/IDynamicsParameter.hpp"
#include "edginc/LazyRiskQuantityFactory.hpp"
#include "edginc/RiskMapping.hpp"
#include <set>
#include <fstream>

DRLIB_BEGIN_NAMESPACE

using std::set;
using std::ofstream;

FORWARD_DECLARE(NaturalDerivRiskQuantity)

/**
 * Index of the first object in @a xs which is IObject::equalTo() @a x, or -1
 * if none
 */

static int findEqualTo(IArrayConstSP xs, IObjectConstSP x) {
    for (int i = 0; i < xs->getLength(); ++i) {
        if (xs->get(i)->equalTo(x.get())) return i;
    }
    return -1;
}

/**
 * A RiskQuantity which computes a one-sided first derivative, expressed in
 * units relative to the nominal shift coefficient ("shift size"), rather than
 * to the actual distance reported by the hypothesis ("divisor").
 *
 * This is like the one-sided derivatives returned by
 * IScalarDerivative::oneSided(), but evalutes to
 *
 * -     (price(altworld.world) - price(world)) / coeff
 *       where altworld = axis.hypothesis(coeff).appliedTo(world).world
 *
 * rather than
 *
 * -     (price(altworld.world) - price(world)) / altworld.distance
 *
 * It also provides a method gearing() which returns
 *
 * -     altworld.distance / coeff
 *
 * The reason we need this is that the distance reported by the hypothesis can
 * be unpredictably different from the nominal coefficient---for instance,
 * Equity::sensShift(const PropertyTweak<Spot>&) treats
 * PropertyTweak<Spot>::coefficient as specifying a proportional spot shift,
 * but reports the absolute change as the TweakOutcome::distance(), so that
 * Delta ends up being expressed in absolute terms.
 *
 * So to reconstruct both a suitable coeff and a suitable distance in
 * RiskMappingMatrix::mapped(), we need an extra degree of freedom per axis.
 * That's RiskMappingMatrix::inputDistanceScaling.
 */

class NaturalDerivRiskQuantity: public RiskQuantity {

    NaturalDerivRiskQuantity(): RiskQuantity(TYPE) {}

    static IObject* emptyOne() {
        return new NaturalDerivRiskQuantity();
    }

    static void load(CClassSP& clazz) {
        REGISTER(NaturalDerivRiskQuantity, clazz);
        SUPERCLASS(RiskQuantity);
        EMPTY_SHELL_METHOD(emptyOne);
        FIELD(axis, "axis");
        FIELD(coeff, "coeff");
		FIELD(order, "order");
    }

    static HypotheticalQuantityArraySP parameters(
            IResultsFunctionConstSP function,
            IRiskAxisConstSP axis,
            double coeff,
			int order) {

		HypotheticalQuantityArraySP hs;
		if (order==2){
			hs = HypotheticalQuantityArraySP(new HypotheticalQuantityArray(3));
			(*hs)[0] = HypotheticalQuantity::SP(axis->hypothesis(0), function);
			(*hs)[1] = HypotheticalQuantity::SP(axis->hypothesis(coeff), function);
			(*hs)[2] = HypotheticalQuantity::SP(axis->hypothesis(-coeff), function);
		} else {
			hs = HypotheticalQuantityArraySP(new HypotheticalQuantityArray(2));
			(*hs)[0] = HypotheticalQuantity::SP(axis->hypothesis(0), function);
			(*hs)[1] = HypotheticalQuantity::SP(axis->hypothesis(coeff), function);
		}
        return hs;
    }

public:

    static CClassConstSP const TYPE;

    IRiskAxisConstSP axis; // $required
    double coeff; // $required
	int order;

    NaturalDerivRiskQuantity(IResultsFunctionConstSP function,
                             IRiskAxisConstSP axis,
                             double coeff,
							 int order):
        RiskQuantity(TYPE, parameters(function, axis, coeff, order)),
        axis(axis),
        coeff(coeff),
		order(order)
    {}

    static NaturalDerivRiskQuantitySP SP(IResultsFunctionConstSP function,
                                         IRiskAxisConstSP axis,
                                         double coeff,
										 int order) {
        return NaturalDerivRiskQuantitySP(new NaturalDerivRiskQuantity(
                                              function, axis, coeff, order));
    }

    double value(const DoubleArray& vals,
                 const DoubleArray& dists) const {
        TRACE_METHOD;
        if (Maths::isZero(coeff)) {
            TRACE("Zero divisor");
            throw ModelException("Zero divisor");
        }
		
		double val;

		if (order==1){
			val = (vals[1] - vals[0])/coeff;
		}else{
			val = (vals[2] + vals[1] - 2*vals[0])/(coeff*coeff);
		}

        return val;
    }

    double gearing(const CDoubleArray& dists) const {
        return Maths::isZero(coeff) ? 1 : dists[1] / coeff;
    }
};

CClassConstSP const NaturalDerivRiskQuantity::TYPE = CClass::registerClassLoadMethod(
    "NaturalDerivRiskQuantity", typeid(NaturalDerivRiskQuantity), load);

/**
 * This is like IScalarDerivative::oneSided() but returns
 * NaturalDerivRiskQuantity (q.v. above) rather than the standard one-sided
 * deriv.
 */

class NaturalDeriv: public CObject,
                    public virtual IScalarDerivative {

    static void load(CClassSP& clazz) {
        REGISTER(NaturalDeriv, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IScalarDerivative);
        EMPTY_SHELL_METHOD(DefaultConstructor<NaturalDeriv>::iObject);
    }

public:

    static CClassConstSP const TYPE;

    NaturalDeriv(): CObject(TYPE) {}

    RiskQuantitySP riskQuantity(IResultsFunctionConstSP function,
                                IRiskAxisConstSP axis,
                                double coeff) const {
        return RiskQuantitySP(
            new NaturalDerivRiskQuantity(function, axis, coeff, order()));
    }

    int order() const {
        return 1;
    }

    static IScalarDerivativeConstSP it() {
        IScalarDerivativeConstSP the(new NaturalDeriv());
        return the;
    }
};

CClassConstSP const NaturalDeriv::TYPE = CClass::registerClassLoadMethod(
    "NaturalDeriv", typeid(NaturalDeriv), load);

FORWARD_DECLARE(AtomicAxesRiskQuantityFactory)

/**
 * Adapter IRiskQuantityFactory, which takes the RiskQuantity's returned from
 * the underlying factories, and substitutes NaturalDerivRiskQuantity's for
 * all the IRiskAxis's they want to make shifts along.
 */

class AtomicAxesRiskQuantityFactory: public CObject,
                                     public virtual IRiskQuantityFactory {

    static IObject* emptyOne() {
        return new AtomicAxesRiskQuantityFactory(
                   IRiskQuantityFactoryArrayConstSP(),1);
    }

    static void load(CClassSP& clazz) {
        REGISTER(AtomicAxesRiskQuantityFactory, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IRiskQuantityFactory);
        FIELD(rqfs, "rqfs");
		FIELD(order, "order");
        EMPTY_SHELL_METHOD(emptyOne);
    }

    IRiskQuantityFactoryArrayConstSP rqfs; // $required
	int order;

public:

    static CClassConstSP const TYPE;

    AtomicAxesRiskQuantityFactory(IRiskQuantityFactoryArrayConstSP rqfs,
		int order):
        CObject(TYPE),
        rqfs(rqfs),
		order(order)
    {}

    static IRiskQuantityFactoryArraySP ArraySP(
           IRiskQuantityFactoryArrayConstSP rqfs,
		   int order) {
        return IRiskQuantityFactoryArraySP(
            new IRiskQuantityFactoryArray(
                1, 
                AtomicAxesRiskQuantityFactorySP(
                    new AtomicAxesRiskQuantityFactory(rqfs, order))));
    }

    NamedRiskQuantityArraySP riskQuantities(
            MultiTweakGroupConstSP world,
            RiskMappingConstSP riskMapping) const {
        TRACE_METHOD;

        IRiskAxisArraySP axes(new IRiskAxisArray());
        vector<double> coeffs;
        vector<IResultsIdentifierConstSP> ids;

        for (int f = 0; f < rqfs->size(); ++f) {
            NamedRiskQuantityArraySP origRqs = (*rqfs)[f]->riskQuantities(
                                                   world, riskMapping);

            for (int q = 0; q < origRqs->size(); ++q) {
                HypotheticalQuantityArrayConstSP params =
                    (*origRqs)[q]->riskQuantity->parameters();
                for (int p = 0; p < params->size(); ++p) {
                    IHypothesisConstSP hyp = (*params)[p]->hypothesis();
                    for (int a = 0; a < hyp->numAtomics(); ++a) {
                        IRiskAxisConstSP axis = hyp->atomic(a)->axis();
                        if (!!axis) {
                            double coeff = hyp->atomic(a)->axisCoefficient();
                            int x = findEqualTo(axes, axis);
                            if (x == -1) {
                                axes->push_back(IRiskAxisSP::constCast(axis));
                                coeffs.push_back(coeff);
                                ids.push_back((*origRqs)[q]->resultsName);
                            }
                            else {
                                if (coeff > coeffs[x]) {
                                    coeffs[x] = coeff;
                                }
                            }
                        }
                    }
                }
            }
        }

        NamedRiskQuantityArraySP rqs(new NamedRiskQuantityArray(axes->size()));

        for (int x = 0; x < axes->size(); ++x) {
            (*rqs)[x] = NamedRiskQuantity::SP(
                NaturalDerivRiskQuantity::SP(IResultsFunction::price(),
                                             (*axes)[x], coeffs[x], order),
                ids[x], 1.);
        }

        return rqs;
    }

    LazyRiskQuantityFactoryArraySP lazies(MultiTweakGroupConstSP world) const {
        // this may need to be implemented for delta-next-day
        return LazyRiskQuantityFactoryArraySP();
    }
};

const CClassConstSP AtomicAxesRiskQuantityFactory::TYPE = CClass::registerClassLoadMethod(
    "AtomicAxesRiskQuantityFactory", typeid(AtomicAxesRiskQuantityFactory), load);

/**
 * Greeks for one-sided tweaks of all the calibratable fields on a list of
 * objects
 *
 * What shift size do we use?  Actually ... it's hardwired to 0.1%.  That's not
 * as stupid as it sounds, though, because the CalibratorFieldTweakHypothesis
 * we use here automatically does exponential tweaks for fields defined with a
 * lower (and/or upper) bound.  Since that's true of most of the fields we're
 * interested in, the tweaks are generally
 *
 * -      field *= exp 0.1%
 *
 * rather than
 *
 * -      field += 0.1%
 */

static IRiskQuantityFactoryArraySP parameterSensitivities(
        ObjectArrayConstSP paramObjs) {

    set<CClassConstSP> paramClasses;

    IRiskQuantityFactoryArraySP sensitivities(
        new IRiskQuantityFactoryArray());

    for (int o = 0; o < paramObjs->size(); ++o) {
        CClassConstSP paramClass = (*paramObjs)[o]->getClass();

        if (paramClasses.find(paramClass) == paramClasses.end()) {
            paramClasses.insert(paramClass);

            CFieldArray fields = Calibrator::IAdjustable::getFields(paramClass);

            double coeff = 0.001;

            for (int f = 0; f < int(fields.size()); ++f) {
                CClassConstSP clazz = fields[f]->getDeclaringClass();

                string na = "price/" + clazz->getName() + "." +
                            fields[f]->getName();

                FieldPathSP field = FieldPath::SP(fields[f]->getName());
 
                if (Calibrator::IAdjustable::hasGetExpiriesMethod(fields[f])) {
                    sensitivities->push_back(
                        PerNameRiskPropertySensitivity<ExpiryWindow>::SP(
                            na,
                            IResultsFunction::price(),
                            fieldRiskProperty::pointwise(
                                clazz, field, FieldPathSP(),
                                FieldPathSP(), FieldPathSP(),
                                0,
                                IFieldTweak::IOperator::numeric(
                                    TweakFunction::adaptive(),
                                    InfiniteRange(), true, true),
                                CDouble::SP(1.),
                                false),
                            NaturalDeriv::it(),
                            1., coeff));
                }
                else {
                    sensitivities->push_back(
                        PerNameRiskPropertySensitivity<Void>::SP(
                            na,
                            IResultsFunction::price(),
                            fieldRiskProperty::scalar(
                                clazz,
                                IFieldTweak::bulk(field, FieldPathSP(),
                                                  CDouble::SP(1.),
                                                  IFieldTweak::IOperator::numeric(
                                                      TweakFunction::adaptive(),
                                                      InfiniteRange(), true, true)),
                                false),
                            NaturalDeriv::it(),
                            1., coeff));
                }
            }
        }
    }

    return sensitivities;
}

/**
 * The set of IRiskAxis's along which @a rqvs were successfully estimated
 */

static IRiskAxisArraySP riskAxesIn(RiskQuantityEvaluator::ValueArraySP rqvs,
                                   const string& description) {
    IRiskAxisArraySP axes(new IRiskAxisArray());
    string irrelevant = ", irrelevant to the calibration vanillas",
            untweakable, notApplicable, zero;

    for (int q = 0; q < rqvs->size(); ++q) {
        RiskQuantityEvaluator::Value& value = *(*rqvs)[q];
        irrelevant = "";

        if (!value.ok()) {
            if (!value.instrumentApplicable || !value.hypothesesApplicable) {
                notApplicable = ", not applicable";
            }
            else if (!!value.oops) {
                untweakable = ", untweakable";
            }
        }
        else {
            const NaturalDerivRiskQuantity* rq =
                dynamic_cast<const NaturalDerivRiskQuantity*>(
                    value.riskQuantity->riskQuantity.get());

            if (rq) {
                if (value.value() == 0) {
                    zero = ", zero deriv";
                }
                else if (findEqualTo(axes, rq->axis) == -1) {
                    axes->push_back(IRiskAxisSP::constCast(rq->axis));
                }
            }
        }
    }

    if (axes->empty()) {
        throw ModelException(
            "All the specified " + description + " greeks are unavailable "
            "(specifically" + irrelevant + untweakable + notApplicable +
            zero + ")");
    }

    return axes;
}

/**
 * Gather values computed for RiskQuantity's into a matrix, with columns
 * corresponding to a vector of IRiskAxis's
 *
 * You can see @a values as a sparse matrix, containing numbers tagged by
 * column = risk axis and row = instrument, and columnBasis as an
 * association of risk axes with actual column indices.
 *
 * Returns a dense matrix (with zeros for instrument/risk axis combinations
 * not mentioned in @a values).
 *
 * Used in CalibratedRiskMappingMatrix::operator()().
 *
 * @a values is a set of RiskQuantityEvaluator::Value's, giving derivatives of
 * the prices of a collection of <I>N</I> instruments along <I>M</I> risk axes.
 *
 *    -  As called, the @a values are generated either from the
 *       IRiskQuantityFactory's given in
 *       CalibratedRiskMappingMatrix::marketGreeks, or from those in
 *       CalibratedRiskMappingMatrix::paramGreeks.  They give the derivatives
 *       of all the vanilla prices to which we're calibrating with respect
 *       to Black-Scholes spot, vols etc. in the former case, and w.r.t.
 *       the parameters of the parametric model in the latter case
 *
 * @a columnBasis is a set of <I>M</I> risk axes extracted from @a values using
 * riskAxesIn().
 *
 *    -  Either the risk axes describing the Black-Scholes quantities w.r.t.
 *       which we've measured derivatives, or those for the model parameters
 *       ditto.
 *
 * @a scaling is optionally an array in which we store the
 * NaturalDerivRiskQuantity::gearing()'s measured for each risk axis.
 */

static CDoubleMatrixSP gathered(int insts,
                                IRiskAxisArraySP columnBasis,
                                RiskQuantityEvaluator::ValueArraySP values,
                                CDoubleArraySP scaling = CDoubleArraySP()) {
    CDoubleMatrixSP it(new DoubleMatrix(columnBasis->size(), insts));

    DoubleMatrix& m = *it;

    if (!!scaling) {
        scaling->resize(columnBasis->size());
        for (int i = 0; i < scaling->size(); ++i) (*scaling)[i] = 0;
    }

    for (int q = 0; q < values->size(); ++q) {
        RiskQuantityEvaluator::Value& qv = *(*values)[q];
        const NaturalDerivRiskQuantity* rq =
            dynamic_cast<const NaturalDerivRiskQuantity*>(
                qv.riskQuantity->riskQuantity.get());

        if (rq && qv.ok()) {
            int a = findEqualTo(columnBasis, rq->axis);
            if (a != -1) {
                m[a][qv.instrument] = qv.value();
                if (!!scaling) (*scaling)[a] = rq->gearing(*qv.dsts);
            }
        }
    }

    return it;
}

/**
 * x + diag(d, d, d, ...)
 */

static CDoubleMatrixSP plusDiag(const DoubleMatrix& x,
                                double d) {
    ASSERT(x.numRows() == x.numCols());

    CDoubleMatrixSP it(new DoubleMatrix(x));

    for (int i = 0; i < x.numRows(); ++i) {
        (*it)[i][i] += d;
    }

    return it;
}

/**
 * Pseudo-inverse with elimination of big singular values.
 *
 * Computes effectively (x^T x)^{-1} x^T but zeroing out all singular
 * values which are bigger than svthreshold.
 */

static CDoubleMatrixSP truncPseudoInverse(CDoubleMatrixConstSP x,
                                          double svthreshold,
                                          double regulariser) {
    try {
        SVDAnalysisSP svd = x->computeSVD();

        // ensure I've got this the right way round ...
        ASSERT(svd->matrixU->numCols() == svd->sValues->size());

        for (int j = 0; j < svd->sValues->size(); ++j) {
            double v = (*svd->sValues)[j];
            double vv = v > 1/svthreshold ? 1 / (v + regulariser) : 0.;
            for (int i = 0; i < svd->matrixU->numRows(); ++i) {
                (*svd->matrixU)[j][i] *= vv;
            }
        }

        CDoubleMatrixSP it(new DoubleMatrix(
            svd->matrixU->mult(*svd->matrixVTranspose)));
        it->transpose();
        return it;
    }
    catch (exception& e) {
        throw ModelException(e, "truncPseudoInverse");
    }
}

/**
 * diag(weights) * x
 */

static CDoubleMatrixSP scaledRows(CDoubleArrayConstSP weights,
                                  const DoubleMatrix& x) {
    ASSERT(weights->size() == x.numRows());

    CDoubleMatrixSP wx(new DoubleMatrix(x));
    for (int i = 0; i < x.numRows(); ++i) {
        for (int j = 0; j < x.numCols(); ++j) {
            (*wx)[j][i] *= (*weights)[i];
        }
    }

    return wx;
}

/**
 * Casts @a senss to array of IRiskQuantityFactory
 */

static IRiskQuantityFactoryArraySP castSensitivitys(SensitivityArraySP senss,
                                                    string name) {
    if (senss->empty()) {
        throw ModelException("'" + name + "' must be a nonempty list of "
                             "sensitivities");
    }

    IRiskQuantityFactoryArraySP them(new IRiskQuantityFactoryArray());

    for (int b = 0; b < senss->size(); ++b) {
        if (!(*senss)[b]) {
            throw ModelException(
                "Entry number " + Format::toString(b) + " in the sensitivities "
                "list '" + name + "' is null");
        }

        if (!IRiskQuantityFactory::TYPE->isInstance(
                                             (*senss)[b])) {
            throw ModelException(
                "'" + name + "' must all be IRiskQuantityFactory's, "
                "(i.e. written in/ported to the 'declarative' framework), but " +
                (*senss)[b]->getClass()->getName() + " isn't");
        }

        them->push_back(
            IRiskQuantityFactorySP::dynamicCast((*senss)[b]));
    }

    return them;
}

/**
 * Addin entry point for generating a RiskMappingMatrix.
 *
 * This is the entry point for the EAS batch process which generates a RiskMappingMatrix
 * to go with a freshly calibrated model parameter (like VolSVJ or SRMEQ::Vol).  See
 * RiskMapping for an overview of the facility.
 *
 * The actual maths is done in operator ()().
 */

class CalibratedRiskMappingMatrix: public CObject,
                                   public virtual ClientRunnable { // $public

protected: // so that Doxygen looks at Arguments

    /**
     * @name Arguments
     */

    //@{

    /**
     * What to call the returned RiskMappingMatrix (its MarketObject::getName())
     *
     * This can be anything the client wants: it doesn't have to be
     * systematically related to any other market data name.
     */

    string name; // $required(Name for returned RiskMappingMatrix)

    /**
     * The Calibrator::ObjFunc which was used to calibrate the dynamics
     * parameter object for which we're making a RiskMappingMatrix
     *
     * Currently must be a VanillaGrid::LeastSquaresSimple, or a
     * Calibrator::ObjFunc::ObjFuncLeastSquareCombo containing a
     * VanillaGrid::LeastSquaresSimple [the former is what we get from Merlin
     * calibration, the latter what we get from the Aladdin calibrator sheet
     * ...].  Used to obtain the model and the calibration instruments
     * (i.e. vanilla grid).
     */

    // $required(Objective function to which parametric model was calibrated)
    Calibrator::ObjFuncSP objFunction;

    /**
     * Threshold constant to help prevent unreasonably large coefficients
     * from appearing in the risk mapping matrix
     *
     * The basic regression calculation we do
     *
     * -     RMM = A! (dV/dBS - dV_param/dBS)
     *
     * where A! is the pseudo-inverse of the matrix of sensitivities of vanilla
     * prices to model parameters.
     *
     * -     A = dV_param/dParam
     * -     A! = (A^T A) \ A^T
     *
     * You can see A! as giving the sensitivities of the model parameters to
     * vanilla prices --- how they would change under recalibration to
     * different market conditions.  So you can see each column of the RMM
     * as encapsulating a "tweak the market and recalibrate" operation.
     *
     * Problems arise, though, around (combinations of) model parameters which
     * have very little influence on any vanilla prices, such as, it turns out,
     * SRMEQ::smileA3.  A simple regression scheme will try to get the
     * prices to come out right even at the cost of making giant changes
     * those parameters (to compensate for their tiny sensitivities).  And
     * that will take them outside any kind of sensible range, or at least
     * outside the zone in which our linearity assumption holds good.
     *
     * To avoid that we explicitly find all A!'s singular values which are
     * bigger than @a svthreshold and set them to zero.  Because they correspond
     * precisely to parameter combinations with little effect on vanilla prices,
     * the infinitesimal sensitivities are still OK and the finite-shift
     * sensitivities are much better.
     */

    // $optional(Threshold on singular values of dParam/dVanillaPrices, defaulting to 10)
    double svthreshold;

    // $transient(Regulariser on singular values of dParam/dVanillaPrices, defaulting to 0)
    double regulariser;

    /**
     * The Black-Scholes greeks which the RiskMappingMatrix should support
     *
     * They must actually be IRiskQuantityFactory's, i.e. Sensitivity's written
     * in/ported to the 'declarative' sensitivities framework.  Currently that
     * includes Delta, VegaPointwise, VegaParallel, ParSpreadRhoParallel, ...
     *
     * The RiskMappingMatrix returned will have columns corresponding to each
     * kind of tweak used in computing the marketGreeks---i.e. to each IRiskAxis
     * along which they request a shift.
     *
     * marketGreeks which are not applicable to the calibration vanillas are
     * silently ignored.  An exception is raised if none of the marketGreeks
     * are applicable (so if a RiskMappingMatrix is successfully returned,
     * it always has at least one column).
     *
     * The marketGreeks are not applied directly to the vanillas over which we
     * calibrate: instead, they are interrogated for the RiskQuantity's they
     * want to compute (via IRiskQuantityFactory::riskQuantities()), and the
     * IRiskAxis's used in constructing those quantities are collated.  Then
     * one sided greeks, normalised by the nominal "tweak coefficient" rather
     * than the effective distance, are measured for each axis.  See
     * AtomicAxesRiskQuantityFactory above.
     *
     * The point of this is
     *
     *    -  we can take out the obfuscating effect of the (in principle
     *       arbitrary) relationship between coefficient and divisor,
     *       without having to introduce a confusing back door into the
     *       generic RiskQuantity;
     *
     *    -  the user can directly specify the greeks they want the matrix to
     *       support (potentially DDeltaDVol), without having to break them
     *       down into first-order single-property greeks
     *       (Delta, VegaParallel).
     */

    SensitivityArraySP marketGreeks; // $required(Greeks to map from)

	SensitivityArraySP marketGreeks2; // (second order Greeks to map from)

    /**
     * The parameter greeks in terms of which the RiskMappingMatrix should
     * support the marketGreeks
     *
     * If absent, defaults to "all fields registered with
     * Calibrator::IAdjustable::registerField() on all the IDynamicsParameter
     * objects which are reachable from the Model and VanillaGrid defined in
     * objFunction".  If you want to specify the fields you want to use
     * manually, you can do so by giving a list of FlexibleSensitivity's.
     *
     * The RiskMappingMatrix returned will have rows corresponding to each
     * field.
     *
     * paramGreeks which are not applicable to the model specified in
     * objFunction are silently ignored.  An exception is raised if none of
     * them are applicable (so if a RiskMappingMatrix is successfully returned,
     * it always has at least one row).
     */

    // $optional(Greeks to map to (determined automatically if absent\))
    SensitivityArraySP paramGreeks;

    /**
     * The same market data as was provided for the Calibrator run, plus the
     * calibrated dynamics parameters object (VolSVJ, SRMEQ::Vol, ...) it
     * produced
     */

    CMarketDataSP market; // $required(Market data)

    //@}

private:

    /**
     * "Global" variables (just for communicating between validatePop2Object(),
     * getMarket() and operator()())
     */

    //@{

    VanillaGrid::LeastSquareSimpleSP vgObjFunc;                // $transient
    IModelSP bsModel, paModel;                                 // $transient
    VanillaGridInstrumentCollectionSP bsVanillas, paVanillas;  // $transient
    IRiskQuantityFactoryArraySP bsTweaks;                      // $transient
	IRiskQuantityFactoryArraySP bsTweaks2;                      // $transient
    IRiskQuantityFactoryArraySP paTweaks;                      // $transient

    //@}

    CalibratedRiskMappingMatrix():
        CObject(TYPE),
        svthreshold(10),
        regulariser(0)
    {}

public:

    static CClassConstSP const TYPE;

    void validatePop2Object() {
        try {
            // This is not great

            vgObjFunc.reset(
                dynamic_cast<VanillaGrid::LeastSquareSimple*>(objFunction.get()));

            if (!vgObjFunc) {
                Calibrator::ObjFuncLeastSquareCombo* lsc =
                    dynamic_cast<Calibrator::ObjFuncLeastSquareCombo*>(
                        objFunction.get());
                if (lsc) {
                    Calibrator::ObjFuncLeastSquareArray ofs =
                        lsc->getObjFuncArray();
                    for (int o = 0; o < ofs.size(); ++o) {
                        vgObjFunc.reset(
                            dynamic_cast<VanillaGrid::LeastSquareSimple*>(
                                ofs[o].get()));
                        if (!!vgObjFunc) break;
                    }
                }

                if (!vgObjFunc) {
                    throw ModelException(
                        "The only kind of objective function which "
                        "CalibratedRiskMappingMatrix understands "
                        "at the moment is VanillaGrid::LeastSquareSimple "
                        "(or ObjFuncLeastSquareCombo containing a "
                        "VanillaGrid::LeastSquareSimple), so you can't "
                        "use " +
                        (!objFunction ? string("NULL") :
                                        objFunction->getClass()->getName()));
                }
            }

            paModel = vgObjFunc->getModel();

            bsTweaks = castSensitivitys(marketGreeks, "marketGreeks");

            if (!!paramGreeks) {
                paTweaks = castSensitivitys(paramGreeks, "paramGreeks");
            }

			if (!!marketGreeks2) {

				for (int i=0; i<marketGreeks2->size(); i++){
					if (findEqualTo(marketGreeks,(*marketGreeks2)[i])==-1){
						 throw ModelException("One of the second order sensitivity is no among the first order sensitivities");
					}
				}
                bsTweaks2 = castSensitivitys(marketGreeks2, "marketGreeks2");
            }
        }
        catch (exception& e) {
            throw ModelException(e, "CalibratedRiskMappingMatrix::validatePop2Object()");
        }
    }

    void getMarket() {
        try {
            vgObjFunc->getMarket(market.get());

            bsVanillas.reset(new VanillaGridInstrumentCollection(
                                 vgObjFunc->getVanillaGrid()));

            if (bsVanillas->size() == 0) {
                throw ModelException(
                    "The calibration grid of vanillas specified by objFunction "
                    "is empty");
            }

            paVanillas.reset(new VanillaGridInstrumentCollection(
                                 vgObjFunc->getVanillaGrid()));

            bsModel.reset(new CClosedFormLN("VolSurface", true));
            bsModel->getInstrumentsAndModelMarket(market, bsVanillas);
            paModel->getInstrumentsAndModelMarket(market, paVanillas);

            if (!paTweaks) {
                ObjectArrayConstSP params = SensMgrConst(
                    MultiTweakGroup::SP(paVanillas, paModel)).
                        all(IDynamicsParameter::TYPE, 0);

                if (params->empty()) {
                    throw ModelException(
                       "The model specified in 'objFunction' is a '" +
                       paModel->getClass()->getName() + "', which is not a "
                       "parametric model --- at least I can't find a "
                       "designated 'IDynamicsParameter' vol object in it");
                }

                if (!MarketObject::TYPE->isInstance((*params)[0])) {
                    throw ModelException(
                        "Found a parametric-model vol object of type " +
                        (*params)[0]->getClass()->getName() +
                        " but it's not a MarketObject");
                }

                paTweaks = parameterSensitivities(params);

                if (paTweaks->size() == 0) {
                    throw ModelException(
                        "No tweakable fields were found on the "
                        "'IDynamicsParameter' vol object(s)");
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, "CalibratedRiskMappingMatrix::getMarket()");
        }
    }

    /**
     * Maths of risk mapping calibration
     *
     * See the presentation in the QLib doc db article "Risk mapping" for
     * an introduction to the methodology.
     */

    RiskMappingMatrixSP operator()() {
        try {
            typedef RiskQuantityEvaluator::ValueArraySP RQVs;

            // Estimate derivs of Black-Scholes vanilla prices wrt
            // Black-Scholes market data

            RQVs vByBS_ = RiskQuantityEvaluator().values(
                AtomicAxesRiskQuantityFactory::ArraySP(bsTweaks,1),
                MultiTweakGroup::SP(bsVanillas, bsModel));
            IRiskAxisArraySP bsBasis = riskAxesIn(vByBS_, "market");
            DoubleArraySP bsScaling(new DoubleArray());
            CDoubleMatrixSP vByBS = gathered(bsVanillas->size(), bsBasis,
                                             vByBS_, bsScaling);

            // Estimate derivs of parametric vanilla prices wrt Black-Scholes
            // market data (with model params held constant)

            RQVs vByBS_param_ = RiskQuantityEvaluator().values(
                AtomicAxesRiskQuantityFactory::ArraySP(bsTweaks,1),
                MultiTweakGroup::SP(paVanillas, paModel));
            CDoubleMatrixSP vByBS_param = gathered(bsVanillas->size(), bsBasis,
                                                   vByBS_param_);
                        
            // Gather the differences dV/dBS - dV_param/dBS

            DoubleMatrix resid = vByBS->minus(*vByBS_param);

            // On the other side, estimate derivs of parametric vanilla prices
            // wrt model parameters

            RQVs vByPa_ = RiskQuantityEvaluator().values(
                AtomicAxesRiskQuantityFactory::ArraySP(paTweaks,1),
                MultiTweakGroup::SP(paVanillas, paModel));
            IRiskAxisArraySP paBasis = riskAxesIn(vByPa_, "parameter");
            CDoubleMatrixSP vByPa = gathered(paVanillas->size(), paBasis,
                                             vByPa_);

            // Construct RMM = dV_param/dParam \ (dV/dBS - dV_param/dBS)
            // or more precisely = (A^T A) \ A^T (dV/dBS - dV_param/dBS)
            // where A = dV_param/dParam, but zeroing out singular values
            // in (A^T A) \ A^T which are above svthreshold.

            DoubleMatrix coeffs = truncPseudoInverse(
                  scaledRows(paVanillas->weights(), *vByPa),
                  svthreshold, regulariser
                )->mult(*scaledRows(paVanillas->weights(), resid));

            // Figure out error bars (experimental)

            DoubleMatrix vByBS_mapped = vByPa->mult(coeffs);

            DoubleArraySP errorBars(new DoubleArray(bsBasis->size()));
            for (int a = 0; a < bsBasis->size(); ++a) {
                double e = 0;
                for (int i = 0; i < paVanillas->size(); ++i) {
                    e = max(e, fabs(vByBS_mapped[a][i] - resid[a][i]));
                }
                (*errorBars)[a] = min(e, 99.e99); // avoid #Inf
            }

			// second order
			if(!!bsTweaks2){
				            
				// Estimate derivs of Black-Scholes vanilla prices wrt
				// Black-Scholes market data (second order)

				RQVs vByBS2_ = RiskQuantityEvaluator().values(
					AtomicAxesRiskQuantityFactory::ArraySP(bsTweaks2,2),
					MultiTweakGroup::SP(bsVanillas, bsModel));
				IRiskAxisArraySP bsBasis2 = riskAxesIn(vByBS2_, "market");
				DoubleArraySP bsScaling2(new DoubleArray());
				CDoubleMatrixSP vByBS2 = gathered(bsVanillas->size(), bsBasis2,
												vByBS2_, bsScaling2);

				// Estimate derivs of parametric vanilla prices wrt Black-Scholes
				// market data (second order with risk mapping)

				RiskMappingMatrixArraySP riskMappingMatrices(new RiskMappingMatrixArray(1));
				(*riskMappingMatrices)[0] = RiskMappingMatrix::SP(name, bsBasis, bsScaling, paBasis,
																	coeffs, errorBars);
				RiskMappingSP riskMapping(new RiskMapping(riskMappingMatrices));

				RQVs vByBS2_param_ = RiskQuantityEvaluator(riskMapping).values(
					AtomicAxesRiskQuantityFactory::ArraySP(bsTweaks2,2),
					MultiTweakGroup::SP(paVanillas, paModel));
				CDoubleMatrixSP vByBS2_param = gathered(bsVanillas->size(), bsBasis2,
													vByBS2_param_);

				// Gather the differences

				DoubleMatrix resid2 = vByBS2->minus(*vByBS2_param);

				// Do linear regression of residuals

				DoubleMatrix coeffs2 = truncPseudoInverse(
					scaledRows(paVanillas->weights(), *vByPa),
					svthreshold, regulariser
					)->mult(*scaledRows(paVanillas->weights(), resid2));

				DoubleMatrix coeffs2Full(coeffs.numCols(),coeffs.numRows());

				for (int i=0; i < coeffs.numCols(); i++){
					int found = findEqualTo(bsBasis2, (*bsBasis)[i]);
					for (int j=0; j < coeffs.numRows(); j++){
						if (found==-1){
							coeffs2Full[i][j] = 0;
						} else {
							coeffs2Full[i][j] = coeffs2[found][j];
						}
					}
				}

				return RiskMappingMatrix::SP(name, bsBasis, bsScaling, paBasis,
					coeffs, coeffs2Full, errorBars);

			}else{
				return RiskMappingMatrix::SP(name, bsBasis, bsScaling, paBasis,
					coeffs, errorBars);
			}
        }
        catch (exception& e) {
            throw ModelException(e, "CalibratedRiskMappingMatrix::()");
        }
    }

    /** The 'addin function' - delegates to object */
    static IObjectSP addinEntryPoint(CalibratedRiskMappingMatrix* params) {
        return params->run();
    }

    /**
     * Entry point for clients
     *
     * Returns a calibrated RiskMappingMatrix.
     */

    IObjectSP run() {
        getMarket();
        return (*this)();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(CalibratedRiskMappingMatrix, clazz);
        clazz->setPublic();
        SUPERCLASS(CObject);
        FIELD(name, "Name for returned RiskMappingMatrix");
        FIELD(objFunction, "Objective function to which parametric model was calibrated");
        FIELD(svthreshold, "Threshold on singular values of dParam/dVanillaPrices, defaulting to 10");
        FIELD_MAKE_OPTIONAL(svthreshold);
        FIELD(regulariser, "Regulariser on singular values of dParam/dVanillaPrices, defaulting to 0");
        FIELD_MAKE_OPTIONAL(regulariser);
        FIELD(marketGreeks, "Greeks to map from");
		FIELD(marketGreeks2, "Greeks2 to map from");
		FIELD_MAKE_OPTIONAL(marketGreeks2);
        FIELD(paramGreeks, "Greeks to map to (determined automatically if absent)");
        FIELD_MAKE_OPTIONAL(paramGreeks);
        FIELD(market, "Market data");
        FIELD(vgObjFunc, "vgObjFunc");
        FIELD_MAKE_TRANSIENT(vgObjFunc);
        FIELD(bsModel, "bsModel");
        FIELD_MAKE_TRANSIENT(bsModel);
        FIELD(paModel, "paModel");
        FIELD_MAKE_TRANSIENT(paModel);
        FIELD(bsVanillas, "bsVanillas");
        FIELD_MAKE_TRANSIENT(bsVanillas);
        FIELD(paVanillas, "paVanillas");
        FIELD_MAKE_TRANSIENT(paVanillas);
        FIELD(bsTweaks, "bsTweaks");
        FIELD_MAKE_TRANSIENT(bsTweaks);
        FIELD(paTweaks, "paTweaks");
        FIELD_MAKE_TRANSIENT(paTweaks);
        EMPTY_SHELL_METHOD(emptyShell);

        Addin::registerClassObjectMethod(
            "CALIBRATED_RISK_MAPPING_MATRIX",
            Addin::RISK,
            "Compute sensitivity of model parameters to BS tweaks",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)addinEntryPoint);
    }

    static IObject* emptyShell() {
        return new CalibratedRiskMappingMatrix();
    }
};

CClassConstSP const CalibratedRiskMappingMatrix::TYPE = CClass::registerClassLoadMethod(
    "CalibratedRiskMappingMatrix", typeid(CalibratedRiskMappingMatrix), load);

bool CalibratedRiskMappingMatrixLoad() {
    return CalibratedRiskMappingMatrix::TYPE != NULL;
}

DRLIB_END_NAMESPACE
