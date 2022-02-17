#include "edginc/config.hpp"
#include "edginc/DiscreteDistribution.hpp"
#include "edginc/Function.hpp"

DRLIB_BEGIN_NAMESPACE

DEFINE_TEMPLATE_TYPE(DiscreteDistributionArray);

/** TYPE (for reflection) */        
CClassConstSP const DiscreteDistribution::TYPE =
CClass::registerClassLoadMethod(
    "DiscreteDistribution",
    typeid(DiscreteDistribution),
    DiscreteDistribution::load);

/** Internal constructor */
DiscreteDistribution::DiscreteDistribution(
    DoubleArraySP values,
    DoubleArraySP probabilities):
        CObject(TYPE), values(values), probabilities(probabilities), delta(1.0)
{
    validatePop2Object();
}

/** Constructor for binary distributions {(0,1-prob); (value, prob)}*/
DiscreteDistribution::DiscreteDistribution(
    double value,
    double probability): CObject(TYPE), delta(1.0)
{
    values.reset(new DoubleArray(2));
    probabilities.reset(new DoubleArray(2));
    
    (*values)[0] = 0.0;
    (*probabilities)[0] = 1.0 - probability;

    (*values)[1] = value;
    (*probabilities)[1] = probability;

    validatePop2Object();
}

/** Virtual destructor */
DiscreteDistribution::~DiscreteDistribution()
{
}

/** Only build instances of that class using reflection */
DiscreteDistribution::DiscreteDistribution():
    CObject(TYPE), values(0), probabilities(0), delta(1.0)
{
}

void DiscreteDistribution::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(DiscreteDistribution, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IDistribution1D);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(values,        "Values in the distribution");
    FIELD(probabilities, "Probabilities assigned to each value");
    FIELD(delta, "minimum difference between 2 consecutive values");
    FIELD_MAKE_TRANSIENT(delta);
}


IObject* DiscreteDistribution::defaultConstructor()
{
    return new DiscreteDistribution();
}


/** Access values */
DoubleArrayConstSP DiscreteDistribution::getValues() const
{
    return values;
}

/** Access probabilities */
DoubleArrayConstSP DiscreteDistribution::getProbabilities() const
{
    return probabilities;
}

/** Update probability */
void DiscreteDistribution::setProbability(int index, double probability)
{
    (*probabilities)[index] = probability;
}

 /**
 * Bisection search for "x" in field "values", using a "tolerance" such that:
 * - return "true" if found
 *   (in that case output "index" is such that |x - (*values)[index])| < tolerance
 * - return "false" if not found
 *   (in that case output "index" is such that (*values)[index-1] < x - tolerance and
 *    x + tolerance < (*values)[index])
 * 
 * NB: this method internally checks that tolerance < delta / 2 where delta is 
 *     the minimum difference between 2 consecutive values
 *     (this ensures we can find at MOST one x such that |x - (*values)[index])| < tolerance)
 * */
bool DiscreteDistribution::search(double x, int& index, double tolerance) const
{
    if (tolerance >= (delta / 2))
    {
        throw ModelException(
            "DiscreteDistribution::search",
            "Tolerance (" +
            Format::toString(tolerance) +
            ") is too big (should be below " + 
            Format::toString(delta / 2) +
            ")");            
    }
    
    // bisection search for x in the (sorted) array "values"
    int s = 0;
    int e = values->size() - 1;
    int mid;

    // check if x is inside range of values
    if (x < (*values)[s] - tolerance)
    {
        // not found
        index = s;
        return false;
    }
    if (x > (*values)[e] + tolerance)
    {
        // not found
        index = e + 1;
        return false;
    }
    
    // bisection loop
    do {
        if (abs(x - (*values)[s]) <= tolerance)
        {
            // found !
            index = s;
            return true;
        }
        if (abs(x - (*values)[e]) <= tolerance)
        {
            // found !
            index = e;
            return true;
        }

       mid = (s + e) / 2;
       if (x <= (*values)[mid])
       {
            e = mid;
       }
       else
       {
            s = mid + 1;
       }
    } while (e > s);

    // not found (NB: here e==s)
    if (x < (*values)[s])
    {
        index = s;
    }
    else
    {
        index = e + 1;
    }
    return false;
}

/**
 * Probability density function
 * [implements IDistribution1D]
 * */
double DiscreteDistribution::pdf(double x) const
{
    int index;
    if (search(x, index))
    {
        // found x in values
        return (*probabilities)[index];
    }
    else
    {
        // did not find x in values
        return 0.0;
    }
}

/**
 * Cumulative density function
 * [implements IDistribution1D]
 * */
double DiscreteDistribution::cdf(double x) const
{
    double res = 0.0;
    for (int i = 0; i < values->size() && (*values)[i] <= x; ++i) {
		res += (*probabilities)[i];
	}
    
    return res;
}

/**
 * nth derivative of "moment generating function" mgf(z)=E[exp(zX)]
 * [implements IDistribution1D]
 * */
double DiscreteDistribution::mgf(double z, int n) const
{
    double x;
    double res = 0.0;
    for (int i = 0; i < values->size(); ++i) {
        x = (*values)[i];
        res += pow(x,n) * exp(z * x) * (*probabilities)[i];
    }
    
    return res;
}

/**
 * Returns the expectation value
 * [implements IDistribution1D]
 * */
double DiscreteDistribution::expectation() const
{
    double res = 0.0;
    for (int i = 0; i < values->size(); ++i) {
        res += (*values)[i] * (*probabilities)[i];
    }
    
    return res;
}


/**
 * Returns the expectation value
 * [implements IDistribution1D]
 * */
double DiscreteDistribution::expectation(const Function1DDouble& payoff) const
{
    double res = 0.0;
    for (int i = 0; i < values->size(); ++i) {
        res += payoff((*values)[i]) * (*probabilities)[i];
    }
    
    return res;
}


/**
 * Returns the variance
 * [implements IDistribution1D]
 * */
double DiscreteDistribution::variance() const
{
    double res = 0.0;
    for (int i = 0; i < values->size(); ++i) {
        res += (*values)[i] * (*values)[i] * (*probabilities)[i];
    }

    double mean = expectation();
    
    return res - mean*mean;
}

/**
 * Returns a "discretised" version of itself
 * [implements IDistribution1D]
 * */
DiscreteDistributionConstSP DiscreteDistribution::discretise() const
{
    return DiscreteDistributionConstSP(this);
}

/** Called immediately after object constructed */
void DiscreteDistribution::validatePop2Object()
{
    static const string method("DiscreteDistribution::validatePop2Object");        
    try
    {
        // checks "values" and "probabilities" have the same length
        int nbValues = values->size();
        if (nbValues != probabilities->size())
        {
            throw ModelException(
                method,
                "Fields 'values' and 'probabilities' don't have "
                "the same number of elements.");            
        }
        
        // checks that values are increasing
        int i;
        for (i = 1; i < nbValues; ++i)
        {
            if ((*values)[i-1] >= (*values)[i]) {
                throw ModelException(
                    method,
                    "Values are not sorted in strictly increasing order.");            
            }
		}
        
        // computes "delta" = minimum difference between 2 consecutive values
        delta = 1.0; // initialisation needed in case of just 1 value (1.0 is arbitrary choice)
        for (i = 1; i < nbValues; ++i)
        {
            delta = min(delta, (*values)[i] - (*values)[i-1]);
        }
        
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/* external symbol to allow class to be forced to be linked in */
bool DiscreteDistributionLoad(){
    return (DiscreteDistribution::TYPE != 0);
}

DRLIB_END_NAMESPACE
