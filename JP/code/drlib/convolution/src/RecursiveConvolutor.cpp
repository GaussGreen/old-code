//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : RecursiveConvolutor.cpp
//
//   Description : Convolution Algorithm
//
//   Date        : June 2006
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/RecursiveConvolutor.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CCMLossUnit.hpp"
#include <map>
#include "edginc/DiscreteDistribution.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/FlatCDO2LossConfig.hpp"
#include "edginc/CDOPortfolio.hpp"
#include "edginc/SingleCreditAsset.hpp"

DRLIB_BEGIN_NAMESPACE

// Class responsible to perform the recursive convolution.
// Would probably internally represent the "DiscreteDistribution"s in a different way
// (discretising "values" according to loss unit).

/** Constructor */
RecursiveConvolutor::RecursiveConvolutor( 
                        int                         maxNbSlice,         // max nb of slices for loss distribution
                        int                         nbSubdivision,      // nb by which to divide the GCD of notionals
                        double                      lossUnitOverride,   // loss unit override
                        int                         discretisationType):// discretisation type
    CObject(TYPE),
    maxNbSlice(maxNbSlice),
    nbSubdivision(nbSubdivision),
    lossUnitOverride(lossUnitOverride),
    discretisationType(discretisationType)
{
}

/** empty Constructor */
RecursiveConvolutor::RecursiveConvolutor():
    CObject(TYPE),
    maxNbSlice(MAX_NB_SLICE),
    nbSubdivision(NB_SUBDIVISION),
    lossUnitOverride(LOSS_UNIT_OVERRIDE),
    discretisationType(DISCRETISATION_TYPE)
{
}

IObject* RecursiveConvolutor::defaultConstructor()
{
    return new RecursiveConvolutor();
}


/* Destructor */
RecursiveConvolutor::~RecursiveConvolutor()
{}

/** TYPE (for reflection) */        
CClassConstSP const RecursiveConvolutor::TYPE =
CClass::registerClassLoadMethod(
    "RecursiveConvolutor",
    typeid(RecursiveConvolutor),
    RecursiveConvolutor::load);

const int RecursiveConvolutor::MAX_NB_SLICE = 5000;
const int RecursiveConvolutor::NB_SUBDIVISION = 1;
const double RecursiveConvolutor::LOSS_UNIT_OVERRIDE = -1;
const int RecursiveConvolutor::DISCRETISATION_TYPE = 3;

void RecursiveConvolutor::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(RecursiveConvolutor, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IConvolutor);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(maxNbSlice,
        "Maximum number of slices for loss distribution "
        "[default value is " +
        Format::toString(MAX_NB_SLICE) +
        "]");
    FIELD_MAKE_OPTIONAL(maxNbSlice);
    FIELD(nbSubdivision,
        "Number by which to divide the GCD of notionals "
        "[default value is " +
        Format::toString(NB_SUBDIVISION) +
        "]");
    FIELD_MAKE_OPTIONAL(nbSubdivision);
    FIELD(lossUnitOverride,
        "Loss unit override "
        "[default value is " +
        Format::toString(LOSS_UNIT_OVERRIDE) +
        "]");
    FIELD_MAKE_OPTIONAL(lossUnitOverride);
    FIELD(discretisationType,
        "Discretisation flag (1 for lower, 2 for upper, 3 for split) "
        "[default value is " +
        Format::toString(DISCRETISATION_TYPE) +
        "]");
    FIELD_MAKE_OPTIONAL(discretisationType);
}

/**
 * "Pure" convolution algorithm that returns the whole convoluted distribution
 * [Implements IConvolutor]
 * */
// Result could be cached if this method is called more than once.
IDistribution1DConstSP RecursiveConvolutor::convolute(IDistribution1DArrayConstSP distributions) const
{
    return convolute(distributions, 0);
}

/**
 * "Pure" convolution algorithm that returns the whole convoluted distribution
 * up to the cutoff loss points after which all losses are aggregated
 * [Implements IConvolutor]
 * */
// Result could be cached if this method is called more than once.
IDistribution1DConstSP RecursiveConvolutor::convolute(IDistribution1DArrayConstSP distributions, double* cutoff) const
{
    int nbDist = (*distributions).size();

    // compute the loss unit or take the loss unit given as an override
    double lu = lossUnitOverride>0?lossUnitOverride:lossUnit(distributions);

    // compute the max and min loss units of each name and of the convoluted distribution
    int maxConvDistLossUnits = 0;   // max nb of loss units in the convolution distribution
    int minConvDistLossUnits = 0;   // min nb of loss units in the convolution distribution
    int maxMaxDistLossUnits = 0;    // max nb of loss units for a single input distribution
    int minMinDistLossUnits = 0;    // min nb of loss units for a single input distribution
    IntArray maxDistLossUnits(nbDist,0); // maxDistLossUnits[i] is the max nb of loss units in the input distribution #i
    IntArray minDistLossUnits(nbDist,0); // minDistLossUnits[i] is the min nb of loss units in the input distribution #i

    // local variables to store the distribution values and probabilities
    DiscreteDistributionConstSP tmpDiscreteDist;
    DoubleArrayConstSP tmpValues;
    DoubleArrayConstSP tmpProbabilities;

    int dcutoff = 0;
    for (int i=0; i<nbDist; ++i)
    {
        tmpValues = (*distributions)[i]->discretise()->getValues();
        double maxDistValue = (*tmpValues)[tmpValues->size() - 1];
        double minDistValue = (*tmpValues)[0];
        int minl2 = 0;
        int minl1 = 0;
        int maxl2 = 0;
        int maxl1 = 0;
        
        // discretisation type defines how the double value in the input distribution
        // is split onto integer amounts of loss units
        switch (discretisationType)
        {
        case 1: // discretisation is done on the closest lower integer amount of loss units.
                // In that case the expected loss might not be conserved
            maxl1 = (int) floor(maxDistValue / lu);
            minl1 = (int) floor(minDistValue / lu);
            maxDistLossUnits[i] = maxl1;
            minDistLossUnits[i] = minl1;
            break;
        case 2: // discretisation is done on the closest upper integer amount of loss units.
                // In that case the expected loss might not be conserved
            maxl2 = (int) ceil(maxDistValue / lu);
            minl2 = (int) ceil(minDistValue / lu);
            maxDistLossUnits[i] = maxl2;
            minDistLossUnits[i] = minl2;
            break;
        case 3: // discretisation is done on the 2 closest integer amount of loss units 
                // such that the expected loss is conserved
            minl1 = (int) floor(minDistValue / lu);
            maxl2 = (int) ceil(maxDistValue / lu);
            maxDistLossUnits[i] = maxl2;
            minDistLossUnits[i] = minl1;
            break;
        default:
            break;
            // throw error
            ;
        }
        // compute the cutoff in terms of loss units
        if (cutoff)
        {
            dcutoff = (int) ceil(*cutoff / lu);
        }

        if (i>0)
        {
            maxMaxDistLossUnits = Maths::max( maxMaxDistLossUnits, maxDistLossUnits[i]);
            minMinDistLossUnits = Maths::min( minMinDistLossUnits, minDistLossUnits[i]);
        }
        else
        {
            maxMaxDistLossUnits = maxDistLossUnits[i];
            minMinDistLossUnits = minDistLossUnits[i];                
        }
        maxConvDistLossUnits = Maths::max(maxConvDistLossUnits, maxDistLossUnits[i] + maxConvDistLossUnits);
        minConvDistLossUnits = Maths::min(minConvDistLossUnits, minDistLossUnits[i] + minConvDistLossUnits);
    }

    // we do a cutoff of the max lu amount of the final loss distribution
    if (cutoff)
    {
        maxConvDistLossUnits = Maths::min(maxConvDistLossUnits, dcutoff - minConvDistLossUnits);
    }

    // allocate array that will contain each input discretised distribution
    // and the array that will contain the discretised convoluted distribution
    DoubleArray dist(maxMaxDistLossUnits  - minMinDistLossUnits  + 1, 0);
    DoubleArray convDist(maxConvDistLossUnits - minConvDistLossUnits + 1, 0);

    // offset is the index in the array that corresponds to 0 value
    int distOffset      = - minMinDistLossUnits;
    int convDistOffset  = - minConvDistLossUnits;

    /**************************/
    /* convolution algorithm  */
    /**************************/

    // initialise the range of the convoluted distribution : initialy a single point is in use
    // proba(convDist = 0) = 1 
    int minConvDistIndex = convDistOffset;
    int maxConvDistIndex = convDistOffset;
    convDist[convDistOffset] = 1.;

    // loop over all distributions
    for (int i=0; i<nbDist; ++i)
    {
        // store the distribution values and probabilities into local variables
        tmpDiscreteDist = (*distributions)[i]->discretise();
        tmpValues = tmpDiscreteDist->getValues();
        tmpProbabilities = tmpDiscreteDist->getProbabilities();
        
        // allign the distribution on integer amounts of loss units
        // and store it in the array dist

        // first reset the range of the dist array that will be used by the distribution #i
        int pmax = maxDistLossUnits[i];
        int pmin = minDistLossUnits[i];
        for (int p=pmin; p<=pmax; ++p)
        {
            dist[p + distOffset] = 0.;
        }

        // We store in nonZeroProbaLossUnit the integer amounts of loss units
        // for which there is a non zero probability.
        IntArray nonZeroProbaLossUnit(maxMaxDistLossUnits - minMinDistLossUnits + 1);
        int nbNonZeroProbaLossUnit = 0;
        int lmax = tmpValues->size();

        // loop over all possible values for the distribution #i
        for (int l=0; l<lmax; ++l)
        {
            int l1, l2;
            double d, p1, p2;

            // discretisation type defines how the double value in the input distribution
            // is split onto integer amounts of loss units
            switch (discretisationType)
            {
            case 1: // discretisation is done on the closest lower integer amount of loss units.
                    // In that case the expected loss might not be conserved
                l1 = (int) floor( (*tmpValues)[l] / lu);
                dist[l1 + distOffset] += (*tmpProbabilities)[l];
                if (nbNonZeroProbaLossUnit == 0 || nonZeroProbaLossUnit[nbNonZeroProbaLossUnit] != l1)
                {
                    nonZeroProbaLossUnit[nbNonZeroProbaLossUnit] = l1;
                    ++nbNonZeroProbaLossUnit;
                }
                break;
            case 2: // discretisation is done on the closest upper integer amount of loss units.
                    // In that case the expected loss might not be conserved
                l2 = (int) ceil( (*tmpValues)[l] / lu);
                dist[l2 + distOffset] += (*tmpProbabilities)[l];                    
                if (nbNonZeroProbaLossUnit == 0 || nonZeroProbaLossUnit[nbNonZeroProbaLossUnit] != l2)
                {
                    nonZeroProbaLossUnit[nbNonZeroProbaLossUnit] = l2;
                    ++nbNonZeroProbaLossUnit;
                }
                break;
            case 3: // discretisation is done on the 2 closest integer amount of loss units 
                    // such that the expected loss is conserved
                // in that case, we solve the system
                // p1 * l1 *lu + p2 * l2 *lu = p * l
                // p1          + p2          = p
                // where the unknows are p1 and p2 and p and l are respectively
                // the input proba and value #l of the input distribution #i
                l1 = (int) floor( (*tmpValues)[l] / lu);
                l2 = (int) ceil(  (*tmpValues)[l] / lu);
                d = l1 * lu - l2 * lu;
                if (Maths::isZero(d))
                {
                    dist[l1 + distOffset] += (*tmpProbabilities)[l];
                    if (nbNonZeroProbaLossUnit == 0 || nonZeroProbaLossUnit[nbNonZeroProbaLossUnit - 1] != l1)
                    {
                        nonZeroProbaLossUnit[nbNonZeroProbaLossUnit] = l1;
                        ++nbNonZeroProbaLossUnit;
                    }
                }
                else
                {
                    p1 = (*tmpProbabilities)[l] * ( (*tmpValues)[l] - l2 * lu ) / d;
                    p2 = (*tmpProbabilities)[l] * ( l1 * lu - (*tmpValues)[l] ) / d;
                    dist[l1 + distOffset] += p1;
                    dist[l2 + distOffset] += p2;
                    if (nbNonZeroProbaLossUnit == 0 || (nbNonZeroProbaLossUnit >= 0 && nonZeroProbaLossUnit[nbNonZeroProbaLossUnit] != l1 && nbNonZeroProbaLossUnit-1 >=0 && nonZeroProbaLossUnit[nbNonZeroProbaLossUnit-1] !=l1))
                    {
                        nonZeroProbaLossUnit[nbNonZeroProbaLossUnit] = l1;
                        ++nbNonZeroProbaLossUnit;
                    }
                    if (nbNonZeroProbaLossUnit == 0 || (nbNonZeroProbaLossUnit >= 0 && nonZeroProbaLossUnit[nbNonZeroProbaLossUnit] != l2 && nbNonZeroProbaLossUnit-1 >=0 && nonZeroProbaLossUnit[nbNonZeroProbaLossUnit-1] !=l2))
                    {
                        nonZeroProbaLossUnit[nbNonZeroProbaLossUnit] = l2;
                        ++nbNonZeroProbaLossUnit;
                    }
                }
                break;
            default:
                break;
                // throw error
                ;
            }
        }
           
        //-----------------
        // In the following:
        //------------------
        // We loop over the relevant points k to update in convDist.
        // The range to update grows with i.
        
        // The range of loss units created by the distribution #i is
        // as computed above minDistLossUnits[i] .. maxDistLossUnits[i]
        // hence the loop on l is done on that range.

        // The convolution algorithm aminly involves the following line
        // convDist[k] = convDist[k - l] * dist[l + distOffset];
        // If for a given distribution all l are >= 0 then 
        // for a given k we only need the elements convDist[k'] with k'<=k.
        // This is why below in that case we loop backward over k.
        // If for a given distribution all l are <= 0 then 
        // for a given k we only need the elements convDist[k'] with k'>=k
        // This is why below in that case we loop upward over k.
        int maxDistLossUnits_i = maxDistLossUnits[i];
        int minDistLossUnits_i = minDistLossUnits[i];
        if (minDistLossUnits_i >= 0) // all the values of distribution #i are >= 0 in that case
        {
            // 2nd version 
            // in this second version, we have stored the values of the distribution for which
            // there is a non 0 probability. We loop over these values only instead of looping 
            // as above over all values including the values that have 0 probability
            int max_k;
            if (cutoff)
            {
                max_k = Maths::min(maxConvDistIndex + maxDistLossUnits_i, dcutoff - minConvDistLossUnits);
            }
            else
            {
                max_k = maxConvDistIndex + maxDistLossUnits_i;
            }

            for (int k = max_k; k >= minConvDistIndex + minDistLossUnits_i ;--k )
            {
                // treat the initial value outside the loop because l can be 0
                // in that case convDist[k] and convDist[k - l] can be the same element
                int l = nonZeroProbaLossUnit[0];
                if (k - l >= 0 && !Maths::isZero(convDist[k - l]))
                {
                    convDist[k] = convDist[k - l] * dist[l + distOffset];
                }
                
                // loop over the remaining values
                for (int p = 1; p < nbNonZeroProbaLossUnit && k - nonZeroProbaLossUnit[p] >= minConvDistIndex; ++p)
                {
                    l = nonZeroProbaLossUnit[p];
                    if (!Maths::isZero(convDist[k - l]))
                    {
                        convDist[k] += convDist[k - l] * dist[l + distOffset];
                    }
                }
                for (int k = minConvDistIndex + minDistLossUnits_i -1 ; k>=0 ;--k )
                {
                    convDist[k] = 0;
                }

            }
        }
        else if (maxDistLossUnits_i <= 0) // all the values of distribution #i are <= 0 in that case
        {
            // 2nd version 
            // in this second version, we have stored the values of the distribution for which
            // there is a non 0 probability. We loop over these values only instead of looping 
            // as above over all values including the values that have 0 probability
            int max_k;
            if (cutoff)
            {
                max_k = Maths::min(maxConvDistIndex + maxDistLossUnits_i, dcutoff - minConvDistLossUnits);
            }
            else
            {
                max_k = maxConvDistIndex + maxDistLossUnits_i;
            }

            for (int k = minConvDistIndex + minDistLossUnits_i; k <= max_k;++k )
            {
                // treat the last value outside the loop because l can be 0
                // in that case convDist[k] and convDist[k - l] can be the same element
                if (nbNonZeroProbaLossUnit > 0)
                {
                    int l = nonZeroProbaLossUnit[nbNonZeroProbaLossUnit - 1];
                    if (k - l >= 0 && !Maths::isZero(convDist[k - l]))
                    {
                        convDist[k] = convDist[k - l] * dist[l + distOffset];
                    }
                }                    
                // loop over the remaining values
                for (int p = nbNonZeroProbaLossUnit - 2; p >= 0 && k - nonZeroProbaLossUnit[p] <= maxConvDistIndex; --p)
                {
                    int l = nonZeroProbaLossUnit[p];
                    if (!Maths::isZero(convDist[k - l]))
                    {
                        convDist[k] += convDist[k - l] * dist[l + distOffset];
                    }
                }
                for (int k = maxConvDistIndex + maxDistLossUnits_i + 1 ; k<=maxConvDistIndex ;++k )
                {
                    convDist[k] = 0;
                }
            }
        }
        else // the distribution #i has values <0 and >0
        {
            // 2nd version 
            // in this second version, we have stored the values of the distribution for which
            // there is a non 0 probability. We loop over these values only instead of looping 
            // as above over all values including the values that have 0 probability

            // else we probably need to copy the convDist array
            DoubleArray convDistCopy(maxConvDistLossUnits - minConvDistLossUnits + 1, 0);

            // TO DO : replace by a more efficient copy
            for (int k=0; k<= convDist.size(); ++k)
            {
                convDistCopy[k] = convDist[k];
            }

            int max_k;
            if (cutoff)
            {
                max_k = Maths::max(maxConvDistIndex + maxDistLossUnits_i, dcutoff - minConvDistLossUnits);
            }
            else
            {
                max_k = maxConvDistIndex + maxDistLossUnits_i;
            }

            for (int k = minConvDistIndex + minDistLossUnits_i; k <= max_k ;++k )
            {
                // loop over all values
                for (int p = nbNonZeroProbaLossUnit; p >= 0 && k - nonZeroProbaLossUnit[p] <= maxConvDistIndex; --p)
                {
                    int l = nonZeroProbaLossUnit[p];
                    if (!Maths::isZero(convDistCopy[k - l]))
                    {
                        convDist[k] += convDistCopy[k - l] * dist[l + distOffset];
                    }
                }
                for (int k = minConvDistIndex + minDistLossUnits_i -1 ; k>=0 ;--k )
                {
                    convDist[k] = 0;
                }
                for (int k = maxConvDistIndex + maxDistLossUnits_i + 1 ; k<=maxConvDistIndex ;++k )
                {
                    convDist[k] = 0;
                }
            }
        }

#ifdef  PRINT_DISTRIB_ITER
            if (debug==0)
            {
            
            FILE *fp3 = fopen("C:\\temp\\lossdistribution_new_iter.out","a");
//            fprintf(fp3,"\n");
//            fprintf(fp3,"%d\t%d\t%.25lf\t%.25lf\n", l1, l2, (1-pm) * weight1 , (1-pm) * (1-weight1));
//            fprintf(fp3,"\n");
            for (int k =0; k<convDist.size(); ++k)
            {
                fprintf(fp3,"%.25lf\t", convDist[k]);
            }
            fprintf(fp3,"\n");
            fclose(fp3);
            }
#endif


        // update the range of the relevant range of the convoluted distribution
        maxConvDistIndex += maxDistLossUnits_i;
        if (cutoff)
        {
            maxConvDistIndex = Maths::min(maxConvDistIndex, maxConvDistLossUnits);
        }

        minConvDistIndex += minDistLossUnits_i;
    }
    /********************************/
    /* end of convolution algorithm */
    /********************************/

    // build the discreteDistribution to return
    double sumProbas = 0.;
    DoubleArraySP values;
    DoubleArraySP probabilities;
    values.reset(new DoubleArray());
    probabilities.reset(new DoubleArray());
    for (int l = minConvDistLossUnits; l < maxConvDistLossUnits ; ++l)
    {
        if ( !Maths::isZero( convDist[l + convDistOffset] ))
        {
            (*values).push_back( l * lu);
            (*probabilities).push_back(convDist[l + convDistOffset]);
            sumProbas += convDist[l + convDistOffset];
        }
    }
    // last loss amounts dealt separately because of the potential cutoff
    if ( !Maths::isZero( 1. - sumProbas ))
    {
        (*values).push_back( maxConvDistLossUnits * lu);
        (*probabilities).push_back(1. - sumProbas);
    }
    
    DiscreteDistributionSP convolutedDiscreteDistribution( new DiscreteDistribution(values, probabilities));
    return convolutedDiscreteDistribution;
}

/** [Implements IConvolutor] */
double RecursiveConvolutor::convoluteAndIntegrate(
    IDistribution1DArrayConstSP distributions,
    ICreditLossConfigConstSP lossConfig,
    const DateTime& timepoint) const
{
    const ITranchesCombinationPayoff* combinationPayoff =
        dynamic_cast<const ITranchesCombinationPayoff*>(lossConfig.get());
    if (combinationPayoff != 0)
    {
        // initialise "payoff" for given timepoint
        DoubleArray baseStrikes;
        DoubleArray baseStrikesWeights;
        double expectedLossOffset;
        combinationPayoff->linearDecomposition(
            timepoint,
            baseStrikes,
            baseStrikesWeights,
            expectedLossOffset);

        int nbStrikes = baseStrikes.size();
        double maxStrike = baseStrikes[nbStrikes-1];

        double cutoff = maxStrike;
        // convolution (optimised to account for maxStrike)
        DiscreteDistributionConstSP lossDistribution =
            convolute(distributions, &cutoff)->discretise();
    
        DoubleArrayConstSP values = lossDistribution->getValues();
        DoubleArrayConstSP probas = lossDistribution->getProbabilities();

        // integration
        double result = expectedLossOffset;
        for (int i = 0; i < values->size(); ++i) {
            double lossLevel = (*values)[i];
            double weight = 0.0;
            for (int k = 0; k < nbStrikes; ++k) {
				if (lossLevel < baseStrikes[k])
                {
                    weight += lossLevel * baseStrikesWeights[k];
                }
                else
                {
                    weight += baseStrikes[k] * baseStrikesWeights[k];
                }
			}
            result += weight * (*probas)[i];
        }
        return result;
    }
    else if (CDOPortfolio::TYPE->isInstance(lossConfig))
    {
        // Don't need a convolution here, just sum the single-name expected losses
        double ptfExpectedLoss = 0.0;
        for (int i = 0; i < distributions->size(); ++i) {
			ptfExpectedLoss += (*distributions)[i]->expectation();
		}
        return ptfExpectedLoss;
    }
    else if ( dynamic_cast<const FlatCDO2LossConfig*>(lossConfig.get()) )
		//(FunctionLossConfig::TYPE->isInstance(lossConfig))
    {
		const FlatCDO2LossConfig* functionLossConfig = 
			dynamic_cast<const FlatCDO2LossConfig*>(lossConfig.get());
			
		DateTimeArrayConstSP timeLine = functionLossConfig->getTimeline();

		int tIdx = functionLossConfig->findProfileIdx(timepoint);

		double tResult = 0.0;
    
		DiscreteDistributionConstSP lossDistribution =
			convolute(distributions)->discretise();
	
		DoubleArrayConstSP values = lossDistribution->getValues();
		DoubleArrayConstSP probas = lossDistribution->getProbabilities();

		// integration
		for (int i = 0; i < values->size(); ++i) 
		{
			tResult += (*functionLossConfig->getLossProfiles())[tIdx]->map( (*values)[i] ) * (*probas)[i];
		}

		return tResult;
	} else
    {
        throw ModelException (
        "RecursiveConvolutor::convoluteAndIntegrate",
        "Payoff with type " +
        lossConfig->getClass()->getName() +
        " is not supported.");
    }
}

// Computes the loss unit
double RecursiveConvolutor::lossUnit(IDistribution1DArrayConstSP distributions) const
{
    // local variables to store the distribution values and probabilities
    DiscreteDistributionConstSP tmpDiscreteDist;
    DoubleArrayConstSP tmpValues;
    DoubleArrayConstSP tmpProbabilities;
    
    int nbDist = (*distributions).size();
    DoubleArraySP amounts(new DoubleArray());
    for (int i = 0; i < nbDist; ++i)
    {
        // discretise...
        tmpDiscreteDist = (*distributions)[i]->discretise();

        // store the distribution values and probabilities into local variables
        tmpValues = tmpDiscreteDist->getValues();
        tmpProbabilities = tmpDiscreteDist->getProbabilities();
        
        for (int l = 0; l < tmpValues->size(); ++l)
        {
            if ( !Maths::isZero( (*tmpValues)[l] ))
            {
                (*amounts).push_back( (*tmpValues)[l]);
            }
        }
    }
    return CCMLossUnit::calcLossUnit(
            (*amounts),
            maxNbSlice,
            nbSubdivision);
}

/* external symbol to allow class to be forced to be linked in */
bool RecursiveConvolutorLoad(){
    return (RecursiveConvolutor::TYPE != 0);
}

DRLIB_END_NAMESPACE

