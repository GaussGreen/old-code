//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : LossDistribution.cpp
//
//   Description : Convolution Algorithm
//
//   Date        : Feb 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Maths.hpp"
#include "edginc/LossDistribution.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE

#define REALLY_TINY 1e-64

CClassConstSP const LossDistribution::TYPE = CClass::registerClassLoadMethod(
    "LossDistribution", typeid(LossDistribution), load);

DEFINE_TEMPLATE_TYPE(LossDistributionArray);

LossDistribution::LossDistribution() :
    CObject(TYPE),
    lossDist(0),
    left(0),
    right(0),
    zero(0),
    lastRunName(0)
{}

LossDistribution::LossDistribution(const int nbLoss, const int maxShortIdx) :
    CObject(TYPE),
    left(0),
    right(0),
    zero(maxShortIdx),
    lastRunName(0)
{
    lossDist.resize(nbLoss);
    //lossDist.clear();
    lossDist[zero] = 1.0;
}

double LossDistribution::maxAtTiny(double inp)
{
    if (inp < REALLY_TINY) {
        return 0.;
    } else {
        return inp;
    }
}

// convolute a distribution
void LossDistribution::convoluteDistribution(int numPoints,
                                             const double * s, // survival probability array 
                                             const double * l) // loss array
{
    DoubleArray::iterator out = lossDist.begin()+zero;

    // very expensive copy operation, this needs to be improved
    DoubleArray workArray = lossDist;
    DoubleArray::iterator work = workArray.begin()+zero;

    int j;

#ifdef DEBUG
    // check inputs
    double totalProb = 0;
    for (j = 0; j < numPoints; ++j)
        totalProb += 1-s[j];
    ASSERT(abs(totalProb-1)<0.01);
#endif

    double bound = right;
    ASSERT(bound < lossDist.size());
    for (j = 0; j <= bound; ++j)
        out[j] = 0;

    for (j=0; j<= bound; ++j) {
        if (work[j] < 3e-16) continue;
        for (int i = 0; i < numPoints; ++i) {
            if (s[i] < 1.0 - 3e-16) {
                // the mass we want to put on out[j+l[i]]
                double delta = (1-s[i])*work[j];
                double tgt = j + l[i];
                int base = tgt;
                // split ratio between two adjacent buckets
                double p = base + 1 - tgt;
                out[base] += p* delta;
                out[base+1] += (1-p)*delta;
                if (base + 1 > right) 
                    right = base + 1;
            }
        }
    }        

#ifdef DEBUG
    // check output
    totalProb = 0;
    for (j=0; j<lossDist.size(); ++j)
        totalProb += out[j];
    ASSERT(abs(totalProb-1)<0.01);
#endif
}

void LossDistribution::convolute(const int    nmIdx,
                                 const double s,
                                 const double weight1,
                                 const int    l1,
                                 const int    l2)
{
    //const char routine[] = "LossDistribution::convolute";
    //try
    //{
    if (s < 1.0 - 3e-16)
    {
        DoubleArray::iterator out = lossDist.begin()+zero;

        double q1 = (1-s) * weight1;
        double q2 = (1-s) * (1-weight1);

        int j;
        // Code commented out is the original algorithm. It is easier to follow
        // but has a very expensive copy operation
        //         /* save the present array */
        //         workArray = outArray;
        //         DoubleArray::iterator work = workArray.begin()+maxShortIdx;
        //
        //         for (j=-left; j<=right; ++j){
        //             out[j] *= s;
        //         }
        //
        //         for (j=-left; j<=right; ++j) {   
        //             /* Pr(q1=0 or q2=0) will be zero in full CCM, so putting 
        //                an "if" won't make it any faster */
        //             out[j+l1] += q1*work[j];
        //             out[j+l2] += q2*work[j];
        //         }

        if (l2 >= 0){ // NB in which case l1 >= l2
            // run loop backwards to avoid overwriting values.
            int endOfOverlap = Maths::max(right -l2 + 1, -left);
            for (j=right; j>=endOfOverlap; --j) {   
                double d1 = out[j];
                out[j+l2] = q2*d1; // NB out[j+l2] will be zero prior to this
                out[j+l1] += q1*d1;
            }
            for (/* from where we were */; j>=-left; --j) {   
                /* Pr(q1=0 or q2=0) will be zero in full CCM, so putting 
                   an "if" won't make it any faster */
                double d1 = out[j];
                double d2 = out[j+l2]; // NB d2 is unchanged at this point
                out[j+l2] = q2*d1 + s*d2;
                out[j+l1] += q1*d1;
            }
            // now multiply the original by s (but skipping what we've done
            // already)
            int startOfOverlap = Maths::min(-left + l2-1, right);
            for (j = startOfOverlap; j>=-left; --j) {   
                out[j] *= s;
            }
        } else { // NB in which case |l1| > |l2|
            // run the loop forwards for this case.
            int endOfOverlap = Maths::min(-left -l2 -1, right);
            for (j=-left; j<=endOfOverlap; ++j) {   
                double d1 = out[j];
                out[j+l2] = q2*d1; // NB out[j+l2] will be zero prior to this
                out[j+l1] += q1*d1;
            }
            for (/* from where we were */; j<=right; ++j) {   
                /* Pr(q1=0 or q2=0) will be zero in full CCM, so putting 
                   an "if" won't make it any faster */
                double d1 = out[j];
                double d2 = out[j+l2]; // NB d2 is unchanged at this point
                out[j+l2] = q2*d1 + s*d2;
                out[j+l1] += q1*d1;
            }
            // now multiply the original by s (but skipping what we've done
            // already)
            int startOfOverlap = Maths::max(right + l2+1, -left);
            for (j = startOfOverlap; j<=right; ++j) {   
                out[j] *= s;
            }
        }

        //adjust the loss bounds
        if (l1>=0)
        {
            right += l1;
            ASSERT(l1 >= l2); /* requires extra precaution when          */
        }                     /* calibrating but avoids using MAX (slow) */
        else                  /* recoveries are such that R2>R1 and      */
        {                     /* therefore loss abs values are in the    */
            left -= l1;       /* opposite order                          */
            ASSERT(-l1 >= -l2);
        }
    }

    //record last run name if appropriate
    if (nmIdx+1 > lastRunName) lastRunName = nmIdx+1;

    //} catch (exception& e){
    //    throw ModelException(e, routine);
    //}
}

void LossDistribution::deconvolute(const int    nmIdx,
                                   const double s,
                                   const double weight1,
                                   const int    l1,
                                   const int    l2)
{
    const char routine[] = "LossDistribution::deconvolute";
    try
    {
        if (s < 1.0 - 3e-16)
        {
            DoubleArray::iterator out = lossDist.begin()+zero;

            double q1 = (1-s) * weight1;
            double q2 = (1-s) * (1-weight1);

            int j;
            // reverse the convolution steps
            // - need to regain "previous state" to deconvolute properly
            // - can do in 5 steps

            // For a positive notional
            if (l1 > 0)
            {
                //for a reasonably large probability of survival
                //(lower values tend to result in numerical instabilities)
                if (s > q1)
                {
                    // - run forwards across the distribution
                    //
                    // ----|------------------0--------|---------------
                    //     l   2      1                r   2      1
                    //     |   |      |                |   |      |
                    //      /s   -q2       -q1          -q1   -q1
                    //           /s        -q2          -q2
                    //                     /s
                    //     L                                      R

                    for (j=-left; j<-left+l2; ++j)
                    {
                        //out[j] /= s;
                        out[j] = maxAtTiny(out[j] / s);
                    }

                    for (j=-left+l2; j<-left+l1; ++j)
                    {
                        //out[j] = (out[j]-(out[j-l2]*q2))/s;
                        out[j] = maxAtTiny((out[j]-(out[j-l2]*q2))/s);
                    }

                    for (j=-left+l1; j<=right-l1; ++j)
                    {
                        //out[j] = (out[j]-((out[j-l2]*q2)+(out[j-l1]*q1)))/s;
                        out[j] = maxAtTiny((out[j]-((out[j-l2]*q2)+(out[j-l1]*q1)))/s);
                    }

                    //the distribution from right-l1+1 to right
                    //should now go to 0

                    for (j=right-l1+1; j<=right; ++j)
                    {
                        out[j] = 0.;
                    }
                }
                else
                {
                    // - run backwards across the distribution
                    // requires a copy
                    DoubleArray workArray = lossDist;
                    DoubleArray::iterator work = workArray.begin()+zero;

                    for (j=right; j>right-l1 && j>=-left; j--)
                    {
                        out[j] = 0.0;
                    }

                    for (j=right-l1; j>right-l1-(l1-l2) && j>=-left; j--)
                    {
                        //out[j] = work[j+l1]/q1;
                        out[j] = maxAtTiny(work[j+l1]/q1);
                    }

                    for (j=right-l1-(l1-l2); j>right-l1-l1  && j>=-left; j--)
                    {
                        //out[j] = (work[j+l1]-(out[j+l2]*q2))/q1;
                        out[j] = maxAtTiny((work[j+l1]-(out[j+l2]*q2))/q1);
                    }

                    for (j=right-l1-l1; j>=-left; j--)
                    {
                        //out[j] = (work[j+l1] - ((out[j+l1-l2]*q2) + (out[j+l1]*s)))/q1;
                        out[j] = maxAtTiny((work[j+l1] - ((out[j+l1-l2]*q2) + (out[j+l1]*s)))/q1);
                    }
                }
            }
            else //(l1 < 0) short name
            {
                //for a reasonably large probability of survival
                //(lower values tend to result in numerical instabilities)
                if (s > q1)
                {
                    // - run backwards across the distribution
                    //
                    // -------------|--------------0-------------|---
                    //  1       2   l                 1      2   r
                    //  |       |   |                 |      |   |
                    //     -q1   -q1       -q1          -q2   /s
                    //           -q2       -q2          /s
                    //                     /s
                    //  L                                        R

                    //NB l1 and l2 are negative in this case
                    for (j=right; j>right+l2; --j)
                    {
                        //out[j] /= s;
                        out[j] = maxAtTiny(out[j]/s);
                    }

                    for (j=right+l2; j>right+l1; --j)
                    {
                        //out[j] = (out[j]-(out[j-l2]*q2))/s;
                        out[j] = maxAtTiny((out[j]-(out[j-l2]*q2))/s);
                    }

                    for (j=right+l1; j>=-left-l1; --j)
                    {
                        //out[j] = (out[j]-((out[j-l2]*q2)+(out[j-l1]*q1)))/s;
                        out[j] = maxAtTiny((out[j]-((out[j-l2]*q2)+(out[j-l1]*q1)))/s);
                    }

                    //the distribution from -left-l1-1 to -left
                    //should now go to 0
                    for (j=-left-l1-1; j>=-left; --j)
                    {
                        out[j] = 0.;
                    }
                }
                else
                {
                    // - run forwards across the distribution
                    // requires a copy
                    DoubleArray workArray = lossDist;
                    DoubleArray::iterator work = workArray.begin()+zero;

                    //for convenience
                    int L1 = -l1;
                    int L2 = -l2;

                    for (j=-left; j<-left+L1; j++)
                    {
                        out[j] = 0.;
                    }

                    for (j=-left+L1; j<-left+L1+L2; j++)
                    {
                        //out[j] = work[j-l1]/q1;
                        out[j] = maxAtTiny(work[j-L1]/q1);
                    }

                    for (j=-left+L1+L2; j<-left+L1+L1; j++)
                    {
                        //out[j] = (work[j-L1] - (out[j-(L1-L2)]*q2))/q1;
                        out[j] = maxAtTiny((work[j-L1] - (out[j-(L1-L2)]*q2))/q1);
                    }


                    for (j=-left+L1+L1; j<=right; j++)
                    {
                        //out[j] = (work[j-L1] - ((out[j-(L1-L2)]*q2) + (out[j-L1]*s)))/q1;
                        out[j] = maxAtTiny((work[j-L1] - ((out[j-(L1-L2)]*q2) + (out[j-L1]*s)))/q1);
                    }
                }
            }

            //reduce left & right bounds
            if (l1 >=0)
            {
                right -= l1;
            }
            else
            {
                left += l1;
            }

        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

int LossDistribution::getLastRunName()
{
    return lastRunName;
}

DoubleArray& LossDistribution::getLossDistribution()
{
    return lossDist;
}

void LossDistribution::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(LossDistribution, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultLossDistribution);

    //transient fields
    FIELD( lossDist,    "the distribution");
    FIELD( left,        "how far we have progressed along the distribution");
    FIELD( right,       "how far we have progressed along the distribution");
    FIELD( zero,        "the position of the 0 loss element cf maxShoftIdx");
    FIELD( lastRunName, "which names contribute to this distribution");

    //FIELD_MAKE_TRANSIENT(lossDist);
    //FIELD_MAKE_TRANSIENT(left);
    //FIELD_MAKE_TRANSIENT(right);
    //FIELD_MAKE_TRANSIENT(zero);
    //FIELD_MAKE_TRANSIENT(lastRunName);
}

IObject* LossDistribution::defaultLossDistribution()
{
    return new LossDistribution();
}


DRLIB_END_NAMESPACE
