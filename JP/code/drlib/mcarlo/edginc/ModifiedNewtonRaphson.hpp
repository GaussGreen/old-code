#ifndef EDR_SRMROOTFINDER_HPP
#define EDR_SRMROOTFINDER_HPP

#include "edginc/SRMConstants.hpp"

DRLIB_BEGIN_NAMESPACE

template<class F>
class NewtonRootFinder
{
    /* The old algorithm did an extra step after computing the difference. 
       We want to allow this behaviour in order to prevent lots of regression test diffs (sorry about that).
       When things 'calm down', it would be ok to remove the option m_bExtraStepMode (so it's always false).
     */
       
    bool m_bExtraStepMode;
public:

    NewtonRootFinder(bool bExtraStepMode = false) : m_bExtraStepMode(bExtraStepMode) {}

    //Returns 'true' if it succeeded:
    //Use this slightly weird construct (rather than throwing an exception on error, for example), 
    //because of backwards compatibility...

    //Note that exceptions thrown by this method are caught and suppressed when called by the 
    //rates calibrator: in that case we just use the previous value...
    void findRoot(double bracket_low, double bracket_high, double start, double tolerance, F &f, double &ret, bool suppress_check = false)
    {
        const double SIGDELTA = 0.0001 ;
        ret = start;
        int iter = 0;

        if (!suppress_check)
        {
            double valLeft = f(bracket_low);
            QLIB_VERIFY(valLeft <= 0, "Left hand evaluation of function must be negative.");
        } 
        do {
           
            double val = f(ret);

            if (!m_bExtraStepMode && fabs(val) < tolerance) break;

            double offset = 0.0;
            double new_val = 0.0;
            double newValOk = true;

            double derivative = f.approxDerivative(ret, val, SIGDELTA);

            if (fabs(derivative ) < SRMConstants::SRM_TINY)
            {
                newValOk = false;
            }
            else
            {
                /* Newton Raphson */
                offset = val /derivative;
                new_val = ret - offset;

                //only want to do this rarely: backwards compatibility
                if (new_val <= bracket_low || new_val >= bracket_high)
                {
                    newValOk = false;
                }
            } 

            //only want to do this rarely: backwards compatibility
            if (!newValOk)
            {
                double valLeft = f(bracket_low);
                QLIB_VERIFY(valLeft <= 0, "Left hand evaluation of function must be negative.");

                double valRight = f(bracket_high);
                int c = 0;
                while(valRight < 0)
                {
                    bracket_high *= 2;
                    valRight = f(bracket_high);
                    c++;
                    QLIB_VERIFY(c<=10, "'bracket_high' has become too large. Maybe 'f' has no root.");
                }

                ret = ZBrent_solve(f, bracket_low, bracket_high, tolerance, ZBrent::default_ITMAX);
            }
            else
            {
                ret = new_val;
            }
            iter++;

            if (m_bExtraStepMode && fabs(val) < tolerance) break;

        } while (iter < SRMConstants::MAX2QITER);   /* Newton-Raphson condition */    

        QLIB_VERIFY(iter < SRMConstants::MAX2QITER, "Newton-Raphson has failed because too many iterations were needed.");
    }
    virtual ~NewtonRootFinder(void) {}
};

DRLIB_END_NAMESPACE

#endif


