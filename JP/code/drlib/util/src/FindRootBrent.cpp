//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : FindRootBrent.hpp - from CMLib FindRoot.h
//
//   Description : Finds roots of functions using Brent methodology
//
//   Author      : CMLib + ALIB
//     The functions brentMethod and secantMethod are taken from the Alib.
//
//   Date        : Dec 17, 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FindRootBrent.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

const double ONE_PERCENT = 0.01;

class FindRootBrent::Finder {
    enum RootSolverStatus{
        ROOT_NOT_FOUND,
        ROOT_FOUND,
        ROOT_BRACKETED
    };

    struct Bound {
        double value;
        bool tested; // attempted to compute f(value)
        bool testPassed; // value is within domain

        Bound() {
            tested = false;
        }
        void setTestResult(bool testResult) {
            tested = true;
            testPassed = testResult;
        }
        bool withinDomain() {
            return tested && testPassed;
        }
    };


    Function&         m_function;

    double            m_xAccuracy;
    double            m_fAccuracy;

    int               m_maxNumOfIterations;

    double            m_xGuess;

    Bound             m_lowerBound;
    Bound             m_upperBound;

    int               m_lastPoint;
    Point             m_points[3];          

    RootSolverStatus  m_status;


    void pushPoint(double x, double y ) {
        ASSERT(m_lastPoint<2);
        m_points[++m_lastPoint] = Point(x,y);
    }
    void popPoint() {
        ASSERT(m_lastPoint>=0);
        --m_lastPoint;
    }
    Point& lastPoint() {
        return m_points[m_lastPoint];
    }

    // returns false if root is not on the boundary 
    // and boundaries do not bracket root
    bool tryBounds(){
        // try lower bound
        if (tryBound(m_lowerBound) ) {
            if (m_status != ROOT_NOT_FOUND ) {
                return true;
            } else {
                // lower bound does not bracket root
                popPoint();
            }
        }
        // try upper bound
        if (tryBound(m_upperBound) ) {
            if (m_status != ROOT_NOT_FOUND ) {
                return true;
            } else {
                // lower bound does not bracket root
                popPoint();
            }
        }
        return false;
    }

    // returns false if boundary is out of the domain or
    // already has been tested
    bool tryBound(Bound& b)
    {
        if (!b.tested) {
            bool result = tryPoint(b.value);
            b.setTestResult(result);
            return result;

        }
        return false;

    }

    // function tests function at point x
    // if x is out of domain, it adjusts bounds
    // if x is within domain it checks wether x is root or brackets root
    // returns false if x is in domain
    bool tryPoint(double x) {
        ASSERT(m_status == ROOT_NOT_FOUND);
        if (x >= m_lowerBound.value && x<=m_upperBound.value) {
            bool withinDomain;
            double y = m_function(x, withinDomain);
            if (withinDomain) {

                // check wether root is bracketed
                if (Maths::isZero(y) || 
                    (fabs(y) <= m_fAccuracy  && 
                      fabs(x-lastPoint().x) <= m_xAccuracy) )
                {
                    m_status = ROOT_FOUND;
                } else if (y>0 != m_points[0].y>0 ) {
                    m_status=ROOT_BRACKETED;
                }
                pushPoint(x,y);
                
                return true;
            } else {
                // x is not in domain - adjust the bounds
                Bound& failedBound = 
                    (x > m_xGuess ) ? m_upperBound : m_lowerBound;

                // x is between m_xGuess and failedBound.value
                if (failedBound.withinDomain() )
                    throw ModelException("FindRootBrent::Finder::tryPoint",
                                         "functions has fragmented domain" );

                failedBound.value = x;
                failedBound.setTestResult(false);
            }
        }
        // x is not in domain - try bounds
        return false;

    }

    /// return a String showing list of points used so far
    string printPoints() const {
#ifdef CM_STORE_SOLVER_POINTS
        size_t n = m_function.m_point.size();
        if (n==0) return string();
        ostringstream os;
        os.precision(16);
        os << "Function to Solve Evaluated " << n << " times\n";
        for (size_t i=0; i<n; ++i)
            os <<  "x=" << m_function.m_point[i].first
               <<" ,y=" << m_function.m_point[i].second << '\n';
        return os.str();
#else
        return string();
#endif
    }
public:
    Finder(Function& function ): m_function(function) {}

    Point operator ()(double guess,
                      double lowerBound,
                      double upperBound,
                      double xAccuracy,
                      double fAccuracy,
                      size_t maxNumOfIterations,
                      double initialXStepSize )
    {
        static const string method("FindRootBrent::Finder()");
#ifdef CM_STORE_SOLVER_POINTS
        m_function.m_point.reserve(maxNumOfIterations); 
#endif
        m_xGuess = guess;

        // set bounds and 
        // check wether guess is equal to one of the bounds
        m_lowerBound.value = lowerBound;
        if (guess == lowerBound)
            m_lowerBound.setTestResult(true);

        m_upperBound.value = upperBound;
        if (guess == upperBound)
            m_upperBound.setTestResult(true);

        m_xAccuracy = xAccuracy;
        m_fAccuracy = fAccuracy;
        m_maxNumOfIterations = maxNumOfIterations;
        m_status = ROOT_NOT_FOUND;

        /* Make sure lower bound is below higher bound.
         */
        if (lowerBound >= upperBound) {
            throw ModelException(method, printPoints() +
                                 Format::toString("lower bound %ld greater or "
                                                  "equal to higher bound %ld",
                                                  lowerBound, upperBound));
        }

        /* Check to make sure the guess is within the bounds 
         */
        if (guess<lowerBound || guess>upperBound) {
            throw ModelException(
                method, printPoints() +
                Format::toString("guess ld% is out of range [%ld, %ld]",
                                 guess, lowerBound, upperBound));
        }

        bool    withinDomain;
        double  yGuess = m_function(guess, withinDomain);
        if (! withinDomain)
            throw ModelException(method, printPoints() + 
                                 "initial guess has to be within domain" );

        /* Check if guess is the root (use bounds for x-accuracy) */
        if (Maths::isZero(yGuess) || 
            (fabs(yGuess) <= fAccuracy &&
             (fabs(lowerBound-guess) <= xAccuracy ||
              fabs(upperBound-guess) <= xAccuracy)))
        {
            return Point(guess, yGuess); 
        }
        m_lastPoint = -1;
        pushPoint(guess, yGuess);


        /* If the initialXStep is 0, set it to ONE_PERCENT of 
         *  of (upperBound - lowerBound). 
         */
        double smallStep  = ONE_PERCENT * (upperBound - lowerBound); 

        if (initialXStepSize == 0.0 )  {
            initialXStepSize = smallStep; 
        }

        if (!tryPoint(guess + initialXStepSize) && 
             !tryPoint(guess - initialXStepSize) && 
             !tryBound(m_lowerBound) &&
             !tryBound(m_upperBound) &&
             !tryPoint(lowerBound + smallStep) && 
             !tryPoint(upperBound - smallStep) ) 
        {
            throw ModelException(method, printPoints() +"domain is too small");
        }

        switch(m_status ) {
        case ROOT_NOT_FOUND:
            /* Call secant method to find the root, or to get a 
               third point, so that two of the three points bracket the root. */
            if (secantMethod() && m_status == ROOT_BRACKETED ) 
            {
                // make sure m_points[1].x is between m_points[0].x
                // and m_points[2].x since m_points[0].x and
                // m_points[1] have the same sign and
                // fabs(m_points[0].y) < fabs(m_point[1].y), we can
                // assume:
                ASSERT((m_points[1].x<m_points[0].x) == 
                       (m_points[0].x<m_points[2].x));
                std::swap(m_points[1], m_points[0]);
                break;
            }
            if (m_status == ROOT_FOUND )
                break;
            // if bracket was found on the domain boundary 
            // take only guess and bracket and add midpoint
            m_points[1] = m_points[2];
            popPoint();

        case ROOT_BRACKETED:
            m_status = ROOT_NOT_FOUND;
            if (!tryPoint((m_points[0].x + m_points[1].x)/2 ) ) 
               throw ModelException(method, printPoints() +
                                    "functions has fragmented domain" );

            // make sure m_point[0] and m_point[2] bracket root
            // and m_points[1].x is between m_points[0].x and m_points[2].x 
            std::swap(m_points[1], m_points[2]); 
            m_status = ROOT_BRACKETED;
        }

        if (m_status == ROOT_FOUND) {
            return lastPoint();
        } else {
            ASSERT(m_status == ROOT_BRACKETED);
            return brentMethod();
        }
    }

    // function returns true if solution/bracket was found in normal iteration
    // false if it was found on the boundary
    bool secantMethod() {
        static const string method("FindRootBrent::Finder::secantMethod");
        ASSERT(m_status == ROOT_NOT_FOUND );
        int i = 0;
        while (++i<=m_maxNumOfIterations )
        {
            ASSERT(m_lastPoint==1);
            /* Swap points so that yPoints[0] is smaller than yPoints[1] */
            if (fabs(m_points[0].y)>fabs(m_points[1].y)) {   
                std::swap(m_points[0],m_points[1]);
            }

            /* Make sure that you do not divide by a very small value */
            double dy = m_points[1].y - m_points[0].y;

            if (fabs(dy) <= m_fAccuracy ) {
                //dy = CopySign(m_fAccuracy, dy );
                dy = dy < 0.0? -fabs(m_fAccuracy): fabs(m_fAccuracy);
            }

            double dx= (m_points[0].x-m_points[1].x)* m_points[0].y / dy ;
            double x = m_points[0].x + dx;

            /* original code was doing something like this
               ...
               if (!tryPoint(x ) ) {
               if (tryBounds() )
               return;
               CM_THROW RuntimeError();
               }
               ...
            */
            if (!tryPoint(x ) ) {
                // secant method failed by trying to sample outside domain
                // check wether boundaries bracket root 
                if (tryBounds() ) 
                    return false;

                // if corresponding boundary is known and it does not
                // bracket root we fail
                Bound& bound = (x < m_xGuess) ? m_lowerBound : m_upperBound;
                if (bound.withinDomain() )
                    throw ModelException(
                        method, printPoints() +
                        Format::toString("function values at bounds %ld, %ld"
                                         " imply no root exists",
                                         m_lowerBound.value,
                                         m_upperBound.value));
                // we don't know the boundary
                // try to find it using bisection
                x = bound.value;
                double lastValid = m_points[1].x;
                do {
                    if (fabs(lastValid - x) <= m_xAccuracy) {
                        throw ModelException(method, printPoints() +
                                             "failed to bracket root with "
                                             "given accuracy");
                    }
                    
                    if (++i<=m_maxNumOfIterations )
                        throw ModelException(
                            printPoints() + 
                            string("Brent root finder aborted after " ) +
                            Format::toString(m_maxNumOfIterations) + 
                            " iterations" );

                    x = (x + lastValid)/2;
                } while(!tryPoint(x) );

                // better approximation of the boundary found -
                // continue secant method
            }

            if (m_status != ROOT_NOT_FOUND)
                return true;

            m_points[1] = m_points[2];
            popPoint();
        }

        // secant method could not bracket the root
        // within m_maxNumOfIterations iterations
        // check wether boundaries bracket root 
        if (tryBounds() ) 
            return false;

        throw ModelException(method, printPoints() +
                             Format::toString("function values at bounds "
                                              "%ld, %ld imply no root exists",
                                              m_lowerBound.value,
                                              m_upperBound.value));
    }

    /* ----------------------------------------------------------------------
       Function    :   brentMethod
       DESCRIPTION :   Finds the root using a combination of inverse quadratic
       method and bisection method.
       ---------------------------------------------------------------------- */
    Point brentMethod() {
        static const string method("FindRootBrent::Finder::brentMethod");
        double      ratio;              /* (x3-x1)/(x2-x1) */
        double      x31;                /* x3-x1*/
        double      x21;                /* x2-x1*/
        double      f21;                /* f2-f1 */
        double      f31;                /* f3-f1 */
        double      f32;                /* f3-f2 */
        double      xm;                 /* New point found using Brent method*/
        double      fm;                 /* f(xm) */

        double      x1 = m_points[0].x; /* xN short hand for m_points[n].x */
        double      f1 = m_points[0].y;

        double      x2 = m_points[1].x;
        double      f2 = m_points[1].y;

        double      x3 = m_points[2].x;
        double      f3 = m_points[2].y;

        bool        withinDomain;
        ASSERT(f1<0 != f3<0 && x1<x2 == x2<x3);
        for (int i = 0; i<m_maxNumOfIterations; ++i) 
        {
            /* Always want to be sure that f1 and f2 have opposite signs,
             * and f2 and f3 have the same sign.
             */
            if (f2*f1>0.0)
            {   
                std::swap(x1,x3);
                std::swap(f1,f3);
            }
            f21 = f2-f1;
            f32 = f3-f2;
            f31 = f3-f1;
            x21 = x2-x1;
            x31 = x3-x1;
            /* Check whether this is suitable for interpolation. When checking
             * for f21,f31,f32 = 0, Alib didn't use IS_ALMOST_ZERO for "efficiency
             * reasons". If the objective Function has been truncated to a 
             * certain number of digits of accuracy (were in one case) 
             * or it is constant, f21,f31,or f32 could be
             * zero. In this case we need to protect against
             * division by zero. So we use bisection instead of brent.
             */
            ratio = (x3-x1)/(x2-x1);
            if (f3*f31<ratio*f2*f21 || Maths::isZero(f21) ||
                Maths::isZero(f31) || Maths::isZero(f32)) 
            {
                /* Do bisection 
                 */
                x3 = x2;
                f3 = f2; 

            } else {
                xm = x1 - (f1/f21)*x21 + ((f1*f2)/(f31*f32))*x31 - 
                    ((f1*f2)/(f21*f32))*x21;
                fm = m_function(xm, withinDomain);
                if (! withinDomain) {
                    throw ModelException(method, printPoints() + 
                                         "function has fragmented domain" );
                }
                if (fm == 0.0 || (fabs(fm) <= m_fAccuracy && 
                                  fabs(xm-x1) <= m_xAccuracy)) {
                    return Point(xm, fm);
                }
                /* If the new point and point1 bracket the root,
                   replace point3 with the new point */
                if (fm*f1<0.0) {
                    x3=xm;
                    f3=fm;
                }
                /* If the root is not bracketed, replace point1 with new point,
                   and point3 with point2 */
                else {
                    x1=xm;
                    f1=fm;
                    x3=x2;
                    f3=f2;
                }
            }             
            x2 = 0.5*(x1+x3); 
            f2=m_function(x2, withinDomain);
            if (!withinDomain)
                throw ModelException(method, printPoints() +
                                     "functions has fragmented domain" );
            if (f2 == 0.0 || (fabs(f2) <= m_fAccuracy &&
                              fabs(x2 - x1) <= m_xAccuracy)) {
                return Point(x2, f2);
            }

        }

        throw ModelException(method, printPoints() +
                             string("Brent root finder aborted after " ) +
                             Format::toString(m_maxNumOfIterations) +
                             " iterations" );

    }    
};

FindRoot::Point FindRootBrent::solveForFunction(Function& function,   
                                                double guess, 
                                                double lowerBound, 
                                                double upperBound, 
                                                double xAccuracy, 
                                                double fAccuracy, 
                                                int    maxNumOfIterations, 
                                                double initialXStepSize) {
    
    Finder solver(function);
    return solver(guess,
                  lowerBound,
                  upperBound,
                  xAccuracy,
                  fAccuracy,
                  maxNumOfIterations,
                  initialXStepSize);
}

DRLIB_END_NAMESPACE

