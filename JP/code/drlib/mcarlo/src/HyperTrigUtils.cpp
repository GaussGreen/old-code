#include "edginc/config.hpp"
#include "edginc/HyperTrigUtils.hpp"
#include <assert.h>


DRLIB_BEGIN_NAMESPACE



//TODO: granularity should depend on C
HyperLocalVolState::HyperLocalVolState(double a, double b, double c, int tableSize) : m_a(a), m_b(b), m_c(c), m_inverseTableSize(tableSize)
{
    QLIB_VERIFY(c >= 0, "c < 0 (" + Format::toString(c) + ")");
    const double inverseInc = 0.01;
    for (double y = -m_inverseTableSize;y<=m_inverseTableSize;y+=inverseInc)
    {
        double x = 0.0;
        double s = -0.999, sdash = -0.999;
        evaluateKExp(y, &x);
        evaluateLocalVol(y, &s, &sdash);
        addInterpolant(x, y, s, sdash);
    }
}

//x is the image of y under KfuncExp
void HyperLocalVolState::addInterpolant(double x, double y, double sm, double smDash)
{
    QuadraticPoly pol;
    pol.m_deg0 = y;
    pol.m_deg1 = sm;
    pol.m_deg2 = 0.5 * smDash * sm;

    FunctionInterpolant interpolant;
    interpolant.m_centre = x;
    interpolant.m_taylor_series = pol;

    m_quadraticInterpolants.push_back(interpolant);        
}
   

void HyperLocalVolState::evaluateLocalVol(double u, double *s, double *sdash) const 
{
    /* smile =  (1 + A * Tanh(C * x) +
    *                   B * ( 1 - 1/Cosh(C * x)) ) */
    if (m_a == 0 && m_b == 0)
    {
        *s = 1.0;
        *sdash = 0.0;
    }
    else 
    {  
        double w = m_c*u;
        const double SMILE_CUTOFF = 100.0;

        if (w < SMILE_CUTOFF && w > - SMILE_CUTOFF) {
            double y = exp(w);
            double z = 1.0/y;
            double t = (y - z) / (y + z);
            double se = 2.0 / (y + z);

            *s =  1.0 + m_a * t
                + m_b * (1.0 - se);
            *sdash = m_c * (m_a * (1 - t * t) + m_b * t * se);
        }
        else if (w <= - SMILE_CUTOFF) {
            /* y = 0 */
            *s = (1.0 - m_a + m_b);
            *sdash = 0.0;
        } else {
            /* z = 0 */
            *s = (1.0 + m_a + m_b);    
            *sdash = 0.0;
        }
    }
    //Should not be too small: should really be close to 1
    *s = Maths::max(*s, 0.01);
}   

/**** Adapted from: Hyb3_Kfunc ***********************************************/
/**                                                         
Calculates the general smile mapping function K(x)    
ie, K(x) where K'(x) = 1/g(x)                         
g(x) is the local volatility mapping function         
Note: if diffusing the log moneyness, rather than just moneyness,
then pass in  exp(log_moneyness)
*********************************************************/

void HyperLocalVolState::evaluateK(
                             double  M,    /**<(I) S/F, moneyness         */
                             double  logM, /**<(I) log moneyness (avoid computing log if logM already known) */
                             double  *val) const /**<(O)                        */
{
    double a1 = m_a;
    double a2 = m_b;
    double a3 = m_c;

    double x, z, alpha, beta, gamma, delta, kappa_up, kappa_low;
    double a, b, c; /*rational fraction decomposition*/

    /*special lognormal case*/
    if ( (Maths::isZero(a1) && Maths::isZero(a2)) || Maths::isZero(a3) ) {
        (*val) = logM;
        return;
    }
    alpha  =  1. + a1 + a2;
    beta   =  -2. * a2;
    gamma  =  1. - a1 + a2; 
    delta  =  beta * beta - 4. * alpha * gamma;

    x   =   logM; /* log moneyness */
    z   =   exp(a3 * x);

    /*Special cases alpha = 0, gamma = 0*/
    /*alpha == 0 case*/
    if(Maths::isZero(alpha))
    {
        /*Two subcases, beta = 0 or gamma = 0*/
        if(Maths::isZero(beta))
        {   
            *val = z*z/2.0 - 0.5 + log(z);
            *val /= gamma;
            *val /= a3;
            return;

        }else if(Maths::isZero(gamma))
        {
            *val = z - 1/z;
            *val /= beta;
            *val /= a3;
            return;
        }else
        {
            *val = (z - 1)/beta;
            *val -= gamma*(log(beta * z + gamma) - log(beta + gamma))/((beta*beta));
            *val += (log(z) - log(beta*z + gamma) + log(beta + gamma))/gamma;
            *val /= a3;
            return;
        }
    }

    /*gamma == 0 case*/
    if(Maths::isZero(gamma))
    {
        if (Maths::isZero(beta))
        {
            *val = log(z) - 1.0/(2*z*z);
            *val += 0.5;
            *val /= (a3 * alpha);
            return;
        }
        else
        {
            a    = -alpha/(beta*beta);
            b    = 1/beta;
            c    = (alpha*alpha)/(beta*beta);
            *val = (log(alpha*z + beta)- log(alpha + beta))/alpha;
            *val += a*log(z) - b/z + b + c*log(alpha*z + beta)/alpha - c*log(alpha + beta)/alpha;
            *val /= a3;
            return;
        }
    }


    /* kappa  =  \int_1^{y^c} 1/(alpha * z^2 + beta * z+ gamma) dz */
    kappa_up  =  (delta > 0)?
        log(  (2. * alpha * z + beta - sqrt(delta) ) / 
        (2. * alpha * z + beta + sqrt(delta) ) ) / sqrt(delta) :
    2*atan( (2. * alpha * z + beta) / sqrt(-delta) ) / sqrt(-delta);
    kappa_low =  (delta > 0)?
        log(  (2. * alpha + beta - sqrt(delta) ) / 
        (2. * alpha + beta + sqrt(delta) ) ) / sqrt(delta) :
    2*atan( (2. * alpha + beta) / sqrt(-delta) ) / sqrt(-delta);

    /* deal with the case delta == 0 (we overwrite kappa_up and kappa_low) */
    if(Maths::isZero(delta))
    {
        kappa_up  = -2/(2.0 * alpha * z + beta);
        kappa_low = -2/(2.0 * alpha     + beta);
    }


    if ( a3*x < 50.0 ) 
    {
        *val      =  a3*x / (a3 * gamma) + 
            0.5 * (1. / alpha - 1. / gamma) * ( log(alpha*z*z+beta*z+gamma) - log(alpha+beta+gamma) ) / a3 -
            0.5 *  beta * (1. / alpha + 1. /gamma )  * ( kappa_up - kappa_low ) / a3; 
    }
    else /* treat this case seperately for speed and stability purpose (AlexK) */
    {
        /* approximate log(alpha*z*z) by log(alpha) + 2.0*(a3*x) + neglect */
        /* log( exp( a3*x) ) = a3*x : more efficient                       */
        *val      =  a3*x / (a3 * gamma) + 
            0.5 * (1. / alpha - 1. / gamma) * ( log(alpha) + 2.0*a3*x - log(alpha+beta+gamma) ) / a3 -
            0.5 *  beta * (1. / alpha + 1. /gamma )  * ( kappa_up - kappa_low ) / a3; 

    }
}

void HyperLocalVolState::evaluateKExp(double Z, double  *val) const
{
    evaluateK(exp(Z), Z, val);
}
  


/** Use interval bisection to locate interpolant with centre closest to Kt */
int HyperLocalVolState::findClosestInterpolantIdx(double Kt) const
{
    if (m_quadraticInterpolants.empty())
    {
        throw ModelException("HyperTrigUtils::HyperLocalVolState::findClosestInterpolantIdx",
            "No interpolants!");
    }
    size_t s = m_quadraticInterpolants.size();
    int min_loc = 0;
    if (Kt >= m_quadraticInterpolants[s-1].m_centre)
        min_loc = m_quadraticInterpolants.size() - 1;
    else if (Kt <= m_quadraticInterpolants[0].m_centre)
        min_loc = 0;
    else
    {
        int low = 0;
        int high = m_quadraticInterpolants.size() - 1;

        while(high - low > 1)
        {
            int midpoint = (low + high) / 2;
            double comp = m_quadraticInterpolants[midpoint].m_centre;
            if (comp <= Kt)
                low = midpoint;
            else
                high = midpoint;
        }

        double dlow = m_quadraticInterpolants[low].m_centre;
        double dhigh = m_quadraticInterpolants[high].m_centre;
        if (Kt - dlow < dhigh - Kt)
            min_loc = low;
        else
            min_loc = high;
    }
    return min_loc;
}
   
/** Locate closest interpolant centre, then plug in the difference to the quadratic */
void HyperLocalVolState::evaluateInverseKExp(double Kt, double *val) const
{
    int min_location = findClosestInterpolantIdx(Kt);    
    double key = m_quadraticInterpolants[min_location].m_centre;
    double diff = Kt - key;
    *val = m_quadraticInterpolants[min_location].m_taylor_series.evaluate(diff);
}  


DRLIB_END_NAMESPACE

