#include <General/General.h>
#include <Magnet/Magnet.h>
extern "C" {
#include "../student.h"
#include "../gaussian.h"
#include "../dependence.h"
#include "../independence.h"
#include "../proba_utils.h"
#include "../payoff_st.h"
}

using namespace CM;
MAGNET_PREFIX("DR_")
MAGNET_CATEGORY("Copula Utils tp")

// --------------------------------------------------------------------------
// StudentCum
//
MAGNET_DESCRIPTION("Returns the cumulative probability of a Student variable")
MAGNET_PARAMETER("",x,"","")
MAGNET_PARAMETER("",freedomDegree,"","degree of freedom")
MAGNET_FUNCTION2(   double, StudentCum,
                    double x,
                    long freedomDegree)
{
    return StudentCum(  x,
                        freedomDegree);
}

// --------------------------------------------------------------------------
// Chi2Density
//
MAGNET_DESCRIPTION("returns the chi2 distribution")
MAGNET_PARAMETER("",x,"","")
MAGNET_PARAMETER("",freedomDegree,"","degree of freedom")
MAGNET_FUNCTION2(   double, Chi2Density,
                    double x,
                    long freedomDegree)
{
    return Chi2Density( x,
                        freedomDegree);
}