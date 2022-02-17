#include "edginc/config.hpp"
#include "edginc/Addin.hpp"
#include "edginc/IRadarRep.hpp"
#include "edginc/RadarRepUtil.hpp"

DRLIB_BEGIN_NAMESPACE

double RadarRep::getValue(const FittingArray& x ) const {
    TransformedArray transformed =  (*m_transform) (x);
    BasisValArray v = (*m_f)(transformed);
    double result = std::inner_product(m_coef.begin(), m_coef.end(), v.begin(), 0.0);
    return result;
}

RegressionCoeff RadarRep::getCoef(void) const {
    return m_coef;
}

DRLIB_END_NAMESPACE
