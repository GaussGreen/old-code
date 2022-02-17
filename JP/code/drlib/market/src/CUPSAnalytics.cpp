//----------------------------------------------------------------------------
//
//   Group       : DR
//
//   Filename    : CUPSAnalytics.cpp
//
//   Description : Code from CMLib library initially used for quanto adjustment
//                 to CDS
//
//   Author      : CMLib (From file with samename)
//
//   Date        : 10 Dec 2004
//
//

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CUPSAnalytics.hpp"

DRLIB_BEGIN_NAMESPACE
BaseCDSFunction::~BaseCDSFunction(){}

CDSFunctionState::CDSFunctionState() {
    I1 = I2 = 0;
    DF = 1.;
    index = 0;
}

ICovarianceFunctions::~ICovarianceFunctions(){}

UnadjustedCDSFunction::UnadjustedCDSFunction(const DoubleArray&    tOffset,
                                             const DoubleArray&    fwdRlDF,
                                             const DoubleArray&    rRate,
                                             double recoveryRate):
    BaseCDSFunction(recoveryRate){
    m_tOffset = tOffset;
    m_fwdRlDF = fwdRlDF;
    m_rRate = rRate;
}

void UnadjustedCDSFunction::updateState(double lambda) {
    double dI1 = 0;
    for (int n = m_state.index; n<m_benchmarkMaturity; n++)
    {
        double r_plus_l = m_rRate[n] + lambda;
        double deltaT = m_tOffset[n+1] - m_tOffset[n];
        if (Maths::isZero(r_plus_l)){
            dI1 += m_state.DF * deltaT;
        } else {
            double fwdDF = m_fwdRlDF[n] * exp(- lambda * deltaT); 
            dI1 += m_state.DF * (1 - fwdDF) / r_plus_l;
            m_state.DF *= fwdDF;
        }
        
    }
    m_state.I1 += dI1;
    m_state.I2 += lambda * dI1;
    m_state.index = m_benchmarkMaturity;
}

void CUPSAnalyticsBase::initialize(const DoubleArray&    tOffset,
                                   const DoubleArray&    fwdRlDF,
                                   const DoubleArray&    rRate,
                                   const DoubleArray&    rBpVol,
                                   const DoubleArray&    lBpVol,
                                   double spreadIRCorr,
                                   double recoveryRate) {
    m_tOffset = tOffset;
    m_fwdRlDF = fwdRlDF;
    m_rRate = rRate;
    m_rBpVol = rBpVol;
    m_lBpVol = lBpVol;
    m_spreadIRCorr = spreadIRCorr;
    m_recoveryRate = recoveryRate;
    finalize();
}

class CUPSAnalyticsR: public CUPSAnalyticsBase  {
protected:
    double SpreadIrCov(int maturity) {
        // initialize with CRK[t_n]
        double result = 0;
        double betaLrec = 1/m_betaL;
        int N = maturity + 1;
        double tn = m_tOffset[N];
        double T2 = m_expBetaL[N] * betaLrec;
        for (int k = 0; k<=maturity; k++)
        {
            int K= k+1;
            double tk = m_tOffset[K];
            double tk1 = m_tOffset[K-1];
            result += m_rBpVol[k] * m_lBpVol[k] * (
                (tk-tk1) * (tn - (tk+tk1)/2) - 
                T2 * (1 / m_expBetaL[K]   * (betaLrec + tn - tk) - 
                       1 / m_expBetaL[K-1] * (betaLrec + tn - tk1)));
        }
        return result * m_spreadIRCorr * betaLrec ;

    }

    double SpreadFxCov(const DoubleArray&    fxVol,
                       double                spreadFXCorr,
                       int                   maturity) {
        // initialize with CRK[t_n]
        double result = 0;
        int N = maturity + 1;
        double T3 = m_expBetaL[N] / m_betaL;
        for (int k = 0; k<=maturity; k++) {
            int K= k+1;
            result += fxVol[k] * m_lBpVol[k] * (
                m_tOffset[K] - m_tOffset[K-1] 
                - T3 *(1/m_expBetaL[K] - 1/m_expBetaL[K-1]));

        }
        return result * spreadFXCorr /m_betaL;
    }

    void finalize() {
        m_expBetaL.resize(m_tOffset.size());
        int n;
        for (n = 0; n<m_tOffset.size(); n++) {
            m_expBetaL[n] = exp(- m_betaL * m_tOffset[n]);
        }

        int numDates = m_rRate.size();

        m_T1N_coef.resize(numDates);
        m_T11N_coef.resize(numDates);
        m_T4N_coef.resize(numDates);
        double betaLrec = 1/m_betaL;
        for (n = 0; n<numDates; n++) {
            int N = n+1;
            double deltaTn = m_tOffset[N] - m_tOffset[N-1];

            // set initial coefficients to J3NN
            double T = m_rBpVol[n] * m_lBpVol[n] * betaLrec;
            m_T11N_coef[n] = T;
            double T1 = T * (deltaTn - betaLrec);
            double T4 = T / m_betaL;
            double T40 = m_expBetaL[N-1] * betaLrec;
            // sum J3NK coefficients
            for (int k = 0; k < n; k++) {
                double T = m_rBpVol[k] * m_lBpVol[k] * betaLrec;
                int K = k+1; // to be used with base-0 arrays
                double deltaTk = m_tOffset[K] - m_tOffset[K-1];

                T1 += T * deltaTk;;
                T4 += - T * T40 * (1/m_expBetaL[K] - 1/m_expBetaL[K-1]);
            }
            m_T1N_coef[n] = T1;
            m_T4N_coef[n] = T4;
        }
    }

    void updateState(double lambda) {
        for (int n = m_state.index; n<m_benchmarkMaturity; n++) {
            double r_plus_l = m_rRate[n] + lambda;
            double deltaT = m_tOffset[n+1] - m_tOffset[n];
            double expBL = m_expBetaL[n+1]/m_expBetaL[n]; // = exp(-betaL*dt)


            double fwdDF;
            double T1N, T11N;
            if (Maths::isZero(r_plus_l)) {
                fwdDF= 1.;
                T1N = deltaT;
                T11N = - deltaT*deltaT/2;
            } else {
                fwdDF= m_fwdRlDF[n] * exp(- lambda * deltaT); 
                T1N = (1 - fwdDF) / r_plus_l;
                T11N = (T1N - deltaT) / r_plus_l;
            }
            double T4N = (1 - fwdDF * expBL)/(r_plus_l+m_betaL);

            m_state.I1 += m_state.DF * T1N;
            m_state.I2 += m_state.DF * 
                (lambda * T1N + m_spreadIRCorr * 
                  (T1N*m_T1N_coef[n] + 
                    T11N*m_T11N_coef[n] + T4N*m_T4N_coef[n]));

            m_state.DF *= fwdDF;
        };
        m_state.index = m_benchmarkMaturity;
    }

    double m_betaL;
    // cached exp(- betaL *(tk-t0))
    // this array is 0-based
    DoubleArray    m_expBetaL;
    DoubleArray    m_T1N_coef;
    DoubleArray    m_T11N_coef;
    DoubleArray    m_T4N_coef;
public:
    CUPSAnalyticsR(
        double betaL) 
    {
        m_betaL = betaL;
    }
};

// CDS pricing that accounts for correlation between spread and interest rate
// betaR != 0 && betaL != 0
// inherits SpreadFxCov from CUPSAnalyticsR

class CUPSAnalytics: public CUPSAnalyticsR {
    
    double SpreadIrCov(int maturity) {
        // initialize with CRK[t_n]
        double result = 0;
        int N = maturity + 1;
        double T1 = m_expBetaR[N]*m_expBetaL[N] / (m_betaR + m_betaL);
        double T2 = m_expBetaR[N] / m_betaR;
        double T3 = m_expBetaL[N] / m_betaL;
        for (int k = 0; k<=maturity; k++)
        {
            int K= k+1;
            result += m_rBpVol[k] * m_lBpVol[k] * (
                m_tOffset[K] - m_tOffset[K-1] 
                + T1 *(1/m_expBetaR[K]/m_expBetaL[K] - 
                       1/m_expBetaR[K-1]/m_expBetaL[K-1])
                - T2 *(1/m_expBetaR[K] - 1/m_expBetaR[K-1])
                - T3 *(1/m_expBetaL[K] - 1/m_expBetaL[K-1]));

        }
        return result * m_spreadIRCorr /m_betaR /m_betaL;
    }

    void finalize() {
        m_expBetaR.resize(m_tOffset.size());
        m_expBetaL.resize(m_tOffset.size());
        int n;
        for (n = 0; n<m_tOffset.size(); n++) {
            m_expBetaR[n] = exp(- m_betaR * m_tOffset[n]);
            m_expBetaL[n] = exp(- m_betaL * m_tOffset[n]);
        }

        int numDates = m_rRate.size();

        m_T1N_coef.resize(numDates);
        m_T2N_coef.resize(numDates);
        m_T3N_coef.resize(numDates);
        for (n = 0; n<numDates; n++) {
            int N = n+1;
            // set initial coefficients to J3NN
            double vol = m_rBpVol[n] * m_lBpVol[n];
            m_T1N_coef[n] = vol / (m_betaR*(m_betaR+m_betaL));
    
            double T2 = - vol / (m_betaR*m_betaL);
            double T3 = vol / ((m_betaR+m_betaL)*m_betaL);
            // sum J3NK coefficients
            for (int k = 0; k < n; k++) {
                double T = m_rBpVol[k] * m_lBpVol[k] / m_betaL;
                int K = k+1; // to be used with base-0 arrays

                T2 += T /m_betaR * m_expBetaR[N-1] *
                    (1/m_expBetaR[K] - 1/m_expBetaR[K-1]);

                T3 += - T /(m_betaR+m_betaL) * 
                    m_expBetaR[N-1] * m_expBetaL[N-1] * 
                    (1/(m_expBetaR[K]*m_expBetaL[K]) -
                     1/(m_expBetaR[K-1]*m_expBetaL[K-1]));
            }
            m_T2N_coef[n] = T2;
            m_T3N_coef[n] = T3;
        }
    }

    void updateState(double lambda) {
        for (int n = m_state.index; n<m_benchmarkMaturity; n++) {
            double r_plus_l = m_rRate[n] + lambda;
            double deltaT = m_tOffset[n+1] - m_tOffset[n];
            double expBR = m_expBetaR[n+1]/m_expBetaR[n]; // = exp(-betaR*dt)
            double expBL = m_expBetaL[n+1]/m_expBetaL[n]; // = exp(-betaL*dt)
            

            double fwdDF;
            double T1N;
            if (Maths::isZero(r_plus_l)) {
                fwdDF= 1.;
                T1N = deltaT;
            } else {
                fwdDF= m_fwdRlDF[n] * exp(- lambda * deltaT); 
                T1N = (1 - fwdDF) / r_plus_l;
            }
            double T2N = (1 - fwdDF * expBR)/(r_plus_l+m_betaR);
            double T3N = (1 - fwdDF * expBR * expBL)/(r_plus_l+m_betaR+m_betaL);

            m_state.I1 += m_state.DF * T1N;
            m_state.I2 += m_state.DF *
                (lambda * T1N + m_spreadIRCorr * 
                  (T1N*m_T1N_coef[n] + T2N*m_T2N_coef[n] + T3N*m_T3N_coef[n]));
            
            m_state.DF *= fwdDF;
        };
        m_state.index = m_benchmarkMaturity;
    }

    double m_betaR;
    // cached exp(- betaR *(tk-t0)) 
    // these arrays are 0-based
    DoubleArray    m_expBetaR;
    DoubleArray    m_T1N_coef;
    DoubleArray    m_T2N_coef;
    DoubleArray    m_T3N_coef;
public:
    CUPSAnalytics(double betaR, double betaL) : CUPSAnalyticsR(betaL) {
        m_betaR = betaR;
    }
};
// betaR == 0 && betaL == 0
class CUPSAnalyticsRL :public CUPSAnalyticsBase {
protected:
    double SpreadIrCov(int maturity) {
        // initialize with CRK[t_n]
        double result = 0;
        int N = maturity + 1;
        double tn = m_tOffset[N];
        for (int k = 0; k<=maturity; k++)
        {
            int K= k+1;
            double tk = m_tOffset[K];
            double tk1 = m_tOffset[K-1];
            double dtk = tk-tk1;
            result += m_rBpVol[k] * m_lBpVol[k] * dtk * 
                ((tn - tk)*(tn - tk1) + dtk*dtk*dtk/3) ;
        }
        return result * m_spreadIRCorr ;
    }

    double SpreadFxCov(const DoubleArray&    fxVol,
                       double                spreadFXCorr,
                       int                   maturity) {
        // initialize with CRK[t_n]
        double result = 0;
        int N = maturity + 1;
        double tn = m_tOffset[N];
        for (int k = 0; k<=maturity; k++) {
            int K= k+1;
            double tk = m_tOffset[K];
            double tk1 = m_tOffset[K-1];
            result += fxVol[k] * m_lBpVol[k] * 
                (tk-tk1) * (tn - (tk+tk1)/2);

        }
        return result * spreadFXCorr;
    }

    void finalize() {
        int numDates = m_rRate.size();

        m_T1N_coef.resize(numDates);
        m_T11N_coef.resize(numDates);
        m_T12N_coef.resize(numDates);
        for (int n = 0; n<numDates; n++) {
            int N = n+1;
            double deltaTn = m_tOffset[N] - m_tOffset[N-1];
            // set initial coefficients to J3NN
            double vol = m_rBpVol[n] * m_lBpVol[n];
            
            m_T12N_coef[n] = vol;
    
            double T11 = vol * deltaTn;
            double T1 = vol * deltaTn * deltaTn / 2;
            // sum J3NK coefficients
            for (int k = 0; k < n; k++) {
                int K = k+1; // to be used with base-0 arrays
                double deltaTk = m_tOffset[K] - m_tOffset[K-1];
                double T = m_rBpVol[k] * m_lBpVol[k] * deltaTk;

                T11 += T;

                T1 += T * (m_tOffset[N] - (m_tOffset[K] + m_tOffset[K-1])/2);
            }
            m_T1N_coef[n] = T1;
            m_T11N_coef[n] = T11;
        }
    }

    void updateState(double lambda) {
        for (int n = m_state.index; n<m_benchmarkMaturity; n++) {
            double r_plus_l = m_rRate[n] + lambda;
            double deltaT = m_tOffset[n+1] - m_tOffset[n];
            double deltaT2 = deltaT * deltaT;


            double fwdDF;
            double T1N, T11N, T12N;
            if (Maths::isZero(r_plus_l)) {
                fwdDF= 1.;
                T1N = deltaT;
                T11N = - deltaT2 / 2;
                T12N = deltaT2 * deltaT / 6;
            } else {
                fwdDF= m_fwdRlDF[n] * exp(- lambda * deltaT); 
                T1N = (1 - fwdDF) / r_plus_l;
                T11N = (T1N - deltaT) / r_plus_l;
                T12N = (T11N + deltaT2/2) / r_plus_l;
            }

            m_state.I1 += m_state.DF * T1N;
            m_state.I2 += m_state.DF * 
                (lambda * T1N + m_spreadIRCorr * 
                  (T1N*m_T1N_coef[n] + T11N*m_T11N_coef[n] +
                   T12N*m_T12N_coef[n]));

            m_state.DF *= fwdDF;
        };
        m_state.index = m_benchmarkMaturity;
    }

    DoubleArray    m_T12N_coef;
    DoubleArray    m_T11N_coef;
    DoubleArray    m_T1N_coef;
public:
    CUPSAnalyticsRL() {}
};

// betaR != 0 && betaL == 0
// inherits SpreadFxCov from CUPSAnalyticsRL
class CUPSAnalyticsL: public CUPSAnalyticsRL {
protected:
    double SpreadIrCov(int maturity) {
        // initialize with CRK[t_n]
        double result = 0;
        double betaRrec = 1/m_betaR;
        int N = maturity + 1;
        double tn = m_tOffset[N];
        double T2 = m_expBetaR[N] * betaRrec;
        for (int k = 0; k<=maturity; k++)
        {
            int K= k+1;
            double tk = m_tOffset[K];
            double tk1 = m_tOffset[K-1];
            result += m_rBpVol[k] * m_lBpVol[k] * (
                (tk-tk1) * (tn - (tk+tk1)/2) - 
                T2 * (1 / m_expBetaR[K]   * (betaRrec + tn - tk) - 
                       1 / m_expBetaR[K-1] * (betaRrec + tn - tk1)));
        }
        return result * m_spreadIRCorr * betaRrec ;

    }

    void finalize() {
        m_expBetaR.resize(m_tOffset.size());
        int n;
        for (n = 0; n<m_tOffset.size(); n++) {
            m_expBetaR[n] = exp(- m_betaR * m_tOffset[n]);
        }

        int numDates = m_rRate.size();
        double betaRrec = 1/m_betaR;

        m_T1N_coef.resize(numDates);
        m_T2N_coef.resize(numDates);
        m_T21N_coef.resize(numDates);
        for (n = 0; n<numDates; n++) {
            int N = n+1;
            double deltaTn = m_tOffset[N] - m_tOffset[N-1];

            // set initial coefficients to J3NN
            double T = m_rBpVol[n] * m_lBpVol[n] * betaRrec;
            m_T1N_coef[n] = T * betaRrec;
            double T2 = - T * (betaRrec + deltaTn);
            double T21 = - T;
            // sum J3NK coefficients
            for (int k = 0; k < n; k++) {
                double T = m_rBpVol[k] * m_lBpVol[k] * m_expBetaR[N-1]/m_betaR;
                int K = k+1; // to be used with base-0 arrays

                T2 += T * 
                    ((betaRrec + m_tOffset[N] - m_tOffset[K])/m_expBetaR[K] - 
                      (betaRrec + m_tOffset[N] - 
                       m_tOffset[K-1])/m_expBetaR[K-1]);

                T21 += T * (1/m_expBetaR[K] - 1/m_expBetaR[K-1]);
            }
            m_T2N_coef[n] = T2;
            m_T21N_coef[n] = T21;
        }
    }

    void updateState(double lambda) {
        for (int n = m_state.index; n<m_benchmarkMaturity; n++)
        {
            double r_plus_l = m_rRate[n] + lambda;
            double deltaT = m_tOffset[n+1] - m_tOffset[n];
            double expBR = m_expBetaR[n+1]/m_expBetaR[n]; // = exp(-betaR*dt)


            double fwdDF;
            double T1N;
            if (Maths::isZero(r_plus_l)) {
                fwdDF= 1.;
                T1N = deltaT;
            } else {
                fwdDF= m_fwdRlDF[n] * exp(- lambda * deltaT); 
                T1N = (1 - fwdDF) / r_plus_l;
            }
            double T2N = (1 - fwdDF * expBR)/(r_plus_l+m_betaR);
            double T21N = (T2N - deltaT)/(r_plus_l+m_betaR);

            m_state.I1 += m_state.DF * T1N;
            m_state.I2 += m_state.DF * (lambda * T1N + m_spreadIRCorr * 
                                         (T1N*m_T1N_coef[n] + 
                                           T2N*m_T2N_coef[n] +
                                           T21N*m_T21N_coef[n]));

            m_state.DF *= fwdDF;
        };
        m_state.index = m_benchmarkMaturity;
    }

    double m_betaR;
    // cached exp(-betaR *(tk-t0))
    // these arrays are 0-based
    DoubleArray    m_expBetaR;
    DoubleArray    m_T21N_coef;
    DoubleArray    m_T2N_coef;
    DoubleArray    m_T1N_coef;
public:
    CUPSAnalyticsL(double betaR) {
        m_betaR = betaR;
    }
};


CUPSAnalyticsBase* CUPSAnalyticsBase::create(double betaR, double betaL) {
    CUPSAnalyticsBase *result;
    if (Maths::isZero(betaR)) {
        if (Maths::isZero(betaL)){
            result = new CUPSAnalyticsRL();
        } else {
            result = new CUPSAnalyticsR(betaL);
        }
    } else {
        if (Maths::isZero(betaL)) {
            result = new CUPSAnalyticsL(betaR);
        } else {
            result = new CUPSAnalytics(betaR, betaL);
        }
    }
    return result;
}

DRLIB_END_NAMESPACE
