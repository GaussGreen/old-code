//----------------------------------------------------------------------------
//
//   Filename    : SRMRatesLiborDiffuse.cpp
//
//   Description : forward libor model path generation
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/QMCRatesDiffuse.hpp"
#include "edginc/SRMRatesLiborDiffuse.hpp"
#include "edginc/SRMRatesUtil.hpp"
#include "edginc/SRMFXDiffuse.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/SRMSwap.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/QMCHelperCachingTimeLogic.hpp"
#include "edginc/IQMCRNGManager.hpp"

#include <cassert>
#include <algorithm>
#include <fstream> // for logging to file

DRLIB_BEGIN_NAMESPACE

// make sure all elementary types are initialized to some unreal values
/** constructor */
SRMRatesLiborDiffuse::SRMRatesLiborDiffuse() :
    QMCRatesDiffuse(),
    liborUtilSP(NULL),
    qLeft(-999),
    qRight(-999),
    m_corrLower(-100),
    m_corrUpper(-100),
    zeroQ(true),
    pivotRatio(-999),
    m_lastSimTime(70.), // in years
    m_liborInterval(6), // in months
    m_predCorrStep(5.), // five year step size as default
    diffYCIdx(-1),
    discYCIdx(-1),
    invSqrtTwelfth( sqrt(1./12.) )
{
}

///////////////////////////////////////////////
// end constructor
///////////////////////////////////////////////

/** destructor */
SRMRatesLiborDiffuse::~SRMRatesLiborDiffuse()
{
}

///////////////////////////////////////////////

/** set the diffusion model */
void SRMRatesLiborDiffuse::setSRMRatesLiborDiffuse(
    int                    _randomIndex,
    const DateTime&        _today,
    SRMRatesLiborUtilSP    _srmRatesLiborUtil,
    bool                   _saveDomLnMoney,
    bool                   _saveSigmaR,
    double                 _NbSigmasMax,
    double                 _NbSigmasMin,
    const vector<double>&  _df,          // for historic dates
    const vector<double>&  _irFxCorr,    // TODO : not used at the moment
    bool                   _calibrateAgainstSwaptionVols,
    const string&          _calibStyle,
    const string&          _calibMaturity,
    const string&          _calibMaturityCMS)
{
    static const string method("SRMRatesLiborDiffuse::setSRMRatesLiborDiffuse");

    randomIndex = _randomIndex;
    saveDomLnMoney = _saveDomLnMoney;
    saveSigmaR = _saveSigmaR;
    NbSigmasMax = _NbSigmasMax;
    NbSigmasMin = _NbSigmasMin;
    today = _today;
    calibrateAgainstSwaptionVols = _calibrateAgainstSwaptionVols;
    liborUtilSP = _srmRatesLiborUtil;

    origSVol=liborUtilSP->getSpotVols();

    //VolProcessedBSIRSP tempVol = liborUtilSP->getProcessedVol();
    if ( !VolProcessedBSIR::TYPE->isAssignableFrom(liborUtilSP->getProcessedVol()->getClass()))
    {
		throw ModelException(method, "Invalid cast (SRMRatesLiborDiffuse::setSRMRatesLiborDiffuse)");

    }
    VolProcessedBSIRSP tempVol = VolProcessedBSIRSP::dynamicCast(liborUtilSP->getProcessedVol());

    DateTimeArray swaptionExpiries;
    IntArray swaptionTenors;
    vector< vector<double> > marketVols;

    MaturityPeriodSP swapFrequency(tempVol->getSwapFrequency());
    DayCountConventionSP swapDCC(tempVol->getSwapDCC());

    int accrInMonths = swapFrequency->toMonths();

    m_nbFactors = liborUtilSP->numFactors();

    m_humpParams.resize(m_nbFactors);
    m_aInit.resize(m_nbFactors); 
    m_bInit.resize(m_nbFactors);
    m_dInit.resize(m_nbFactors); 
    m_a.resize(m_nbFactors); 
    m_b.resize(m_nbFactors);
    m_d.resize(m_nbFactors);

    CholMatrix_11.resize(m_nbFactors);  
    CholMatrix_21.resize(m_nbFactors); 
    CholMatrix_31.resize(m_nbFactors);  
    CholMatrix_22.resize(m_nbFactors);  
    CholMatrix_32.resize(m_nbFactors);  
    CholMatrix_33.resize(m_nbFactors); 

    /** Get benchmark details */
    tempVol->getSwaptionGrid(swaptionExpiries,  // in years, offset from today
                             swaptionTenors,    // in months, offset from swap start
                             marketVols);

    int nbExpiries = swaptionExpiries.size();
    int nbTenors = swaptionTenors.size();
    int gridSize = nbExpiries * nbTenors;
    m_gridSwap.resize(gridSize);

    int nbSwap = 0;
    m_marketExpiries.resize(nbExpiries);
    m_marketTenors.resize(nbTenors);
    m_marketSwaptionVols.resize(nbExpiries);
	m_marketExpiriesDate.resize(nbExpiries);

    const IYieldCurve *zc( (liborUtilSP->getDiffYC()).get() );

    for(int i = 0; i < nbExpiries; i++)
    {
        m_marketExpiries[i] = today.yearFrac(swaptionExpiries[i]);
		m_marketExpiriesDate[i] = swaptionExpiries[i];
        (m_marketSwaptionVols[i]).resize(nbTenors);
        for(int j = 0; j < nbTenors; j++)
        {
            if (i == 0)
            {
                m_marketTenors[j] = swaptionExpiries[i].yearFrac(
                     swaptionExpiries[i].rollDateInMonths(swaptionTenors[j]));
            }

            // modify
            m_marketSwaptionVols[i][j] = marketVols[i][j];

            (m_gridSwap[nbSwap]).reset( new SRMSwapClass(
                                             swaptionExpiries[i], // (I) Date instrument begins at
                                             true,  // stub
                                             *swapDCC,
                                             accrInMonths,
                                             swaptionTenors[j]) );

            (m_gridSwap[nbSwap])->setWeights( *zc );

            nbSwap++;
        }
    }

    MaturityPeriod calibMat(_calibMaturity);
    MaturityPeriod calibMatCMS(_calibMaturityCMS);

    m_calibMat = calibMat.toYears();
    m_calibMatCMS = calibMatCMS.toYears();

    setLiborModel(
        m_nbFactors,
        predictorCorrector);  

    setCalibDetails(
        _calibStyle,
        m_calibMat,
        m_calibMatCMS);

    double calibAccuracy = 0.001;
    
    analyticCalib(calibAccuracy);

    return;
}

//****************************************************************************//
//  end : setSRMRatesLiborDiffuse                                             //
//****************************************************************************//

/** setParameters */
// this function checked 1 August
void SRMRatesLiborDiffuse::setParameters(
                                    const vector<double> &a, 
                                    const vector<double> &b,
                                    const vector<double> &d,
                                    const double &lambda,
                                    const vector<double> &h)
{
    static const string method = "SRMRatesLiborDiffuse::setParameters";

    int libor = -1,
        factor = -1;

    double temp = 0., temp2 = 0. , tau = 0.,
           term = -1.,
           time2Reset = -1.,
           fudgeA = -1., fudgeD = -1.,
           time2  = -1., time1  = -1.,
           exp2   = -1., exp1   = -1.,
           i0     = -1., i1     = -1., i2  = -1.,
           u1     = -1., u2     = -1., u3  = -1.,
           c11    = -1., c12    = -1., c13 = -1.,
           c22    = -1., c23    = -1., c33 = -1.,
           sqrtLiborExpiry = -1., len = -1., capletVol = -1.;

	QLIB_VERIFY(m_nbLibors > 0,
		"failed because m_nbLibors < 1");

    QLIB_VERIFY(m_nbFactors > 0,
                "failed because m_nbFactors < 1");

	vector<double> norm(m_nbLibors);
    vector<double> aRotate(m_nbFactors);
    vector<double> bRotate(m_nbFactors);
    vector<double> dRotate(m_nbFactors);

    if ( (a.size()            != m_nbFactors) ||
         (b.size()            != m_nbFactors) ||
         (d.size()            != m_nbFactors) ||
         (m_humpParams.size() != m_nbFactors) )
    {
        throw ModelException(method, "Not enough parameters (SRMRatesLiborDiffuse::setParameters)");
    }

    for (factor = 0; factor < m_nbFactors; factor++)
    {
        m_a[factor]          = a[factor];
        m_b[factor]          = b[factor];
        m_d[factor]          = d[factor];
        m_humpParams[factor] = h[factor];
    }

    // normalize to remain consistent with caplet calibration
    for (libor = 0; libor < m_nbLibors; libor++)
    {
        norm[libor] = 0.0;

        // set exponential matrix and rotations
        tau = m_resetDates[libor];
        // the max ensures effect only active after h[0] years
        temp = Maths::max(0., (tau - h[0]) * h[1]);

        temp2 = temp;
        // x/(1+x^2) has maximum value 0.5 attained at x = 1
        temp2 *= 2. / (1. + temp * temp);

        double sinTheta = temp2;
        double cosTheta = sqrt(1. - temp2 * temp2);

        temp = m_lambda * m_resetDates[libor];
        m_expReset[libor] = exp(-temp);

        // note : rotates only first and second factor
        if (m_nbFactors > 1)
        {
            aRotate[0] = a[0] * cosTheta - a[1] * sinTheta;
            dRotate[0] = d[0] * cosTheta - d[1] * sinTheta;
            aRotate[1] = a[0] * sinTheta + a[1] * cosTheta;
            dRotate[1] = d[0] * sinTheta + d[1] * cosTheta;

            // d rotation only
            aRotate[0] = a[0] + d[0] - dRotate[0];
            aRotate[1] = a[1] + d[1] - dRotate[1];

            // c rotation only
            // aRotate[0] += dRotate[0] - d[0];
            // aRotate[1] += dRotate[1] - d[1];
        }
        else
        {
            aRotate[0] = a[0];
            dRotate[0] = d[0];
        }

        for(factor = 2; factor < m_nbFactors; factor++)
        {
            aRotate[factor] = a[factor];
            dRotate[factor] = d[factor];
        }

        for (factor=0; factor<m_nbFactors; factor++)
        {
            time2 = m_resetDates[libor];
            time1 = 0.0;

            if ( m_lambda > TOO_SMALL )
            {
                exp2 = exp(m_lambda * time2);
                exp1 = exp(m_lambda * time1);

                i0 = (exp2 - exp1)/m_lambda;

                i1 = (time2 * exp2 - i0)/m_lambda;

                c11 = time2;
                c12 = i0;
                c13 = i1;

                exp2 = exp(2.0 * m_lambda * time2);
                exp1 = exp(2.0 * m_lambda * time1);

                i0 = (exp2 - exp1) / (2.0*m_lambda);
                i1 = (time2 * exp2 - i0) / (2.0*m_lambda);
                i2 = (time2 * time2 * exp2 - 2.0 * i1) / (2.0*m_lambda);

                c22 = i0;
                c23 = i1;
                c33 = i2;
            }
            else
            {
                c11 = time2;
                c22 = time2;
                c33 = (time2 * time2 * time2)/3.0;
                c12 = time2;
                c13 = 0.5 * time2 * time2;
                c23 = 0.5 * time2 * time2;
            }

            fudgeA = (factor == 0) ? 1. + m_fudgeA[libor] : 1.;

            u1 = dRotate[factor];
            u2 = (aRotate[factor] * fudgeA + b[factor] * time2) *
                  m_expReset[libor];
            u3 = - b[factor] * m_expReset[libor];

            norm[libor] += (u1 * c11 * u1 + u2 * c22 * u2 + u3 * c33 * u3 +
                            2.0 * c12 * u1 * u2+
                            2.0 * c13 * u1 * u3 +
                            2.0 * c23 * u2 * u3);
        } // for factor

        time2 = m_resetDates[libor];
        sqrtLiborExpiry = sqrt(time2);

        if ( sqrtLiborExpiry > TOO_SMALL )
        {
            len = sqrt(norm[libor]) / sqrtLiborExpiry;
            capletVol = m_liborVols[libor];

            for (factor = 0; factor < m_nbFactors; factor++)
            {
                fudgeA = (factor == 0) ? 1. + m_fudgeA[libor] : 1.;

                m_A[libor][factor] = capletVol
                    * fudgeA
                    * aRotate[factor] / len;
                m_B[libor][factor] = capletVol
                    * b[factor] / len;
                m_D[libor][factor] = capletVol
                    * dRotate[factor] / len;
            }
        }

    }// for libor

    return;
}

//****************************************************************************//
//  end : setParameters                                                       //
//****************************************************************************//

/** getLiborVols */
void SRMRatesLiborDiffuse::getLiborVols(
                                vector< vector<double> > &instLiborVols,
                                double time,
                                long firstLibor,
                                long lastLibor) const
{
    double timeToExpiry = -1., expApprox = -1., tmp = -1., tmpVol = -1.;

    int libor = -1,
        factor = -1;

    const int nbLibors = instLiborVols.size(); 

    QLIB_VERIFY( ( firstLibor < 0 ) && ( firstLibor >= m_nbLibors ),
        "Libor index out of range in SRMRatesLiborDiffuse::getLiborVols");

    QLIB_VERIFY( (lastLibor < 0 ) && ( lastLibor >= m_nbLibors ),
        "Libor index out of range in SRMRatesLiborDiffuse::getLiborVols");

    const int nbFactors = instLiborVols[0].size();

    QLIB_VERIFY( (nbFactors != m_nbFactors),
        "invalid factor number in SRMRatesLiborDiffuse::getLiborVols");

    tmp = m_lambda * time;
    expApprox = exp(tmp);

    for (libor = firstLibor; libor <= lastLibor; libor++)
    {
        timeToExpiry = m_resetDates[libor] - time;

        if (timeToExpiry > 0.)
        {
            for(factor = 0; factor < m_nbFactors; factor++)
            {
                tmpVol = (m_A[libor][factor] + m_B[libor][factor] * timeToExpiry) *
                    m_expReset[libor] * expApprox + m_D[libor][factor];
                instLiborVols[libor][factor] = tmpVol;
            }
        }
        else // beyond reset date so volatility = 0
        {
            for(factor = 0; factor < m_nbFactors; factor++)
            {
                instLiborVols[libor][factor] = 0.;
            }
        }
    }


    return;
}
//****************************************************************************//
//  end : getLiborVols                                                        //
//****************************************************************************//

/** getLiborSkew */
const vector<double> & SRMRatesLiborDiffuse::getLiborSkew(void) const
{
    return (m_liborSkew);
}


//****************************************************************************//
//  end : getLiborSkew                                                        //
//****************************************************************************//

/** "expands" skew from size m_nbParams to m_nbLibors */
void SRMRatesLiborDiffuse::setLiborSkew(vector<double> &skew)
{
    static const string method = "SRMRatesLiborDiffuse::setLiborSkew";

    int nbParams = m_paramDates.size();

    QLIB_VERIFY(skew.size() == nbParams,
       "Error : skew.size() != nbParams in setLiborSkew"); 

    long i, j = 0;

    for(i = 0; i < m_nbLibors; i++)
    {
        // first on or after
        while((j < nbParams - 1) &&
              (m_paramDates[j+1] <= m_resetDates[i]))
        {
           j++;
        }

        m_liborSkew[i] = Maths::max(1. - skew[j],1.e-6); // following convention of HJM

        //QLIB_VERIFY(skew[j] > TOO_SMALL,
        //  "BGM skew must be strictly positive"); 

        QLIB_VERIFY(skew[j] < 1. + TOO_SMALL,
            "BGM skew must be less than 1"); 

    }

    return;
}

//****************************************************************************//
//  end : setLiborSkew                                                        //
//****************************************************************************//

/** getFactorVols */
void SRMRatesLiborDiffuse::getFactorVols(
                                        vector< vector<double> > &facLdArray,
                                        long   firstLibor,
                                        long   lastLibor,
                                        double timeFirst,
                                        double timeLast) 
{
    static const string method = "SRMRatesLiborDiffuse::getFactorVols";

    long factor = -1,
         libor  = -1,
         firstFullLibor = -1;

    double u1 = -1.,
           u2 = -1.,
           u3 = -1.;

    if ( ( timeFirst < 0 ) || ( timeLast < 0 )  )
    {
        throw ModelException(method, "Invalid input for SRMRatesLiborDiffuse::getFactorVols ( timeFirst < 0  )");
    }

    if ( timeLast < timeFirst )
    {
        throw ModelException(method, "Invalid input for SRMRatesLiborDiffuse::getFactorVols ( timeLast < timeFirst )");
    }

    if ( ( firstLibor < 0) || ( firstLibor >= m_nbLibors) )
    {
        throw ModelException(method, "libor index out of range in SRMRatesLiborDiffuse::getFactorVols");
    }

    if ( ( lastLibor < 0) || ( lastLibor >= m_nbLibors) )
    {
        throw ModelException(method, "libor index out of range in SRMRatesLiborDiffuse::getFactorVols");
    }

    if ( firstLibor > lastLibor )
    {
        throw ModelException(method, " firstLibor must be < lastLibor in SRMRatesLiborDiffuse::getFactorVols");
    }

    // NOTE : for predictor-corrector in the risk-neutral measure, we have to allow for 
    //        many libors with resets < timeLast

    firstFullLibor = firstLibor;
    while (m_resetDates[firstFullLibor] < timeLast)
    {
        firstFullLibor++;
    }

    // set begin tail to zero
    for (libor = 0; libor < firstLibor; libor++)
    {
        for (factor = 0; factor < m_nbFactors; factor++)
        {
            facLdArray[libor][3*factor]   = 0.;
            facLdArray[libor][3*factor+1] = 0.;
            facLdArray[libor][3*factor+2] = 0.;
        }
    }

    for(libor = firstLibor; libor < firstFullLibor; libor++)
    {
        choleskyFactor(timeFirst, m_resetDates[libor]);

        for (factor = 0; factor < m_nbFactors; factor++)
        {
            // get factor loadings from model's parameters
            u1  =  m_D[libor][factor];
            u2  =  ( m_A[libor][factor] +
                m_B[libor][factor] * m_resetDates[libor] )
                * m_expReset[libor];
            u3 = - m_B[libor][factor] * m_expReset[libor];

            facLdArray[libor][3*factor]   = (CholMatrix_11[factor] * u1) +
                (CholMatrix_21[factor] * u2) +
                (CholMatrix_31[factor] * u3);
            facLdArray[libor][3*factor+1] = (CholMatrix_22[factor] * u2) +
                (CholMatrix_32[factor] * u3);
            facLdArray[libor][3*factor+2] =  CholMatrix_33[factor] * u3;
        }
    }

    // normal case
    choleskyFactor(timeFirst,timeLast);

    // when firstFullLibor = firstLibor + 1
    // the factor loadings have already been set
    for (libor = firstFullLibor; libor <= lastLibor; libor++)
    {
        for (factor = 0; factor < m_nbFactors; factor++)
        {
            // get factor loadings from model's parameters
            u1  =  m_D[libor][factor];
            u2  =  ( m_A[libor][factor] +
                m_B[libor][factor] * m_resetDates[libor] )
                * m_expReset[libor];
            u3 = - m_B[libor][factor] * m_expReset[libor];

            facLdArray[libor][3*factor]   = (CholMatrix_11[factor] * u1) +
                (CholMatrix_21[factor] * u2) +
                (CholMatrix_31[factor] * u3);
            facLdArray[libor][3*factor+1] = (CholMatrix_22[factor] * u2) +
                (CholMatrix_32[factor] * u3);
            facLdArray[libor][3*factor+2] =  CholMatrix_33[factor] * u3;
        }
    }

    // set tail to zero
    for (libor = lastLibor + 1; libor < m_nbLibors; libor++)
    {
        for (factor = 0; factor < m_nbFactors; factor++)
        {
            facLdArray[libor][3*factor]   = 0.;
            facLdArray[libor][3*factor+1] = 0.;
            facLdArray[libor][3*factor+2] = 0.;
        }
    }

    return;
}

//****************************************************************************//
//  end: getFactorVols                                                    //
//****************************************************************************//

// Cholesky factorization of three correlated BM
// see documentation
void SRMRatesLiborDiffuse::choleskyFactor(double timeFirst, double timeLast) 
{

    double i0, i1, i2, sqrtArg;

	double deltaT = timeLast - timeFirst;
	double sqrtDeltaT = sqrt(deltaT);

	double tempExpLast  = exp(m_lambda * timeLast);
	double tempExpFirst = exp(m_lambda * timeFirst);

	for (int factor = 0; factor < m_nbFactors; factor++)
    {
        // initialize Cholesky matrix
        CholMatrix_11[factor] = 0.; CholMatrix_21[factor] = 0.; CholMatrix_31[factor] = 0.;
        CholMatrix_22[factor] = 0.; CholMatrix_32[factor] = 0.;
        CholMatrix_33[factor] = 0.;

        // for the given factor construct Cholesky matrix
        if ( m_lambda > TOO_SMALL )
        {
			double invLambda = 1. / m_lambda;
            // non zero elements of first column of Cholesky matrix

            if ( deltaT > TOO_SMALL )
            {
                CholMatrix_11[factor] = sqrtDeltaT;

                i0 = invLambda * (tempExpLast - tempExpFirst);

                i1 = invLambda * (timeLast * tempExpLast -
                                  timeFirst * tempExpFirst - i0);

                CholMatrix_21[factor] = i0/CholMatrix_11[factor];
                CholMatrix_31[factor] = i1/CholMatrix_11[factor];
            }
            // else singular matrix  a11 = a21 = a31 = 0
            // AND ...? TODO : add error message ?

            // non zero elements of second column of Cholesky matrix
            tempExpLast  *= tempExpLast;
            tempExpFirst *= tempExpFirst;

            i0 = 0.5 * invLambda * (tempExpLast - tempExpFirst);

            i1 = 0.5 * invLambda * (timeLast * tempExpLast -
                                    timeFirst * tempExpFirst - i0);

            i2 = 0.5 * invLambda * (timeLast * timeLast * tempExpLast -
                                    timeFirst * timeFirst * tempExpFirst - 2.*i1);

            sqrtArg = i0 - (CholMatrix_21[factor] * CholMatrix_21[factor]);
            if ( sqrtArg > 0. )
            {
                CholMatrix_22[factor] = sqrt( sqrtArg );
                CholMatrix_32[factor] = (i1 - (CholMatrix_21[factor] * CholMatrix_31[factor]))/CholMatrix_22[factor];
            }
            // else singular matrix a22 = a32 = 0
            // TODO : add error message ?

            // non zero element of third column of Cholesky marix
            sqrtArg = i2 - (CholMatrix_31[factor] * CholMatrix_31[factor])
                         - (CholMatrix_32[factor] * CholMatrix_32[factor]);
            if ( sqrtArg > 0. )
            {
                CholMatrix_33[factor] = sqrt( sqrtArg );
            }
            // else singular matrix a33 = 0
        }
        else
        {
            // Cholesky factors for zero lambda
            CholMatrix_11[factor] = CholMatrix_21[factor] = sqrtDeltaT;
            CholMatrix_31[factor] = 0.5 * sqrtDeltaT * (timeLast + timeFirst);
            CholMatrix_33[factor] = invSqrtTwelfth * sqrtDeltaT * deltaT;
        }

        } // factor

   return;
}

//****************************************************************************//
//  end: choleskyFactor                                                       //
//****************************************************************************//

/** setLiborModel */
void SRMRatesLiborDiffuse::setLiborModel(
                   long nbFactors,
                   simulationType simType)
{
    static const string method = "SRMRatesLiborDiffuse::setLiborModel";

    int
         period = -1,
         libor = -1,
         factor = -1,
         i = -1,
         j  = -1;

    double
         zcbReset = -1.,
         zcbPay = -1.,
         volTime = -1.,
         instVol = -1.,
         invNbSteps = -1.,
         temp = -1.;

    m_nbFactors = nbFactors;
    m_simulationType = simType;

    m_gaussianVector.resize(3 * m_nbFactors);

    // paranoid programming
    m_resetDates.clear();
    m_payDates.clear();

    m_resetDates.reserve(500); 
    m_payDates.reserve(500);

    // now set up the schedule of libor reset and payment dates
	// the schedule has three parts : 
	// 1) an irregular schedule to ensure resolution at short end
	// 2) a regular schedule, with first reset at 12 months
	// 3) benchmark swaption expiries, ensuring that each swap start is a libor reset 
    const DateTime& spotDate = (liborUtilSP->getDiffYC())->getSpotDate();

	// to ensure high resolution in short end
	DateTimeArray irregularSchedule(5);
	irregularSchedule[0] = spotDate.rollDate(7);
	irregularSchedule[1] = spotDate.rollDate(14);
	irregularSchedule[2] = spotDate.rollDateInMonths(1);
	irregularSchedule[3] = spotDate.rollDateInMonths(3);
	irregularSchedule[4] = spotDate.rollDateInMonths(6);

    // regular schedule starts after one year
    DateTime startRegular = spotDate.rollDateInMonths(12);

    DateTimeArray* regularSchedule(
        SwapTool::dateArray(
                            startRegular,    // start here
                            m_liborInterval, // interval = count periods
                            "M",             // e.g. Y, M, W, D
                            0,               // 0=start @ basedate, 1=start @ baseDate + interval
                            1,               // arrayIncrement, usually +1 or -1
                            m_lastSimTime * 12 / m_liborInterval)); // how many dates, going out to 70 years

    vector<const DateTimeArray*> tempArray(3);
	tempArray[0] = &m_marketExpiriesDate;
    tempArray[1] = &irregularSchedule;
    tempArray[2] = regularSchedule;

    DateTimeArray TPDate(DateTime::merge(tempArray));
    // sort dates and remove duplicates
    DateTime::doSortUniq(TPDate);

    m_nbLibors = TPDate.size() - 1;

    // allocate memory for arrays of base class
    m_libor.resize(m_nbLibors);
    m_df.resize(m_nbLibors+1);
    m_liborVols.resize(m_nbLibors);
    m_liborSkew.resize(m_nbLibors);
    m_initialLibor.resize(m_nbLibors);
    m_accruals.resize(m_nbLibors);
    m_invAccruals.resize(m_nbLibors);

    m_fudgeA.resize(m_nbLibors); 
    m_fudgeD.resize(m_nbLibors);
    m_cosTheta.resize(m_nbLibors); 
    m_sinTheta.resize(m_nbLibors);

    m_A.resize(m_nbLibors);
    m_B.resize(m_nbLibors);
    m_D.resize(m_nbLibors);
    m_expReset.resize(m_nbLibors);
    m_facLdArray.resize(m_nbLibors);
    m_instLiborVol.resize(m_nbLibors);

    m_modelCorrelations.resize(m_nbLibors);

    double liborAccr = 0.,
           tmpTime = 0.;

    // returns Libor * accrual
    m_initialLibor = simpleFwdCurve(liborUtilSP->getDiffYC(),TPDate);

    for(i = 0; i < m_nbLibors; i++)
    {
       tmpTime = spotDate.yearFrac(TPDate[i]);
       liborAccr = (TPDate[i]).yearFrac(TPDate[i+1]);

       if (!Maths::isPositive(liborAccr))
       {
           throw ModelException(method, "zero accrual in libor model");
       }

       m_resetDates[i] = tmpTime;
       m_payDates[i] = tmpTime + liborAccr;
       m_accruals[i] = liborAccr;
       m_invAccruals[i] = 1. / liborAccr;
       m_initialLibor[i] /= liborAccr;

       m_A[i].resize(m_nbFactors);
       m_B[i].resize(m_nbFactors);
       m_D[i].resize(m_nbFactors);
       m_facLdArray[i].resize(3*m_nbFactors);
       m_instLiborVol[i].resize(m_nbFactors);
       m_modelCorrelations[i].resize(m_nbLibors);
    }

	m_a[0] = m_aInit[0] = 1.;      
	m_a[1] = m_aInit[1] = 0.;
	m_b[0] = m_bInit[0] = 0.3;      
	m_b[1] = m_bInit[1] = 0.1;
	m_d[0] = m_dInit[0] = 1.;      
	m_d[1] = m_dInit[1] = 1.;
	m_lambda = m_lambdaInit = 0.6; 

    const DateTime lastSimDate = liborUtilSP->getSimDates().back();

	// for now, we need to do all the libors (can't use the lastsimdate as done for HJM)
	m_lastLiborUsed = m_nbLibors - 1;

    // allocate memory for arrays of derived class
    if ( m_simulationType == predictorCorrector )
    {
        // for predictor-corrector extra arrays are needed
        m_previousLibor.resize(m_nbLibors);
        m_drift.resize(m_nbLibors);
        m_driftNew.resize(m_nbLibors);
        m_shock.resize(m_nbLibors);

        // we need 3*m_nbFactors random numbers per draw to
        // implement correlation structure with m_nbFactors
        m_v.resize(3 * m_nbFactors);

        m_factorVol.resize(m_nbLibors);
        for (libor = 0; libor < m_nbLibors; libor++ )
        {
            m_factorVol[libor].resize(3*m_nbFactors);
        }
    }
    else
    {
        m_v.resize(m_nbFactors);
    }

    // construct random lite generator
    if (m_simulationType == predictorCorrector)
    {
        const size_t rngSeed = 12345657UL;
        ranNumLiteSP.reset(new RngGaussianNRClass(randomIndex + rngSeed));
    }

    // set the skew
    const int nbParamDates = 1;
    m_paramDates.resize(nbParamDates,0.);

    qLeft=liborUtilSP->getQLeft();
    qRight=liborUtilSP->getQRight();
    zeroQ = Maths::isZero(qLeft) && Maths::isZero(qRight);

    QLIB_VERIFY(qLeft==qRight,
        "failed because qLeft!=qRight");
    // for now, same skew for all libors
    vector<double> skew(nbParamDates, qRight);
    setLiborSkew(skew);

    // 2 target correlations input 
    int nbCorr = 2L; 
    m_targetCorr.resize(nbCorr);
    m_corrSwapIndices.resize(2 * nbCorr);

    m_corrLower = 0.87;
    m_corrUpper = 0.95;

    QLIB_VERIFY(
        fabs(m_corrLower) <= 1., 
        "Error : m_corrLower is not a valid correlation");

    QLIB_VERIFY(fabs(m_corrUpper) <= 1., 
        "Error : m_corrUpper is not a valid correlation");

    // m_corrLower should not be larger than m_corrUpper
    if ( m_corrLower > m_corrUpper )
    {
       throw ModelException(method, "first correlation greater than second");
    }

    m_targetCorr[0] = m_corrLower;
    m_targetCorr[1] = m_corrUpper;

    return;
}
///////////////////////////////////////////////////////////
// end setLiborModel
///////////////////////////////////////////////////////////

/** finalize */
void SRMRatesLiborDiffuse::finalize(DateTimeArrayConstSP allDatesSP)
{
    static const string method("SRMRatesLiborDiffuse::finalize");

    DateTimeArray dfRequestedDates  = getSpotDates();
    DateTimeArray edfRequestedDates = getForwardDates();
    DateTimeArray edfForwardDates   = getForwardForwardDates();

	QLIB_VERIFY(liborUtilSP.get()!=NULL,
        "Internal error : liborUtil void"); 

    DateTimeArray diffusionDates = SRMUtil::calcDiffusionDates(today, *allDatesSP);

    const size_t Nall =  SRMUtil::getNumSimDates(today, *allDatesSP); // rename later
    const size_t Nedf =  edfRequestedDates.size(); // rename later
    const size_t Ndf  =  dfRequestedDates.size();

    const size_t numSimDates = SRMUtil::getNumSimDates(today, *allDatesSP);

    QLIB_VERIFY(numSimDates == diffusionDates.size(),
        "Internal error : problem in timeline"); 

    domLnMONEY = vector<double>(saveDomLnMoney? numSimDates-1: 0);
    sigmaR = vector<double>(saveSigmaR? numSimDates-1: 0);

    effRSVol = vector<double>(numSimDates-1); // for zero q, effR * svol[Nall]
    QPivot = vector<double>(numSimDates-1); // [Nall]
    svol = vector<double>(numSimDates-1); /* contains SpotVol initially. If zeroQ is false
                                             then processed in calcSigmaRParams */ // [Nall]

    // cached variables
    // TODO getSubIndexes should throw an error if input date is not found in allDates
    QLIB_VERIFY(DateTime::isSubset((*allDatesSP), dfRequestedDates),
        "Internal error : dfRequestDates not in AllDates"); 

    dfIndexes = DateTime::getIndexes((*allDatesSP),  dfRequestedDates); // indexes for when we save discount factors [Ndf]

    QLIB_VERIFY(DateTime::isSubset((*allDatesSP), edfRequestedDates),
        "Internal error : edfRequestDates not in AllDates"); 

    expDFIndexes = DateTime::getIndexes((*allDatesSP), edfRequestedDates); // indexes for when we compute expected  [Nedf]

    if (!fx && expDFIndexes.empty() && dfIndexes.empty())
    {
        throw ModelException(method, "Internal error - no "
                                     "results required from diffused path");
    }

    todayIdx = today.find((*allDatesSP));
    // df contains history disc factors only. So must resize // i.e. increase the size
    
    //Initialise to zero: some products (e.g. TARN) use zero discount factor in order to ignore coupon payments for past dates.
    df.resize(Ndf, 0.0);
    getDiscYCIdx();
    getDiffYCIdx();

    todayIndex = today.findUpper(dfRequestedDates);
    if (todayIndex >= 0 && todayIndex != dfRequestedDates.size() && dfRequestedDates[todayIndex] == today)
        ++todayIndex; // skip today as diffusion results are always in futureDates

    QLIB_VERIFY( todayIndex == dfRequestedDates.size() || dfRequestedDates[todayIndex] > today,
        "Internal error : problem in dfRequestedDates"); 

    // now make life easier - add request for index off the end
    size_t errIndex = (*allDatesSP).size()+numSimDates+1;
    dfIndexes.push_back(errIndex);
    expDFIndexes.push_back(errIndex);

	QLIB_VERIFY( expDFIndexes.size() > 0,
		"Internal error : not enough expDF dates"); 

	VectorDiffusedStates.clear();
    VectorDiffusedStates.resize(expDFIndexes.size());

    // assert that the diffusion dates have no past dates
    // TODO relax this requirement
    QLIB_VERIFY(diffusionDates[0] >= today,
        "Internal error : diffusion can't have past dates"); 

    liborUtilSP->computeLogDiscFactor(logDiscFactor); // [Nall]
    sqrtYearFrac = SRMUtil::computeSqrtYearFrac(diffusionDates);       // [Nall]                       // discount factors
    m_diffuseTime = SRMUtil::computeYearFrac(today, diffusionDates); 

    for(size_t i=0; i < ycForwardsDB.size(); ++i)
    {
        IYieldCurveConstSP  yc = ycForwardsDB[i].first;
        vector<double> &   fwd = ycForwardsDB[i].second;
        liborUtilSP->computeLogDiscFactor(edfForwardDates, fwd, yc); // initialize fwd between t0 and T
    }

    if (liborUtilSP->getMomentMatchingFlag())
    {
        liborUtilSP->computeLogDiscFactor(dfRequestedDates, originalDFs, liborUtilSP->getDiscYC());
    }


    assetDatesFixed = true; // no more requested dates can be added
    initialized = true; // fully initialized
}
//////////////////////////////////////////////////////////
// end finalize
//////////////////////////////////////////////////////////

/** generate path across all dates. */
void SRMRatesLiborDiffuse::generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/)
{
    static const string method = "SRMRatesLiborDiffuse::generatePath";

    // random numbers for this path
    switch ( m_simulationType )
    {
    case smallStep:
        {
            generatePathSmallStep(rngMgr);
        }
        break;
    case predictorCorrector:
        {
            generatePathPredCorr(rngMgr);
        }
        break;
    default :
        {
            throw ModelException(method, "Invalid simulation method for the Libor model");
        }
    }

    return;
}

////////////////////////////////////////////////////////////
// generatePath 
////////////////////////////////////////////////////////////

// generatePathSmallStep
/** generate path across all dates. */
// NOT WORKING !!!
void SRMRatesLiborDiffuse::generatePathSmallStep(IQMCRNGManagerSP rngMgr)
{
    static const string method = "SRMRatesLiborDiffuse::generatePath";

    int numDates = logDiscFactor.size();

    int dfDatePos = todayIndex; // position in dfIndexes/df.
    int expDatePos = expDFIndexes.front() == 0? 1: 0;

    int dfDateIdx = dfIndexes[dfDatePos];
    int expDateIdx = expDFIndexes[expDatePos];
    int stopIdx = Maths::min(dfDateIdx, expDateIdx); // when to do something

    const vector<double> liborSkew = getLiborSkew();

    double LnMONEY = 0.,
        shortRate = 0.;

    double sigmaFX = 0; // FX data
    if (fx.get())
    {
        // storing random numbers and returning vol
        sigmaFX = fx->begin(rngMgr);
    }

    // initialize libors
    for (int libor = 0; libor < m_nbLibors; libor++)
    {
        double tempLibor = m_initialLibor[libor];
        m_libor[libor] = tempLibor;
    }

    shortRate = log(1. + m_accruals[0] * m_libor[0]) / m_accruals[0];

    // random numbers for this path
    for(int factor = 0; factor < m_nbFactors; factor++)
    {
        m_gaussianVector[factor] = rngMgr->getCorrelatedRandoms(randomIndex, factor);
    }

    double TNow = 0., TNext = 0.;
    int firstLiborAlive = 0;

    // libors already initialized

    // loop over all dates
    for (int i = 0; i < numDates; /* increment in body */)
    {
        double sqrtDeltaT = sqrtYearFrac[i];
        double deltaT = sqrtDeltaT * sqrtDeltaT;

        TNow = TNext;
        TNext += deltaT;

        for(int factor = 0; factor < m_nbFactors; factor++)
        {
            getLiborVols(
                m_instLiborVol,
                TNow,
                firstLiborAlive,
                m_lastLiborUsed);
        }

        /* Step2: update LnMONEY */
        /*-----------------------*/
        LnMONEY += shortRate * deltaT;

        // update libors
        // 1. initialize vectors for Libors update
        //    to brownian motion increments
        // 2. for each factor get libor instantaneous volatilities
        for (int factor = 0; factor < m_nbFactors; factor++)
        {
            m_v[factor] = sqrtDeltaT * m_gaussianVector[factor][i];
        }

        // initialize update term
        double update = 0.0;

        // update Libors which have not reset yet
        // simulating in risk-neutral measure
        for (int libor = firstLiborAlive; libor <= m_lastLiborUsed;
            libor++)
        {
            double deltL = m_accruals[libor] * m_libor[libor];

            double inv1deltL = 1.;

            if ( deltL > -1.0 )
            {
                inv1deltL = 1.0 / (1.0 + deltL);
            }
            else
            {
                throw ModelException("forwardMC", "Error in Monte Carlo: division by zero!");
            }

            double skew = liborSkew[libor];

            // now take into account log-norm shift
            double shift = m_accruals[libor] * (1. - skew)
                * m_initialLibor[libor];
            deltL *= skew;
            deltL += shift;

            // this is (L_{n}+a)delta_{n}/(1 + L_{n}delta_{n})
            double driftFactor = deltaT * deltL * inv1deltL;

            for (int factor = 0; factor < m_nbFactors; factor++)
            {
                double tempVol = m_instLiborVol[libor][factor];

                // add to the multiplicative update factor
                // the inner product of the vectors for update
                // with the n-th volatility vector
                update += m_v[factor] * tempVol;

                // update the vectors for next run of the loop
                m_v[factor] = driftFactor * tempVol;
            }

            m_libor[libor] = m_libor[libor] + update *
                ( skew * m_libor[libor] +
                (1. - skew) * m_initialLibor[libor]);
        } // n - end libors update loop

        // tracking which libors remain alive
        while (TNext > m_resetDates[firstLiborAlive])
        {
            firstLiborAlive++;
        }

        int liborIndex = Maths::max(0, firstLiborAlive - 1);

        shortRate = log(1. + m_accruals[liborIndex] 
        * m_libor[liborIndex]) /
            m_accruals[liborIndex];

        /* Step5: add the CUPS component if necessary */
        if (fx.get())
        {
            //double F0 = IrFxCorr0 * sigmaFX * RootDelT;
            sigmaFX = fx->moveToNextDate(LnMONEY, i); // get new sigmaFX
        }

        // We need this for each equity ccy
        if (!domLnMONEY.empty())
        {
            domLnMONEY[i] = LnMONEY; // needed by fx's linked to dom
        }

        i++; /* up the loop counter. This comes back to the fact with eg
             5 dates there are 4 diffusion steps. It's easier to
             identify the index of the date you want rather than the
             index of the following one */

        if (i + todayIdx == stopIdx)
        {
            // TODO double check later with past dates
            // hit an 'event' - record something
            if (stopIdx == dfDateIdx)
            {
                // save the value
                df[dfDatePos] = LnMONEY;
                dfDatePos++;
                dfDateIdx = dfIndexes[dfDatePos];
            }

            if (stopIdx == expDateIdx)
            {
                // save 
                // index of last "frozen" libor
                VectorDiffusedStates[expDatePos].freezeIdx = 
                    Maths::max(0, firstLiborAlive - 1);
                // offset from today
                VectorDiffusedStates[expDatePos].obsTime = TNext; 
                double dt = m_resetDates[firstLiborAlive] - TNext;
                VectorDiffusedStates[expDatePos].dt = dt; 
                VectorDiffusedStates[expDatePos].shortRate = shortRate; 

                // first, fill in discount factors maturing in the past (<= today)
                int libor = 0;
                while( libor < firstLiborAlive )
                {
                    m_df[libor] = 1.;
                    libor++;
                }

                m_df[firstLiborAlive] = exp( - shortRate * dt);

                double previousDf = 1.;
                // forward ZCB, with forward date m_resetDates[firstLiborAlive]
                while( libor <= m_lastLiborUsed )
                {
                    double deltL = m_accruals[libor] * m_libor[libor];
                    double inv1deltL = 1.;

                    if (deltL > -1.)
                    {
                        inv1deltL = 1. / (1. + deltL);
                    }

                    previousDf *= inv1deltL;
                    m_df[libor] = previousDf;
                    libor++;
                }

                // could / should drop next loop
                // libor rates constant from index = m_lastLiborUsed + 1 onwards
                while( libor < m_nbLibors )
                {
                    double deltL = m_accruals[libor] * m_libor[m_lastLiborUsed];
                    double inv1deltL = 1.;

                    if (deltL > -1.)
                    {
                        inv1deltL = 1. / (1. + deltL);
                    }

                    previousDf *= inv1deltL;
                    m_df[libor] = previousDf;
                    libor++;
                }

                VectorDiffusedStates[expDatePos].stateDf = m_df; // engine discount factors 
                expDatePos++;
                // that's why we pushed an extra element to expDFIndexes;
                expDateIdx = expDFIndexes[expDatePos];

            }

            // figure out when we need to stop next time
            stopIdx = Maths::min(dfDateIdx, expDateIdx); // refresh
        }

    } // end iterating forward smallStep

    // then do exp on df vector
    for (size_t j = todayIndex; j < df.size(); j++){
        df[j] = exp(-df[j]);
    }

    if (fx.get())
    {
        fx->end();
    }
}
/////////////////////////////////////////////////////////
// end generatePathSmallStep
/////////////////////////////////////////////////////////

// generatePathPredCorr
/** generate path across all dates. */
void SRMRatesLiborDiffuse::generatePathPredCorr(IQMCRNGManagerSP rngMgr)
{
    static const string method = "SRMRatesLiborDiffuse::generatePath";

    int numDates = logDiscFactor.size();

    int dfDatePos = todayIndex; // position in dfIndexes/df.
    int expDatePos = expDFIndexes.front() == 0? 1: 0;
    
    int dfDateIdx = dfIndexes[dfDatePos];
    int expDateIdx = expDFIndexes[expDatePos];
    int stopIdx = Maths::min(dfDateIdx, expDateIdx); // when to do something

    const vector<double> liborSkew = getLiborSkew();

    double LnMONEY = 0.,
           shortRate = 0.;

    double sigmaFX = 0; // FX data

    if (fx.get())
    {
        // storing random numbers and returning vol
        sigmaFX = fx->begin(rngMgr);
    }

    // initialize libors
    for (int libor = 0; libor < m_nbLibors; libor++)
    {
        double tempLibor = m_initialLibor[libor];
        m_libor[libor] = tempLibor;
        m_previousLibor[libor] = tempLibor;
    }

    shortRate = log(1. + m_accruals[0] * m_libor[0]) / m_accruals[0];

    // random numbers for this path
    // Pack random numbers supplied from outside into first m_nbFactor arrays
    // put internally generated random numbers into next 2 * m_nbFactor arrays
    // as a consequence, only first array is correlated with other assets

    for(int factor = 0; factor < m_nbFactors; factor++)
    {
        m_gaussianVector[factor] = rngMgr->getCorrelatedRandoms(randomIndex, factor);
    }

    // TODO - FIXME - IE, NEED TO CHANGE IQMCRNGManager TO AVOID HACK USING ranNumLite
    // get vector of independent gaussian random numbers
    int dim = 2 * m_nbFactors * numDates;
    double *dummy = new double[dim];
    ranNumLiteSP->getRandVector(dummy, dim);
    for(int factor = 0; factor < 2 * m_nbFactors; factor++)
    {
        m_gaussianVector[m_nbFactors + factor] = dummy + factor * numDates;
    }

    // start simulation
    double TNow = 0., TNext = 0.;
    int firstLiborAlive = 0;

    int currentIdx = todayIdx, 
        nextIdx = todayIdx;

    while (currentIdx < todayIdx + numDates)
    {
        // we jump to smaller of :
        //     a) m_predCorrStep
        //     b) next "stopIdx" time - i.e., next time mmkt or zcb needed
        //     c) next libor reset
        TNow = TNext;

        nextIdx = stopIdx;
        // if the jump is too large, reduce index
        while ( m_diffuseTime[nextIdx] > 
            Maths::max(m_diffuseTime[currentIdx+1], TNow + m_predCorrStep) )
        {
            nextIdx = Maths::max(currentIdx + 1, (int) ( (currentIdx + nextIdx)/2 ));
        }

        TNext = m_diffuseTime[nextIdx];

        double deltaT = TNext - TNow;

        // Step2: update LnMONEY 
        //-----------------------
        LnMONEY += shortRate * deltaT;

        // update firstLiborAlive to ensure that Libors that have already
        // fixed are not evolved by the simulation
        // Note that we may have TNext > m_resetDates[firstLibor]
        // - see below
        if (firstLiborAlive < m_lastLiborUsed)
        {
            while( 
                  (firstLiborAlive < m_lastLiborUsed) &&
				  (TNow > m_resetDates[firstLiborAlive])
                 )
            {
                firstLiborAlive++;
            }
        }

        // initialize vectors for libor drift updates
        for (int factor = 0; factor < m_nbFactors; factor++)
        {
            m_v[3 * factor]     = 0.;
            m_v[3 * factor + 1] = 0.;
            m_v[3 * factor + 2] = 0.;
        }

        getFactorVols(
            m_factorVol, 
            firstLiborAlive, 
            m_lastLiborUsed, 
            TNow,
            TNext); 

        // simulation in risk-neutral measure!! 
        for (int libor = firstLiborAlive; 
             libor <= m_lastLiborUsed;
             libor++
            )
        {
            m_shock[libor] = 0.;
            m_drift[libor] = 0.;

            double deltL = m_accruals[libor] * m_previousLibor[libor];
            double inv1deltL = 1.;

            if (deltL > -1.0)
            {
                inv1deltL = 1.0 / (1.0 + deltL);
            }
            else
            {
                throw ModelException("forwardMC", "Error in Monte Carlo: Libor rates negative and too large"); 
            }

            double skew = liborSkew[libor];

            double shift = m_accruals[libor] * (1. - skew) * m_initialLibor[libor];
            deltL *= skew;
            deltL += shift;

            // i.e delta_{n}(L_{n}+a)/(1+delta_{n}L_{n})
            double driftFactor = deltL * inv1deltL;

            for (int factor = 0; factor < m_nbFactors; factor++)
            {
                double volFactor0 = m_factorVol[libor][3*factor]; 
                double volFactor1 = m_factorVol[libor][3*factor+1];
                double volFactor2 = m_factorVol[libor][3*factor+2];

                double skewVolFactor0 = skew * volFactor0;
                double skewVolFactor1 = skew * volFactor1;
                double skewVolFactor2 = skew * volFactor2;

                m_v[3*factor]   += driftFactor * volFactor0;
                m_v[3*factor+1] += driftFactor * volFactor1;
                m_v[3*factor+2] += driftFactor * volFactor2;

                m_shock[libor] +=   
                    ( 
                    skewVolFactor0 * m_gaussianVector[3*factor][currentIdx]   +
                    skewVolFactor1 * m_gaussianVector[3*factor+1][currentIdx] +
                    skewVolFactor2 * m_gaussianVector[3*factor+2][currentIdx] );

                m_drift[libor] += 
                    ( 
                    skewVolFactor0 * (m_v[3*factor]   - 0.5 * skewVolFactor0) +
                    skewVolFactor1 * (m_v[3*factor+1] - 0.5 * skewVolFactor1) +
                    skewVolFactor2 * (m_v[3*factor+2] - 0.5 * skewVolFactor2) );

            } // for loop: factor=0 to m_nbFactors-1    
        } // for loop: n=m_lastLiborUsed downto m_firstLiborAlive

        for (int libor = firstLiborAlive; 
             libor <= m_lastLiborUsed;
             libor++
            )
        {
            double skew = liborSkew[libor];
            skew = Maths::max(skew,1.e-6); // TODO : properly

            m_libor[libor] = 
                ( ( skew * m_previousLibor[libor] + 
                    (1. - skew) * m_initialLibor[libor] 
                  ) * exp( m_drift[libor] + m_shock[libor] ) 
                      - (1. - skew) * m_initialLibor[libor] ) / skew; 
        }

        // SECOND STEP: new drift                         //
        // initialize vectors for libor drift updates
        for (int factor = 0; factor < m_nbFactors; factor++)
        {
            m_v[3 * factor]     = 0.;
            m_v[3 * factor + 1] = 0.;
            m_v[3 * factor + 2] = 0.;
        }

        for (int libor = firstLiborAlive; 
             libor <= m_lastLiborUsed;
             libor++
            )
        {
            m_driftNew[libor] = 0.;

            double deltL = m_accruals[libor] * m_libor[libor];
            double inv1deltL = 1.;

            if (deltL > -1.)
            {
                inv1deltL = 1. / (1. + deltL);
            }
            else
            {
                //throw ModelException( "Error in Monte Carlo: Libor rates negative and too large"); 
            }

            double skew = liborSkew[libor];

            double shift = m_accruals[libor] * (1. - skew) * m_initialLibor[libor];
            deltL *= skew;
            deltL += shift;

            // i.e delta_{n}(L_{n}+a)/(1+delta_{n}L_{n})
            double driftFactor = deltL * inv1deltL;

            for (int factor = 0; factor < m_nbFactors; factor++)
            {
                double volFactor0 = m_factorVol[libor][3 * factor];
                double volFactor1 = m_factorVol[libor][3 * factor+1];
                double volFactor2 = m_factorVol[libor][3 * factor+2];

                double skewVolFactor0 =  skew * volFactor0;
                double skewVolFactor1 =  skew * volFactor1;
                double skewVolFactor2 =  skew * volFactor2;

                m_v[3*factor]   += driftFactor * volFactor0;
                m_v[3*factor+1] += driftFactor * volFactor1;
                m_v[3*factor+2] += driftFactor * volFactor2;

                m_driftNew[libor] += 
                    ( 
                    skewVolFactor0 * (m_v[3*factor]   - 0.5 * skewVolFactor0) +
                    skewVolFactor1 * (m_v[3*factor+1] - 0.5 * skewVolFactor1) +
                    skewVolFactor2 * (m_v[3*factor+2] - 0.5 * skewVolFactor2) );

            } // for loop: factor=0 to m_nbFactors-1

            m_drift[libor] = 
                0.5 * (m_drift[libor] + m_driftNew[libor]);

        } // for loop: n=m_lastLiborUsed down to m_firstLiborAlive

        for (int libor = firstLiborAlive; 
             libor <= m_lastLiborUsed;
             libor++
            )
        {
            double skew = liborSkew[libor];
            skew = Maths::max(skew,1.e-6); // TODO : properly

            m_libor[libor] = 
                ( ( skew * m_previousLibor[libor] + 
                   (1. - skew) * m_initialLibor[libor] 
                  ) * exp( m_drift[libor] + m_shock[libor] ) 
                      - (1. - skew) * m_initialLibor[libor] 
                ) / skew; 

            m_previousLibor[libor] = m_libor[libor]; 
        }

        /////////////////////////////////////////////////

        // tracking which libors remain alive
		while ((firstLiborAlive < m_lastLiborUsed) &&
			   (TNext >= m_resetDates[firstLiborAlive]))
        {
            firstLiborAlive++;
        }

        int liborIndex = Maths::max(0, firstLiborAlive - 1);

        shortRate = log(1. + m_accruals[liborIndex] * m_libor[liborIndex]) /
            m_accruals[liborIndex];

        // Step5: add the CUPS component if necessary 
        if (fx.get())
        {
            // TODO - add CUPS adjustment
            //double F0 = IrFxCorr0 * sigmaFX * RootDelT;
            // get new sigmaFX   
            sigmaFX = fx->moveToNextDate(LnMONEY, currentIdx);      
        }      

        // We need this for each equity ccy 
        if (!domLnMONEY.empty()) 
        {
            domLnMONEY[currentIdx] = LnMONEY; // needed by fx's linked to dom
        }

        ////////////////////////////
        // record 
        ///////////////////////////
        if (nextIdx == stopIdx)
        { 
            // TODO double check later with past dates
            // hit an 'event' - record something
            if (stopIdx == dfDateIdx)
            {
                // save the value
                df[dfDatePos] = LnMONEY;

                // find index of next dfDate
                dfDatePos++;
                dfDateIdx = dfIndexes[dfDatePos];
            }

            if (stopIdx == expDateIdx)
            {
                // save 
                // index of last "frozen" libor
                VectorDiffusedStates[expDatePos].freezeIdx = 
                    Maths::max(0, firstLiborAlive - 1);
                // offset from today
                VectorDiffusedStates[expDatePos].obsTime = TNext; 
                double dt = m_resetDates[firstLiborAlive] - TNext;
                VectorDiffusedStates[expDatePos].dt = dt; 
                VectorDiffusedStates[expDatePos].shortRate = shortRate; 

                // first, fill in discount factors maturing in the past (<= today)
                int libor = 0;
                while( libor < firstLiborAlive )
                {
                    m_df[libor] = 1.;
                    libor++;
                }

                m_df[firstLiborAlive] = exp( - shortRate * dt);

                double previousDf = 1.;
                // forward ZCB, with forward date m_resetDates[firstLiborAlive]
                while( libor <= m_lastLiborUsed )
                {
                    double deltL = m_accruals[libor] * m_libor[libor];
                    double inv1deltL = 1.;

                    if (deltL > -1.) 
                    {
                        inv1deltL = 1. / (1. + deltL);
                    } 
                    previousDf *= inv1deltL;
                    m_df[libor] = previousDf;
                    libor++;
                }

                // could / should drop next loop 
                // libor rates constant from index = m_lastLiborUsed + 1 onwards
                while( libor < m_nbLibors )
                {
                    double deltL = m_accruals[libor] * m_libor[m_lastLiborUsed];
                    double inv1deltL = 1.;

                    if (deltL > -1.)
                    {
                        inv1deltL = 1. / (1. + deltL);
                    }
                    previousDf *= inv1deltL;
                    m_df[libor] = previousDf;
                    libor++;
                }

                VectorDiffusedStates[expDatePos].stateDf = m_df; // engine discount factors 

                // find index of next expDate
                expDatePos++;
                // that's why we pushed an extra element to expDFIndexes;
                expDateIdx = expDFIndexes[expDatePos]; 

            }

            // figure out when we need to stop next time
            stopIdx = Maths::min(dfDateIdx, expDateIdx); // refresh
        }               

        currentIdx = nextIdx;

    } // dates

    // then do exp on df vector
    for (size_t j = todayIndex; j < df.size(); j++){
        df[j] = exp(-df[j]);
    }

    if (fx.get())
    {
        fx->end();
    }

	delete dummy;

}
/////////////////////////////////////////////////////////
// end generatePathPredCorr
/////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// Notes :
// 1) at the moment, this function can return only the diffused curve
// n) needs to be optimized
///////////////////////////////////////////////////////////////////////////
/** getExpectedDiscFactor */
double SRMRatesLiborDiffuse::getExpectedDiscFactor(size_t idx, FwdIdx i, FwdIdx j)
{
    static const string method = "SRMRatesLiborDiffuse::getExpectedDiscFactor";

    const int iEDF = getFwdIdx2EdfIdx(i);

    const DiffusedState & state = VectorDiffusedStates[iEDF];

    // wasteful: FIXME: create yearFrac[] array for all ForwardForward dates
    // so that T-t == yearFrac[j] - yearFrac[i]

    DateTimeArray edfDates = getForwardForwardDates();
    DateTime firstMatDate = edfDates[i];

    double stubDf = 1.;

    double t = state.obsTime;
    double T = t + firstMatDate.yearFrac(edfDates[j]) ;

    int indexStart = state.freezeIdx, // index of last frozen libor
        maxIndex = (state.stateDf).size() - 1;

    int index = indexStart;

    if (index > maxIndex)
    {
        throw ModelException(method, "Indexing error");   
    }

    // index of last libor resetting before maturity of interest
    while (index < maxIndex - 1)
    {
        if ( m_resetDates[index+1] > T ) break;
        index++;
    }

    double df = (state.stateDf)[index]; 
    double dt = T - m_resetDates[index]; 

    double libor = 0.;

    // now calculate the libor corresponding to last reset before
    // interesting maturity
    if (dt > 0.)
    {
        if (index < maxIndex - 1) 
        {
            double dfAfter = (state.stateDf)[index+1]; 
            double accr = m_resetDates[index+1] - m_resetDates[index];
            libor = (df/dfAfter - 1.) / accr;  
            stubDf = 1./(1. + dt * libor);
        }
        else if (index > indexStart) 
        {
            double dfBefore = (state.stateDf)[index]; 
            double accr = m_resetDates[index] - m_resetDates[index-1];
            libor = (dfBefore/df - 1.) / accr;  
            stubDf = 1./(1. + dt * libor);
        }
        else
        {
            // use last fwd rate to extrapolate
            double dfAfter = (state.stateDf)[maxIndex - 1]; 
            double accr = m_resetDates[index+1] - m_resetDates[index];
            libor = (df/dfAfter - 1.) / accr;  
            double rate = log(df/dfAfter) / accr;
            stubDf = exp(-accr * rate);
        }

        df *= stubDf;
    }

    return df;
}

////////////////////////////////////////////////////
// end getExpectedDiscFactor 
////////////////////////////////////////////////////

void SRMRatesLiborDiffuse::printModel(void) const
{
    ofstream outFileParams("c:/temp/BGM_Params.txt");
    ofstream outFileCalib("c:/temp/BGM_Calibration.txt");
    ofstream outFileCor("c:/temp/BGM_Correlation.txt");

    int i = -1, 
        libor = -1, 
        libor2 = -1, 
        factor = -1, 
        row = -1, 
        column = -1;


    if ( (outFileParams != NULL) &&
         (outFileCalib != NULL) &&
         (outFileCor != NULL)
        )
    {
        outFileParams.precision(6);
        outFileCalib.precision(2);
        outFileCor.precision(3);

        outFileParams << "** number of libors in model = " << m_nbLibors << std::endl;
        outFileParams << std::endl;

        outFileParams << "** libor expiry dates and caplet implied vols" << std::endl;
        for (libor = 0; libor < m_nbLibors; libor++)
        {
            outFileParams << m_resetDates[libor] << "   " 
                << (100.*m_liborVols[libor]) << std::endl;
        }
        outFileParams << std::endl;

        outFileParams << std::endl;
        outFileParams << "** Parameter vector lambda:\n\n";
        
        outFileParams.width(12);
        outFileParams.setf(ios::left | ios::fixed);
        outFileParams  << std::fixed << m_lambda;

        outFileParams << std::endl;
        outFileParams << "** Parameter matrix A:\n\n";
        for (libor = 0; libor < m_nbLibors; libor++)
        {
            for (factor = 0; factor < m_nbFactors; factor++)
            {
                outFileParams.width(12);
                outFileParams.setf(ios::left | ios::fixed);
                outFileParams << std::fixed  << m_A[libor][factor];
            }
            outFileParams << std::endl;
        }

        outFileParams << std::endl;
        outFileParams << "** Parameter matrix B:\n\n";
        for (libor = 0; libor < m_nbLibors; libor++)
        {
            for (factor = 0; factor < m_nbFactors; factor++)
            {
                outFileParams.width(12);
                outFileParams.setf(ios::left | ios::fixed);
                outFileParams << std::fixed  << m_B[libor][factor];
            }
            outFileParams << std::endl;
        }

        outFileParams << std::endl;
        outFileParams << "** Parameter matrix D:\n\n";
        for (libor = 0; libor < m_nbLibors; libor++)
        {
            for (factor = 0; factor < m_nbFactors; factor++)
            {
                outFileParams.width(12);
                outFileParams.setf(ios::left | ios::fixed);
                outFileParams << std::fixed << m_D[libor][factor];
            }
            outFileParams << std::endl;
        }

        ///////////////////////////////////////////////////////////
        //  swaption vol grid
        ///////////////////////////////////////////////////////////

        if ( m_modelSwaptionVols.size() > 1 )
        {
            outFileCalib << "** market swaption matrix" << std::endl;

			for (column = 0; column < (int)m_marketTenors.size(); column++)
            {
                outFileCalib.width(8);
                outFileCalib.setf(ios::left | ios::fixed);
                outFileCalib << std::fixed 
                    << "________" ;
            }
            outFileCalib << std::endl;
            ///////////////////////////////////////////////////

			for (column = 0; column < (int)m_marketTenors.size(); column++)
            {
                outFileCalib.width(8);
                outFileCalib.setf(ios::left | ios::fixed);
                outFileCalib << std::fixed 
                    << m_marketTenors[column];
            }
            outFileCalib.width(8);
            outFileCalib.setf(ios::left | ios::fixed);
            outFileCalib << std::fixed 
                << "|" ;
            outFileCalib << std::endl;

			for (column = 0; column < (int)m_marketTenors.size(); column++)
            {
                outFileCalib.width(8);
                outFileCalib.setf(ios::left | ios::fixed);
                outFileCalib << std::fixed 
                    << "________" ;
            }
            outFileCalib << std::endl;
            ////////////////////////////////////////////////////

			for (row = 0; row < (int)m_marketExpiries.size(); row++)
            {
				for (column = 0; column < (int)m_marketTenors.size(); column++)
                {
                    outFileCalib.width(8);
                    outFileCalib.setf(ios::left | ios::fixed);
                    outFileCalib << std::fixed 
                        << (100. * m_marketSwaptionVols[row][column]);
                }
                outFileCalib << std::endl;
            }

            outFileCalib << std::endl;
            outFileCalib << "** analytic swaption matrix" << std::endl;

			for (column = 0; column < (int)m_marketTenors.size(); column++)
            {
                outFileCalib.width(8);
                outFileCalib.setf(ios::left | ios::fixed);
                outFileCalib << std::fixed 
                    << "________" ;
            }
            outFileCalib << std::endl;
            /////////////////////////////////////////////////////////

			for (column = 0; column < (int)m_marketTenors.size(); column++)
            {
                outFileCalib.width(8);
                outFileCalib.setf(ios::left | ios::fixed);
                outFileCalib << std::fixed 
                    << m_marketTenors[column];
            }
            outFileCalib.width(8);
            outFileCalib.setf(ios::left | ios::fixed);
            outFileCalib << std::fixed 
                << "|" ;
            outFileCalib << std::endl;
            
			for (column = 0; column < (int)m_marketTenors.size(); column++)
            {
                outFileCalib.width(8);
                outFileCalib.setf(ios::left | ios::fixed);
                outFileCalib << std::fixed 
                    << "________" ;
            }
            outFileCalib << std::endl;
            /////////////////////////////////////////////////////////

			for (row = 0; row < (int)m_marketExpiries.size(); row++)
            {
				for (column = 0; column < (int)m_marketTenors.size(); column++)
                {
                    outFileCalib.width(8);
                    outFileCalib.setf(ios::left | ios::fixed);
                    outFileCalib << std::fixed 
                        << (100. * m_modelSwaptionVols[row][column]);
                }
                outFileCalib << std::endl;
            }

            outFileCalib << std::endl;
            outFileCalib << "** analytic calibration errors " << std::endl;

			for (column = 0; column < (int)m_marketTenors.size(); column++)
            {
                outFileCalib.width(8);
                outFileCalib.setf(ios::left | ios::fixed);
                outFileCalib << std::fixed 
                    << "________" ;
            }
            outFileCalib << std::endl;
            //////////////////////////////////////////////////////////

			for (column = 0; column < (int)m_marketTenors.size(); column++)
            {
                outFileCalib.width(8);
                outFileCalib.setf(ios::left | ios::fixed);
                outFileCalib << std::fixed 
                    << m_marketTenors[column];
            }
            outFileCalib.width(8);
            outFileCalib.setf(ios::left | ios::fixed);
            outFileCalib << std::fixed 
                << "|" ;
            outFileCalib << std::endl;
            
			for (column = 0; column < (int)m_marketTenors.size(); column++)
            {
                outFileCalib.width(8);
                outFileCalib.setf(ios::left | ios::fixed);
                outFileCalib << std::fixed 
                    << "________" ;
            }
            outFileCalib << std::endl;
            //////////////////////////////////////////////////////////

            int calibIndex = 0;
			for (row = 0; row < (int)m_marketExpiries.size(); row++)
            {
				for (column = 0; column < (int)m_marketTenors.size(); column++)
                {
                    outFileCalib.width(8);
                    outFileCalib.setf(ios::left | ios::fixed);
                    if ((column == 0) || (SRMRound(m_marketExpiries[row]) > 0))
                    {
                        outFileCalib 
                            << std::fixed
                            << 100. * (m_modelSwaptionVols[row][column] - m_marketSwaptionVols[row][column]);
                    }
                    else
                    {
                        outFileCalib << std::fixed << " -- ";
                    }
                }
                outFileCalib << std::endl;
            }
        }

        // Print correlations to different file
        // note : next line must be consistent with the one in setModelCorrelationMatrix above 
        const int maxNbCorr = Maths::min(100L,m_nbLibors-1);

        outFileCor << "\n** Initial correlation matrix" << std::endl;

        for (libor = 0; libor < maxNbCorr; libor++)
        {
            for (libor2 = 0; libor2 < maxNbCorr; libor2++)
            {
                outFileCor.width(8);
                outFileCor.setf(ios::left | ios::fixed);
                outFileCor << std::fixed << m_modelCorrelations[libor][libor2] ;
            }
            outFileCor << std::endl;
        }

        outFileParams.close();
        outFileCalib.close();
        outFileCor.close();
    } // if 

    return;
}

//****************************************************************************//
//  end : printModel                                                          //
//****************************************************************************//

DRLIB_END_NAMESPACE


