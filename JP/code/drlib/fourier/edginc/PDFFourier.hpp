//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PDFFourier.hpp
//
//   Description : Implementation of PDFFourier for SVJJ etc.
//
//   Date        : 4 June 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_PDFFOURIER_HPP
#define EDR_PDFFOURIER_HPP

#include "edginc/VolBase.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolParam.hpp"
#include "edginc/PDFCalculator.hpp"

DRLIB_BEGIN_NAMESPACE

class FOURIER_DLL PDFFourier : public PDFCalculator, public FourierProduct, public StFourierProductLogRtn
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz){
        REGISTER(PDFFourier, clazz);
        SUPERCLASS(PDFCalculator);
    }

    static PDFCalculator* pdfCreator(const DateTime& valueDate, 
                                   const IModel* model,
                                   const Asset* asset, 
                                   const CVolProcessed* vol);

    virtual ~PDFFourier(){};

    /** constructor */
    PDFFourier(const DateTime&				valueDate,        // value date
               const CAsset*				asset,        // single factor
               const FourierEngine*         model,
               const FFTIntegrator1D*		fftintegrator);

    /** constructor */
    PDFFourier(const DateTime&				valueDate,        // value date
               const CAsset*				asset,        // single factor
               const FourierEngine*			model,
               const Integrator1D*			integrator);

	/** Calculate the cumulative probability at each strike.
        Each prob is estimated as a "small" call spread (i.e., by tweaking the strike)
        Returns false if if the difference between 2 consecutive probs is negative.
    */
    virtual void probabilities(const DoubleArray& strikes,
                               const DateTime&    maturity,
                               DoubleArray&       probs) const;

    virtual void probabilities(const CLatticeDouble&   strikes,
                               const DateTimeArray&    maturities,
                               CLatticeDouble&         probs) const
    {
        // do nothing for now
    }

    /** Calculate the "local" density at each strike.
        Each density is estimated as a "small" butterfly (using same strikes 
        as above together with the center strike).
        Returns false if any of the densities is negative
    */
    virtual void localDensity(const DoubleArray& strikes,
                              const DateTime&    maturity,
                              DoubleArray&       density) const
    {
        // do nothing for now
    }

    /** Calculate the "integrated" density at each strike.
        Each density is computed by taking the difference between 2 consecutive
        probabilities divided by the difference between the 2 consecutive strikes.
        false should be returned if any of the so-computed densities is negative;
    */        
    virtual void integratedDensity(const DoubleArray& strikes,
                                   const DateTime&    maturity,
                                   DoubleArray&       density) const
    {
        // do nothing for now
    }


    Function1DComplexArrayConstSP FFTIntegrand(
										    const DateTime& maturity,
										    const DoubleArray& strikes)const;

    Function1DDoubleArrayConstSP Integrand(
										    const DateTime& maturity,
										    const DoubleArray& strikes)const;

    static void setUseFFT(bool use){useFFT = use;}
private:
    PDFFourier();

	const DateTime				today; // $unregistered
	const CAsset*               asset;        // single factor $unregistered
    const FourierEngine*        model; // $unregistered
    const FFTIntegrator1D*		fftintegrator; // $unregistered
    const Integrator1D*			integrator; // $unregistered
	static bool					useFFT;
};

DRLIB_END_NAMESPACE
#endif
