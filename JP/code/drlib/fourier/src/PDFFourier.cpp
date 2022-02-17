//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PDFFourier.cpp
//
//   Description : Implementation of PDFFourier for SVJJ etc.
//
//   Date        : 4 June 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FourierEngine.hpp"
#include "edginc/PDFFourier.hpp"
#include "edginc/PDFCalculatorMaker.hpp"

DRLIB_BEGIN_NAMESPACE

/* local helper */
/** Algorithm that forces the implied probabilities to be strictly decreasing */
void massageDistribution2(
    const CSliceDouble& sliceStrikes,  // the strikes
    CSliceDouble&       sliceDist,     // the distributions
    bool                midPercentile  /* whether to start from the midStrike or
                                          the midProbability */ ) {
    static const string method = "CDFMapping::massageDistribution2";

    try {
        int nStrikes = sliceDist.size();

        // Start from left to right and Max, Min each point to be in
        // the interval [0.0, 1.0]
        int iStrike;
        for(iStrike = 0; iStrike < nStrikes; iStrike++) {
            sliceDist[iStrike] = Maths::min(Maths::max(
                                                sliceDist[iStrike],0.0),1.0);
        }

        // Now we have a distribution that is in [0.0, 1.0] but it may
        // not be decreasing We will start massaging from left to
        // middle and from right to middle If the two points we end up
        // are inconsistent then we can't do it NOTE: middle is either
        // the midStrike or the midProb
        int position;
        if(midPercentile) {
            // We want to pick the strike that is nearer to 50% probaility
            double deviation = 1.1; // that is safe
            double target = 0.5;    // that is our target probability
            for(iStrike = 0, position = 0; iStrike < nStrikes; iStrike++) {
                if(fabs(sliceDist[iStrike] - target) < fabs(deviation)) {
                    position = iStrike;
                    deviation = fabs(sliceDist[iStrike] - target);
                }
            }
        } else {
            // We pick the middle strike
            position = nStrikes / 2;
        }

        // From left to middle
        for(iStrike = 1 ; iStrike < position; iStrike++) {
            sliceDist[iStrike] = Maths::min(sliceDist[iStrike-1], sliceDist[iStrike]);
        }

        // From right to middle
        for(iStrike = nStrikes - 2 ; iStrike >= position; iStrike--) {
            sliceDist[iStrike] = Maths::max(sliceDist[iStrike], sliceDist[iStrike + 1]);
        }

        // if in the middle the two massages are inconsistent fail
        string message = "Fixing distribution from left is inconsistent with fixing it"
                         " from right around strike " + Format::toString(sliceStrikes[position]);
        if(position > 0 && Maths::isNegative(sliceDist[position-1] - sliceDist[position])) {
            throw ModelException(message);
        }

        if(position < nStrikes - 1 && Maths::isNegative(sliceDist[position] - sliceDist[position+1])) {
            throw ModelException(message);
        }

        // So now we have a distribution that is in [0.0, 1.0] and is decreasing
        // But we need to make it strictly decreasing

        // We test if we are dealing with a very pathological case
        if(Maths::equals(sliceDist[0], sliceDist[nStrikes-1])) {
            throw ModelException("Lowest and highest points have equal mass");
        } else if(Maths::equals(sliceDist[0], 0.0)) {
            throw ModelException("Distribution at lowest point is 0.0");
        } else if(Maths::equals(sliceDist[nStrikes-1], 1.0)) {
            throw ModelException("Distribution at highest point is 1.0");
        } else if(Maths::isPositive(sliceDist[nStrikes-1]- sliceDist[0])) {
            throw ModelException("Highest point has more mass than lowest");
        }

        // We can now assume that we are not dealing with a totally paranoid distribution.
        // We make it strictly decreasing
        for(iStrike = 0; iStrike < nStrikes-1; iStrike++) {

            // for each strike: if the next one has equal distribution with the current,
            // move rightwards until we find the next strike where it is strictly less
            if(Maths::equals(sliceDist[iStrike], sliceDist[iStrike+1])) {
                int j = iStrike+1;
                while(j < nStrikes - 1 && Maths::equals(sliceDist[iStrike], sliceDist[j]) ) {
                    j+=1;
                }

                // 2 cases: i) such a strike exists or ii) it does not exist but there must be
                //          a strike leftwards that has strictly greater probability
                double slope;
                if(!Maths::equals(sliceDist[iStrike], sliceDist[j])) {
                    // we found the next decreasing point
                    slope = (sliceDist[iStrike]    - sliceDist[j]) /    // from current point
                            (sliceStrikes[iStrike] - sliceStrikes[j]);  // to next decreasing point
                    int i;
                    for(i = iStrike + 1; i <= j; i++){
                        // We replace the in between strikes by a straight line
                        sliceDist[i] = sliceDist[iStrike] +
                                       slope * (sliceStrikes[i] - sliceStrikes[iStrike]);
                    }
                } else {
                    // We ran out of strikes on the right, but there must exist a strike on the left
                    // (otherwise would have failed in the paranoid tests)
                    slope = (sliceDist[iStrike-1]    - sliceDist[j]) /   // from previous point
                            (sliceStrikes[iStrike-1] - sliceStrikes[j]); // to the end
                    int i;
                    for(i = iStrike; i <= j; i++){
                        // We replace the current and subsequent strikes
                        sliceDist[i] = sliceDist[iStrike-1] +
                                       slope * (sliceStrikes[i] - sliceStrikes[iStrike-1]);
                    }
                }
            }
        }
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}
/* end local helper */

/** Double integrand - Default integrator case */
template <class Process, class Product>
class VanillaImpliedGridDoubleIntegrand: public Function1DDouble {
public:
    VanillaImpliedGridDoubleIntegrand(const Process&   process,
                    const Product&   product,
                    double           strike,
                    double           fwd,
                    const DateTime&  maturity,
                    int              j):    // 0 or 1
    Function1DDouble(Range(OpenBoundary(0.0), Infinity(Infinity::Plus))),
    logMoneyness(log(strike/fwd)),
    matdate(maturity),
    j(j),
    process(process),
    product(product){}

    virtual double operator()(double  u) const {    // u == frequency
        Complex  z(0.0, u); // = i * u
        Complex  Laplace = exp(process.scalelessCumulant(product, z + j, matdate) - z * logMoneyness) / z;
        return Laplace.real();
    }

private:
    double logMoneyness;
    const DateTime& matdate;
    double j;
    const Process& process;
    const Product& product;
};
/** Double integrand - FFT integrator case */
template <class Process, class Product>
class VanillaImpliedGridDoubleIntegrandFFT: public Function1DComplex {
public:
    VanillaImpliedGridDoubleIntegrandFFT(const Process&   process,
                       const Product&   product,
                       const DateTime&  maturity,
                       int              j):    // 0 or 1
    // Infinite range by default
    matdate(maturity),
    j(j),
    process(process),
    product(product){}

    virtual Complex operator()(double  u) const {    // u == frequency
        if (Maths::isZero(u)){
            throw ModelException("VanillaImpliedGridDoubleIntegrandFFT::operator(double)",
                                 "Zero frequency is not supported yet");
        }
        Complex  z(0.0, u); // = i * u
        Complex  Laplace = exp(process.scalelessCumulant(product, z + j, matdate)) / z;
        return Laplace;
    }

private:
    const DateTime& matdate;
    double j;
    const Process& process;
    const Product& product;
};

CClassConstSP const PDFFourier::TYPE =
        CClass::registerClassLoadMethod("PDFFourier", typeid(PDFFourier), load);

PDFFourier::PDFFourier(const DateTime&             valueDate,        // value date
                       const CAsset*               asset,        // single factor
                       const FourierEngine*         model,
                       const FFTIntegrator1D*      fftintegrator)
                       : PDFCalculator(TYPE),
                       FourierProduct(asset, valueDate, 0, 0),
                       model(model),
					   asset(asset),
                       fftintegrator(fftintegrator)
					   {
                            setUseFFT(true);
                            FourierProcess& process = const_cast<FourierProcess&>(this->model->getProcess());
                            process.validate(this);
                       }

PDFFourier::PDFFourier(const DateTime&             valueDate,        // value date
                       const CAsset*               asset,        // single factor
                       const FourierEngine*         model,
                       const Integrator1D*			integrator)
                       : PDFCalculator(TYPE),
                       FourierProduct(asset, valueDate, 0, 0),
                       model(model),
					   asset(asset),
                       integrator(integrator)
					   {
                            setUseFFT(false);
                            FourierProcess& process = const_cast<FourierProcess&>(this->model->getProcess());
                            process.validate(this);
                       }

/** create integrand */
Function1DComplexArrayConstSP PDFFourier::FFTIntegrand(
										const DateTime& maturity,
										const DoubleArray& strikesinp)const{
	static const string method = "PDFFourier::Integrand";

    try{
        int nbStrikes = strikesinp.size();
        Function1DComplexArraySP functions(new Function1DComplexArray(0));
        functions->reserve(nbStrikes);

        DoubleArray fwd(1);
		fwd[0]=asset->fwdValue(maturity);


        const StFourierProductLogRtn& thisProd = *this;
        const StFourierProcessLogRtn* thisProc = dynamic_cast<const StFourierProcessLogRtn*>(&model->getProcess());
        if(!thisProc) {
            throw ModelException(method, "Process does not support StFourierProcessLogRtn interface.");
        }

        for(int iStrike = 0; iStrike < nbStrikes; ++iStrike){
            double strike = strikesinp[iStrike];
            functions->push_back(Function1DComplexSP(
                new VanillaImpliedGridDoubleIntegrandFFT<StFourierProcessLogRtn, StFourierProductLogRtn>
                                                                    (*thisProc,
                                                                     thisProd,
                                                                     maturity,
                                                                     0)));
        }
        return functions;
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

Function1DDoubleArrayConstSP PDFFourier::Integrand(
										const DateTime& maturity,
										const DoubleArray& strikesinp)const{
	static const string method = "PDFFourier::Integrand";

    try{
        int nbStrikes = strikesinp.size();
        Function1DDoubleArraySP functions(new Function1DDoubleArray(0));
        functions->reserve(nbStrikes);

        DoubleArray fwd(1);
		fwd[0]=asset->fwdValue(maturity);


        const StFourierProductLogRtn& thisProd = *this;
        const StFourierProcessLogRtn* thisProc = dynamic_cast<const StFourierProcessLogRtn*>(&model->getProcess());
        if(!thisProc) {
            throw ModelException(method, "Process does not support StFourierProcessLogRtn interface.");
        }

        for(int iStrike = 0; iStrike < nbStrikes; ++iStrike){
            double strike = strikesinp[iStrike];
            functions->push_back(Function1DDoubleSP(
                new VanillaImpliedGridDoubleIntegrand<StFourierProcessLogRtn, StFourierProductLogRtn>
                                                                    (*thisProc,
                                                                     thisProd,
																	 strike,
																	 fwd[0],
                                                                     maturity,
                                                                     0)));
        }
        return functions;
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}


/** Calculate the cumulative probability at each strike.
    prob is estimated using Fourier engine */
void PDFFourier::probabilities(const DoubleArray& strikesinp,
                           const DateTime&    maturity,
                           DoubleArray&       probs) const
{
    static const string method = "PDFFourier::probabilities";
    try {

		if(useFFT){
			int nbStrikes = strikesinp.size();

			DoubleArray fwd(1);
			fwd[0]=asset->fwdValue(maturity);

			DoubleArray outputCDF(nbStrikes);

			Function1DComplexArrayConstSP myfuncs = FFTIntegrand(maturity,strikesinp);

			FourierProductFFTIntegrator1D::IntegralArray integrals(nbStrikes);
			int iIntegral;
			for(iIntegral = 0; iIntegral < nbStrikes; iIntegral++) {
				try {
					integrals[iIntegral] = fftintegrator->integrate(*(*myfuncs)[iIntegral]);
				} catch(exception&) {};
			}

			iIntegral = 0;
			int iStrike;
			for (iStrike=0; iStrike < nbStrikes; ++iStrike){
				double strike = strikesinp[iStrike];
				double logMoneyness = log(strike / fwd[0]);
				FourierProductFFTIntegrator1D::Integral integral0;
				integral0 = integrals[iIntegral++];
				outputCDF[iStrike]=integral0->getValue(logMoneyness) / Maths::PI + 0.5;
			}

			DoubleArray strikesCopy(strikesinp); // to avoid const problems
			CSliceDouble sliceStrikes (&strikesCopy[0],strikesinp.size());

			DoubleArray outputCDFCopy(outputCDF);
			CSliceDouble sliceoutputCDF (&outputCDFCopy[0],outputCDF.size());

			try{
				massageDistribution2(sliceStrikes, sliceoutputCDF, true);
			}catch(exception& ){};

			for (iStrike=0; iStrike < nbStrikes; ++iStrike){
				probs[iStrike]=sliceoutputCDF[iStrike];
			}
		}else if(!useFFT){

			int nbStrikes = strikesinp.size();

			DoubleArray fwd(1);
			fwd[0]=asset->fwdValue(maturity);

			DoubleArray outputCDF(nbStrikes);

			Function1DDoubleArrayConstSP myfuncs = Integrand(
															maturity,
															strikesinp);
			FourierProductIntegrator1D::IntegralArray integrals(nbStrikes);
			int iIntegral;
			for(iIntegral = 0; iIntegral < nbStrikes; iIntegral++) {
				try {
					integrals[iIntegral] = integrator->integrate(*(*myfuncs)[iIntegral]);
				} catch(exception& ) {};
			}

			iIntegral = 0;
			int iStrike;
			for (iStrike=0; iStrike < nbStrikes; ++iStrike){
				double integral0 = integrals[iIntegral++];
				outputCDF[iStrike]=integral0 / Maths::PI + 0.5;
			}
			DoubleArray strikesCopy(strikesinp); // to avoid const problems
			CSliceDouble sliceStrikes (&strikesCopy[0],strikesinp.size());

			DoubleArray outputCDFCopy(outputCDF);
			CSliceDouble sliceoutputCDF (&outputCDFCopy[0],outputCDF.size());

			try{
				massageDistribution2(sliceStrikes, sliceoutputCDF, true);
			}catch(exception& ){};

			for (iStrike=0; iStrike < nbStrikes; ++iStrike){
				probs[iStrike]=sliceoutputCDF[iStrike];
			}
		}
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}
// pdfCreator
PDFCalculator* PDFFourier::pdfCreator(const DateTime& valueDate,
                                   const IModel* model,
                                   const Asset* asset,
                                   const CVolProcessed* vol)
{
    const FourierEngine* fourier = dynamic_cast<const FourierEngine*>(model);
    if(!fourier)
        throw ModelException("pdfCreator", "invalid model to create PDFFourier");

    if (!asset || !model)
        throw ModelException("pdfCreator", "failed to create pdf: no asset or model inputs");

    if (useFFT)
        return new PDFFourier(valueDate,
				            asset,
				            fourier,
				            fourier->fftIntegrator.get());
    else
        return new PDFFourier(valueDate,
				            asset,
				            fourier,
				            fourier->integrator.get());

}

bool PDFFourier::useFFT = false;

/* external symbol to allow class to be forced to be linked in */
bool PDFFourierLoad(){
    PDFCalculatorMaker::setCreator(PDFFourier::pdfCreator, "PDFFourier");
    return (PDFFourier::TYPE != 0);
}

DRLIB_END_NAMESPACE


