//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MertonDependence.cpp
//
//   Description : Holds Merton dependence parameters and methods
//
//   Author      : Oliver Brockhaus
//
//   Date        : 21 Jan 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Dependence.hpp"
#include "edginc/VolMerton.hpp"
#include "edginc/Random.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/Integrator.hpp"
#include "edginc/MultiFactors.hpp"

DRLIB_BEGIN_NAMESPACE

///////////////////////////////////////////////////////////////////////////////////////////////
// DEPENDENCE MAKER MERTON
///////////////////////////////////////////////////////////////////////////////////////////////

class DependenceMakerMerton : public DependenceMaker{
public:
    static CClassConstSP const TYPE;

    class Support : virtual public DependenceMaker::ISupport {
    public:
        virtual void getMertonData(DateTimeArray& toDates,
                                   IMultiFactors& mAsset) const  = 0; 
    };

    // interface that the instrument must implement
    virtual DependenceSP createDependence(const DependenceMaker::ISupport* support) const;

    void getData( const MarketData* market, const IModel* model, const string& name );
    
    void checkDependenceMaker(const IObject* obj); 

    /** Overrides clone method to copy the map */
    IObject* clone() const {
        static const string routine = "DependenceMakerMerton::clone";        
        try {
            IObject* Copy = CObject::clone(); // call parent
            DependenceMakerMerton* copy = dynamic_cast<DependenceMakerMerton*>(Copy);
            if(!copy) {
                throw ModelException("Clone method failed");
            }
            copy->volMap = volMap;

            return copy;

        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    // Invoked when Class is 'loaded'
    static void load(CClassSP& clazz);

    DependenceMakerMerton();
    static IObject* defaultDependenceMakerMerton();

private:

//    CStringArraySP                   assetNames;
	bool							doCumulativesFlag; // $unregistered
    int                             nbAssets; // $unregistered
    mutable vector<VolMertonSP>     vol; // $unregistered
    map<string,VolMertonSP>         volMap; // $unregistered
};

CClassConstSP const DependenceMakerMerton::TYPE = CClass::registerClassLoadMethod(
    "DependenceMakerMerton", typeid(DependenceMakerMerton), DependenceMakerMerton::load);

///////////////////////////////////////////////////////////////////////////////////////////////
// MERTON: Poisson jump times, normal jumps plus diffusion with mean reverting vol
///////////////////////////////////////////////////////////////////////////////////////////////
class Merton : public Dependence
{
public:

    Merton();
    ~Merton();

    // validation
    virtual void validatePop2Object();

    // constructor
    Merton( DoubleMatrix&               randoms,
            const vector<VolMertonSP>&  vols,
            const DoubleMatrix&         correlation,
            const DateTimeArray&        dates,
            const TimeMetricConstSP&    timeMetric,
            const DoubleMatrix&         driverSpotVars,
            const bool                  useEffectiveCorrelation);

    // correlate a matrix of random numbers: [2*nbAssets+1][nbSteps]
    virtual void correlateSeries( DoubleMatrix& noise, int pathIdx );

    // generate a matrix of random numbers and from there a copula
    void generateCumulatives(     
            const IRandomSP&    rand,
            DoubleMatrix&       cumulatives );
    
    // calculate effective correlation
    double effectiveCorrelation( 
        const double    correlation, 
        const DateTime  &date1,
        const DateTime  &date2,
        const int       iAsset, 
        const int       jAsset ) const;

    DoubleMatrix getCorrelations(int index){
        throw ModelException("MertonDependence::getCorrelations",
                             "Not done yet");
    }
	double getCorrelationIJ(int i, int j, int index){
		throw ModelException("MertonDependence::getCorrelationIJ", "Not implemented yet");
	}
	int getNumAssets() {return nbAssets;}

private:
    DoubleMatrix                    correlation;
    int                             nbAssets;
    int                             nbSteps;
    DoubleMatrixArraySP             correlCoeffsTerm;
    PoissonSP                       poisson;
    CStringArraySP                  assetNames;
    vector<VolProcessedMertonSP>    volProcessed;
    TimeMetricConstSP               timeMetric;
    bool                            useEffectiveCorrelation;
    DoubleMatrix                    driverSpotVars;
};

// integrand used in effectiveCorrelation
class ProductOfVols : public Function1DDouble
{
public:
    ProductOfVols(
        const double            lower,
        const double            upper,
        VolProcessedMertonSP    vol1,
        VolProcessedMertonSP    vol2) 
        : 
    Function1DDouble(Range(ClosedBoundary(lower), ClosedBoundary(upper))), v1(vol1), v2(vol2) {}

    double operator()(double  x) const;

private:
    VolProcessedMertonSP    v1;
    VolProcessedMertonSP    v2;
};

///////////////////////////////////////////////////////////////////////////////////////////////
// MERTON IMPLEMENTATION
///////////////////////////////////////////////////////////////////////////////////////////////

/** Destructor */
Merton::~Merton(){ /*empty*/ }

/** Validation */
void Merton::validatePop2Object(){  /*empty*/ }

/** Default Constructor */
Merton::Merton(){}

/** Constructor */
Merton::Merton(
    DoubleMatrix&               randoms,
    const vector<VolMertonSP>&  vols,
    const DoubleMatrix&         correlation,
    const DateTimeArray&        dates,
    const TimeMetricConstSP&    timeMetric,
    const DoubleMatrix&         driverSpotVars,
    const bool                  useEffectiveCorrelation) :
    timeMetric(timeMetric),
    useEffectiveCorrelation(useEffectiveCorrelation),
    driverSpotVars(driverSpotVars)
{
    static const string method("Merton::Merton");

    try{
        // dimension checking ...
        nbAssets = correlation.numCols();
        nbSteps = randoms.numRows();
        if( correlation.numRows() != nbAssets ) {
            throw ModelException(method, 
                                 "Have correlation.numRows = " +
                                 Format::toString(correlation.numRows()) + 
                                 " and correlation.numCols = " +
                                 Format::toString(nbAssets) );
        }
        if( randoms.numCols() != 2*nbAssets+1 ) {
            throw ModelException(method, 
                                 "Have randoms.numCols = " +
                                 Format::toString(randoms.numCols()) + 
                                 ": should be 2*nbAssets+1 = " +
                                 Format::toString(2*nbAssets+1) );
        }
        /*
        if( vols.size() != nbAssets ) {
            throw ModelException(method, 
                                 "Have vols.size = " +
                                 Format::toString(vols.size()) + 
                                 ": should be nbAssets = " +
                                 Format::toString(nbAssets) );
        }
        */

        // need to set up bigger matrix for jumps...
        int iAsset, jAsset;
        DoubleMatrix correlationTotal = DoubleMatrix(2*nbAssets+1,2*nbAssets+1);

        for( iAsset=0; iAsset<nbAssets; iAsset++ ) {
            correlationTotal[iAsset][iAsset] = 1;
            correlationTotal[iAsset+nbAssets][iAsset+nbAssets] = 1;
            double beta_i = vols[iAsset]->getBeta();
            for( jAsset=iAsset+1; jAsset<nbAssets; jAsset++ ) {
                double beta_j = vols[jAsset]->getBeta();
                correlationTotal[iAsset][jAsset] 
                    = correlation[iAsset][jAsset];
                correlationTotal[jAsset][iAsset] 
                    = correlation[jAsset][iAsset];
                correlationTotal[iAsset+nbAssets][jAsset+nbAssets]
                    = beta_i * beta_j;
                correlationTotal[jAsset+nbAssets][iAsset+nbAssets]
                    = beta_i * beta_j;
            }
        }
        correlationTotal[2*nbAssets][2*nbAssets] = 1;

        // store vols
        DateTimeArray toDates(dates.size()-1);
        int iDate;
        for( iDate=0; iDate<dates.size()-1; iDate++ ) {
            toDates[iDate] = dates[iDate+1];
        }
        volProcessed.resize(nbAssets);

        for( iAsset=0; iAsset<nbAssets; iAsset++ ) {
            volProcessed[iAsset] = VolProcessedMertonSP(new VolProcessedMerton(vols[iAsset]));
            volProcessed[iAsset]->setupCumulatives(dates[0], toDates, driverSpotVars[iAsset]);
        }

        // for number of jumps simulation...
        double crashRate = volProcessed[0]->getCrashRate();
        poisson = PoissonSP( new Poisson(crashRate));

        // time as double...
        DoubleArray timeLine(dates.size()-1);
        for( iDate=0; iDate<dates.size()-1; iDate++ ) {
            timeLine[iDate] = timeMetric->yearFrac(dates[iDate],dates[iDate+1]);
        }

        // cumulative distribution of Poisson...
        poisson->numJumpsPreprocess(timeLine,100);

        // effective correlation

        // allocations for speed in correlateSeries
        correlCoeffsTerm = DoubleMatrixArraySP(new DoubleMatrixArray(nbSteps));
        if( useEffectiveCorrelation ) {

            double targetCorr;
            double x1, x2, dx, f, fmid, xmid, rtb;
            int jAsset, j;
            for( iDate=0; iDate<dates.size()-1; iDate++ ) {
                for( iAsset=0; iAsset<nbAssets; iAsset++ ) {
                    for( jAsset=iAsset+1; jAsset<nbAssets; jAsset++ ) {
                        // bisection
                        targetCorr  = correlation[iAsset][jAsset];
                        x1          = -1;
                        x2          = 1;
                        f           = effectiveCorrelation(x1,dates[iDate],dates[iDate+1],iAsset,jAsset) - targetCorr;
                        fmid        = effectiveCorrelation(x2,dates[iDate],dates[iDate+1],iAsset,jAsset) - targetCorr;
                        if( f*fmid >= 0 ) {
                            throw ModelException(method, 
                                                 "No Gauss correlation in [-1,1] gives target correlation " +
                                                 Format::toString(targetCorr) + " between asset " +
                                                 Format::toString(iAsset) + " and " +
                                                 Format::toString(jAsset) + " dates " +
                                                 dates[iDate].toString() + " and " +
                                                 dates[iDate+1].toString() );
                        }
                        rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
                        for(  j=1;j<=100;j++) {
                            fmid = effectiveCorrelation(xmid=rtb+(dx *= 0.5),dates[iDate],dates[iDate+1],iAsset,jAsset) - targetCorr;
                            if (fmid <= 0.0 ) rtb = xmid;
                            if (fabs(dx) < 1.0e-10 || fmid == 0.0) break;
                        }
                        correlationTotal[iAsset][jAsset] = rtb;
                        correlationTotal[jAsset][iAsset] = rtb;
                    } // iAsset
                } // jAsset
                CDoubleMatrixSP correlCoeffsTemp(
                    new CDoubleMatrix(correlationTotal.computeSquareRoot()));
                // upper triangular matrix
                (*correlCoeffsTerm)[iDate] = correlCoeffsTemp;
            } // iDate
        } // useEffectiveCorrelation
        else {
            DoubleMatrix correlCoeffsTemp=correlationTotal.computeSquareRoot();
            // MAR: Do we really need to copy the matrix at each point? - seems
            // very wasteful.
            for (iDate=0; iDate<dates.size()-1; iDate++) {
                 // upper triangular matrix
                (*correlCoeffsTerm)[iDate].reset(
                    new CDoubleMatrix(correlCoeffsTemp));
            }
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }

}

void Merton::correlateSeries( DoubleMatrix& noise, int pathIdx ) {

    static const string method("Merton::correlateSeries");

    try{
        // pathIdx not needed here
        // Start from last dimension as we are changing them, first one unchanged
        // Correlates asset and jump noise, number of jumps per step
        int iAsset, jAsset, iStep;
        for (iAsset = noise.numCols() - 1; iAsset > 0; iAsset --) {
            double* mainSeries = noise[iAsset];
            for (iStep = 0; iStep < nbSteps; iStep ++){
                mainSeries[iStep] *= (*(*correlCoeffsTerm)[iStep])[iAsset][iAsset];
            }
            for (jAsset = iAsset - 1; jAsset > -1; jAsset --) {
                double* otherSeries = noise[jAsset];
                for (iStep = 0; iStep < nbSteps; iStep ++) {
                    mainSeries[iStep] += otherSeries[iStep] * (*(*correlCoeffsTerm)[iStep])[iAsset][jAsset];
                }
            }
        }
    } catch (exception& e){
        throw ModelException(e,method);
    }
}

// generate a matrix of random numbers and from there a copula
void Merton::generateCumulatives(     
        const IRandomSP&    rand,
        DoubleMatrix&       cumulatives ) {

    static const string method("Merton::generateCumulatives");

    try{

        int iAsset,iStep;
        DoubleMatrix randoms = DoubleMatrix(2*nbAssets+1, nbSteps);

        // fetch random numbers
        if (true/*isCarefulRandoms*/) {
            // Generate by asset then by date from latest, so passing sample dates
            // drops a contiguous block of randoms "from the back" and can maintain
            // same randoms on same sample dates more easily
            for(iStep=nbSteps-1; iStep>=0; iStep--){
                for(iAsset=0; iAsset<2*nbAssets+1; iAsset++){
                    rand->fetch(1, &randoms[iAsset][iStep]);
                }
            }
        } else {
            for(iAsset=0; iAsset<2*nbAssets+1; iAsset++){
                rand->fetch(randoms.numRows(), randoms[iAsset]);
            }
        }

        // correlate (asset and jump) noises between assets
        correlateSeries(randoms, 0); // dummy pathIdx

        // generate number of jumps per period
        poisson->numJumps(randoms[2*nbAssets],nbSteps);

        // do Merton copula
        for(iAsset=0; iAsset<nbAssets; iAsset++) {
            (volProcessed[iAsset])->generateCumulatives(
                randoms[iAsset],
                randoms[iAsset+nbAssets],
                randoms[2*nbAssets],
                cumulatives[iAsset]);
        }

    } catch (exception& e){
        throw ModelException(e,method);
    }
}

// calculate effective correlation
double Merton::effectiveCorrelation( 
    const double    correlation, 
    const DateTime  &date1,
    const DateTime  &date2,
    const int       iAsset, 
    const int       jAsset ) const {

    static const string method("Merton::effectiveCorrelation");

    try{
        
        // variance
        double varGauss_i = volProcessed[iAsset]->CalcVar(date1,date2);
        double varGauss_j = volProcessed[jAsset]->CalcVar(date1,date2);

        // covariance
        double t1 = volProcessed[iAsset]->calcTradingTime(volProcessed[iAsset]->getBaseDate(),date1);
        double t2 = volProcessed[iAsset]->calcTradingTime(volProcessed[iAsset]->getBaseDate(),date2);
        ProductOfVols   productOfVols(t1,t2,volProcessed[iAsset],volProcessed[jAsset]);
        ClosedRomberg1D integrator(1.0e-6,20,5); // accuracy, nbIterMax, RichardsonPolyOrder
        double covGauss = correlation * integrator.integrate(productOfVols);

        // tau
        int iJump;
        double covJump 
            = volProcessed[iAsset]->getBeta() * volProcessed[iAsset]->getCrashSizeUncertainty()
            * volProcessed[jAsset]->getBeta() * volProcessed[jAsset]->getCrashSizeUncertainty();
        double varJump_i
            = volProcessed[iAsset]->getCrashSizeUncertainty() * volProcessed[iAsset]->getCrashSizeUncertainty();
        double varJump_j
            = volProcessed[jAsset]->getCrashSizeUncertainty() * volProcessed[jAsset]->getCrashSizeUncertainty();
        double tau = 0;
        double proba = exp((t1-t2)*volProcessed[iAsset]->getCrashRate());

        int nbJumps = 0;
        double quantile = 0;
        volProcessed[iAsset]->Quantile(date1,date2,.00001,100,&quantile,&nbJumps );
        for( iJump=0; iJump<nbJumps; iJump++ ) {
            double rho_i 
                = (covGauss + iJump*covJump) 
                / sqrt( (varGauss_i + iJump*varJump_i) * (varGauss_j + iJump*varJump_j) );
            tau += proba * asin(rho_i);
            proba *= volProcessed[iAsset]->getCrashRate() * (t2-t1) / (double)(iJump+1);
        }
        return sin(tau);

    } catch (exception& e){
        throw ModelException(e,method);
    }
}

double ProductOfVols::operator()(double  x) const { 
    double var1 = v1->CalcVar(x);
    double var2 = v2->CalcVar(x);
    return sqrt( var1*var2 );
}

///////////////////////////////////////////////////////////////////////////////////////////////
// MERTON DEPENDENCE MAKER IMPLEMENTATION
///////////////////////////////////////////////////////////////////////////////////////////////
DependenceSP DependenceMakerMerton::createDependence(const DependenceMaker::ISupport* support) const {

    static const string method("DependenceMakerMerton::createDependence");

    try{
        const DependenceMakerMerton::Support* supportDmm = 
            dynamic_cast<const DependenceMakerMerton::Support*>(support);
        if (!supportDmm) {
            throw ModelException(method, "DependenceMakerMerton is not a valid dependence maker.");
        }

        DateTimeArray dates;
        IMultiFactors* mAsset = NULL;
        supportDmm->getMertonData(dates, *mAsset);

        DateTimeArray toDates(dates.size()-1);
        int iDate;
        for( iDate=0; iDate<dates.size()-1; iDate++ ) {
            toDates[iDate] = dates[iDate+1];
        }

        int                   numAssets   = mAsset->NbAssets();
        CDoubleMatrixConstSP   correl(mAsset->factorsCorrelationMatrix());
        vol.resize(numAssets);
        CStringArraySP assetNames(new CStringArray(numAssets));

        // holds random numbers for process dW^n_i + #J_i * J^n_i: number of jumps #J same for all assets
        DoubleMatrix randoms = DoubleMatrix(2*numAssets+1, dates.size()-1);
        TimeMetricConstSP timeMetric;
        DoubleMatrix driverSpotVars = DoubleMatrix(numAssets, dates.size()-1);

        int iAsset;
        for( iAsset=0; iAsset<numAssets; iAsset++ ) {
            string name = mAsset->factorGetName(iAsset);
            map<string, VolMertonSP>::const_iterator iter = volMap.find(name);
            if (iter != volMap.end()){
                vol[iAsset] = iter->second;
            }
            else{
                throw ModelException(method, 
                                     "Vol for " + name + " not found." );
            }

            // Driver vols

            // override with an ATM vol
            ATMVolRequestSP requestATM(new ATMVolRequest());            
            requestATM->allowNegativeFwdVar(false);
            CVolProcessedSP volTemp(mAsset->factorGetProcessedVol(0,requestATM.get())); 
            CVolProcessedBSSP volTempBS(CVolProcessedBSSP::dynamicCast(volTemp));

            // get Total and Fwd variance
            DoubleArray totalVar(dates.size());

            volTempBS->CalcVar(dates[0], toDates, volTempBS->fromFirst, totalVar);
            int iStep;
            for( iStep=0; iStep<dates.size()-1; iStep++ ) {
                driverSpotVars[iAsset][iStep] = totalVar[iStep];
            }

            // time metric
            if( iAsset==0 ) {
                timeMetric = TimeMetricConstSP( volTempBS->GetTimeMetric().get() );
            }

        }

        return DependenceSP(new Merton( randoms, 
                                        vol, 
                                        *correl, 
                                        dates, 
                                        timeMetric, 
                                        driverSpotVars,
                                        true /*useEffectiveCorrelation*/ ));
    } catch (exception& e){
        throw ModelException(e,method);
    }
}

// creates vector of required vol objects if required for generateCumulatives
void DependenceMakerMerton::getData(const MarketData* market, const IModel* model, const string& name ) {

    static const string method("Merton::getData");
    CClassConstSP volClass(CClass::forName("VolMerton"));

    try{
        MarketObjectSP data(market->GetData(name, volClass));
        if (data.get()){
            map<string, VolMertonSP>::const_iterator iter = volMap.find(name);
            if (iter == volMap.end()){
                VolMertonSP temp(VolMertonSP::dynamicCast(data));
                volMap[name] = VolMertonSP(copy(temp.get()));
                volMap[name]->getMarket(model, market);
            }
        }
    } catch (exception& e){
        throw ModelException(e,method);
    }
}

void DependenceMakerMerton::checkDependenceMaker(const IObject* obj) {}

// Invoked when Class is 'loaded'
void DependenceMakerMerton::load(CClassSP& clazz){
    // clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(DependenceMakerMerton, clazz);
    SUPERCLASS(DependenceMaker);
    EMPTY_SHELL_METHOD(defaultDependenceMakerMerton);
}

IObject* DependenceMakerMerton::defaultDependenceMakerMerton(){
    return new DependenceMakerMerton();
}

DependenceMakerMerton::DependenceMakerMerton():
DependenceMaker(TYPE){};

// external symbol to allow class to be forced to be linked in
bool MertonDependenceLoad(){
    return (DependenceMakerMerton::TYPE != 0 && DependenceMaker::TYPE != 0);
}

DRLIB_END_NAMESPACE

