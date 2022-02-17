//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CountingTree.cpp
//
//   Author      : Charles Morcom
//
//   Date        : June 26, 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CountingTree.hpp"
#include <math.h>

#define MESH_SCALE 1.5
#define SQRT_MESH_SCALE 1.224744871391589

/* d*->d' transition types */
#define DTRANS_DU 0
#define DTRANS_D0 1
#define DTRANS_DD 2
#define DTRANS_JU 3
#define DTRANS_JD 4

/* Node iteration types */
#define INNER_ITER 0
#define MID_ITER 1
#define OUTER_ITER 2

DRLIB_BEGIN_NAMESPACE

CClassConstSP const CountingTree::TYPE = CClass::registerClassLoadMethod(
    "CountingTree", typeid(CountingTree), load);

void CountingTree::load(CClassSP& clazz){
    REGISTER(CountingTree, clazz);
    SUPERCLASS(FDModel);
}


CountingTree::CountingTree(const CClassConstSP &type) : FDModel(type), 
        F(0), xForFStart(), dfSpotVolatility(), N(0), 
        localJumpHazard(), localJumpMeanSize(), localJumpWidth(), 
        dfCorrelation(), 
        hasCState(), cMax(), cTailProbCutOff(1.0),
        T(0), length(), lengthD(), fForX(), sqrtDFCorrelation(), D(), DInverse(),
        marketFactorValues(), xIsOU(), xMeanReversion(), xCEVParameter(), hDriftAdjustment(),
        transitionHazardFactorIndex(), hasLocalJumps(false), localJumpFactorIndex(),
        localJumpUProb(), localJumpDProb(), localJumpNJProb(), localJumpUFactor(),
        localJumpDFactor(), range(), lastProbabilities(),
        cTransitionProbs(), dTransitionProbs(), shift(), localJumpUDestination(),
        localJumpDDestination() {};

CountingTree::~CountingTree()
{}

void CountingTree::addMarketFactor(bool addF, double mr, bool isCEV, double beta, int f, int n) {
    const static char* method = "CountingTree::addMarketFactor";
    if (mr<0) throw ModelException(method, "Mean-reversion parameters may not be negative.");
    
    fForX.push_back(f);
    if (addF) {
        F++;
        xForFStart.push_back(n);
    }
    N++;

    localJumpHazard.push_back(0.0);
    localJumpMeanSize.push_back(0.0);
    localJumpWidth.push_back(0.0);

    hasCState.push_back(false);
    cMax.push_back(0);

    xIsOU.push_back(!isCEV);
    xMeanReversion.push_back(mr);
    xCEVParameter.push_back(beta);

    transitionHazardFactorIndex.push_back(-1); 
}

/** Add a CEV market factor: the 'x' variable modeled is x = y^(1-beta)/(1-beta),
        where y is CEV with parameter beta, given mean-reversion, and time-dependent
        mean-reversion level determined during calibration. */
void CountingTree::addCEVMarketFactor(double meanReversionRate, double cevBetaParameter) {
    addMarketFactor(true, meanReversionRate, true, cevBetaParameter, F, N);
}

    /** Add a single dimension mapped OU market factor (usu. 2Q) */
void CountingTree::addOUMarketFactor(double meanReversionRate) {
    addMarketFactor(true, meanReversionRate, false, -1.0, F, N);
}
    /** Add a multi-dimensional mapped OU market factor */
void CountingTree::addOUMarketFactor(const vector<double>& meanReversionRates) {
    int f = F;
    int n = N;
    for (int i=0; i<(int)meanReversionRates.size(); n++) {
        this->addMarketFactor((i==0 ? true : false), meanReversionRates[n], false, -1.0, f, n);
    }
}

/** Adds a non-trivial counting-state to the nth dimension driver whose transitions
    depend on the value of market-factor f. The states can vary from 0 (intially) up
    to maxCState. */
void CountingTree::setCountingState(int n, int f, int maxCState) {
    static const char* method = "CountingTree::setCountingState";
    if (n>=N) throw ModelException(method, "Illegal attempt to add c-state to a (c,d) dimension which does not exist (n>=N)!");
    if (f>=F) throw ModelException(method, "Illegal attempt to use market factor which doesn't exist (f>=F)!");

    if (maxCState<0) throw ModelException(method, "maxCState may not be negative");
    if (maxCState==0) return; // this is not really a counting state at all.

    hasCState[n] = true;
    cMax[n] = maxCState;
    transitionHazardFactorIndex[n] = f;
}

/** Adds local jumps to market factor f (where f has a multi-dimensional driver,
    the jumps are added to the _first_ driver). */
void CountingTree::setLocalJumps(int f, double jumpHazard, double jumpMeanSize, double jumpWidth) {
    const static char* method = "CountingTree::setLocalJumps";
    if (f>=F) throw ModelException(method, "Attempt to add local jumps to a market factor which does not exist (f>=F).");
    if (jumpHazard<0) throw ModelException(method, "The local jump hazard may not be negative");
    if (jumpWidth<1.0) throw ModelException(method, "The jump width must be>=1.0 (i.e. the U jump must be bigger than the D jump).");
    

    if (jumpHazard==0.0) return; // no jumps after all
    if (jumpWidth==1.0 && jumpMeanSize==1.0) return; // no jumps after all

    int idx = xForFStart[f];
    localJumpHazard[idx] = jumpHazard;
    localJumpMeanSize[idx] = jumpMeanSize;
    localJumpWidth[idx] = jumpWidth;
}

/** Set the local diffusion correlations for the drivers (x). These are constant
    for all time. */
void CountingTree::setDriverCorrelation(const DoubleMatrix& corr) {
    const static char* method = "CountingTree::setDriverCorrelation";
    if (corr.numRows()!=corr.numCols()) throw ModelException(method, "The driver correlation matrix must be square!");
    dfCorrelation = corr;
}



/** Set the (time-dependent) driver diffusion volatilities. vol[t][n] is the
    volatility of driver n from t to (t+1). */
void CountingTree::setFactorVolatility(const DoubleMatrix& vol) {
    const static char* method = "CountingTree::setFactorVolatility";
    vol.checkNonNegative();
    dfSpotVolatility = vol;
}

/** Checks to see that market factors, drivers, counting-states, local jumps
    and dependencies are a AOK, once they have all been added, and the time-line has been
    set-up. */
void CountingTree::checkFactorSetUpOK() {

    const static char* method = "CountingTree::checkFactorSetUpOK";

    /* Dimensions, factor relations and ordering. */
    if (F<=0) throw ModelException(method, "No market factors (F<=0)!");
    if (N<F) throw ModelException(method, "There must be at least as many drivers as market factors (N>=F)");
    if (xForFStart.size()!=F) throw ModelException(method, "xForFStart should have size F.");
    if (xForFStart[0]<=0) throw ModelException(method, "xForFStart[0] should be 0: it is not.");
    for (int f=1; f<F; f++) {
        int n = xForFStart[f];
        if (n<=xForFStart[f-1])
            throw ModelException(method, "elements of xForFStart must be strictly increasing");
        if (n>=N) throw ModelException(method, "All elements of xForFStart must point to valid drivers.");
    }
    if (fForX.size()!=N) throw ModelException(method, "fForX must have size N.");
    if (localJumpHazard.size()!=N || localJumpMeanSize.size()!=N || localJumpWidth.size()!=N)
        throw ModelException(method, "local jump structures should all have size N.");
    if (hasCState.size()!=N) throw ModelException(method, "hasCState must have size N.");
    if (cMax.size()!=N) throw ModelException(method, "cMax must have size N.");
    if (xIsOU.size()!=N) throw ModelException(method, "xIsOU must have size N.");
    if (xMeanReversion.size()!=N) throw ModelException(method, "xMeanReversion must have size N.");
    if (xCEVParameter.size()!=N) throw ModelException(method, "xCEVParameter must have size N.");
    if (transitionHazardFactorIndex.size()!=N) throw ModelException(method, "transitionHazardFactorIndex must have size N."); 
    for (int n=0; n<N; n++) {
        if (fForX[n]<0 || fForX[n]>=F || (n>0 && fForX[n]<fForX[n-1]))
            throw ModelException(method, "fForX must index values in [0, ..., F-1], and must be increasing.");
        if (hasCState[n]) {
            if (cMax[n]<=0) throw ModelException(method, "If a factor has a c-state, cMax must be strictly positive.");
            if (transitionHazardFactorIndex[n]<0 || transitionHazardFactorIndex[n]>=F)
                throw ModelException(method, "If a factor has a c-state, its transition hazard factor index must have value in (0, ..., F-1)");
        } else {
            if (cMax[n]!=0) throw ModelException(method, "If a factor has no c-state, cMax should be zero!");
        }
    }

    /* Vols and correlations */
    if (dfSpotVolatility.numCols()<T) 
        throw ModelException(method, 
            "You must provide T spot volatilities for each factor.");
    if (dfSpotVolatility.numRows()!=N) 
        throw ModelException(method, 
            "You must provide spot vols for every factor.");
    dfSpotVolatility.checkPositive();
    if (dfCorrelation.numRows()!=dfCorrelation.numCols()) 
        throw ModelException(method, 
            "The correlation matrix must be square.");
    if (dfCorrelation.numRows()!=N) 
        throw ModelException(method, 
            "The correlation matrix must have the same size as there are driver factors.");
    sqrtDFCorrelation = dfCorrelation.computeSquareRoot();// this also checks positivity, etc.
}

/* Turns x and f vols and weights into D matrices and their inverses. */
void CountingTree::generateDMatrices() {
    for (int t=0; t<T; t++) {
        D[t] = sqrtDFCorrelation;
        double deltaN = sqrt(lengthD[t]) * SQRT_MESH_SCALE;
        // multiply each row by the factor vol
        for (int n=0; n<N; n++) {
            D[t].rowMultiply(n, deltaN * dfSpotVolatility[t][n] );
            DInverse[t] = D[t].computeInverse();
        }
        D[T] = D[T-1];
        DInverse[T] = DInverse[T-1];
    }
}

/** collect model initialisation data, set up timeline */
void CountingTree::initModel(){
    try{
		
        /* TIME LINE: SET UP BY CALLING PARENT initModel */
		FDModel::initModel(); // setup timeline
        T = timeLine->NumOfStep+1;

        /* CHECK THAT THE MARKET FACTORS ARE SET-UP CORRECTLY */
        checkFactorSetUpOK();

        /* ALLOCATE SIZES/MEMORY FOR PRECOMPUTED VARIABLES */
        length.resize(T+1);              // length[t] = timeLine->TradeYrFrac[t+1]
	    lengthD.resize(T+1);             // L^D asset-space time-steps
        D.resize(T+1);            // x = D.d matrices at each time-step
        DInverse.resize(T+1);            // DInverse.x = d matrices at each time-step
		
        /* SETUP TIME-LINE EXTRAS */
        double stepSize = 1.0/stepsPerYearFD;
        /* Note that tIdx refers to what happens between tIdx and tIdx+1. This is not the
           same as TradeYrFrac which refers to what happens up to
           tIdx - i.e. between tIdx-1 and tIdx. */
		for (int t=0; t<=T; t++) {
            if (t==T) {
                length[t] = length[t-1];
                lengthD[t] = lengthD[t-1];
            } else {
                double dt = timeLine->TradeYrFrac[t+1];
                length[t] = dt;
                double dtD = stepSize;
                lengthD[t] = dtD;
            }
		}

        /* COMPUTE CORRELATION ORTHOGONALIZATION AND LOCAL VOL MATRICES */
        
        generateDMatrices();

    }
    catch(exception& e){
        throw ModelException(e, "CountingTree::initModel()");
    }
}


double CountingTree::getPrice0( const TreeSlice& price ) const
{
    //!!! TEMPORARY FIX - REVIEW THIS
    return price.getCentre();
}

// type loading
bool CountingTreeLoad(){
    return CountingTree::TYPE != 0;
}

/* On forward induction, this calculates transition probabilities, state-probabilities
   and limits. This is done inductively, so assume know Inner and Outer limits for time tIdx,
   and tIdx state-probabilities. */
void CountingTree::computeTransitionAndLimits() {

    int t = range.tIdx;
    /*========================================================================
     * YOU ARE AT THE END OF THE TREE. NO NEXT SLICE, BUT IT'S TIME TO
     * RECOMPUTE GLOBAL LIMITS BASED ON MID AND INNER LIMITS, RATHER THAN
     * OUTER LIMITS
     *=======================================================================*/
    if (t==T) {
        resizeGlobalLimits();
        return;
    }

    double dt = length[t];
        double dtD = lengthD[t];
        /* If mu (drift error from discretization) is zero, this the 
           diffusion transition probability. If meshScale=3/2 and dt=dtD, then
           PU=PD=P0 = 1/3 */
        double baseDTransP = dt / (dtD * MESH_SCALE) ;
        /* You need to calculate D which is the matrix that expresses the
           relation between x (correlated drivers) and d (uncorrelated drivers): x = Dd */
        double deltaN = SQRT_MESH_SCALE * sqrt(dtD);

        /* You need the D mesh-transition matrix and its inverse */
        const DoubleMatrix& locD = D[t];
        const DoubleMatrix& locDPlus = D[t+1];
        const DoubleMatrix& locDInverse = DInverse[t];

    /* Check to see what Mid limits should be: make sure that
       outer limits are >= mid limits. 
       This is rather inefficient - you should be able to do this as
       you go though the value calculation, with a little more care. I have also
       drastically simplified things in the interests of getting things done 
       before I leave - NO TRANSIENT CONTAGION */
    if (range.hasCStates) {
        // the current implementation has NO TRANSIENT CONTAGION. This means
        // that d*=d, so the d component of the mid limits is always the same
        // as the inner limits. For the purposes of the forward induction,
        // for now I always assume that all possible c-transitions can exist
        // so that the mid upper limits are always going to be cMax.
        // Since this is always the case for the outer limits, there is
        // never any need to resize anything. This is not very efficient, and
        // WILL BREAK WITH TRANSIENT CONTAGION.
        range.cLMidLimits[t] = range.cLInnerLimits[t];
        range.cUMidLimits[t] = cMax;
        range.dLMidLimits[t] = range.dLInnerLimits[t];
        range.dUMidLimits[t] = range.dUInnerLimits[t];
    } else {
        // there are no (c,d)->(c',d*) transitions, so mid limits
        // are the same as current inner limits, outer limits
        // are OK as they are, and there's nothing much to do.
        range.cLMidLimits[t] = range.cLInnerLimits[t];
        range.cUMidLimits[t] = range.cUInnerLimits[t];
        range.dLMidLimits[t] = range.dLInnerLimits[t];
        range.dUMidLimits[t] = range.dUInnerLimits[t];
    }

    vector<int> c;
    vector<int> d;
    if (range.hasCStates) {

        /* Get next time state probability slice and initialize to zero */
        TreeSliceNDimCounting tempP(range);
        tempP.allocDim(N);
        /* Iterate through tIdx outer limits to compute mid limits
        and (c,d)->(c',d*) transitions: only if have cstates. */
        for (nodeIterBegin(OUTER_ITER, range, N, c, d); 
             !nodeIterEnd(OUTER_ITER, range, N, c, d);
             nodeIterNext(OUTER_ITER, range, N, c, d)) {

            int offset = sliceOffset(range, N, c, d);

            /* Calculate market factor values (not just inner, since
               may need for local jumps later) */
            
            for (int f=0; f<F; f++) {
                marketFactorValues[f][offset] = factorMapping(f, vMapping(f,d), c);
            }

            /* Only calculate c->c' transitions within inner limits */
            if (withinInnerLimits(range, N, c, d)) {
                // compute probabilities of transitions and mid limits
                /* You need to calculate c->c' transition probabilities at this node */
                for (int n=0; n<N; n++) {
                    if (hasCState[n] && c[n]<cMax[n]) {
                        // compute transition probability: NB no need
                        // if c=cMax, since no further transitions possible
                        double q = exp(- dt * marketFactorValues[transitionHazardFactorIndex[n]][offset]);
                        double p = 1-q;
                        double cumP = 0.0;
                        int K = (cMax[n]-c[n]);
                        double lastTrans = ::pow(q, K); // = q^K.p^0.C(K,0) = P(no transitions)
                        cTransitionProbs[n][0][offset] = lastTrans;
                        for (int dc=1; cumP<cTailProbCutOff && dc<range.cUMidLimits[t][n]-c[n]; dc++) {
                            lastTrans *= p*(K-dc+1)/(q*dc); // q^dc.p^(K-c').C(K,dc) = P(dc transitions)
                            cumP += lastTrans;
                            if (cumP>=cTailProbCutOff) {
                                // you want to allocate all remaining probability to this state
                                // transition
                                cTransitionProbs[n][dc][offset] = lastTrans + 1.0-cumP;
                                tempP.values[offset] += lastProbabilities->values[offset]*(lastTrans+1.0-cumP);
                            } else {
                                cTransitionProbs[n][dc][offset] = lastTrans;
                                tempP.values[offset] += lastProbabilities->values[offset]*lastTrans;
                            }

                        }
                    } else {
                        // no transition possible - all probability goes
                        // to next unchanged
                        tempP.values[offset] += lastProbabilities->values[offset];
                    }
                }
            } else {
                // no transitions - all probability stays where it is
                tempP.values[offset] += lastProbabilities->values[offset];
            }
        }

        /* Overwrite state probabilities with the c-transition-adjusted ones
           you have just calculated */
         lastProbabilities->takeFrom(tempP);
        
    } else {
        /* No C-States, so
           no calculations to do yet for state-probabilities: leave lastProbabilities alone */
    }


    /* Now iterate again and calculate diffusive transitions and shifts,
       tIdx+1 outer limits are just the outermost shifted destinations. */
    range.cLOuterLimits[t+1] = range.cLOuterLimits[t];
    range.cUOuterLimits[t+1] = range.cUOuterLimits[t];
    range.dLOuterLimits[t+1] = range.dLOuterLimits[t];
    range.dUOuterLimits[t+1] = range.dUOuterLimits[t];
    for (nodeIterBegin(OUTER_ITER, range, N, c, d); 
         !nodeIterEnd(OUTER_ITER, range, N, c, d);
         nodeIterNext(OUTER_ITER, range, N, c, d)) {

     int offset = sliceOffset(range, N, c, d);

        if (hasLocalJumps) {
            /* You need to compute the jump destination nodes - the jump probabilities
                are node-independent and have been calculated already. */
            for (int n=0; n<N; n++) {
                if (!localJumpHazard[n]) continue; // no jump = no calculations
                
                int fi = localJumpFactorIndex[n]; // which factor drives the jumps?
                /* Is is the value of the market factor at (c,d) */
                double baseF = marketFactorValues[fi][offset];
                /* This is probably not efficient: you know baseV already just from
                    d, since x = Dd, and v is a sum of a subset of the x_i. */
                double baseV = inverseFactorMapping(fi, baseF, c);
                
                double targetUV = inverseFactorMapping(fi, baseF*localJumpUFactor[n], c);
                double targetDV = inverseFactorMapping(fi, baseF*localJumpDFactor[n], c);

                double upD = singleDInverseMapping(fi, targetUV-baseV, n);
                double downD = singleDInverseMapping(fi, targetDV-baseV, n);

                // local jump destinations are rounded to the nearest integer
                // point in the d-grid.
                int juDest = (int)floor(upD+0.5);
                int jdDest = (int)floor(downD+0.5);
                localJumpUDestination[n][offset] = juDest;
                localJumpDDestination[n][offset] = jdDest;
                range.dLOuterLimits[t+1][n] = min(jdDest, range.dLOuterLimits[t+1][n]);
                range.dUOuterLimits[t+1][n] = max(juDest, range.dUOuterLimits[t+1][n]);
                
            }
        }

        /* Whether you have local jumps or not, you still need to calculate the
            DU, D0, DD probabilities and the shift offset slices. You should probably
            combine this loop with the one that computes the jump-transitions, by the
            way. */
        DoubleArray x_now = locD.mult( d ); // matrix multiplication
        DoubleArray x_plus = x_now;
        for (int n=0; n<N; n++) {
            /* First, you need to compute where (c,d) would evolve (assuming no
                jumps) by (t+1), in expectation, using the Euler approximation
                of the diffusion SDE. This will either be 
                E[x_+] = x + dt.(-mean reversion).x, for a simple OU x, or a more
                complex expression for a CEV process. */
            if (xIsOU[n]) {
                // if O-U, dx = -beta.x.dt + diffusion
                x_plus[n] *= (1-xMeanReversion[n]*dt);
            } else {
                double beta = xCEVParameter[n];
                double oneMinusBeta = 1.0 - beta;
                double betaOverOneMinusBeta = beta/oneMinusBeta;
                double x = x_now[n];
                x_plus[n] = x ; //+ 
                    //dt * (pow(x*oneMinusBeta,-betaOverOneMinusBeta)*xCEVMRLevel[n]
                    //     - x*oneMinusBeta*xMeanReversion[n]
                    //     - betaOverOneMinusBeta * sigmaSq / x);
            }
        }
        DoubleArray dDoublePlus = locDInverse.mult( x_plus); // matrix multiplication
        vector<int> dPlus;
        for (int n=0; n<N; n++) {
            // this is the nearest grid coordinate
            dPlus[n] = (int)floor(dDoublePlus[n] + 0.5);
            
            // shift[n] is an (n+1) dimensional slice - so you have
            // to recompute this for each dimension.
            // note that you can only do this here becuse you know that DPlusInv is
            // lower triangular, so that shiftIndex can only depend on
            // dPlus for dimensions up to and including n
            int so = sliceOffset(range, n+1, c, d);
            shift[n][so] = dPlus[n];
            /* Adjust next stage limits */
            range.dLOuterLimits[t+1][n] = min(dPlus[n]-1, range.dLOuterLimits[t+1][n]);
            range.dUOuterLimits[t+1][n] = max(dPlus[n]+1, range.dUOuterLimits[t+1][n]);

            // this is the bias you need in the trinomial random variable
            double mu = dDoublePlus[n] - dPlus[n];
            double P = mu/deltaN;
            double PSquared = P*P;
            double noJumpProb = localJumpNJProb[t][n];
            /* Note the use of so - each dimension of transition probabilities is a slice
                with only n dimensions, not N. */
            dTransitionProbs[n][DTRANS_DU][so] = 0.5 * noJumpProb * (baseDTransP + P + PSquared);
            dTransitionProbs[n][DTRANS_DD][so] = 0.5 * noJumpProb * (baseDTransP - P + PSquared);
            dTransitionProbs[n][DTRANS_D0][so] = noJumpProb * (1 - baseDTransP - PSquared);
        }
        
    }

    /*========================================================================
     * NOW YOU CAN TELL IF YOU NEED TO RESIZE THE GLOBAL LIMITS
     *=======================================================================*/
     resizeGlobalLimits();

    /*========================================================================
     * COMPUTE (t+1) STATE-PROBABILITIES
     * REMEMBER THAT YOU ALREADY CALCULATED THE EFFECTS OF C-TRANSITIONS,
     * SO YOU ONLY NEED TO DO THE DIFFUSION/LOCAL-JUMP CALCULATIONS
     *=======================================================================*/
     TreeSliceNDimCounting tempP(range);
     tempP.allocDim(N);
     for (nodeIterBegin(OUTER_ITER, range, N, c, d); 
         !nodeIterEnd(OUTER_ITER, range, N, c, d);
         nodeIterNext(OUTER_ITER, range, N, c, d)) {

             int offset = sliceOffset(range, N, c, d);
             
             
                vector<int> dPrime;
                vector<int> dTrans;
                double pDTrans;
                int offsetDTrans; // offset into (t+1) probability slice
				for (dPrimeIterBegin(N, c, d, offset, dTrans, dPrime, &pDTrans, &offsetDTrans);
                    dPrime[0]>=0;
                    dPrimeIterNext(N, c, d, offset, dTrans, dPrime, &pDTrans, &offsetDTrans)) {
    				if (pDTrans) tempP.values[offsetDTrans] += lastProbabilities->values[offset]*pDTrans;
			    }
                
                
     }
             /* Overwrite state probabilities with the d-transition-adjusted ones
           you have just calculated */
         lastProbabilities->takeFrom(tempP);

         /* NOW CLIP TO COMPUTE (t+1) INNER LIMITS */
    computeInnerLimits();
}

/* Computes probabilities etc. required for transition between tIdx and tIdx+1. 
   If computeLimits is true, also calculates MidLimits[tIdx] and OuterLimits[tIdx], assuming
   that InnerLimits[tIdx] have been done. Finally, clips OuterLimits[tIdx] to find
   InnerLimits[tIdx+1] */
void CountingTree::computeTransition() {
    try {
        int t = range.tIdx;
        // if at end of tree, nothing to do
		if (t==T) return;

        vector<int> c;
        vector<int> d;
        double dt = length[t];
        double dtD = lengthD[t];
        /* If mu (drift error from discretization) is zero, this the 
           diffusion transition probability. If meshScale=3/2 and dt=dtD, then
           PU=PD=P0 = 1/3 */
        double baseDTransP = dt / (dtD * MESH_SCALE) ;
        /* You need to calculate D which is the matrix that expresses the
           relation between x (correlated drivers) and d (uncorrelated drivers): x = Dd */
        double deltaN = SQRT_MESH_SCALE * sqrt(dtD);

        /* You need the D mesh-transition matrix and its inverse */
        const DoubleMatrix& locD = D[t];
        const DoubleMatrix& locDPlus = D[t+1];
        const DoubleMatrix& locDInverse = DInverse[t];

        // make sure you have calculated all the mapped factors that you need
        // to compute c-transition probabilities and local jump destinations
        // it is ugly to have to do this, but for local jumps you have to
        // calculate these in a separate loop in case. For now, I am (rather inefficiently)
        // calculating them anyway, even if I don't need them for sure.
        /* NOT NEEDED, I THINK: CAN COMPUTE JUMP DESTINATION WITHOUT KNOWING VALUES IN ADVANCE 
           SEE BELOW */
        for (nodeIterBegin(MID_ITER,range,N,c,d); 
             !nodeIterEnd(MID_ITER,range,N,c,d); 
             nodeIterNext(MID_ITER,range,N,c,d)) {
            int offset = sliceOffset(range,N,c,d);
            for (int f=0; f<F; f++) {
                marketFactorValues[f][offset] = factorMapping(f, vMapping(f,d), c);
            }
        }

        /* Iterate across MidLimits (since need dTransitionProbs for all possible (c', d*) nodes,
           as well as cTransitionProbabilities for (c, d) nodes. */
        for (nodeIterBegin(MID_ITER,range,N,c,d); 
             !nodeIterEnd(MID_ITER,range,N,c,d); 
             nodeIterNext(MID_ITER,range,N,c,d)) {

            // this is the full N-dimension slice offset - relative to _mid_ limits
            int offset = sliceOffset(range,N,c,d);

            if (range.hasCStates && withinInnerLimits(range,N,c,d)) {
                /* You need to calculate c->c' transition probabilities at this node */
                for (int n=0; n<N; n++) {
                    if (hasCState[n] && c[n]<cMax[n]) {
                        // compute transition probability: NB no need
                        // if c=cMax, since no further transitions possible
                        double q = exp(- dt * marketFactorValues[transitionHazardFactorIndex[n]][offset]);
                        double p = 1-q;
                        double cumP = 0.0;
                        int K = (cMax[n]-c[n]);
                        double lastTrans = ::pow(q, K); // = q^K.p^0.C(K,0) = P(no transitions)
                        cTransitionProbs[n][0][offset] = lastTrans;
                        for (int dc=1; cumP<cTailProbCutOff && dc<range.cUMidLimits[t][n]-c[n]; dc++) {
                            lastTrans *= p*(K-dc+1)/(q*dc); // q^dc.p^(K-c').C(K,dc) = P(dc transitions)
                            cumP += lastTrans;
                            if (cumP>=cTailProbCutOff) {
                                // you want to allocate all remaining probability to this state
                                // transition
                                cTransitionProbs[n][dc][offset] = lastTrans + 1.0-cumP;
                            } else {
                                cTransitionProbs[n][dc][offset] = lastTrans;
                            }
                        }
                    }
                }

                /* This would probably be the right moment to calculate dStar values
                   as well, if you are going to have transient contagion, unless you
                   do it on the fly in the DEV method: for now, not done. */
                /* TODO: INTRODUCE TRANSIENT CONTAGION */
            }

            /* All nodes within mid limits need to have d*->d' transition probabilities
               calculated (i.e. not just inner limits) 
               These are of the form pTransitionProbs[n][DTRANS_TYPE][offset],
               where DTRANS_TYPE is _DU, _D0, _DD for the diffusions and _JU and
               _JD for the local jumps. */
            if (hasLocalJumps) {
                /* You need to compute the jump destination nodes - the jump probabilities
                   are node-independent and have been calculated already. */
                for (int n=0; n<N; n++) {
                    if (!localJumpHazard[n]) continue; // no jump = no calculations
                    
                    int fi = localJumpFactorIndex[n]; // which factor drives the jumps?
                    /* Is is the value of the market factor at (c,d) */
                    double baseF = marketFactorValues[fi][offset];
                    /* This is probably not efficient: you know baseV already just from
                       d, since x = Dd, and v is a sum of a subset of the x_i. */
                    double baseV = inverseFactorMapping(fi, baseF, c);
                    
                    double targetUV = inverseFactorMapping(fi, baseF*localJumpUFactor[n], c);
                    double targetDV = inverseFactorMapping(fi, baseF*localJumpDFactor[n], c);

                    double upD = singleDInverseMapping(fi, targetUV-baseV, n);
                    double downD = singleDInverseMapping(fi, targetDV-baseV, n);

                    // local jump destinations are rounded to the nearest integer
                    // point in the d-grid.
                    int juDest = (int)floor(upD+0.5);
                    int jdDest = (int)floor(downD+0.5);
                    localJumpUDestination[n][offset] = juDest;
                    localJumpDDestination[n][offset] = jdDest;
                    /* NOTE: these don't worry about the (t+1) inner limits at
                       all. Whether or not one should actually use these destinations
                       in EV computations is entirely up to the (D)EV method. */

                    
                }
            }

            /* Whether you have local jumps or not, you still need to calculate the
               DU, D0, DD probabilities and the shift offset slices. You should probably
               combine this loop with the one that computes the jump-transitions, by the
               way. */
            DoubleArray x_now = locD.mult( d ); // matrix multiplication
            DoubleArray x_plus = x_now;
            for (int n=0; n<N; n++) {
                /* First, you need to compute where (c,d) would evolve (assuming no
                   jumps) by (t+1), in expectation, using the Euler approximation
                   of the diffusion SDE. This will either be 
                   E[x_+] = x + dt.(-mean reversion).x, for a simple OU x, or a more
                   complex expression for a CEV process. */
                if (xIsOU[n]) {
                    // if O-U, dx = -beta.x.dt + diffusion
                    x_plus[n] *= (1-xMeanReversion[n]*dt);
                } else {
                    double beta = xCEVParameter[n];
                    double oneMinusBeta = 1.0 - beta;
                    double betaOverOneMinusBeta = beta/oneMinusBeta;
                    double x = x_now[n];
                    x_plus[n] = x ; //+ 
                        //dt * (pow(x*oneMinusBeta,-betaOverOneMinusBeta)*xCEVMRLevel[n]
                        //     - x*oneMinusBeta*xMeanReversion[n]
                        //     - betaOverOneMinusBeta * sigmaSq / x);
                }
            }
            DoubleArray dDoublePlus = locDInverse.mult( x_plus); // matrix multiplication
            vector<int> dPlus;
            for (int n=0; n<N; n++) {
                // this is the nearest grid coordinate
                dPlus[n] = (int)floor(dDoublePlus[n] + 0.5);
                
                // shift[n] is an (n+1) dimensional slice - so you have
                // to recompute this for each dimension.
                // note that you can only do this here becuse you know that DPlusInv is
                // lower triangular, so that shiftIndex can only depend on
                // dPlus for dimensions up to and including n
                int so = sliceOffset(range, n+1, c, d);
                shift[n][so] = dPlus[n];

                // this is the bias you need in the trinomial random variable
                double mu = dDoublePlus[n] - dPlus[n];
                double P = mu/deltaN;
                double PSquared = P*P;
                double noJumpProb = localJumpNJProb[t][n];
                /* Note the use of so - each dimension of transition probabilities is a slice
                   with only n dimensions, not N. */
                dTransitionProbs[n][DTRANS_DU][so] = 0.5 * noJumpProb * (baseDTransP + P + PSquared);
                dTransitionProbs[n][DTRANS_DD][so] = 0.5 * noJumpProb * (baseDTransP - P + PSquared);
                dTransitionProbs[n][DTRANS_D0][so] = noJumpProb * (1 - baseDTransP - PSquared);
            }
     
        }

    } catch (exception& e) {
        throw ModelException(e, "CountingTree::computeTransition");
    }
}

/**(D)EV of a valueSlice to the current time point from the next. The operation is EV if the 
	   discountSlice pointer is null, else DEV by multiplying by the discountSlice values. */
void CountingTree::sliceEV(
    TreeSliceNDimCounting& valueSlice,      /**<(t+1) slice to (D)EV to (t) slice */
    TreeSliceNDimCounting* discountSlice    /**<(t) discount factor slice, if want DEV. If null, EV. */
        ) const {

	try {

        int t = range.tIdx;

		// if at end of tree, nothing to do
		if (t==T) return;

        int sliceDim = valueSlice.dim;
        double* sliceVal = valueSlice.values;
        TreeSliceNDimCounting tempSlice(range);
        tempSlice.allocDim(sliceDim);
        double* tempSliceVal = tempSlice.values;

		// need to think about how to make sure that valueSlice values are sane outside
		// where they were calculated at the last EV. How do we compute whether they are 
		// needed or not?

		// loop over all states in this (tIdx) time slice relevant to this sliceDim.
        // these are only the inner limits
		vector<int> d; // diffusion grid states
		vector<int> c; // counting states
		for (nodeIterBegin(INNER_ITER,range,sliceDim,c,d); 
             !nodeIterEnd(INNER_ITER,range,sliceDim,c,d); 
             nodeIterNext(INNER_ITER,range,sliceDim,c,d)) {

            // this is the index into the tIdx slice
            int offsetCD = sliceOffset(range, sliceDim, c, d);
            double* localTempSliceVal = tempSliceVal + offsetCD;
            localTempSliceVal[0] = 0.0;
            
            // you always want the discount from the current node, _not_
            // the node offset for cPrime, dStar.
            double disc = (discountSlice 
                ? (discountSlice->dim==sliceDim 
                        ? discountSlice->values[offsetCD] 
                        : discountSlice->val(c,d)) 
                : 1.0);

			// loop over all possible c-state changes from (c,d)->(c',d*)
            vector<int> cPrime; // possible c-state transitions
			double pCPrime; // probability of cPrime transition
			vector<int> dStar; // diffusion grid transition, given cPrime
			for (cPrimeDStarIterBegin(sliceDim, c, d, offsetCD, cPrime, dStar, &pCPrime);
                 cPrime[0]>=0;
                 cPrimeDStarIterNext(sliceDim, c, d, offsetCD, cPrime, dStar, &pCPrime)) {

                int offsetCPrimeDStar = tempSlice.offset(cPrime, dStar);

				// do all the possible diffusion transitions from (c',d*) to the
                // +1/0/-1/JU/JD nodes in each dimension
                double ev=0;
                vector<int> dPrime;
                vector<int> dTrans;
                double pDTrans;
                int offsetDTrans; // offset into (t+1) (cPrime, dPrime)
				for (dPrimeIterBegin(sliceDim, cPrime, dStar, offsetCPrimeDStar, dTrans, dPrime, &pDTrans, &offsetDTrans);
                    dPrime[0]>=0;
                    dPrimeIterNext(sliceDim, cPrime, dStar, offsetCPrimeDStar, dTrans, dPrime, &pDTrans, &offsetDTrans)) {
    				
                    ev += pDTrans * sliceVal[offsetDTrans];
			    }
                if (discountSlice) {
                    ev *= pCPrime*disc;
                } else {
                    ev *= pCPrime;
                }
                localTempSliceVal[0] += ev;
			}

        }
            /* Copy temp slice values back to valueSlice */
        /* WRITE THIS, CHARLES! */

	} catch (exception& e) {
		throw ModelException(e, "CountingTree::sliceEV");
	}


}

        void CountingTree::nodeIterBegin(int limitType, const TreeSliceNDimCounting::Range& range, 
                                int sliceDim, vector<int>& c, vector<int>& d) {

    int t = range.tIdx;
    switch (limitType) {
        case INNER_ITER:
            c = range.cLInnerLimits[t];
            d = range.dLInnerLimits[t];
            return;
        case MID_ITER:
            c = range.cLMidLimits[t];
            d = range.dLMidLimits[t];
            return;
        case OUTER_ITER:
            c = range.cLOuterLimits[t];
            d = range.dLOuterLimits[t];
            return;
        default:
            throw ModelException("Unrecognised node iteration type.");
    }
}
void CountingTree::nodeIterNext(int limitType, const TreeSliceNDimCounting::Range& range, 
                                int sliceDim, vector<int>& c, vector<int>& d) {
    int t = range.tIdx;
    const vector<int>* cLP;
    const vector<int>* cUP;
    const vector<int>* dLP;
    const vector<int>* dUP;
    switch (limitType) {
        case INNER_ITER:
            cLP = &range.cLInnerLimits[t];
            cUP = &range.cUInnerLimits[t];
            dLP = &range.dLInnerLimits[t];
            dUP = &range.dUInnerLimits[t];
            break;
        case MID_ITER:
            cLP = &range.cLMidLimits[t];
            cUP = &range.cUMidLimits[t];
            dLP = &range.dLMidLimits[t];
            dUP = &range.dUMidLimits[t];
            break;
        case OUTER_ITER:
            cLP = &range.cLOuterLimits[t];
            cUP = &range.cUOuterLimits[t];
            dLP = &range.dLOuterLimits[t];
            dUP = &range.dUOuterLimits[t];
            break;
        default:
            throw ModelException("Unrecognised node iteration type.");
    }
    const vector<int>& cL = *cLP;
    const vector<int>& cU = *cUP;
    const vector<int>& dL = *dLP;
    const vector<int>& dU = *dUP;

    int n=0;
    while (n<sliceDim) {
        if (d[n]<dU[n]) {
            // just increment this (innermost) dimension since still 
            // less than max
            d[n]++;
            return;
        } else if (range.hasCStates) {
            if (c[n]<cU[n]) {
                // increment c state and reset d state
                c[n]++;
                d[n] = dL[n];
                return;
            } else {
                // reset c and d counters and loop to increment higher dimension
                d[n] = dL[n];
                c[n] = cL[n];
                n++;
            }
        } else {
            // reset to lower limit and try to increment next one up
            d[n] = dL[n];
            n++;
        }
    }
    
    // if you are here, you have run out of (d,c)-states to increment,
    // so you are finished: signal this to nodeIterEnd()
    c[0]=-1;
}




bool CountingTree::withinInnerLimits(const TreeSliceNDimCounting::Range& range,
                                     int sliceDim, const vector<int>& c, const vector<int>& d) {
    int t = range.tIdx;
    for (int n=0; n<sliceDim; n++) {
        int dn = d[n];
        if (dn<range.dLInnerLimits[t][n] || dn>range.dUInnerLimits[t][n]) return false;
        if (range.cSize[n]>0) {
            int cn = c[n];
            if (cn<range.cLInnerLimits[t][n] || cn>range.cUInnerLimits[t][n]) return false;
        }
    }
    return true;
}

/*=========================================================================
* METHODS TO LOOP THROUGH POSSIBLE C-STATE TRANSITIONS AT A GIVEN NODE
*=======================================================================*/
/* The start point for this is the zero-change state, so cPrime=c and dStar=d.
   the probability of this is computed.*/
void CountingTree::cPrimeDStarIterBegin(
    int sliceDim,
    const vector<int>& c, const vector<int>& d, int cdOffset,
    vector<int>& cPrime, vector<int>& dStar, double* transProb
    ) const {

    cPrime = c; // base case is no transitions at all
    dStar = d; // no transitions means no dStar jumps

    // just remains to compute probability of transition
    if (!range.hasCStates) {
        // if no c-states, this is trivial
        // the probability is 1.0 for sure
        *transProb = 1.0;
    } else {
        double pLoc = 1.0;
        int n=0;
        while (n<sliceDim) {
            // the probability of being in the base state is only <1 if
            // there are possible transitions
            if (hasCState[n]) pLoc *= cTransitionProbs[n][0][cdOffset];
            n++;
        }
        *transProb = pLoc;
    }
    return;
}
/* Increments to next possible (c',d*) given (c,d) and computes transition
    probability. When gets to end, sets c'[0]=-1. */
void CountingTree::cPrimeDStarIterNext(
    int sliceDim, 
    const vector<int>& c, const vector<int>& d, int cdOffset,
    vector<int>& cPrime, vector<int>& dStar, double* transProb
    ) const {

    if (!range.hasCStates) {
        // if there are no c-transitions, then the only possible (c',d*)
        // is (c,d), so any iteration reaches the end since begin
        // was the only state.
        cPrime[0] = -1;
        return;
    }

    int n=0;
    bool exit = false; // loop until one cPrime entry changed
    const vector<int>& cPrimeULimits = range.cUMidLimits[range.tIdx];
    while (!exit && n<sliceDim) {
        if (hasCState[n]) {
            if (cPrime[n]<cPrimeULimits[n]) {
                // cPrime[n] still has headroom: increment
                // and signal that done
                cPrime[n]++;
                exit = true;
            } else {
                // got to max cPrime[n] - reset, and
                // loop to n+1
                cPrime[n] = c[n];
            }
        }
        n++;
    }

    // If exit was not set to true, _all_ cPrime[n]
    // are at the maximum, so you're done
    if (!exit) {
        cPrime[0] = -1;
        return;
    }

    // now you need to calculate the transient contagion state, dStar
    // and the probability of the transition, transProb
    double pLoc = 1.0;
    /* I have removed all transient contagion entirely, so d*=d always */
    dStar = d;
    for (n=0; n<sliceDim; n++) {
        // update transition probability: it's a slice
        // depending on the change in the c-state
        if (hasCState[n]) {
            pLoc *= cTransitionProbs[n][cPrime[n]-c[n]][cdOffset];
        }

        /* TRANSIENT CONTAGION REMOVED ENTIRELY: NEED TO THINK MORE */
        // update d* if this d-factor has transient contagion from 
        // this or another c-factor
        //int driverN = transientContagionDriver[n];
        //if (driverN>=0) {
            // this d-state has transient contagion from c-state
            // changes in factor driverN.
        //    int cStateChange = cPrime[driverN] - c[driverN];
        //    dStar[n] = 
        //        (cStateChange==0 
        //            ? d[n] // in ground state, so d*=d
        //            : dStarMapping[driverN][cStateChange][cdOffset]);
        //} else {
        //    // no transient contagion for this d-factor, so d*=d
        //    // and no extra probability computation needed
        //    dStar[n] = d[n];
        //}
    }
    *transProb = pLoc;
}



    /*=========================================================================
     * METHODS TO LOOP THROUGH POSSIBLE d*->d' TRANSITIONS AT A GIVEN NODE
     *=======================================================================*/
    /* Given (cPrime, dStar), find the starting point of the accessible set 
       of dPrime in the t+1 slice. Also compute the transition probability and
       the shiftedOffset into the t+1 slice. */
    void CountingTree::dPrimeIterBegin(
        int sliceDim,
        const vector<int>& cPrime, 
        const vector<int>& dStar, 
        int offset,
        vector<int>& dTransState, // have values DTRANS_DU etc.
        vector<int>& dPrime, // d-values in t+1 slice
        double* p, // probability of transition
        int* shiftedOffset // offset into t+1 value slice
        ) const {

            double pLoc = 1.0;
            for (int i=0; i<sliceDim; i++) {
                dTransState[i] = DTRANS_DU;
                // the ith shift is an (i+1) dimensional slice and DTRANS_DU
                // is always +1
                dPrime[i] = shift[i][sliceOffset(range,i+1, cPrime, dStar)] + 1;
                pLoc *= dTransitionProbs[i][DTRANS_DU][offset];
            }
            *p = pLoc;
            *shiftedOffset = sliceOffset(range,sliceDim, cPrime, dPrime);
    }
    void CountingTree::dPrimeIterNext(
        int sliceDim,
        const vector<int>& cPrime, 
        const vector<int>& dStar, 
        int offset, // offset into this (t) slice DO YOU NEED THIS?
        vector<int>& dTransState, // have values DTRANS_DU etc.
        vector<int>& dPrime, // d-values in (t+1) slice
        double* p, // probability of transition
        int* shiftedOffset // offset into (t+1) value slice
        ) const {
        
        double pLoc = 1.0;
        int n=0;
        bool exit = false;
        while (!exit && n<sliceDim) {
            if (localJumpHazard[n] && dTransState[n]<DTRANS_JD) {
                // dimension with local jumps: OK to increment
                dTransState[n]++;
                exit=true;
            } else if (dTransState[n]<DTRANS_DD) {
                // dimension with no local jumps: OK to increment
                dTransState[n]++;
                exit=true;
            } else {
                // need to reset this dimension and try to increment the next one
                dTransState[n] = DTRANS_DU;
            }
        }
        if (!exit) {
            dPrime[0] = -1; // you have reached the end
            return;
        }

        for (int i=0; i<sliceDim; i++) {
            // the ith shift is an (i+1) dimensional slice and DTRANS_DU
            // is always +1
            int tState = dTransState[i];
            int so = sliceOffset(range,i+1, cPrime, dStar);
            switch (tState) {
                case DTRANS_DU:
                    dPrime[i] = shift[i][so] + 1;
                    pLoc *= dTransitionProbs[i][tState][so];
                    break;
                case DTRANS_D0:
                    dPrime[i] = shift[i][so];
                    pLoc *= dTransitionProbs[i][tState][so];
                    break;
                case DTRANS_DD:
                    dPrime[i] = shift[i][so] - 1;
                    pLoc *= dTransitionProbs[i][tState][so];
                    break;
                case DTRANS_JU:
                    dPrime[i] = localJumpUDestination[i][so];
                    pLoc *= localJumpUProb[range.tIdx][i];
                    break;
                case DTRANS_JD:
                    dPrime[i] = localJumpDDestination[i][so];
                    pLoc *= localJumpDProb[range.tIdx][i];
                    break;
            }
            if (i==sliceDim-1) *shiftedOffset = so;
        }
        *p = pLoc;
        return;
    }

    /* Computes the slice value array offset, given a set of (c,d) coordinates
       and a slice dimension. c and d coordinates are interleaved into
       p = (d0, c0, d1, c1, ... ) and each are converted into [0, Si] values
       where Si is upperlimit - lowerlimit + 1, by subtracting the lower limits. 
       Note that slice array offsets are controlled by the _global_ limits, not
       those for this particular time-point. */
        int CountingTree::sliceOffset(const TreeSliceNDimCounting::Range& range,
            int sliceDim, const vector<int>& c, const vector<int>& d) {    
	    int os = d[0] - range.dLGlobalLimits[0];
        bool hasC = range.hasCStates;
        int i = 0;
	    while (i++ && i<sliceDim) {
		    os  = os*range.dSize[i-1] + d[i] - range.dLGlobalLimits[i];
            if (hasC && range.cSize[i-1]>1) {
                os = os*range.cSize[i-1] + c[i] - range.cLGlobalLimits[i];
            }
	    }
	    return os;
    }

    /* Maps a state (c,v=sum(x)) to the value of the fth factor. To get v, you should
       use vMapping(f, d). */
    double CountingTree::factorMapping(int f, double v, const vector<int>& c) const {
        // TODO: you need to write this and decide on the factor smile representation
        // structures.
        return v;
    }
    /* Returns the v = sum(x) value that corresponds to a given market factorVal
       and a given c-state. */
    double CountingTree::inverseFactorMapping(int f, double factorVal, const vector<int>& c) const {
        // TODO: write this.
        return factorVal;
    }
    /* Given v_f = sum(x), what single d-value in dimension dDimension maps to this? 
       This is used to compute the change in d required to produce a required
       change in v. */
    double CountingTree::singleDInverseMapping(int f, double v, int dDimension) const {
        double dsum = 0.0;
        const DoubleMatrix& locD = D[range.tIdx];
        if (F!=N) {
            int startN = xForFStart[f];
            int endN = (f<F-1 ? xForFStart[f+1] : N);
            for (int n=startN; n<endN; n++) {
                dsum += locD[dDimension][n];
            }
            return v/dsum;
        } else {
            return v/locD[dDimension][f];
        }
    }
    /* Maps a diffusion grid vector, d, to v_f = sum(x), where x = Dd */
    double CountingTree::vMapping(int f, const vector<int>& d) const {
        DoubleArray x = D[range.tIdx].mult(d);
        if (F!=N) {
            double vsum = 0.0;
            int startN = xForFStart[f];
            int endN = (f<F-1 ? xForFStart[f+1] : N);
            for (int n=startN; n<endN; n++) {
                vsum += x[n];
            }
            return vsum;
        } else {
            return x[f];
        }
    }

    void CountingTree::resizeGlobalLimits() {
        throw ModelException("CountingTree::resizeGlobalLimits: this hasn't been implemented, yet!");
        if (range.tIdx==T) {
            // this is a final resizing to inner and mid limits.
        } else {
            // check that (t+1) outer limits are within global limits
            // if not, you need to resize
        }
    }

     void CountingTree::computeInnerLimits() {
        throw ModelException("CountingTree::computeInnerLimits: this hasn't been implemented, yet!");
        
    }

DRLIB_END_NAMESPACE
