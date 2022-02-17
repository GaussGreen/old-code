
// Calibration logic
#include <time.h>

#include "Calibrator.h"
#include "Quadrature.h"
#include "Convolution.h"
#include "Tranche.h"
#include "CorrSkew.h"

#include "imsl.h"
#include "imsls.h"

//#include "gsl/gsl_multimin.h"
//#include "gsl/gsl_multifit_nlin.h"

#define ERR 999999999.99

namespace CM {

SymbolMap<CalibrationMode> CalibrationModeSymbolMap(
    QUASI_NEWTON, "quasi newton unconstrained",
    LM, "lm unconstrained",
    NON_LINEAR, "non linear constrained"
);

}

struct Param_2fact{
    int                 n;
    Array<double>       seeds;
    Array<int>          masks;
    Array<double>       targetSpreads; 
    Array<int>          loss;
    Array<double>       defaultProb;
    int                 notional;
    Array<double>       attachment;
    double              noIntervals;
    double              lbound;
    double              ubound;
    double              mbound;
} param_2fact;

struct Param_3factRR {
   int					n;
   Array<double>		seeds;
   Array<int>           masks;
   Array<double>        targetSpreads;
   double               r2;
   Array<double>		Recovery;
   int					indNotional;
   double				DefaultProb;
   Array<double>		Attachment;
   double				NoIntervals;
   double				lbound;							
   double				ubound;
   double				mbound;
} param_3factRR;

//////////////////////////////////////////////////////////////////////
//
// Calibration drivers
//
//////////////////////////////////////////////////////////////////////

// IMSL driver

double Calibrate2Fact_Driver_QuasiNewton(
        int n,
        double params[])
{
    double beta1;
    double beta2;
    double theta;

    int i = 0;
    if (param_2fact.masks[0])
        beta1 = params[i++];
    else beta1 = param_2fact.seeds[0];
    if (param_2fact.masks[1])
        beta2 = params[i++];
    else beta2 = param_2fact.seeds[1];
    if (param_2fact.masks[2])
        theta = params[i++];
    else theta = param_2fact.seeds[2];

    int noTranches = param_2fact.attachment.size() -1;

    Array<double> Beta(2);
    Beta[0] = beta1; Beta[1] = beta2;
    Array<double> Theta(1);
    Theta[0] = theta;

    // step 1, calculates expected losses
	Array<int> LossParam=CalcLossUnit(param_2fact.n,param_2fact.loss);
    Array<double> lossDist = LossDistRandomNFactorsSimp(
                                param_2fact.n,
                                Beta,
                                Theta,
                                param_2fact.loss, 
                                param_2fact.defaultProb, 
                                param_2fact.noIntervals, 
                                param_2fact.lbound, 
                                param_2fact.ubound, 
                                param_2fact.mbound);
    if (lossDist[0] < 0)
        return ERR; // xxx

	Array<double> TrancheLoss = ExpectedLossTranche(param_2fact.notional,
                                                    LossParam[0],
                                                    param_2fact.attachment,
                                                    lossDist);

    // step 2, calculate expected spreads using quick approximation
    Array<double> spreads(noTranches);
    spreads[0] = TrancheLoss[0] - (1 - exp(-500.*5/10000));
    for (i = 1; i < noTranches; ++i)
        spreads[i] = -log(1.-TrancheLoss[i])*10000/5;

    // step 3, calculate SSRE
    double error(0);
    double off(0);
    for (i = 0; i < noTranches; ++i) {
        off = (spreads[i] - param_2fact.targetSpreads[i])/param_2fact.targetSpreads[i];
        error += off*off;
    }

    return error;
}

void Calibrate2Fact_Driver_LM(
        int m,
        int n,
        double params[],
        double distance[])
{
    double beta1;
    double beta2;
    double theta;

    int i = 0;
    if (param_2fact.masks[0])
        beta1 = params[i++];
    else beta1 = param_2fact.seeds[0];
    if (param_2fact.masks[1])
        beta2 = params[i++];
    else beta2 = param_2fact.seeds[1];
    if (param_2fact.masks[2])
        theta = params[i++];
    else theta = param_2fact.seeds[2];

    int noTranches = param_2fact.attachment.size() -1;

    Array<double> Beta(2);
    Beta[0] = beta1; Beta[1] = beta2;
    Array<double> Theta(1);
    Theta[0] = theta;

    // step 1, calculates expected losses
	Array<int> LossParam=CalcLossUnit(param_2fact.n,param_2fact.loss);
    Array<double> lossDist = LossDistRandomNFactorsSimp(
                                param_2fact.n,
                                Beta,
                                Theta,
                                param_2fact.loss, 
                                param_2fact.defaultProb, 
                                param_2fact.noIntervals, 
                                param_2fact.lbound, 
                                param_2fact.ubound, 
                                param_2fact.mbound);

    if (lossDist[0] < 0) {
        for (i = 0; i < noTranches; ++i)
            distance[i] = ERR;
        return;
    }

	Array<double> TrancheLoss = ExpectedLossTranche(param_2fact.notional,
                                                    LossParam[0],
                                                    param_2fact.attachment,
                                                    lossDist);

    // step 2, calculate expected spreads using quick approximation
    Array<double> spreads(noTranches);
    spreads[0] = TrancheLoss[0] - (1 - exp(-500.*5/10000));
    for (i = 1; i < noTranches; ++i)
        spreads[i] = -log(1.-TrancheLoss[i])*10000/5;

    // step 3, calculate distance
    for (i = 0; i < noTranches; ++i)
        distance[i] = (spreads[i] - param_2fact.targetSpreads[i])/param_2fact.targetSpreads[i];
}

void Calibrate2Fact_Driver_NONLINEAR(
        int m,
        int meq,
        int n,
        double params[],
        int active[],
        double *f,
        double g[])
{
    double beta1;
    double beta2;
    double theta;

    int i = 0;
    if (param_2fact.masks[0])
        beta1 = params[i++];
    else beta1 = param_2fact.seeds[0];
    if (param_2fact.masks[1])
        beta2 = params[i++];
    else beta2 = param_2fact.seeds[1];
    if (param_2fact.masks[2])
        theta = params[i++];
    else theta = param_2fact.seeds[2];

    int noTranches = param_2fact.attachment.size() -1;

    Array<double> Beta(2);
    Beta[0] = beta1; Beta[1] = beta2;
    Array<double> Theta(1);
    Theta[0] = theta;

    // step 1, calculates expected losses
	Array<int> LossParam=CalcLossUnit(param_2fact.n,param_2fact.loss);
    Array<double> lossDist = LossDistRandomNFactorsSimp(
                                param_2fact.n,
                                Beta,
                                Theta,
                                param_2fact.loss, 
                                param_2fact.defaultProb, 
                                param_2fact.noIntervals, 
                                param_2fact.lbound, 
                                param_2fact.ubound, 
                                param_2fact.mbound);

    if (lossDist[0] < 0) {
        *f = ERR;
        g[0] = beta1 - beta2;
        return;
    }

	Array<double> TrancheLoss = ExpectedLossTranche(param_2fact.notional,
                                                    LossParam[0],
                                                    param_2fact.attachment,
                                                    lossDist);

    // step 2, calculate expected spreads using quick approximation
    Array<double> spreads(noTranches);
    spreads[0] = TrancheLoss[0] - (1 - exp(-500.*5/10000));
    for (i = 1; i < noTranches; ++i)
        spreads[i] = -log(1.-TrancheLoss[i])*10000/5;

    // step 3, calculate SSRE
    double error(0);
    double off(0);
    for (i = 0; i < noTranches; ++i) {
        off = (spreads[i] - param_2fact.targetSpreads[i])/param_2fact.targetSpreads[i];
        error += off*off;
    }
    *f = error;

    // calculate constraint
    if (active[0])
        g[0] = beta1 - beta2;

}

/*GSL driver, to be implemented
double Calibrate2Fact_Driver_GSL(
         const gsl_vector *v,
         void* params)
{

    return 0;
}*/


///////////////////////////////////////////////////////////////////
//
// Main calibration entry point
//
///////////////////////////////////////////////////////////////////

// IMSL version
Array<double> Calibrate2Fact(
        CalibrationMode     mode,
        int                 n,
        Array<double>&      seeds, // beta1, beta2, theta
        Array<int>&         masks,
        Array<double>&      targetSpreads, 
        Array<int>&         loss,
        Array<double>&      defaultProb,
        int                 notional,
        Array<double>&      attachment,
        double              noIntervals,
        double              lbound,
        double              ubound,
        double              mbound)
{
    // set up static paramsters
    param_2fact.n = n;
    param_2fact.seeds = seeds;
    param_2fact.masks = masks;
    param_2fact.targetSpreads = targetSpreads;
    param_2fact.loss = loss;
    param_2fact.defaultProb = defaultProb;
    param_2fact.notional = notional;
    param_2fact.attachment = attachment;
    param_2fact.noIntervals = noIntervals;
    param_2fact.lbound = lbound;
    param_2fact.ubound = ubound;
    param_2fact.mbound = mbound;

    int i;
    int numOfVars(0);
    for (i = 0; i < 3; ++i) {
        if (masks[i]) 
            ++numOfVars;
    }

    // set up initial guess & default answer
    i = 0;
    Array<double> rslt (5);
    double xguess[3];
    if (masks[0]) {
        xguess[i++] = seeds[0];
        rslt[0] = ERR;
    } else rslt[0] = seeds[0];
    if (masks[1]) {
        xguess[i++] = seeds[1];
        rslt[1] = ERR;
    } else rslt[1] = seeds[1];
    if (masks[2]) {
        xguess[i++] = seeds[2];
        rslt[2] = ERR;
    } else rslt[2] = seeds[2];
    rslt[3] = rslt[4] = ERR;

    double xscale[3]; xscale[0] = xscale[1] = xscale[2] = 0.00001;

    clock_t start = clock();
    double* sol;
    // call IMSL quasi-newton calibrate beta1, beta2, theta to target spread
    if (mode == QUASI_NEWTON)
        sol = imsl_d_min_uncon_multivar(Calibrate2Fact_Driver_QuasiNewton, numOfVars, 
                                            IMSL_XGUESS, xguess,
                                            IMSL_XSCALE, xscale,
                                            IMSL_REL_FCN_TOL, 0.1,
                                            IMSL_MAX_FCN, 300,
                                            0); 

    // call IMSL Levenberg-Marquardt calibrate
    if (mode == LM)
        sol = imsl_d_nonlin_least_squares(Calibrate2Fact_Driver_LM, 
                                              attachment.size()-1, 
                                              numOfVars,
                                              IMSL_XGUESS, xguess,
                                              // IMSL_REL_FCN_TOL, 0.1,
                                              IMSL_XSCALE, xscale,
                                              IMSL_MAX_FCN, 300,
                                              0);

    // call IMSL general nonlinear solver with constraint beta0 > beta1
    if (mode == NON_LINEAR) {
        double xlb[3]; xlb[0] = xlb[1] = xlb[2] = -1.0e6;
        double xub[3]; xub[0] = xub[1] = xub[2] = 1.0e6; 
        sol = imsl_d_min_con_nonlin(Calibrate2Fact_Driver_NONLINEAR,
                                    1, // # constraints
                                    0, // # equality constraints
                                    numOfVars, // # variables
                                    0, // user supply bounds
                                    xlb, // low bound
                                    xub, // up bound
                                    IMSL_XGUESS, xguess,
                                    IMSL_XSCALE, xscale,
                                    IMSL_ERR_REL, 0.1,
                                    IMSL_ITMAX, 300,
                                    0);
    }

    if (sol == 0) {
        long err = imsl_error_code();
        rslt[3] = err;
        return rslt;
    }

    clock_t end = clock();
    double error = Calibrate2Fact_Driver_QuasiNewton(numOfVars, sol);

    i = 0;
    if (masks[0])
        rslt[0] = sol[i++];
    if (masks[1])
        rslt[1] = sol[i++];
    if (masks[2])
        rslt[2] = sol[i++];
    rslt[3] = error;
    rslt[4] = static_cast<double>(end - start) / CLOCKS_PER_SEC;

    free(sol);

    return rslt;
}

/* GSL Version, to be implemented
Array<double> Calibrate2Fact_GSL(
        CalibrationMode     mode,
        int                 n,
        double              beta1Seed, // initial guesses
        double              beta2Seed,
        double              thetaSeed,
        Array<double>&      targetSpreads, 
        Array<int>&         loss,
        Array<double>&      defaultProb,
        int                 notional,
        Array<double>&      attachment,
        double              noIntervals,
        double              lbound,
        double              ubound,
        double              mbound)
{
    // initialize minimizer state
    const gsl_multimin_fdfminimizer_type *T;
    T = gsl_multimin_fdfminimizer_conjugate_fr; // Fletcher-Reeves conjugate
    // T = gsl_multimin_fdfminimizer_conjugate_pr; // Polak-Ribiere conjugate
    // T = gsl_multimin_fdfminimizer_vector_bfgs; // vector BFGS conjugate
    // T = gsl_multimin_fdfminimizer_steepest_descent; // 
    // T = gsl_multimin_fminimizer_nmsimplex; // NM simplex

    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(T, 3);
    
    gsl_multimin_function_fdf driver;
    // driver.f = 

    // run iteration

    Array<double> rslt(3);
    return rslt;
} */


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

double Calibrate3FactRR_Driver_QuasiNewton(
        int n,
        double params[])
{    
    double beta1;
    double beta2;
    double beta3;
    double theta1;
    double theta2;
    double r1;

    int i = 0;
    if (param_3factRR.masks[0])
        beta1 = params[i++];
    else beta1 = param_3factRR.seeds[0];
    if (param_3factRR.masks[1])
        beta2 = params[i++];
    else beta2 = param_3factRR.seeds[1];
    if (param_3factRR.masks[2])
        beta3 = params[i++];
    else beta3 = param_3factRR.seeds[2];
    if (param_3factRR.masks[3])
        theta1 = params[i++];
    else theta1 = param_3factRR.seeds[3];
    if (param_3factRR.masks[4])
        theta2 = params[i++];
    else theta2 = param_3factRR.seeds[4];
    if (param_3factRR.masks[5])
        r1 = params[i++];
    else r1 = param_3factRR.seeds[5];

    int noTranches = param_3factRR.Attachment.size() -1;

    // start computation

    Array<double> Beta(3);
    Beta[0] = beta1; Beta[1] = beta2; Beta[2] = beta3;
    Array<double> Theta(2);
    Theta[0] = theta1; Theta[1] = theta2;
    Array<double> ThresholdR(2);
    ThresholdR[0] = r1;
    ThresholdR[1] = param_3factRR.r2;

    // step 1, calculates expected losses

	double loss0=100*(1-param_3factRR.Recovery[0]);
	double loss1=100*(1-param_3factRR.Recovery[1]);
	double loss2=100*(1-param_3factRR.Recovery[2]);

	int minloss=GCD(int(loss0),GCD(int(loss1),int(loss2)));

	int rat0=loss0/minloss;
	int rat1=loss1/minloss;
	int rat2=loss2/minloss;
	
	double Notional=param_3factRR.indNotional*param_3factRR.n;
	double lossUnit=(1-param_3factRR.Recovery[0])*param_3factRR.indNotional/(rat0);
	
	Array<double> lossDist = 
        LossDistRandomNFactorsHomogenSimpRR(param_3factRR.n, Beta, Theta, 0, 
                                            param_3factRR.DefaultProb, 
                                            param_3factRR.NoIntervals, 
                                            param_3factRR.Recovery, 
                                            ThresholdR,
                                            param_3factRR.lbound, 
                                            param_3factRR.ubound, 
                                            param_3factRR.mbound);

    if (lossDist[0] < 0)
        return ERR;

	Array<double> TrancheLoss = ExpectedLossTranche(Notional,
                                                    lossUnit,
                                                    param_3factRR.Attachment,
                                                    lossDist);

    // step 2, calculate expected spreads using quick approximation
    Array<double> spreads(noTranches);
    spreads[0] = TrancheLoss[0] - (1 - exp(-500.*5/10000));
    for (i = 1; i < noTranches; ++i) {
        spreads[i] = -log(1.-TrancheLoss[i])*10000/5;
    }

    // step 3, calculate SSRE
    double error(0);
    double off;
    for (i = 0; i < noTranches; ++i) {
        off = (spreads[i] - param_3factRR.targetSpreads[i])/param_3factRR.targetSpreads[i];
        error += off*off;
    }

    return error;
}


void Calibrate3FactRR_Driver_LM(
        int m,
        int n,
        double params[],
        double distance[])
{
    double beta1;
    double beta2;
    double beta3;
    double theta1;
    double theta2;
    double r1;

    int i = 0;
    if (param_3factRR.masks[0])
        beta1 = params[i++];
    else beta1 = param_3factRR.seeds[0];
    if (param_3factRR.masks[1])
        beta2 = params[i++];
    else beta2 = param_3factRR.seeds[1];
    if (param_3factRR.masks[2])
        beta3 = params[i++];
    else beta3 = param_3factRR.seeds[2];
    if (param_3factRR.masks[3])
        theta1 = params[i++];
    else theta1 = param_3factRR.seeds[3];
    if (param_3factRR.masks[4])
        theta2 = params[i++];
    else theta2 = param_3factRR.seeds[4];
    if (param_3factRR.masks[5])
        r1 = params[i++];
    else r1 = param_3factRR.seeds[5];

    int noTranches = param_3factRR.Attachment.size() -1;

    // start computation

    Array<double> Beta(3);
    Beta[0] = beta1; Beta[1] = beta2; Beta[2] = beta3;
    Array<double> Theta(2);
    Theta[0] = theta1; Theta[1] = theta2;
    Array<double> ThresholdR(2);
    ThresholdR[0] = r1;
    ThresholdR[1] = param_3factRR.r2;

    // step 1, calculates expected losses

	double loss0=100*(1-param_3factRR.Recovery[0]);
	double loss1=100*(1-param_3factRR.Recovery[1]);
	double loss2=100*(1-param_3factRR.Recovery[2]);

	int minloss=GCD(int(loss0),GCD(int(loss1),int(loss2)));

	int rat0=loss0/minloss;
	int rat1=loss1/minloss;
	int rat2=loss2/minloss;
	
	double Notional=param_3factRR.indNotional*param_3factRR.n;
	double lossUnit=(1-param_3factRR.Recovery[0])*param_3factRR.indNotional/(rat0);
	
	Array<double> lossDist = 
        LossDistRandomNFactorsHomogenSimpRR(param_3factRR.n, Beta, Theta, 0, 
                                            param_3factRR.DefaultProb, 
                                            param_3factRR.NoIntervals, 
                                            param_3factRR.Recovery, 
                                            ThresholdR,
                                            param_3factRR.lbound, 
                                            param_3factRR.ubound, 
                                            param_3factRR.mbound);

    if (lossDist[0] < 0) {
        for (int i = 0; i < noTranches; ++i)
            distance[i] = ERR;
        return;
    }

	Array<double> TrancheLoss = ExpectedLossTranche(Notional,
                                                    lossUnit,
                                                    param_3factRR.Attachment,
                                                    lossDist);

    // step 2, calculate expected spreads using quick approximation
    Array<double> spreads(noTranches);
    spreads[0] = TrancheLoss[0] - (1 - exp(-500.*5/10000));
    for (i = 1; i < noTranches; ++i) {
        spreads[i] = -log(1.-TrancheLoss[i])*10000/5;
    }

    // step 3, calculate distance
    for (i = 0; i < noTranches; ++i) {
        distance[i] = (spreads[i] - param_3factRR.targetSpreads[i])/param_3factRR.targetSpreads[i];
    }
}

void Calibrate3FactRR_Driver_NONLINEAR(
        int m,
        int meq,
        int n,
        double params[],
        int active[],
        double *f,
        double g[])
{
    double beta1;
    double beta2;
    double beta3;
    double theta1;
    double theta2;
    double r1;

    int i = 0;
    if (param_3factRR.masks[0])
        beta1 = params[i++];
    else beta1 = param_3factRR.seeds[0];
    if (param_3factRR.masks[1])
        beta2 = params[i++];
    else beta2 = param_3factRR.seeds[1];
    if (param_3factRR.masks[2])
        beta3 = params[i++];
    else beta3 = param_3factRR.seeds[2];
    if (param_3factRR.masks[3])
        theta1 = params[i++];
    else theta1 = param_3factRR.seeds[3];
    if (param_3factRR.masks[4])
        theta2 = params[i++];
    else theta2 = param_3factRR.seeds[4];
    if (param_3factRR.masks[5])
        r1 = params[i++];
    else r1 = param_3factRR.seeds[5];

    int noTranches = param_3factRR.Attachment.size() -1;

    // start computation

    Array<double> Beta(3);
    Beta[0] = beta1; Beta[1] = beta2; Beta[2] = beta3;
    Array<double> Theta(2);
    Theta[0] = theta1; Theta[1] = theta2;
    Array<double> ThresholdR(2);
    ThresholdR[0] = r1;
    ThresholdR[1] = param_3factRR.r2;

    // step 1, calculates expected losses

	double loss0=100*(1-param_3factRR.Recovery[0]);
	double loss1=100*(1-param_3factRR.Recovery[1]);
	double loss2=100*(1-param_3factRR.Recovery[2]);

	int minloss=GCD(int(loss0),GCD(int(loss1),int(loss2)));

	int rat0=loss0/minloss;
	int rat1=loss1/minloss;
	int rat2=loss2/minloss;
	
	double Notional=param_3factRR.indNotional*param_3factRR.n;
	double lossUnit=(1-param_3factRR.Recovery[0])*param_3factRR.indNotional/(rat0);
	
	Array<double> lossDist = 
        LossDistRandomNFactorsHomogenSimpRR(param_3factRR.n, Beta, Theta, 0, 
                                            param_3factRR.DefaultProb, 
                                            param_3factRR.NoIntervals, 
                                            param_3factRR.Recovery, 
                                            ThresholdR,
                                            param_3factRR.lbound, 
                                            param_3factRR.ubound, 
                                            param_3factRR.mbound);

    if (lossDist[0] < 0) {
        *f = ERR;
        if (active[0])
            g[0] = beta1 - beta2;
        if (active[1])
            g[1] = beta2 - beta3;
        if (active[2])
            g[2] = theta2 - theta1;
        return;
    }

	Array<double> TrancheLoss = ExpectedLossTranche(Notional,
                                                    lossUnit,
                                                    param_3factRR.Attachment,
                                                    lossDist);

    // step 2, calculate expected spreads using quick approximation
    Array<double> spreads(noTranches);
    spreads[0] = TrancheLoss[0] - (1 - exp(-500.*5/10000));
    for (i = 1; i < noTranches; ++i) {
        spreads[i] = -log(1.-TrancheLoss[i])*10000/5;
    }

    // step 3, calculate SSRE
    double error(0);
    double off;
    for (i = 0; i < noTranches; ++i) {
        off = (spreads[i] - param_3factRR.targetSpreads[i])/param_3factRR.targetSpreads[i];
        error += off*off;
    }

    *f = error;

    if (active[0])
        g[0] = beta1 - beta2;
    if (active[1])
        g[1] = beta2 - beta3;
    if (active[2])
        g[2] = theta2 - theta1;
    
}



Array<double> Calibrate3FactHomogRR(
       CalibrationMode      mode,
	   int					n,
	   Array<double>		&seeds,
       Array<int>           &masks,
       Array<double>        &targetSpreads,
       double               r2,
	   Array<double>		&Recovery,
	   int					indNotional,
	   double				DefaultProb,
	   Array<double>		&Attachment,
	   double				NoIntervals,
	   double				lbound,							
	   double				ubound,
	   double				mbound)
{
    // set up static parameters
    param_3factRR.n = n;
    param_3factRR.seeds = seeds;
    param_3factRR.masks = masks;
    param_3factRR.targetSpreads = targetSpreads;
    param_3factRR.r2 = r2;
    param_3factRR.Recovery = Recovery;
    param_3factRR.indNotional = indNotional;
    param_3factRR.DefaultProb = DefaultProb;
    param_3factRR.Attachment = Attachment;
    param_3factRR.NoIntervals = NoIntervals;
    param_3factRR.lbound = lbound;
    param_3factRR.ubound = ubound;
    param_3factRR.mbound = mbound;

    int i;
    int numOfVars(0);
    for (i = 0; i < 6; ++i) {
        if (masks[i]) 
            ++numOfVars;
    }

    // set up initial guess & default answer
    i = 0;
    Array<double> rslt (8);
    double xguess[6];
    if (masks[0]) {
        xguess[i++] = seeds[0];
        rslt[0] = ERR;
    } else rslt[0] = seeds[0];
    if (masks[1]) {
        xguess[i++] = seeds[1];
        rslt[1] = ERR;
    } else rslt[1] = seeds[1];
    if (masks[2]) {
        xguess[i++] = seeds[2];
        rslt[2] = ERR;
    } else rslt[2] = seeds[2];
    if (masks[3]) {
        xguess[i++] = seeds[3];
        rslt[3] = ERR;
    } else rslt[3] = seeds[3];
    if (masks[4]) {
        xguess[i++] = seeds[4];
        rslt[4] = ERR;
    } else rslt[4] = seeds[4];
    if (masks[5]) {
        xguess[i++] = seeds[5];
        rslt[5] = ERR;
    } else rslt[5] = seeds[5];
    rslt[6] = rslt[7] = -1;

    double xscale[6]; 
    xscale[0] = xscale[1] = xscale[2] = 0.0001;
    xscale[3] = xscale[4] = xscale[5] = 0.0001;

    clock_t start = clock();
    double* sol;
    
    if (mode == QUASI_NEWTON)
        sol = imsl_d_min_uncon_multivar(Calibrate3FactRR_Driver_QuasiNewton, numOfVars, 
                                            IMSL_XGUESS, xguess,
                                            IMSL_XSCALE, xscale,
                                            IMSL_MAX_FCN, 1600,
                                            0);
    if (mode == LM)
        sol = imsl_d_nonlin_least_squares(Calibrate3FactRR_Driver_LM,
                                            Attachment.size()-1,
                                            numOfVars,
                                            IMSL_XGUESS, xguess,
                                            IMSL_XSCALE, xscale,
                                            IMSL_MAX_FCN, 1600,
                                            0);

    if (mode == NON_LINEAR) {
        double xlb[6]; 
        xlb[0] = xlb[1] = xlb[2] = xlb[3] = xlb[4] = xlb[5] = -1.0e6;
        double xub[6]; 
        xub[0] = xub[1] = xub[2] = xub[3] = xub[4] = xub[5] = 1.0e6; 
        sol = imsl_d_min_con_nonlin(Calibrate3FactRR_Driver_NONLINEAR,
                                    3, // # constraints
                                    0, // # equality constraints
                                    numOfVars, // # variables
                                    0, // user supply bounds
                                    xlb, // low bound
                                    xub, // up bound
                                    IMSL_XGUESS, xguess,
                                    IMSL_XSCALE, xscale,
                                    IMSL_ERR_REL, 0.1,
                                    IMSL_ITMAX, 1600,
                                    0);
    }

    if (sol == 0) {
        long err = imsl_error_code();
        rslt[6] = err;
        return rslt;
    }

    clock_t end = clock();

    i = 0;
    if (masks[0])
        rslt[0] = sol[i++];
    if (masks[1])
        rslt[1] = sol[i++];
    if (masks[2])
        rslt[2] = sol[i++];
    if (masks[3])
        rslt[3] = sol[i++];
    if (masks[4])
        rslt[4] = sol[i++];
    if (masks[5])
        rslt[5] = sol[i++];
    rslt[6] = Calibrate3FactRR_Driver_QuasiNewton(numOfVars, sol);
    rslt[7] = static_cast<double>(end - start) / CLOCKS_PER_SEC;

    free(sol);

    return rslt;
}
