//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMCorrelation.cpp
//
//   Description : Helper for SRM - used for splitting up correlation between
//                 the different factors
//
//   Author      : Mark A Robson
//
//   Date        : 12 July 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/SRMCorrelation.hpp"
#include "edginc/SRMConstants.hpp"
#include "edginc/SRMGenDiffusibleAsset.hpp"

DRLIB_BEGIN_NAMESPACE


// Various static functions to expand corrs for different asset combinations:
// IR vs IR corr split
vector<double> SRMCorrelation::ExpandAssetAssetCorrs(
    SRMGenDiffusibleIR* irAsset1, 
    SRMGenDiffusibleIR* irAsset2, 
    double corr, 
    OutputType choice)
{
    return irToIr(
        corr, 
        *irAsset1->getSRMRatesUtil().get(), 
        *irAsset2->getSRMRatesUtil().get(), 
        choice);
}

// IR vs Energy corr split
vector<double> SRMCorrelation::ExpandAssetAssetCorrs(
    SRMGenDiffusibleIR* irAsset, 
    SRMGenDiffusibleEnergy* enAsset,
    double corr, 
    OutputType irChoice, 
    OutputType enChoice)
{
    return enrgToIr(
        corr,
        enChoice,
        irChoice,
        *enAsset->getSRMEnrgUtil().get(),
        *irAsset->getSRMRatesUtil().get(),
        false // IR is the primary asset in this case
        );
}

// 1-factor vs IR corr split
vector<double> SRMCorrelation::ExpandAssetAssetCorrs(
    SRMGenDiffusibleIR* irAsset, 
    double corr, 
    OutputType choice)
{
    return assetIr(corr, *irAsset->getSRMRatesUtil().get(), choice);
}

// 1-factor vs Energy corr split
vector<double> SRMCorrelation::ExpandAssetAssetCorrs(
    SRMGenDiffusibleEnergy* enAsset, 
    double corr, 
    OutputType choice)
{
    return enrgToAsset(corr, choice, *enAsset->getSRMEnrgUtil().get());
}

// Energy vs Energy corr split
vector<double> SRMCorrelation::ExpandAssetAssetCorrs(
    SRMGenDiffusibleEnergy* enAsset1, 
    SRMGenDiffusibleEnergy* enAsset2, 
    double corr, 
    OutputType choice)
{
    return enrgToEnrg(
        corr, 
        choice, 
        *enAsset1->getSRMEnrgUtil().get(), 
        *enAsset2->getSRMEnrgUtil().get());
}

// Energy vs IR corr split
vector<double> SRMCorrelation::ExpandAssetAssetCorrs(
    SRMGenDiffusibleEnergy* enAsset, 
    SRMGenDiffusibleIR* irAsset, 
    double corr, 
    OutputType enChoice, 
    OutputType irChoice)
{
    return enrgToIr(
        corr, 
        enChoice, 
        irChoice, 
        *enAsset->getSRMEnrgUtil().get(), 
        *irAsset->getSRMRatesUtil().get(),
        true // Energy is primary asset
        );
}


/***** From util_s.c:IR_IR_Corr ******************************************
 *
 *      Produces the implied correlation structure between YC driving factors
 *      ie the correlation between "dom" and "for" factors  
 *
 *      choice for output :
 *
 *      1 = EXPONENTIAL FACTORS
 *      2 = TRIANGULAR
 *      3 = PC
 * 
 *      Note: in case dimension of yield curve is 1, DomRho and ForRho
 *            can be both input as NULL.
 *
 *      This function is Flat Forward ready
 *  Really should return a DoubleMatrix but that requires too many changes
 */
vector<double> SRMCorrelation::irToIr( 
    double         RhoDswapFswap,  /*(I) corr. between "dom" and "for" rates */
    const SRMRatesUtil& dom,     //(I) 'domestic' IR data
    const SRMRatesUtil& fgn,     //(I) 'foreign' IR data
    OutputType     choice){      /*(I) specifies factors wrt which correl is
                                   output: 1=expon.,2=triang;3=PC      */
    static const string method("irToIr");
    try{
        /* MAR: In the SRM3 code the zero curve has today equal to the spot date
           so need to do all yearFracs etc from spot date */
        const DateTime& domSpotDate = dom.getDiffYC()->getSpotDate();
        /* DomB and ForB are B factors in Christian's paper */
        vector<double> DomB(dom.bFactor());
        double length = domSpotDate.yearFrac(dom.getCorrSwapStart());
        vector<double> DomBeta(DomB.size());
        dom.getBeta(DomBeta);
        vector<double> DomAlpha(DomB.size());
        dom.getAlpha(DomAlpha);
        const vector<double>& DomRho = dom.getRho();
        int DomNbFactor = DomBeta.size(); // for ease
        int i;
        for (i = 0; i < DomNbFactor; i++) {
            DomB[i] *= exp(-DomBeta[i] * length);
        }
        // do the same for foreign
        /* DomB and ForB are B factors in Christian's paper */
        vector<double> ForB(fgn.bFactor());
        vector<double> ForAlpha(ForB.size());
        fgn.getAlpha(ForAlpha);
        vector<double> ForBeta(DomB.size());
        fgn.getBeta(ForBeta);
        const vector<double>& ForRho = fgn.getRho();
        const DateTime& fgnSpotDate = fgn.getDiffYC()->getSpotDate();
        length = fgnSpotDate.yearFrac(fgn.getCorrSwapStart());
        int ForNbFactor = ForBeta.size();
        for (i = 0; i < ForNbFactor; i++) {
            ForB[i] *= exp(-ForBeta[i] * length);
        }
        /* calculate now normalised weights of orthogonal factors*/
        vector<double> DomA(DomNbFactor);
        vector<double> ForA(ForNbFactor);
        DoubleMatrix DomTransform;
        DoubleMatrix ForTransform;
        switch (choice) {
        case EXPONENTIAL_FACTORS:
        case TRIANGULAR:
            DomTransform = dom.triangulation();
            ForTransform = fgn.triangulation();
            if (DomNbFactor > 0) {
                DomA[0] = DomAlpha[0] * DomB[0];
            }
            if (ForNbFactor > 0) {
                ForA[0] = ForAlpha[0] * ForB[0];
            }
            if (DomBeta.size() > 1) {                
                /* domestic weights */
                DomA[0] += DomAlpha[1] * DomB[1] * DomTransform[1][0];
                DomA[1]  = DomAlpha[1] * DomB[1] * DomTransform[1][1];
            }/* end of "DomNbFactor > 1"*/
            if (DomNbFactor > 2) {                
                /* domestic weights */
                DomA[0] += DomAlpha[2] * DomB[2] * DomTransform[2][0];
                DomA[1] += DomAlpha[2] * DomB[2] * DomTransform[2][1];
                DomA[2]  = DomAlpha[2] * DomB[2] * DomTransform[2][2];
            }/* end of DomNbFactor>2*/
            if (ForNbFactor > 1) {                
                /* foreign weights */
                ForA[0] += ForAlpha[1] * ForB[1] * ForTransform[1][0];
                ForA[1]  = ForAlpha[1] * ForB[1] * ForTransform[1][1];
            }/* end of ForNbFactor>1*/
            if (ForNbFactor > 2) {
                /* foreign weights */
                ForA[0] += ForAlpha[2] * ForB[2] * ForTransform[2][0];
                ForA[1] += ForAlpha[2] * ForB[2] * ForTransform[2][1];
                ForA[2]  = ForAlpha[2] * ForB[2] * ForTransform[2][2];
            }/* end of ForNbFactor>2*/
            break;
        case PC:
            /* PC for Domestic curve*/
            DoubleMatrix DomOmtx(dom.get3A());
            for (i = 0; i < DomNbFactor; i++) {
                for (int j = 0; j < DomNbFactor; j++) {
                    DomA[i] += DomAlpha[j] * DomB[j] * DomOmtx[j][i];
                } /* for j*/
            }/* for i*/

            /* PC for Foreign curve */
            DoubleMatrix ForOmtx(fgn.get3A());
            for (i = 0; i < ForNbFactor; i++) {
                for (int j = 0; j < ForNbFactor; j++) {
                    ForA[i] += ForAlpha[j] * ForB[j] * ForOmtx[j][i];
                } /* for j*/
            
            }/* for i*/
            break;
        }/*end of switch on choice*/

        /* normalizing domestic weights*/
        double sum = 0.;
        for (i = 0; i < DomNbFactor; i++){
            sum += DomA[i] * DomA[i];
        }
        sum = sqrt(sum);
        if (DomNbFactor > 0 && sum < SRMConstants::SRM_TINY) {
            throw ModelException(method, "Problem in Customised_Correl_Mtx: "
                                 "sum of the domestic weights is zero, unable"
                                 " to normalise");
        }
        for (i = 0; i < DomNbFactor; i++){
            DomA[i] /= sum;
        }
        /* normalizing foreign weights*/
        sum = 0.;
        for (i = 0; i < ForNbFactor; i++){
            sum += ForA[i] * ForA[i];
        }
        sum = sqrt(sum);
        if (ForNbFactor > 0 && sum < SRMConstants::SRM_TINY) {
            throw ModelException(method, "Problem in Customised_Correl_Mtx: "
                                 "sum of the foreign weights is zero, unable"
                                 " to normalise");
        }
        for (i = 0; i < ForNbFactor; i++){
            ForA[i] /= sum;
        }
        // create space for intermediate results
        vector<double> DummyRhoDomFor(DomNbFactor * ForNbFactor);
    
        for (i = 0; i< DomNbFactor; i++) {
            for (int j = 0; j < ForNbFactor; j++) {
                DummyRhoDomFor[j + i * ForNbFactor] = 
                    DomA[i] * ForA[j] * RhoDswapFswap;
            }/* for j */
        } /* for i */
    
        // create space for results
        vector<double> RhoDomFor(DomNbFactor * ForNbFactor);
        if (choice == EXPONENTIAL_FACTORS) {   
            /* correlations between domestic and foreign factors*/
            for (i = 0; i < DomNbFactor; i++) {
                for (int j = 0; j < ForNbFactor; j++) {
                    sum = 0.0;
                    for (int k = 0; k < DomNbFactor; k++) {
                        double cum = 0.0;
                        for (int l = 0; l < ForNbFactor; l++) {
                            cum += ForTransform[j][l] * 
                                DummyRhoDomFor[l + k * ForNbFactor];
                        } /* for l */
                        sum += DomTransform[i][k] * cum;
                    } /* for k */
                    RhoDomFor[j + i * ForNbFactor] = sum;
                } /* for j */
            } /* for i */

            /* now test whether we recover input correlation */
            double test_DVol = 0.;
            for (i = 0; i < DomNbFactor; i++) {
                test_DVol += DomAlpha[i] * DomAlpha[i]* DomB[i]* DomB[i];
            }
            if (DomNbFactor  > 1) {
                test_DVol += 2.0 * DomAlpha[0] * DomAlpha[1] * 
                    DomB[0] * DomB[1]* DomRho[0];
            }
            if (DomNbFactor > 2) {
                test_DVol += 2.0 * DomAlpha[0] * DomAlpha[2] * 
                    DomB[0] * DomB[2] * DomRho[1];
                test_DVol += 2.0 * DomAlpha[1] * DomAlpha[2] * 
                    DomB[1] * DomB[2] * DomRho[2];
            }
            if (DomNbFactor > 0 && test_DVol < SRMConstants::SRM_TINY) {
                throw ModelException(method, "Problem calculating test Vol "
                                     "for domestic rate");
            }
            double test_FVol = 0;
            for (i = 0; i < ForNbFactor; i++) {
                test_FVol += ForAlpha[i] * ForAlpha[i] * ForB[i] * ForB[i];
            }
            if (ForNbFactor  > 1) {
                test_FVol += 2.0 * ForAlpha[0] * ForAlpha[1] * 
                    ForB[0] * ForB[1] * ForRho[0];
            }
            if (ForNbFactor > 2) {
                test_FVol += 2.0 * ForAlpha[0] * ForAlpha[2] * 
                    ForB[0] * ForB[2] * ForRho[1];
                test_FVol += 2.0 * ForAlpha[1] * ForAlpha[2] * 
                    ForB[1] * ForB[2] * ForRho[2];
            }
            if (ForNbFactor > 0 && test_FVol <SRMConstants::SRM_TINY) {
                throw ModelException(method,
                                     "Problem calculating test Vol for foreign "
                                     "rate");
            }
            /* Corr For versus Dom*/
            double test_RhoDswapFswap = 0.;
            for (i = 0; i < DomNbFactor; i++) {
                for (int j = 0; j < ForNbFactor; j++) {                
                    test_RhoDswapFswap += DomAlpha[i] * DomB[i] *
                        ForAlpha[j] * ForB[j] * 
                        RhoDomFor[j + i * ForNbFactor];
                }/* for j*/
            }/* for i*/
            test_RhoDswapFswap /= sqrt(test_FVol);
            test_RhoDswapFswap /= sqrt(test_DVol);
            if ((fabs(test_RhoDswapFswap - RhoDswapFswap) > SRMConstants::SRM_TINY)) {
                throw ModelException(method, "Failed to reproduce"
                                     "input correlation coefficient");
            }
        }/* end of choice=1*/
        else if (choice == TRIANGULAR || choice == PC) {   
            /* the correlation in this case have already been obtained */
            return DummyRhoDomFor;
        }
        return RhoDomFor;
    } catch (exception& e){
        throw ModelException(e, method, "Failed for currencies "+
                             dom.getDiscYC()->getCcy()+" [domestic] and "+
                             fgn.getDiscYC()->getCcy()+" [foreign]");
    }
}/* IR_IR_Corr */


/***** From util_s.c::AssetIrCorr ************************************/
//
//      Produces the implied corr structure between
//      Asset and the factors driving the IR curve.
// 
//      the factors can be:
//      choice 1 = EXPONENTIAL 
//      choice 2 = TRIANGULAR
//      choice 3 = PC
//
//      In case the YC dim is 1, Rho can be input as NULL.
//      Rk: Asset is the equivalent of the FX rate in the doc
//
//      This function is Flat Forward ready
//
vector<double> SRMCorrelation::assetIr(
    double         RhoAssetSwap,  /*(I) corr between Asset and swap rate */
    const SRMRatesUtil& ir,            //(I) IR data
    OutputType     choice){       /*(I) specifies factors wrt which correl is
                                    output: 1=expon.,2=triang;3=PC      */
    static const string method("SRMCorrelation::assetIR");
    try{
        /* MAR: In the SRM3 code the zero curve has today equal to the spot date
           so need to do all yearFracs etc from spot date */
        const DateTime& spotDate = ir.getDiffYC()->getSpotDate();
        int NbFactor = ir.numFactors(); // for ease

        //Zero factor: deterministic rates: 
        //size of return vector is the number of factors, i.e. zero.
        if (NbFactor == 0)
        {
            return vector<double>();
        }

        // get the 'B factor'
        vector<double> B(ir.bFactor());
        double length = spotDate.yearFrac(ir.getCorrSwapStart());
        vector<double> Beta;
        ir.getBeta(Beta);
        vector<double> Alpha;
        ir.getAlpha(Alpha);
        int i;
        for (i = 0; i < NbFactor; i++) {          
            B[i] *= exp(- Beta[i] * length);
        }
        vector<double> A(NbFactor);
        DoubleMatrix Transform;
        /* calculate now normalised weights */
        switch (choice) {
        case EXPONENTIAL_FACTORS:
        case TRIANGULAR:
        {
            Transform = ir.triangulation();
            A[0] = Alpha[0] * B[0];
            if (NbFactor > 1) {                
                A[0] += Alpha[1] * B[1] * Transform[1][0];
                A[1]  = Alpha[1] * B[1] * Transform[1][1];
            }/* end of "NbFactor > 1"*/
            if (NbFactor > 2) {                
                A[0] += Alpha[2] * B[2] * Transform[2][0];
                A[1] += Alpha[2] * B[2] * Transform[2][1];
                A[2]  = Alpha[2] * B[2] * Transform[2][2];
            
            }/* end of NbFactor>2*/
            break;
        }
        case PC:
        {
            DoubleMatrix OMtx(ir.get3A());
            for (int i = 0; i < NbFactor; i++) {
                for (int j = 0; j < NbFactor; j++) {
                    A[i] += Alpha[j] * B[j] * OMtx[j][i];
                } /* for j*/
            }/* for i*/
            break;
        }
        }/*end of switch on choice*/

        /* normalizing domestic weights*/
        double sum = 0.0;
        for (i = 0; i < NbFactor; i++){
            sum += A[i] * A[i];
        }
        sum = sqrt(sum);
        if (sum < SRMConstants::SRM_TINY) {
            throw ModelException(method, "Sum of the normalising weights is "
                                 "zero, unable to normalise");
        }

        for (i = 0; i < NbFactor; i++){
            A[i] /= sum;
        }
        /* DummyRhoAssetIr are the corr mentioned in doc*/ 
        vector<double> DummyRhoAssetIr(NbFactor);
        for (i = 0; i < NbFactor; i++) {
            DummyRhoAssetIr[i] = A[i] * RhoAssetSwap;
        }
    
        if (choice == EXPONENTIAL_FACTORS) {   
            vector<double> RhoAssetIr(NbFactor); // return value
            for (i = 0; i < NbFactor; i++) {
                sum = 0.;
                for (int j = 0; j < NbFactor; j++) {
                    sum += Transform[i][j] * DummyRhoAssetIr[j];
                } /* for j */
                RhoAssetIr[i] = sum;
            }/* for i */
                
            /* now test whether we recover input correlations*/
        
            /* Corr Fx versus Dom*/
            double test_RhoAssetSwap = 0.0;
            double test_Vol = 0.0;
            vector<double> Rho(ir.getRho());
            for (i = 0; i < NbFactor; i++) {
                test_Vol += Alpha[i] * Alpha[i]* B[i]* B[i];
            }

            if (NbFactor  > 1) {
                test_Vol += 2.0 * Alpha[0] * Alpha[1] * B[0] * B[1]* Rho[0];
            }

            if (NbFactor > 2) {
                test_Vol += 2.0 * Alpha[0] * Alpha[2] * B[0] * B[2] * Rho[1];
                test_Vol += 2.0 * Alpha[1] * Alpha[2] * B[1] * B[2] * Rho[2];
            }

            if (test_Vol < SRMConstants::SRM_TINY) {
                throw ModelException(method,
                                     "Problem calculating test Vol for IR");
            }

            for (i = 0 ; i < NbFactor; i++) {
                test_RhoAssetSwap += Alpha[i] * B[i] * RhoAssetIr[i];
            }
            test_RhoAssetSwap /= sqrt(test_Vol);

            if (fabs(test_RhoAssetSwap - RhoAssetSwap) > SRMConstants::SRM_TINY) {
                throw ModelException(method, "Failed to reproduce"
                                     "input correlation coefficient");
            }
            return RhoAssetIr;
        }/* end of choice=1*/
        else {
            /* the correlations in this case have already been obtained */
            return DummyRhoAssetIr;
        }
    } catch (exception& e){
        throw ModelException(e, method, "Failed for ccy "+ir.getCcy());
    }
}


/***** From util_s.c:enrgToEnrg ******************************************
*
*      Produces the implied correlation structure between driving factors
*      of two energy products (e.g. 2 energy forward prices)
*
*      choice for output :
*
*      1 = EXPONENTIAL FACTORS
*      2 = TRIANGULAR
*      3 = PC
* 
*      Note: in case dimension of energy curve is 1, rho1 and rho2
*            can be both input as NULL.
*
*      This function is Flat Forward ready
*  Really should return a DoubleMatrix but that requires too many changes
*/
vector<double> SRMCorrelation::enrgToEnrg(
	double rhoEnrgEnrg,   /*(I) corr between two energy instruments */
	OutputType choice,         /*(I) specifies factors wrt which correl is */
	const SRMEnergyUtilBase& enrgUtil1,
	const SRMEnergyUtilBase& enrgUtil2
	/*const DateTime & today,
	const DateTime &maturity1,
	int nbFactor1,
	double *rho1,
	double x1,
	double y1,
	double z1,
	double *beta1,          //(I) mean reversion for energy curve         
	const DateTime &maturity2,
	int nbFactor2,
	double *rho2,
	double x2,
	double y2,
	double z2,
	double *beta2*/
	)
{
	static const string method = "SRMCorrelation::enrgToEnrg";

	int     i, j;
	/*
	vector<double> A1(nbFactor1);	// the epsilon in Dickinson's paper
	vector<double> A2(nbFactor2);
	vector<double> gamma1(nbFactor1);
	vector<double> gamma2(nbFactor2);
	DoubleMatrix OMtx1(nbFactor1, nbFactor1);	// the A in Dickinson's paper
	DoubleMatrix OMtx2(nbFactor2, nbFactor2);

	if (((nbFactor1 != 1) && (nbFactor1 != 2)) ||
		((nbFactor2 != 1) && (nbFactor2 != 2)) )
		throw ModelException(method, "Invalid number of energy curve "
		"factors. Number of energy factor can either be 1 or 2" );

	double length1 = today.yearFrac(maturity1);
	double length2 = today.yearFrac(maturity2);
	*/
	
	int nbFactor1 = enrgUtil1.numFactors();
	int nbFactor2 = enrgUtil2.numFactors();
	vector<double> A1(nbFactor1);
	vector<double> A2(nbFactor2);
	vector<double> gamma1;
	vector<double> gamma2;
	DoubleMatrix OMtx1;
	DoubleMatrix OMtx2;

	// compute the pre-normalized epsilon vectors
	gamma1 = enrgUtil1.getVectorGamma();
	gamma2 = enrgUtil2.getVectorGamma();
	switch (choice)	{
	case PC:
		/*
		if (Get2A(nbFactor1, x1, y1, z1, beta1, rho1, OMtx1) == FAILURE) {
			DR_Error("Get2A failed in %s\n", routine);
			goto RETURN;
		}

		switch (nbFactor1) {
		case 1:
			A1[0] = y1  + x1 * exp(-beta1[0] * length1);
			break;
		case 2:
			gamma1[0] = (y1  + x1 * exp(-beta1[0] * length1));
			gamma1[1] = z1 * exp(-beta1[1] * length1);
			A1[0] = gamma1[0] * OMtx1[0][0] + gamma1[1] * OMtx1[1][0];
			A1[1] = gamma1[0] * OMtx1[0][1] + gamma1[1] * OMtx1[1][1];
			break;
		}

		if (Get2A(nbFactor2, x2, y2, z2, beta2,	rho2, OMtx2) == FAILURE) {
		DR_Error("Get2A failed in %s\n", routine);
		goto RETURN;
		}

		switch (nbFactor2) {
		case 1:
			A2[0] = y2  + x2 * exp(-beta2[0] * length2);
			break;
		case 2:
			gamma2[0] = (y2  + x2 * exp(-beta2[0] * length2));
			gamma2[1] = z2 * exp(-beta2[1] * length2);
			A2[0] = gamma2[0] * OMtx2[0][0] + gamma2[1] * OMtx2[1][0];
			A2[1] = gamma2[0] * OMtx2[0][1] + gamma2[1] * OMtx2[1][1];
			break;
		}*/

		OMtx1 = enrgUtil1.getMatrixAByPCA();
		OMtx2 = enrgUtil2.getMatrixAByPCA();

		for (i = 0; i < nbFactor1; ++i) {
			for (j = 0; j < nbFactor1; ++j) {
				A1[i] += gamma1[j] * OMtx1[j][i];
			}
		}

		for (i = 0; i < nbFactor2; ++i) {
			for (j = 0; j < nbFactor2; ++j) {
				A2[i] += gamma2[j] * OMtx2[j][i];
			}
		}
		break;

	case EXPONENTIAL_FACTORS: 
	case TRIANGULAR:
		/*
		switch (nbFactor1) {
		case 1:
			A1[0] = y1 + x1 * exp(-beta1[0] * length1);
			break;
		case 2:
			A1[0] = y1 + x1 * exp(-beta1[0] * length1);
			A1[1] = z1;
			break;
		}

		switch (nbFactor2) {
		case 1:
			A2[0] = y2 + x2 * exp(-beta2[0] * length2);
			break;
		case 2:
			A2[0] = y2 + x2 * exp(-beta2[0] * length2);
			A2[1] = z2;
			break;
		}*/

		for (i = 0; i < nbFactor1; ++i) {
			A1[i] = gamma1[i]; 
		}

		for (i = 0; i < nbFactor2; ++i) {
			A2[i] = gamma2[i];
		}

		break;
	}     

	// normalize the epsilon vectors
	double sum = 0.0;
	for (i = 0; i < nbFactor1; i++)
		sum += A1[i] * A1[i];

	sum = sqrt(sum);
	if (sum < SRMConstants::SRM_TINY)
		throw ModelException(method, "Sum of the normalising weights is zero, unable to normalise!");

	for (i = 0; i < nbFactor1; i++)
		A1[i] /= sum;

	sum = 0.0;
	for (i = 0; i < nbFactor2; i++)
		sum += A2[i] * A2[i];

	sum = sqrt(sum);
	if (sum < SRMConstants::SRM_TINY)
		throw ModelException(method, "Sum of the normalising weights is zero, unable to normalise!");

	for (i = 0; i < nbFactor2; i++)
		A2[i] /= sum;

	// prepare to output results:
	vector<double> corrEnrgEnrg(nbFactor1 * nbFactor2);
	for (i = 0; i < nbFactor1; ++i) {
		for (j = 0; j < nbFactor2; ++j) {
			corrEnrgEnrg[i*nbFactor2 + j] = A1[i] * A2[j] * rhoEnrgEnrg;
		}
	}

	return (corrEnrgEnrg);
} // end enrgToEnrg


/***** From util_s.c::enrgToAsset ************************************/
//
//      Produces the implied corr structure between
//      Asset and the factors driving the energy curve.
// 
//      the factors can be:
//      choice 1 = EXPONENTIAL 
//      choice 2 = TRIANGULAR
//      choice 3 = PC
//
//      In case the energy dim is 1, rho can be input as NULL.
//      Rk: Asset is the equivalent of the FX rate in the doc
//
//      This function is Flat Forward ready
//
vector<double> SRMCorrelation::enrgToAsset(
	double rhoEnrgAsset,   /*(I) corr between Asset and swap rate          */
	OutputType choice,         /*(I) specifies factors wrt which correl is */
	const SRMEnergyUtilBase& enrgUtil
	/*const DateTime & today,           
	const DateTime & maturity,
	int nbFactor,
	double *rho,
	double x,
	double y,
	double z,
	double *beta*/          /*(I) mean reversion for energy curve           */
	)
{
	static const string method = "SRMCorrelation::enrgToAsset";

	int     i, j;
	/*
	vector<double> A(nbFactor);	// the epsilon in Dickinson's paper
	vector<double> gamma(nbFactor);
	DoubleMatrix OMtx(nbFactor, nbFactor); // the A in Dickinson's paper

	if ((nbFactor != 1) && (nbFactor != 2))
		throw ModelException(method, "Invalid number of energy curve "
										"factors: " + nbFactor);

	double length = today.yearFrac(maturity);
	*/

	int nbFactor = enrgUtil.numFactors();
	vector<double> A(nbFactor); // the epsilon vector in Dickinson's paper
	vector<double> gamma;	// gamma's in Dickinson's paper
	DoubleMatrix OMtx; // the A matrix in Dickinson's paper

	// compute the pre-normalized epsilon vector
	gamma = enrgUtil.getVectorGamma();
	switch (choice)	{
	case PC:
		/*
		//if (Get2A(NbFactor, x, y, z, beta, rho, OMtx) == FAILURE)
		//{
		//DR_Error("Get2A failed in %s\n", routine);
		//goto RETURN;
		//}
		switch (nbFactor) {
		case 1:
			A[0] = y  + x * exp(-beta[0] * length);
			break;
		case 2:
			gamma[0] = (y  + x * exp(-beta[0] * length));
			gamma[1] = z * exp(-beta[1] * length);
			A[0] = gamma[0] * OMtx[0][0] + gamma[1] * OMtx[1][0];
			A[1] = gamma[0] * OMtx[0][1] + gamma[1] * OMtx[1][1];
			break;
		}*/

		OMtx = enrgUtil.getMatrixAByPCA();
		for (i = 0; i < nbFactor; ++i) {
			for (j = 0; j < nbFactor; ++j) {
				A[i] += gamma[j] * OMtx[j][i];
			}
		}
		break;

	case EXPONENTIAL_FACTORS: 
	case TRIANGULAR:
		/*
		switch (nbFactor) {
		case 1:
			A[0] = y + x * exp(-beta[0] * length);
			break;
		case 2:
			A[0] = y + x * exp(-beta[0] * length);
			A[1] = z;
			break;
		}*/

		gamma = enrgUtil.getVectorGamma();
		for (i = 0; i < nbFactor; ++i) {
			A[i] = gamma[i];
		}

		break;
	}     

	// normalize the epsilon vector
	double sum = 0.0;
	for (i = 0; i < nbFactor; i++)
		sum += A[i] * A[i];

	sum = sqrt(sum);

	if (sum < SRMConstants::SRM_TINY)
		throw ModelException(method, "Sum of the normalising weights is zero, unable to normalise!");

	for (i = 0; i < nbFactor; i++)
		A[i] /= sum;

	// prepare to output results:
	vector<double> corrEnrgAsset(nbFactor);
	for (i = 0; i < nbFactor; i++)
		corrEnrgAsset[i] = A[i] * rhoEnrgAsset;

	return (corrEnrgAsset);
} // end enrgToAsset 

/***** From util_s.c::enrgToIr ************************************/
//
//      Produces the implied corr structure between
//      driving factors of IR and energy instruments.
// 
//      the factors can be:
//      choice 1 = EXPONENTIAL 
//      choice 2 = TRIANGULAR
//      choice 3 = PC
//
//      In case the energy dim is 1, enrgRho can be input as NULL.
//
//      This function is Flat Forward ready
//
vector<double> SRMCorrelation::enrgToIr(
	double rhoEnrgIr,   /*(I) corr between energy and Ir          */
	OutputType enChoice,
	OutputType irChoice,         /*(I) specifies factors wrt which correl is */
	/*const DateTime & enMaturity, // note: we use IR spot day as today 
	int nbEnFactor,
	double *enRho,
	double x,
	double y,
	double z,
	double *enBeta,*/          //(I) mean reversion for energy curve     
	const SRMEnergyUtilBase& enrgUtil,
	const SRMRatesUtil& ir,           //(I) IR util data
    bool enrgIsPrimaryAsset // true if Energy is the primary asset
                            // e.g. if enrgIsPrimaryAsset = true, then the 
                            // resulting correlation matrix will be:
                            //              IR:
                            //              factor1 factor2 factor3
                            // EN: factor1     x       x       x
                            //     factor2     x       x       x
                            //
                            // which will be flattened row-by-row into a 1D
                            // double vector for output.
                            //
                            // if enrgPrimary = false, then the resulting 
                            // matrix will be:
                            //              EN:
                            //              factor1 factor2
                            // IR: factor1     x       x
                            //     factor2     x       x
                            //     factor3     x       x
	)
{
	static const string method("SRMCorrelation::enrgToIr");
	try{
		/* MAR: In the SRM3 code the zero curve has today equal to the spot date
		so need to do all yearFracs etc from spot date */
		const DateTime& spotDate = ir.getDiffYC()->getSpotDate();
		int nbIrFactor = ir.numFactors(); // for ease
		// get the Ir rho's
		const vector<double>& irRho = ir.getRho();
		// get the 'B factor'
		vector<double> B(ir.bFactor());
		double irLength = spotDate.yearFrac(ir.getCorrSwapStart());
		double sum, cum;
		vector<double> Beta;
		ir.getBeta(Beta);
		vector<double> Alpha;
		ir.getAlpha(Alpha);
		int i, j, k, l;
		for (i = 0; i < nbIrFactor; i++) {          
			B[i] *= exp(- Beta[i] * irLength);
		}
		vector<double> A(nbIrFactor);
		DoubleMatrix Transform;
		/* calculate now normalised weights */
		switch (irChoice) {
		case EXPONENTIAL_FACTORS:
		case TRIANGULAR:
			{
				Transform = ir.triangulation();
				A[0] = Alpha[0] * B[0];
				if (nbIrFactor > 1) {                
					A[0] += Alpha[1] * B[1] * Transform[1][0];
					A[1]  = Alpha[1] * B[1] * Transform[1][1];
				}/* end of "nbIrFactor > 1"*/
				if (nbIrFactor > 2) {                
					A[0] += Alpha[2] * B[2] * Transform[2][0];
					A[1] += Alpha[2] * B[2] * Transform[2][1];
					A[2]  = Alpha[2] * B[2] * Transform[2][2];

				}/* end of nbIrFactor>2*/
				break;
			}
		case PC:
			{
				DoubleMatrix OMtx(ir.get3A());
				for (i = 0; i < nbIrFactor; i++) {
					for (j = 0; j < nbIrFactor; j++) {
						A[i] += Alpha[j] * B[j] * OMtx[j][i];
					} /* for j*/
				}/* for i*/
				break;
			}
		}/*end of switch on choice*/

		/* normalizing domestic weights*/
		sum = 0.0;
		for (i = 0; i < nbIrFactor; i++){
			sum += A[i] * A[i];
		}
		sum = sqrt(sum);
		if (sum < SRMConstants::SRM_TINY) {
			throw ModelException(method, "Sum of the normalising weights for "
				"IR is zero, unable to normalize!");
		}

		for (i = 0; i < nbIrFactor; i++)
			A[i] /= sum;


		// compute the pre-normalized epsilon vector for energy
		/*
		vector<double> EnA(nbEnFactor);	// the epsilon in Dickinson's paper
		vector<double> EnGamma(nbEnFactor);
		DoubleMatrix EnOMtx(nbEnFactor, nbEnFactor); // the A in Dickinson's paper

		if ((nbEnFactor != 1) && (nbEnFactor != 2))
			throw ModelException(method, "Invalid number of energy curve "
			"factors: " + nbEnFactor);

		double enLength = spotDate.yearFrac(enMaturity);
		*/

		int nbEnFactor = enrgUtil.numFactors();
		vector<double> EnA(nbEnFactor);
		vector<double> EnGamma;
		DoubleMatrix EnOMtx;

		// compute the pre-normalized epsilon vector
		EnGamma = enrgUtil.getVectorGamma();
		switch (enChoice)	{
		case PC:
			/*
			if (Get2A(nbEnFactor, x, y,	z, enBeta, enRho, EnOMtx) == FAILURE) {
				DR_Error("Get2A failed in %s\n", routine);
				goto RETURN;
			}
			
			switch (nbEnFactor) {
			case 1:
				EnA[0] = y  + x * exp(-enBeta[0] * enLength);
				break;
			case 2:
				EnGamma[0] = (y  + x * exp(-enBeta[0] * enLength));
				EnGamma[1] = z * exp(-enBeta[1] * enLength);
				EnA[0] = EnGamma[0] * EnOMtx[0][0] + EnGamma[1] * EnOMtx[1][0];
				EnA[1] = EnGamma[0] * EnOMtx[0][1] + EnGamma[1] * EnOMtx[1][1];
				break;
			}*/

			EnOMtx = enrgUtil.getMatrixAByPCA();
			for (i = 0; i < nbEnFactor; ++i) {
				for (j = 0; j < nbEnFactor; ++j) {
					EnA[i] += EnGamma[j] * EnOMtx[j][i];
				}
			}
			break;

		case EXPONENTIAL_FACTORS: 
		case TRIANGULAR:
			/*
			switch (nbEnFactor) {
			case 1:
				EnA[0] = y + x * exp(-enBeta[0] * enLength);
				break;
			case 2:
				EnA[0] = y + x * exp(-enBeta[0] * enLength);
				EnA[1] = z;
				break;
			}*/

			for (i = 0; i < nbEnFactor; ++i) {
				EnA[i] = EnGamma[i];
			}
			break;
		}     

		// normalize the epsilon vector
		sum = 0.0;
		for (i = 0; i < nbEnFactor; i++)
			sum += EnA[i] * EnA[i];

		sum = sqrt(sum);

		if (sum < SRMConstants::SRM_TINY)
			throw ModelException(method, "Sum of the normalising weights for "
				"energy is zero, unable to normalize!");

		for (i = 0; i < nbEnFactor; i++)
			EnA[i] /= sum;
		
        /*
        // dummyCorrIrEnrg are the corr mentioned in doc
		vector<double> dummyCorrEnrgIr(nbEnFactor * nbIrFactor);
		for (i = 0; i < nbEnFactor; ++i) {
			for (j = 0; j < nbIrFactor; ++j) {
				dummyCorrEnrgIr[i*nbIrFactor + j] = EnA[i] * A[j] * rhoEnrgIr;
			}
		}
        */

        // dummyCorrIrEnrg are the corr mentioned in doc
        vector<double> dummyCorrEnrgIr;
        DoubleMatrix dummyCorrMatrix(nbIrFactor/*nbColumn*/, nbEnFactor/*nbRow*/);

        for (i = 0; i < nbEnFactor; ++i) {
            for (j = 0; j < nbIrFactor; ++j) {
                dummyCorrMatrix[j][i] = EnA[i] * A[j] * rhoEnrgIr;
            }
        }

        // check if EN is the primary asset
        if (!enrgIsPrimaryAsset)
            dummyCorrMatrix.transpose(); // transpose the matrix if IR is primary

        // flatten the corr matrix row-by-row into a 1D vector 
        for (i = 0; i < dummyCorrMatrix.numRows(); ++i) {
            for (j = 0; j < dummyCorrMatrix.numCols(); ++j) {
                dummyCorrEnrgIr.push_back(dummyCorrMatrix[j][i]);
            }
        }

		if (irChoice == EXPONENTIAL_FACTORS) 
		{   
            /*
			vector<double> corrEnrgIr(nbEnFactor * nbIrFactor); // return value
			for (i = 0; i < nbEnFactor; i++) {
				for (j = 0; j < nbIrFactor; j++) {
					sum = 0.;
					for (k = 0; k < nbEnFactor; k++) {
						cum = 0.;
						for (l = 0; l < nbIrFactor; l++) {
							cum += Transform[j][l] * 
								dummyCorrEnrgIr[l + k * nbIrFactor];
						} // for l 
						sum += cum;
					} // for k 
					corrEnrgIr[j + i * nbIrFactor] = sum;
				} // for j 
			} // for i 
            */

            vector<double> corrEnrgIr; // return value
            DoubleMatrix corrMatrix(nbIrFactor/*nbColumns*/, nbEnFactor/*nbRows*/);
            for (i = 0; i < nbEnFactor; i++) {
                for (j = 0; j < nbIrFactor; j++) {
                    sum = 0.;
                    for (k = 0; k < nbEnFactor; k++) {
                        cum = 0.;
                        for (l = 0; l < nbIrFactor; l++) {
                            cum += Transform[j][l] * 
                                dummyCorrEnrgIr[l + k * nbIrFactor];
                        } // for l 
                        sum += cum;
                    } // for k 
                    //corrEnrgIr[j + i * nbIrFactor] = sum;
                    corrMatrix[j][i] = sum;
                } // for j 
            } // for i 

			/* now test whether we recover input correlation */
			double test_IrVol = 0.;

			for (i = 0; i < nbIrFactor; i++) {
				test_IrVol += Alpha[i] * Alpha[i]* B[i]* B[i];
			}

			if (nbIrFactor  > 1) {
				test_IrVol += 2.0 * Alpha[0] * Alpha[1] * 
					B[0] * B[1]* irRho[0];
			}

			if (nbIrFactor > 2)	{
				test_IrVol += 2.0 * Alpha[0] * Alpha[2] * 
					B[0] * B[2] * irRho[1];
				test_IrVol += 2.0 * Alpha[1] * Alpha[2] * 
					B[1] * B[2] * irRho[2];
			}

			if (test_IrVol < SRMConstants::SRM_TINY)
				throw ModelException(method, "Problem calculating IR test Vol "
					"for an energy-Ir pair ");

			/* Corr For versus Dom*/
			double test_rhoEnrgIr = 0.;
			for (i = 0; i < nbIrFactor; i++)
			{                
				/*
				// corresponds to the 1st energy factor
				test_rhoEnrgIr += Alpha[i] * B[i] * 
					(x + y * exp(-enBeta[0] * enLength)) * corrEnrgIr[i]; 
				if (nbEnFactor > 1)
					// corresponds to the 2nd energy factor
					test_rhoEnrgIr += Alpha[i] * B[i] * z *
						exp(-enBeta[1] * enLength) * corrEnrgIr[1 * nbIrFactor + i];
						*/
                // corresponds to the 1st energy factor
                test_rhoEnrgIr += Alpha[i] * B[i] * EnGamma[0] * corrMatrix[i][0]; //test_rhoEnrgIr += Alpha[i] * B[i] * EnGamma[0] * corrEnrgIr[i]; 
                if (nbEnFactor > 1)
                    // corresponds to the 2nd energy factor
                    test_rhoEnrgIr += Alpha[i] * B[i] * EnGamma[1] * corrMatrix[i][1]; //test_rhoEnrgIr += Alpha[i] * B[i] * EnGamma[1] * corrEnrgIr[1*nbIrFactor + i];

			}/* for i*/

			test_rhoEnrgIr /= sqrt(test_IrVol);

			if ((fabs (test_rhoEnrgIr - rhoEnrgIr) > SRMConstants::SRM_TINY))
				throw ModelException(method, "Failed to reproduce input "
					"correlation coefficient for an energy-ir pair");

            // check if EN is the primary asset
            if (!enrgIsPrimaryAsset)
                corrMatrix.transpose(); // transpose the matrix if IR is primary

            // flatten the corr matrix to a 1D vector 
            for (i = 0; i < corrMatrix.numRows(); ++i) {
                for (j = 0; j < corrMatrix.numCols(); ++j) {
                    corrEnrgIr.push_back(corrMatrix[j][i]);
                }
            }

			return corrEnrgIr;
		}/* end of choice=1*/
		else {
			/* the correlations in this case have already been obtained */
			return dummyCorrEnrgIr;
		}
	} catch (exception& e){
		throw ModelException(e, method, "Failed for energy-IR pair: ccy "+
			ir.getCcy());
	}

} // end enrgToIr

DRLIB_END_NAMESPACE
