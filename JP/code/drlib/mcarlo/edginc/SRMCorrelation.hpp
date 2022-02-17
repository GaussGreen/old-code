//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMCorrelation.hpp
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
#ifndef EDR_SRMCORRELATION_HPP
#define EDR_SRMCORRELATION_HPP

#include "edginc/SRMRatesUtil.hpp"
#include "edginc/SRMEnergyUtil.hpp"

DRLIB_BEGIN_NAMESPACE

class SRMGenDiffusibleIR;
class SRMGenDiffusibleFX;
class SRMGenDiffusibleEquity;
class SRMGenDiffusibleCredit;
class SRMGenDiffusibleBasisSpread;
class SRMGenDiffusibleEnergy;

/** Helper for SRM - used for splitting up correlation between
    the different factors */
class MCARLO_DLL SRMCorrelation{
public:
    //// How irToIr returns its results
    enum OutputType {
        EXPONENTIAL_FACTORS = 1,
        TRIANGULAR,
        PC
    };

    // IR vs IR corr split
    static vector<double> ExpandAssetAssetCorrs(
        SRMGenDiffusibleIR* irAsset1, 
        SRMGenDiffusibleIR* irAsset2, 
        double corr, 
        OutputType choice);

    // IR vs Energy corr split
    static vector<double> ExpandAssetAssetCorrs(
        SRMGenDiffusibleIR* irAsset, 
        SRMGenDiffusibleEnergy* enAsset,
        double corr, 
        OutputType irChoice, 
        OutputType enChoice);

    // 1-factor vs IR corr split
    static vector<double> ExpandAssetAssetCorrs(
        SRMGenDiffusibleIR* irAsset, 
        double corr, 
        OutputType choice);

    // 1-factor vs Energy corr split
    static vector<double> ExpandAssetAssetCorrs(
        SRMGenDiffusibleEnergy* enAsset, 
        double corr, 
        OutputType choice);

    // Energy vs Energy corr split
    static vector<double> ExpandAssetAssetCorrs(
        SRMGenDiffusibleEnergy* enAsset1, 
        SRMGenDiffusibleEnergy* enAsset2, 
        double corr, 
        OutputType choice);

    // Energy vs IR corr split
    static vector<double> ExpandAssetAssetCorrs(
        SRMGenDiffusibleEnergy* enAsset, 
        SRMGenDiffusibleIR* irAsset, 
        double corr, 
        OutputType enChoice, 
        OutputType irChoice);


    /***** From util_s.c:IR_IR_Corr ******************************************
     *
     *    Produces the implied correlation structure between YC driving factors
     *    ie the correlation between "dom" and "for" factors  
     *
     *   This function is Flat Forward ready
     *   Really should return a DoubleMatrix but that requires too many changes
     */
    static vector<double> irToIr( 
        double         RhoDswapFswap,  //(I) corr. between "dom" and "for" rates
        const SRMRatesUtil& dom,            //(I) 'domestic' IR data
        const SRMRatesUtil& fgn,            //(I) 'foreign' IR data
        OutputType     choice);      /*(I) specifies factors wrt which correl is
                                       output: 1=expon.,2=triang;3=PC      */

    /***** From util_s.c::AssetIrCorr ************************************/
    //
    // Produces the implied corr structure between
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
    static vector<double> assetIr(
        double         RhoAssetSwap, /*(I) corr between Asset and swap rate */
        const SRMRatesUtil& ir,           //(I) IR data
        OutputType     choice);      /*(I) specifies factors wrt which correl is
                                       output: 1=expon.,2=triang;3=PC      */

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
	static vector<double> enrgToEnrg(
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
		);

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
	static vector<double> enrgToIr(
		double rhoEnrgIr, /*(I) corr b/w energy and ir */
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
        bool enrgIsPrimaryAsset // true if energy is the primary asset 
                                // (more comments on this flag in cpp)
		);

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
	static vector<double> enrgToAsset(
		double rhoEnrgAsset,   /*(I) corr between energy and asset */
		OutputType choice,         /*(I) specifies factors wrt which correl is */
		const SRMEnergyUtilBase& enrgUtil
		/*
		const DateTime & today,           
		const DateTime & maturity,
		int nbFactor,
		double *rho,
		double x,
		double y,
		double z,
		double *beta*/          /*(I) mean reversion for energy curve           */
		);
};

DRLIB_END_NAMESPACE
#endif
