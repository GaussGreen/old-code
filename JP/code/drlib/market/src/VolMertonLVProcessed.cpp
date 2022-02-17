//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolMertonLVProcessed.cpp
//
//   Description : implementation for VolMertonLVProcessed object
//
//   Date        : 30 Oct 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolMertonLVProcessed.hpp"
#include "edginc/Black.hpp"

DRLIB_BEGIN_NAMESPACE

// debug testing
//#define PRINT_DEBUG
#ifdef PRINT_DEBUG
    FILE * file1 = fopen("C:\\temp\\file1.txt","w");
    FILE * file2 = fopen("C:\\temp\\file2.txt","w");
    FILE * file3 = fopen("C:\\temp\\file3.txt","w");
    FILE * file4 = fopen("C:\\temp\\file4.txt","w");
    FILE * file5 = fopen("C:\\temp\\file5.txt","w");
    FILE * file6 = fopen("C:\\temp\\file6.txt","w");
    FILE * file7 = fopen("C:\\temp\\file7.txt","w");
#endif
// debug testing

// util functions

static double priceMerton (double yrsToMat,
                           bool isCall,
                           double strike,
                           double forward,
                           double pv,
                           double ATMVol,
                           double JumpRate,
                           double JumpMean,
                           double JumpWidth) {
    const double tolerance=0.001;///if the fraction of price added is less than tolerance time the already estimated price, we stop. 
    const int maxIteration=50;///Maximum of iterations. 
    double MertonPrice=0;
    double increment;
    double factor=1;
    double varianceModification = JumpWidth*JumpWidth;
    double fwrdModification = JumpRate*yrsToMat*(1-exp(JumpMean+0.5*varianceModification));
    int i=0;
    do{
        increment=factor*
            Black::price(isCall, 
                         forward*exp(fwrdModification+i*(JumpMean+0.5*varianceModification)),
                         strike, pv, ATMVol*ATMVol*yrsToMat+i*varianceModification);	
			
        MertonPrice+=increment;
        i++;
        factor*=JumpRate*yrsToMat/i;
    } while(increment>tolerance*MertonPrice && (i<maxIteration));

    return exp(-JumpRate*yrsToMat)*MertonPrice;
}

double N1DensityModified(double  StandardVariate, double Substraction)
{
//    result = exp(-(StandardVariate*StandardVariate-Substraction*Substraction)/2.0);
    double  result=Substraction*Substraction-StandardVariate*StandardVariate;
	result *=0.5;
	result=exp(result);
	return result;
}

/////////// for VolProcessedMertonLVParam class

VolMertonLVProcessed::VolMertonLVProcessed(): CVolProcessedDVFParam(TYPE){}


VolMertonLVProcessed::VolMertonLVProcessed(const CVolBase*         vol,
                          const CVolParamConstSP&  volParam,
                          const CVolRequestDVF*    volRequest, 
                          const CAsset*            asset,
                          const VolSurface*        VolSurf,
						  VolMertonLVSP			   mertonLV) :

	CVolProcessedDVFParam (vol,
                      volParam,
                      volRequest, 
                      asset,
                      VolSurf),
					myVolMertonLV(mertonLV){}
/*
bool VolMertonLVProcessed::computeLocV2(double  yrsToMat,
                                          double  strike,
                                          double  forward,
                                          double  growthRate,
                                          SImpV&  impV,
                                          double  &v2) const
{
    static const string routine = "computeLocV2";

    if (yrsToMat < SMALL_YEARS_TO_MATURITY){
   
           v2  = impV.vol;
           v2 *= v2;
            return true;
    }
    if (impV.vol<0.0005){
        v2 = 0.0;
//        return true;
    }

  	static double aTMVol = getMertonParam("ATMVol");
    static double jumpRate =getMertonParam("JumpRate");
    static double jumpMean =getMertonParam("JumpMean");
    static double jumpWidth =getMertonParam("JumpWidth");
	static double jumpWidthSquare=jumpWidth*jumpWidth;
	static double jumpRateprime=jumpRate*exp(jumpMean+0.5*jumpWidthSquare);
	
	double JumpRateMat=jumpRate*yrsToMat;
	double JumpRatePrimeMat=jumpRateprime*yrsToMat;

    double rootYTM  = sqrt(yrsToMat);
    double stdDev   = impV.vol * rootYTM;
	double impVd_dStrike_Square=impV.d_dStrike*impV.d_dStrike;
	double notUsed;
    DateTime maturity= GetTimeMetric()->impliedTime(startDate, yrsToMat, notUsed);

	double numerator = impV.d_dT + (jumpRate-jumpRateprime+growthRate) * strike *  impV.d_dStrike;
    numerator *= 2.0 * yrsToMat;
    numerator += impV.vol;
    numerator *= impV.vol;
	
	double precision=0.01;
	int precisionMax = 15;
	
	vector<double> An(precisionMax);
	if (JumpRateMat!=0)	An[0]=exp(-JumpRatePrimeMat);
	else An[0]=1.0;
	int i=0;
	do{		
		i++;
		An[i]=An[i-1]*JumpRatePrimeMat/((double)i);
	} while((An[i]>precision*An[i-1]) && (i<precisionMax));
	int numSteps=i+1;

	double d0, d01, vn, dn, temp, numeratorLoc, denominatorLoc, tempNumeratorDenominator, alphan, tempepsilon;

	d0=log(forward / strike)+(JumpRateMat-JumpRatePrimeMat);
	d01=d0/stdDev+0.5*stdDev;

	numeratorLoc=An[0]*N1Density(d01);

	temp=d01*impV.d_dStrike+1/(strike*rootYTM);
	denominatorLoc=numeratorLoc*(impV.d2_dStrike2-impVd_dStrike_Square*d01*rootYTM+temp*temp/impV.vol);

	for(i = 1; i <numSteps; i++){
		vn=sqrt(stdDev*stdDev+i*jumpWidthSquare);
		dn=(d0+i*(jumpMean+0.5*jumpWidthSquare))/vn+0.5*vn;
		alphan=1/(1+i*jumpWidthSquare/(stdDev*stdDev));
		tempNumeratorDenominator=sqrt(alphan)*An[i]*N1Density(dn);
		numeratorLoc+=tempNumeratorDenominator;

		temp=(alphan*dn*impV.d_dStrike+1/(strike*rootYTM));
		temp=impV.d2_dStrike2+impVd_dStrike_Square*((1-alphan)/impV.vol-dn*sqrt(alphan)*rootYTM)+temp*temp/impV.vol;
  		denominatorLoc+=tempNumeratorDenominator*temp;
	}
	denominatorLoc*=strike*strike*impV.vol*yrsToMat;

	double g=aTMVol;
	double pv=equityYC->pv(maturity);
	tempepsilon= priceMerton(yrsToMat,
							true,
							strike,
							forward,
							pv,
							impV.vol,
							jumpRate,
							jumpMean,
							jumpWidth);
	/////////////////////////////////////CONVOLUTION////////////////////////////////////////////////////////
	double strikemin=-3.0;
	double strikemax=-strikemin;
	int precisionMax2=100;
	vector<double> wStrike(precisionMax2);
	double deltastrikeconv=(strikemax-strikemin)/((double)precisionMax2);
	double adjustment=0;

	for(i=0;i<precisionMax2;i++){
		wStrike[i]=strikemin+i*deltastrikeconv;
		wStrike[i]=exp(-0.5*wStrike[i]*wStrike[i])/Maths::ROOT_TWO_PI;
		adjustment+=wStrike[i];
	}

	double epsilon;
	double count=0;
	double impliedVolConvolution;
	double strikeConvolution;
	double strikeCorrection;
	double strikeCorrection0=exp(-jumpMean-jumpWidth*strikemin);

	temp=0;
	for(i=0;i<precisionMax2;i++){
		strikeCorrection=strikeCorrection0*exp(-i*jumpWidth*deltastrikeconv);

		strikeConvolution=strike*strikeCorrection;
		impliedVolConvolution = computeImpVol(maturity, strikeConvolution);

		epsilon=priceMerton(yrsToMat,
							true,
							strikeConvolution,
							forward,
							pv,
							impliedVolConvolution,
							jumpRate,
							jumpMean,
							jumpWidth);
		temp+=wStrike[i]*epsilon/strikeCorrection;
	}

	numeratorLoc=numeratorLoc*numerator
+2*impV.vol*jumpRate*rootYTM*(tempepsilon*adjustment-temp)*deltastrikeconv/forward;

	double v = numeratorLoc / denominatorLoc;
	static double matref=yrsToMat;

#ifdef PRINT_DEBUG
        if (matref==yrsToMat)
        {
            fprintf(file1, "%f,", impV.vol);
            fprintf(file2, "%f,", impV.d_dT);
            fprintf(file3, "%f,", impV.d_dStrike);
            fprintf(file4, "%f,", impV.d2_dStrike2);
            fprintf(file5, "%f,", numerator);
            fprintf(file6, "%f,", numeratorLoc);
            fprintf(file7, "%f,", denominatorLoc);
        }
		else
		{
			matref=yrsToMat;
            fprintf(file1, "\n");
            fprintf(file2, "\n");
            fprintf(file3, "\n");
            fprintf(file4, "\n");
            fprintf(file5, "\n");
            fprintf(file6, "\n");
            fprintf(file7, "\n");
            fprintf(file1, "%f,", impV.vol);
            fprintf(file2, "%f,", impV.d_dT);
            fprintf(file3, "%f,", impV.d_dStrike);
            fprintf(file4, "%f,", impV.d2_dStrike2);
            fprintf(file5, "%f,", numerator);
            fprintf(file6, "%f,", numeratorLoc);
            fprintf(file7, "%f,", denominatorLoc);

		}
#endif

	if (!Maths::isPositive(v)){
	return false;
	}
	v2 = v;
    return true;
}*/
/*give correct prices*///////////////////

bool VolMertonLVProcessed::computeLocV2(double  yrsToMat,
                                          double  strike,
                                          double  forward,
                                          double  growthRate,
                                          SImpV&  impV,
                                          double  &v2) const
{
    static const string routine = "computeLocV2";

    if (yrsToMat < SMALL_YEARS_TO_MATURITY){
   
        v2  = impV.vol;
        v2 *= v2;
        return true;
    }
    if (Maths::isZero(impV.vol)){
        v2 = 0.0;
        return true;
    }

    double aTMVol = getMertonParam("ATMVol");
    double jumpRate =getMertonParam("JumpRate");
    double jumpMean =getMertonParam("JumpMean");
    double jumpWidth =getMertonParam("JumpWidth");
    double jumpWidthSquare=jumpWidth*jumpWidth;
    double jumpRateprime=jumpRate*exp(jumpMean+0.5*jumpWidthSquare);
	
    double JumpRateMat=jumpRate*yrsToMat;
    double JumpRatePrimeMat=jumpRateprime*yrsToMat;

    double rootYTM  = sqrt(yrsToMat);
    double stdDev   = impV.vol * rootYTM;
    double impVd_dStrike_Square=impV.d_dStrike*impV.d_dStrike;
    double notUsed;
    DateTime maturity= GetTimeMetric()->impliedTime(startDate, yrsToMat, notUsed);

    double numerator = impV.d_dT + (jumpRate-jumpRateprime+growthRate) * strike *  impV.d_dStrike;
    numerator *= 2.0 * yrsToMat;
    numerator += impV.vol;
    numerator *= impV.vol;
	
//	int precisionMax = 15;
    int precisionMax = 10;
	
    vector<double> An(precisionMax);
    if (JumpRateMat!=0)	An[0]=exp(-JumpRatePrimeMat);
    else An[0]=1.0;
    int i=0;
/*	do{		
		i++;
		An[i]=An[i-1]*JumpRatePrimeMat/((double)i);
	} while((An[i]>precision*An[i-1]) && (i<precisionMax));
	int numSteps=i+1;*/
    int numSteps=precisionMax;
    for(i=1;i<precisionMax;i++) An[i]= An[i-1]*JumpRatePrimeMat/((double)i);

    double d0, d01, vn, dn, temp, numeratorLoc, denominatorLoc, tempNumeratorDenominator, alphan, tempepsilon;

    d0=log(forward / strike)+(JumpRateMat-JumpRatePrimeMat);
    d01=d0/stdDev+0.5*stdDev;

    numeratorLoc=An[0]*N1Density(d01);

    temp=d01*impV.d_dStrike+1/(strike*rootYTM);
    denominatorLoc=numeratorLoc*(impV.d2_dStrike2-impVd_dStrike_Square*d01*rootYTM+temp*temp/impV.vol);

    for(i = 1; i <numSteps; i++){
        vn=sqrt(stdDev*stdDev+i*jumpWidthSquare);
        dn=(d0+i*(jumpMean+0.5*jumpWidthSquare))/vn+0.5*vn;
        alphan=1/(1+i*jumpWidthSquare/(stdDev*stdDev));
        tempNumeratorDenominator=sqrt(alphan)*An[i]*N1Density(dn);
        numeratorLoc+=tempNumeratorDenominator;

        temp=(alphan*dn*impV.d_dStrike+1/(strike*rootYTM));
        temp=impV.d2_dStrike2+impVd_dStrike_Square*((1-alphan)/impV.vol-dn*sqrt(alphan)*rootYTM)+temp*temp/impV.vol;
        denominatorLoc+=tempNumeratorDenominator*temp;
    }
    denominatorLoc*=strike*strike*impV.vol*yrsToMat;

    double g=aTMVol;
    double pv=equityYC->pv(maturity);
    tempepsilon= priceMerton(yrsToMat,
                             true,
                             strike*exp(-jumpMean-0.5*jumpWidthSquare),
                             forward,
                             pv,
                             sqrt(impV.vol*impV.vol+jumpWidthSquare/yrsToMat),
                             jumpRate,
                             jumpMean,
                             jumpWidth);
    tempepsilon-= priceMerton(yrsToMat,
                              true,
                              strike*exp(-jumpMean-0.5*jumpWidthSquare),
                              forward,
                              pv,
                              sqrt(g*g+jumpWidthSquare/yrsToMat),
                              jumpRate,
                              jumpMean,
                              jumpWidth);
    /////////////////////////////////////CONVOLUTION////////////////////////////////////////////////////////
    double strikemin=-4.0;
    double strikemax=-strikemin;
    int precisionMax2=10;
    vector<double> wStrike(precisionMax2);
    double deltastrikeconv=(strikemax-strikemin)/((double)precisionMax2);
    double adjustment=0;

    for(i=0;i<precisionMax2;i++){
        wStrike[i]=strikemin+i*deltastrikeconv;
        wStrike[i]=exp(-0.5*wStrike[i]*wStrike[i])/Maths::ROOT_TWO_PI;
        adjustment+=wStrike[i];
    }
    ////I am using Closed Newton Cotes formula or order 10.//////
    vector<double> weight(precisionMax2);
    double ratio=9.0/89600.0;
    weight[0]=2857.0*ratio;
    weight[1]=15741.0*ratio;
    weight[2]=1080.0*ratio;
    weight[3]=19344.0*ratio;
    weight[4]=5788.0*ratio;
    weight[5]=5788.0*ratio;
    weight[6]=19344.0*ratio;
    weight[7]=1080.0*ratio;
    weight[8]=15741.0*ratio;
    weight[9]=2857.0*ratio;

    double epsilon;
    double position;
	double impliedVolConvolution;
	double strikeConvolution;

	position=strike*exp(-jumpMean
		-jumpWidth*(jumpWidth+strikemin));

	temp=0;
	for(i=0;i<precisionMax2;i++){

		strikeConvolution=position*exp(-i*jumpWidth*deltastrikeconv);
        impliedVolConvolution = computeImpVol(maturity, strikeConvolution);
		if (Maths::isZero(impliedVolConvolution )){
			impliedVolConvolution = 0.0;
		}

		epsilon=priceMerton(yrsToMat,
							true,
							strikeConvolution,
							forward,
							pv,
							impliedVolConvolution,
							jumpRate,
							jumpMean,
							jumpWidth);
		epsilon-=priceMerton(yrsToMat,
							true,
							strikeConvolution,
							forward,
							pv,
							g,
							jumpRate,
							jumpMean,
							jumpWidth);
		temp+=wStrike[i]*epsilon*weight[i];
	}

	numeratorLoc=numeratorLoc*numerator
+2*impV.vol*JumpRatePrimeMat*(tempepsilon*adjustment-temp)*deltastrikeconv/rootYTM/(forward*pv);

	double v = numeratorLoc / denominatorLoc;
	if (!Maths::isPositive(v)){
	return false;
	}
	v2 = v;
    return true;
}


/*
bool VolMertonLVProcessed::computeLocV2(double  yrsToMat,
                                          double  strike,
                                          double  forward,
                                          double  growthRate,
                                          SImpV&  impV,
                                          double  &v2) const
{
    static const string routine = "computeLocV2";

    if (yrsToMat < SMALL_YEARS_TO_MATURITY){
   
           v2  = impV.vol;
           v2 *= v2;
            return true;
    }
    if (Maths::isZero(impV.vol)){
        v2 = 0.0;
        return true;
    }

    static double jumpRate =getMertonParam("JumpRate");
    static double jumpMean =getMertonParam("JumpMean");
    static double jumpWidth =getMertonParam("JumpWidth");
	static double jumpWidthSquare=jumpWidth*jumpWidth;
	static double jumpRateprime=jumpRate*exp(jumpMean+0.5*jumpWidthSquare);
	
    double rootYTM  = sqrt(yrsToMat);
    double stdDev   = impV.vol * rootYTM;
	double impVd_dStrike_Square=impV.d_dStrike*impV.d_dStrike;
	
    double x = log(forward / strike);//+(jumpRate-jumpRateprime)*yrsToMat;
	x /= stdDev;
	x += 0.5 * stdDev;

    double probDensRatio = impV.vol * impV.d2_dStrike2;
    probDensRatio += x * (x - stdDev) * impV.d_dStrike * impV.d_dStrike;
    probDensRatio *= strike * yrsToMat;
    probDensRatio += 2.0 * x * rootYTM * impV.d_dStrike;
    probDensRatio *= strike;
    probDensRatio += 1.0;

	if (!Maths::isPositive(probDensRatio)){
		return false;
	}

    double numerator = impV.d_dT + (jumpRate-jumpRateprime+growthRate) * strike *  impV.d_dStrike;
    numerator *= 2.0 * yrsToMat;
    numerator += impV.vol;
    numerator *= impV.vol;

    double v = numerator / probDensRatio;

	if (!Maths::isPositive(v)){
	return false;
	}
	v2 = v;
    return true;
}
*/
/*
bool VolMertonLVProcessed::computeLocV2(double  yrsToMat,
                                          double  strike,
                                          double  forward,
                                          double  growthRate,
                                          SImpV&  impV,
                                          double  &v2) const
{
    static const string routine = "computeLocV2";
    if (yrsToMat < SMALL_YEARS_TO_MATURITY){
           v2  = impV.vol;
           v2 *= v2;
            return true;
    
    }
    if (!Maths::isPositive(impV.vol)){
        v2 = 0.0;
        return true;
    }
    double rootYTM  = sqrt(yrsToMat);
    double stdDev   = impV.vol * rootYTM;
    static double jumpRate =getMertonParam("JumpRate");
    static double jumpMean =getMertonParam("JumpMean");
    static double jumpWidth =getMertonParam("JumpWidth");

	static double jumpWidthSquare=jumpWidth*jumpWidth;
	static double jumpRateprime = jumpRate*exp(jumpMean+0.5*jumpWidthSquare);

    if (!Maths::isPositive(stdDev)){
		return false;
    }

    double x = log(forward / strike) / stdDev + 0.5 * stdDev;

    double probDensRatio = impV.vol * impV.d2_dStrike2;
    probDensRatio += x * (x - stdDev) * impV.d_dStrike * impV.d_dStrike;
    probDensRatio *= strike * yrsToMat;
    probDensRatio += 2.0 * x * rootYTM * impV.d_dStrike;
    probDensRatio *= strike;
    probDensRatio += 1.0;
	if (!Maths::isPositive(probDensRatio)){
		return false;
	}
      
    double numerator = impV.d_dT + (jumpRate-jumpRateprime+growthRate) * strike *  impV.d_dStrike;
    numerator *= 2.0 * yrsToMat;
    numerator += impV.vol;
    numerator *= impV.vol;

    double v = numerator / probDensRatio;
    
    if (!Maths::isPositive(v)){
        return false;
    }

    v2 = v;
    return true;
}*/
//////////////// end of VolProcessedMertonLVParam class

CClassConstSP const VolMertonLVProcessed::TYPE = CClass::registerClassLoadMethod(
    "VolMertonLVProcessed", typeid(VolMertonLVProcessed), load);


void VolMertonLVProcessed::load(CClassSP& clazz){
        REGISTER(VolMertonLVProcessed, clazz);
        SUPERCLASS(CVolProcessedDVFParam);
        EMPTY_SHELL_METHOD(defaultVolMertonLVProcessed);
        FIELD(myVolMertonLV, "");
        FIELD_MAKE_TRANSIENT(myVolMertonLV);
        FIELD(asset, "");
        FIELD_MAKE_TRANSIENT(asset);
        FIELD(myParamVol, "");
        FIELD_MAKE_TRANSIENT(myParamVol);
        FIELD(equityYC, "");
        FIELD_MAKE_TRANSIENT(equityYC);
        FIELD(startDate, "");
        FIELD_MAKE_TRANSIENT(startDate);
}

double VolMertonLVProcessed::getMertonParam(const string& param_name) const
{
    static char routine[] = "VolMertonLVProcessed::InterpParam";
    if (param_name == "ATMVol"){
        return myVolMertonLV->ATMVol;
    }
    else if (param_name == "JumpRate"){
        return myVolMertonLV->JumpRate;
    }
    else if (param_name == "JumpMean"){
        return myVolMertonLV->JumpMean;
    }
    else if (param_name == "JumpWidth"){
        return myVolMertonLV->JumpWidth;
    }
    else
        throw ModelException(routine, "Unknown parameter name");
}

VolSurfaceSP VolMertonLVProcessed::getMertonImpliedSurface() const
{
    return myVolMertonLV->mertonSurface;
}

///////////from Olivier Brockhaus Merton Model///
void VolMertonLVProcessed::Quantile(
    const DateTime &date1,
    const DateTime &date2,
    double         epsilon,
    int            maxJumps,
    double*        quantile,
    int*           nbJumps ) const
{

    double dt       = myVolMertonLV->timeMetric->yearFrac(date1,date2);
    double proba    = exp(-dt* myVolMertonLV->JumpRate);
    *quantile = 0;
    int    iStep;

    for( iStep=1; iStep<=maxJumps; iStep++ )
    {
        *quantile += proba;
        proba *= myVolMertonLV->JumpRate * dt / (double)iStep;
        if( *quantile>1-epsilon ) break;
    }
    *nbJumps = iStep;
}

/** Calculates jump factor. */
double VolMertonLVProcessed::CalcJump(	double noise,
										int    numJumps) const
{
    double amplitudeJump = 0;
    if( numJumps>0 ) {
        amplitudeJump +=(double)numJumps * myVolMertonLV->JumpMean
			+ sqrt((double)numJumps)*(myVolMertonLV->JumpWidth)
			* noise;
    }
    return amplitudeJump;
}


DRLIB_END_NAMESPACE
