#include "edginc/config.hpp"
#include "edginc/VnfmTestTime.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/IrConverter.hpp"


DRLIB_BEGIN_NAMESPACE

// static variables
static const double M       = 30.;//30.;      // Principal factors are orthogonal over [0,M]
static const int NUM_TO_M   = 600;//600;      // Compute equidistantly the P.C's at this numer of points  
static const double TINY    = 1.e-6;    // Cut off to avoid numerical 0/0 calculations


/***** OPTIMIZER *********************************************************************************/
// Class for objective function //
class ObjFuncT: public MFunctionND{
public:

    virtual void operator()(const CDoubleArray&  x,
                            CDoubleArray&        f) const
    { f = vnfmPlus.rmse_swap_vol(x); }

    ObjFuncT(const RangeArray&   defRanges,
            int                 nbVars,
            int                 nbFuncs,
            VnfmTPlus&           vnfmPlus,
            string              model):
    MFunctionND(nbVars, nbFuncs, defRanges),vnfmPlus(vnfmPlus),model(model){}

protected:
    VnfmTPlus&   vnfmPlus;                   // access to vnfm class
    string      model;                      // which model do we use
};

/*************************************************************************************************/


/***** NATIVE TYPES *******************************************************************************/
// We must re-construct native types for the reflection to work
class VnfmTDoubleArray : public CObject {
public:
    static CClassConstSP const TYPE;
    VnfmTDoubleArray() :  CObject(TYPE) {}

    CDoubleMatrixSP	param;
private:
    IObjectSP myCalc_PC3toFix3(){ return IObjectSP(VnfmT::pC3toFix3(*param)); }
    IObjectSP myCalc_Fix3toPC3(){ return IObjectSP(VnfmT::fix3toPC3(*param)); }
    IObjectSP myCalc_PC()		{ return IObjectSP(VnfmT::pComponents(*param)); }
    IObjectSP myCalc_Hump()		{ return IObjectSP(VnfmT::hump(*param)); }

    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new VnfmTDoubleArray(); }
};

CClassConstSP const VnfmTDoubleArray::TYPE = CClass::registerClassLoadMethod(
    "VnfmTDoubleArray", typeid(VnfmTDoubleArray), load);

void VnfmTDoubleArray::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("My description");
    REGISTER(VnfmTDoubleArray, clazz);
    SUPERCLASS(CObject);
    FIELD(param,"");
    EMPTY_SHELL_METHOD(defaultConstructor);

    Addin::registerObjectMethod(
            "PC3_TO_Fix3_T", Addin::UTILITIES,
            "Help message",
            false, Addin::expandSimple,
            &VnfmTDoubleArray::myCalc_PC3toFix3);

    Addin::registerObjectMethod(
            "Fix3_TO_PC3_T", Addin::UTILITIES,
            "Help message",
            false, Addin::expandSimple,
            &VnfmTDoubleArray::myCalc_Fix3toPC3);

    Addin::registerObjectMethod(
            "P_COMPONENTS_T", Addin::UTILITIES,
            "Help message",
            false, Addin::expandSimple,
            &VnfmTDoubleArray::myCalc_PC); 

    Addin::registerObjectMethod(
            "TIME_HUMP_T", Addin::UTILITIES,
            "Help message",
            false, Addin::expandSimple,
            &VnfmTDoubleArray::myCalc_Hump); 
}

/*************************************************************************************************/


/*************************************  VNFM  ****************************************************/
// Excel interface.
void VnfmT::load(CClassSP& clazz) {
    clazz->setPublic();                         // make visible to EAS/spreadsheet
    clazz->setDescription("My description");
    REGISTER(VnfmT, clazz);
    SUPERCLASS(CObject);                        // Base class
    FIELD(pC3param,"PC3 parameters");
    FIELD(daycount,"Day count convention");
    FIELD(frequency,"Swap frequency");
    FIELD(reference_tenor,"bootstrap tenor (in years)"); FIELD_MAKE_OPTIONAL(reference_tenor);
    FIELD(irvol,"Swap vol object");
    FIELD(yc,"Yield curve object");
//	FIELD(sv_data, "sv_data");
    FIELD(yc_data, "yc_data");
    FIELD(spotvol, "spotvol");
//	FIELD_MAKE_TRANSIENT(sv_data);
    FIELD_MAKE_TRANSIENT(yc_data);
    FIELD_MAKE_TRANSIENT(spotvol);
    EMPTY_SHELL_METHOD(defaultConstructor);

    Addin::registerObjectMethod(
            "SWAP_VOL_T", Addin::UTILITIES,
            "Help message",
            false, Addin::expandSimple,
            &VnfmT::getSwapVol);    

    Addin::registerObjectMethod(
            "SPOT_VOL_T", Addin::UTILITIES,
            "Help message",
            false, Addin::expandSimple,
            &VnfmT::getSpotVol);    
}

CClassConstSP const VnfmT::TYPE = CClass::registerClassLoadMethod("VnfmT", typeid(VnfmT), load);

// constructor
VnfmT::VnfmT(CDoubleMatrix &param, DayCountConventionSP dcc, int freq, CIntSP ref_tenor, 
           IRVolSP ir, YieldCurveSP yc, CClassConstSP const &type):
CObject(type),
pC3param(param),
daycount(dcc),
frequency(freq),
reference_tenor(ref_tenor)
{
    validatePop2Object();
}

// Need to use this function for the constructor to work both with reflection and normal construction 
void VnfmT::validatePop2Object() 
{
	try
	{
		// local variables for function call
		ExpiryArray selectedTenors;
		ExpiryArray selectedExpiries; 

		// create user friendly market data container
		IrConverter::to_SWAPVOL_DATA(sv_data, selectedTenors, selectedExpiries, irvol.get());   
		set_YIELDCURVE_DATA(yc.getSP());

		spotvol.reset(new DoubleArray(sv_data.NbSwaptionExpiries, 1.));

		if(reference_tenor.get()!=NULL) { set_spotvol(reference_tenor->intValue());}
	}
    catch(exception& e) { throw ModelException(e, "VnfmT::validatePop2Object"); }	
}

// populate simple yield curve container
void VnfmT::set_YIELDCURVE_DATA(YieldCurveSP yc)
{
	try
	{
		CashFlowArraySP ratesAndDates = yc->getRatesAndDates();     

		int row_size = ratesAndDates->size(); 
		yc_data.reset(new CDoubleMatrix(3, row_size));

		// need to investigate how to properly use day count conventions
		DateTime today   = yc->getToday();
		(*yc_data)[0][0] = daycount->years(today,(*ratesAndDates)[0].date);   
		(*yc_data)[1][0] = yc->pv((*ratesAndDates)[0].date);
		(*yc_data)[2][0] = -log ( (*yc_data)[1][0] ) / daycount->years( today , (*ratesAndDates)[0].date );

		for(int i = 1; i < row_size; i++)
		{
			(*yc_data)[0][i] = daycount->years(today,(*ratesAndDates)[i].date);   
			(*yc_data)[1][i] = yc->pv((*ratesAndDates)[i].date);
			(*yc_data)[2][i] = -log ( (*yc_data)[1][i] / (*yc_data)[1][i-1] ) /
								daycount->years( (*ratesAndDates)[i-1].date , (*ratesAndDates)[i].date );
		}

	}
    catch(exception& e) { throw ModelException(e, "VnfmT::set_YIELDCURVE_DATA"); }	
}

// see tecnical document for description of this factor.
double VnfmT::b_factor(int expiry, int tenor, int index) // expiry and tenor are assumed to be in month
{
	// need to put in proper day count fractions for fixed leg
	// however right now we don't know the actual pay dates
	try
	{
		bool not_ready = true;
		int counter = 0;

		while(not_ready) 
		{
			++counter;
			not_ready = 12 * (*yc_data)[0][counter] < expiry;
		}
	    
		// expiry lies between counter-1 and counter
		double DF1 = (*yc_data)[1][counter]; double DF2 = (*yc_data)[1][counter-1];
		double TF1 = (*yc_data)[0][counter]; double TF2 = (*yc_data)[0][counter-1];

		// interpolate in log space
		double kk  = log(DF1 / DF2) / (TF1 - TF2);	
		double DF = DF2 * exp( kk * (expiry/12. - TF2) );

		int k=-1; double a=0.; double nom=0.; double denom=0.; 
		double DF_i; double t2 = expiry/12.;

		// mean reversions are not allowed to be time dependent
		double mr=pC3param[index][0]; 

		for(int i=0; i<frequency*tenor/12; i++)
		{
			double maturity = expiry/12. + (i+1.)/frequency;
			not_ready = true;

			while(not_ready)
			{
				k++;
				double t1 = min(maturity , (*yc_data)[0][counter+k]);
				if(fabs(mr) > TINY)
					a -= (*yc_data)[2][counter+k]*( exp(-mr*( t1 - expiry/12. ))
													-exp(-mr*( t2 - expiry/12. )) ) / mr;
				else
					a += (*yc_data)[2][counter+k]*(t1-t2);

				not_ready = (*yc_data)[0][counter+k] < maturity;
				t2=t1;
			}       
			
			//reset k for next iteration
			k=k-1;

			// expiry lies between counter+k-1 and counter+k
			DF1 = (*yc_data)[1][counter+k]; DF2 = (*yc_data)[1][counter+k-1];
			TF1 = (*yc_data)[0][counter+k]; TF2 = (*yc_data)[0][counter+k-1];

			// interpolate in log space
			kk  = log(DF1 / DF2) / (TF1 - TF2);	
			DF_i = DF2 * exp( kk * (maturity - TF2) );

			nom     += DF_i / DF * a;
			denom   += DF_i / DF;
		}

		return (nom/denom + a * DF_i / DF / (1. - DF_i / DF));
    }
    catch(exception& e) { throw ModelException(e, "VnfmT::b_factor"); }
}

// calculates model swaption volatilities
double VnfmT::swap_vol(int expiry, int tenor)
{
	try{

    double b1 = b_factor(expiry, tenor, 0);
    double b2 = b_factor(expiry, tenor, 1);

    CDoubleMatrixSP fix3param = pC3toFix3(pC3param);

	// mean reversions are not allowed to be time dependent
    double mr1  = (*fix3param)[0][0]; double mr2   = (*fix3param)[1][0]; double out=0.;

    double j11, j12, j22;
    double t1 = 0; double t2 = sv_data.SwaptionExpiries[0]/12.; double t3 = expiry/12.; 
    for(int i=0; sv_data.SwaptionExpiries[i]<=expiry && i<sv_data.NbSwaptionExpiries; i++)
    {
		double w1 = (*fix3param)[2][i]; double w2 = (*fix3param)[3][i]; double rho = (*fix3param)[4][i];
        j11 = ( (fabs(2 * mr1) > TINY) ? (exp(-2*mr1*(t3-t2)) - exp(-2*mr1*(t3-t1)) ) / (2 * mr1) : t2-t1 )*w1*w1;
        j22 = ( (fabs(2 * mr2) > TINY) ? (exp(-2*mr2*(t3-t2)) - exp(-2*mr2*(t3-t1)) ) / (2 * mr2) : t2-t1 )*w2*w2;
        j12 = ( (fabs(mr1+mr2) > TINY) ? (exp(-(mr1+mr2)*(t3-t2)) - exp(-(mr1+mr2)*(t3-t1)) ) / (mr1+mr2) : t2-t1 )*w1*w2*rho;

        out += ::pow( (*spotvol)[i] , 2 ) * (b1 * b1 * j11 + 2 * b1 * b2 * j12 + b2 * b2 * j22);

        t1 = t2;
        t2 = sv_data.SwaptionExpiries[i+1]/12.;
    }

    return (sqrt(out / expiry * 12) );

    }
    catch(exception& e) { throw ModelException(e, "VnfmT::swap_vol"); }
}

// Calculates the spot volatility so that we exactly hit the market swaption vols for a certain tenor
void VnfmT::set_spotvol(int ref_tenor) // ref_tenor is supposed to be in years
{
	try{

    CDoubleMatrixSP   fix3param = pC3toFix3(pC3param);

    int ref=0;
    for(int i=0; i<sv_data.NbSwapTenors; i++)
    {
        if(sv_data.SwapTenors[i] < 12*ref_tenor) { ref++; }
    }

    for(int i=0; i<sv_data.NbSwaptionExpiries; i++)
    {
        int expiry  = sv_data.SwaptionExpiries[i];
        double b1 = b_factor(expiry, 12*ref_tenor, 0);
        double b2 = b_factor(expiry, 12*ref_tenor, 1);
        double j11, j12, j22;

        double t1=0.; double k1=0.;
        double t2=sv_data.SwaptionExpiries[0]/12.;
        double t3=sv_data.SwaptionExpiries[i]/12.;
        for(int j = 0; j<=i; j++)
        {
			double mr1  = (*fix3param)[0][j]; double mr2   = (*fix3param)[1][j]; 
			double w1   = (*fix3param)[2][j]; double w2    = (*fix3param)[3][j];
			double rho  = (*fix3param)[4][j];

            j11 = ( (fabs(2 * mr1) > TINY) ? (exp(-2*mr1*(t3-t2)) - exp(-2*mr1*(t3-t1)) ) / (2 * mr1) : t2-t1 )*w1*w1;
            j22 = ( (fabs(2 * mr2) > TINY) ? (exp(-2*mr2*(t3-t2)) - exp(-2*mr2*(t3-t1)) ) / (2 * mr2) : t2-t1 )*w2*w2;
            j12 = ( (fabs(mr1+mr2) > TINY) ? (exp(-(mr1+mr2)*(t3-t2)) - exp(-(mr1+mr2)*(t3-t1)) ) / (mr1+mr2) : t2-t1 )*w1*w2*rho;

            if(j<i) { k1 += (*spotvol)[j]*(*spotvol)[j]*(b1*b1*j11 + 2*b1*b2*j12 + b2*b2*j22); }

            t1 = t2;
            t2 = sv_data.SwaptionExpiries[j+1]/12.;
        }

        double k2 = (b1*b1*j11 + 2*b1*b2*j12 + b2*b2*j22);
        if( ::pow( sv_data.VolMatrix[i][ref] , 2) * t3 - k1 > 0 )
            (*spotvol)[i] = sqrt( ( ::pow( sv_data.VolMatrix[i][ref] , 2 ) * t3 - k1 ) / k2 );
        else
            (*spotvol)[i] = 0.;

    }

    }
    catch(exception& e) { throw ModelException(e, "VnfmT::set_spotvol"); }
}

// Populate the model swaption vol grid. This is a public function.
CDoubleMatrixSP VnfmT::getSwapVol()
{
    try
    {
		// local variables for function call
		ExpiryArray selectedTenors, selectedExpiries; 
		IrConverter::to_SWAPVOL_DATA(sv_data, selectedTenors, selectedExpiries, irvol.get());  

        int row_size = sv_data.NbSwaptionExpiries; 
        int col_size = sv_data.NbSwapTenors; 

        CDoubleMatrixSP swapVol(new CDoubleMatrix(col_size, row_size));

        for(int i = 0; i < row_size; i++)
            for(int j = 0; j < col_size; j++)
                (*swapVol)[j][i] = swap_vol(sv_data.SwaptionExpiries[i],sv_data.SwapTenors[j]);

        return swapVol;
    }
    catch(exception& e) { throw ModelException(e, "VnfmT::getSwapVol"); }
}

// Returns the calibrated spot volatility. This is a public function.
CDoubleArraySP VnfmT::getSpotVol() { return spotvol; }

// Calculate the principal components. Public function
CDoubleMatrixSP VnfmT::pComponents(CDoubleMatrix &pC3param)    // (I) vector of PC3 parameters
{
    try
    {
        double inc = M / NUM_TO_M;
        // Given a smart pointer to a matrix, then for some reason Excel outputs
        // the transpose of the matrix. Hence if we want to output a nxm matrix 
        // we need to deliver a smart pointer to a mxn matrix.
        CDoubleMatrixSP pC(new CDoubleMatrix(3,NUM_TO_M));

        double a1       = pC3param[0][0];
        double a2       = pC3param[1][0];
        double angle    = pC3param[4][0];

        CDoubleMatrixSP rK = setrKmatrix(a1, a2, angle);

        for(int i = 0; i<NUM_TO_M; i++)
        {
            (*pC)[0][i] = inc*(i+1);
            (*pC)[1][i] = ( (*rK)[0][0] * exp( -a1*inc*(i+1) ) + (*rK)[0][1] * exp( -a2*inc*(i+1) ) );
            (*pC)[2][i] = ( (*rK)[1][0] * exp( -a1*inc*(i+1) ) + (*rK)[1][1] * exp( -a2*inc*(i+1) ) );
        }

        return pC;
    }
    catch(exception& e) { throw ModelException(e, "VnfmT::pComponents"); }
}

// Calculate the hump location for the 1st principal components. Public function
CDoubleArraySP VnfmT::hump(CDoubleMatrix &pC3param)    // (I) vector of PC3 parameters
{
    try
    {
        CDoubleArraySP time(new CDoubleArray(1));

        double a1       = pC3param[0][0];
        double a2       = pC3param[1][0];
        double angle    = pC3param[4][0];

        CDoubleMatrixSP rK = setrKmatrix(a1, a2, angle);

		double r_00 = (*rK)[0][0];
		double r_01 = (*rK)[0][1];

		double t=0.0;
		if(fabs(a1*r_00)>TINY && fabs(a2*r_01)>TINY && (a1*r_00*a2*r_01)<TINY)
		{
			t = log(( -a2 * r_01 ) / ( a1 * r_00 )) / ( a2 - a1 ) ;
		}
		else
		{
			if(r_00*exp(-a1*M)+r_01*exp(-a2*M)>r_00+r_01) { t = M; }	
		}

		(*time)[0] = t;

        return time;
    }
    catch(exception& e) { throw ModelException(e, "VnfmT::hump"); }
}

// set Fix3 parameters from PC3 parameters. Public function
CDoubleMatrixSP VnfmT::pC3toFix3(CDoubleMatrix &pC3param) // (I) vector of PC3 parameters 
{
    try 
    {
        if (pC3param.numCols()!=5) throw ModelException("PC3param.numCols() should be 5");

		int rows = pC3param.numRows();
        CDoubleMatrixSP fix3param(new CDoubleMatrix(5, rows));
        
		for(int i=0; i<rows; i++)
		{
			double a1       = pC3param[0][i];
			double a2       = pC3param[1][i];
			double c1		= pC3param[2][i];
			double c2		= pC3param[3][i];
			double angle    = pC3param[4][i];

	        CDoubleMatrixSP rK = setrKmatrix(a1, a2, angle);

			(*fix3param)[0][i] = a1;
			(*fix3param)[1][i] = a2;
			(*fix3param)[2][i] = sqrt( ::pow( c1*(*rK)[0][0], 2) + ::pow( c2*(*rK)[1][0], 2) );
			(*fix3param)[3][i] = sqrt( ::pow( c1*(*rK)[0][1], 2) + ::pow( c2*(*rK)[1][1], 2) );

			// Need to fix this line so that we never divide by zero
			(*fix3param)[4][i] = ( c1*c1 * (*rK)[0][1] * (*rK)[0][0] + c2*c2 * (*rK)[1][0] * (*rK)[1][1] ) 
								/ (*fix3param)[2][i] / (*fix3param)[3][i];
		}

        return fix3param;
    }
    catch(exception& e) { throw ModelException(e, "VnfmT::pC3toFix3"); }
}

// set Fix3 parameters from PC3 parameters. Public function.
CDoubleMatrixSP VnfmT::fix3toPC3(CDoubleMatrix &fix3param2) // (I) vector of Fix3 parameters
{
    try 
    {
        if (fix3param2.numCols()!=5) throw ModelException("fix3param2.numCols() should be 5");

		DoubleArray fix3param(5); 

		int rows = fix3param2.numRows();
        CDoubleMatrixSP pC3param(new CDoubleMatrix(5,rows));


		for(int i=0; i<rows; i++)
		{
			fix3param[0]=fix3param2[0][i]; fix3param[1]=fix3param2[1][i];
			fix3param[2]=fix3param2[2][i]; fix3param[3]=fix3param2[3][i]; fix3param[4]=fix3param2[4][i];
			// setup numerical equation solver
			double angle; double low=-0.1; double high=0.1; bool not_ready=true; int count=0;
			if(i>0) {low = (*pC3param)[4][i-1]-0.1; high = (*pC3param)[4][i-1]+0.1;}
			
			while(not_ready)
			{
				count++;
				not_ready = ((findAngle(low,&fix3param)*findAngle(high,&fix3param) >= 0.) && (count < 1000));
				if(not_ready)
				{
					low     -= 0.1;
					high    += 0.1;
				}
			}
	    
			try 
			{
				angle = zbrentUseful(
							findAngle,              // (I) The function to find the root of
							&fix3param,             // (I) class input
							low,                    // (I) Low value for x
							high,                   // (I) High value for x
							1.e-10);                // (I) Tolerance
			} 
			catch (exception&) 
			{
				throw ModelException(   // e, // hide "zbrent" message which confuses users
										"Brent algorithm failed");
			}       

			double beta1   = fix3param[0];
			double beta2   = fix3param[1];
			double alpha1  = fix3param[2];
			double alpha2  = fix3param[3];
			CDoubleMatrixSP rK = setrKmatrix(beta1, beta2, angle);

			(*pC3param)[0][i] = beta1;
			(*pC3param)[1][i] = beta2;
			(*pC3param)[2][i] = sqrt(( ::pow((*rK)[1][1] * alpha1 , 2) - ::pow((*rK)[1][0] * alpha2 , 2) ) 
								/ ( ::pow((*rK)[0][0] * (*rK)[1][1] , 2) - ::pow((*rK)[1][0] * (*rK)[0][1] , 2) ));
			(*pC3param)[3][i] = sqrt(( ::pow((*rK)[0][0] * alpha2 , 2) - ::pow((*rK)[0][1] * alpha1 , 2) ) 
								/ ( ::pow((*rK)[0][0] * (*rK)[1][1] , 2) - ::pow((*rK)[1][0] * (*rK)[0][1] , 2) ));
			(*pC3param)[4][i] = angle;
		}

        return pC3param;
    }
    catch(exception& e) { throw ModelException(e, "VnfmT::fix3toPC3"); }
}

// Creates the rK Matrix (see technical document). 
CDoubleMatrixSP VnfmT::setrKmatrix(double a1,        // (I) mean reversion for first factor 
                                    double a2,      // (I) mean reversion for second factor
                                    double angle)   // (I) angle in rotation matrix
{
    try
    {       
        CDoubleMatrixSP rK(new DoubleMatrix(2,2));

        // This is the covariance matrix
        double c_00 = (a1 + a1 > TINY) ? (1 - exp(-(a1 + a1) * M)) / (a1 + a1) : M;
        double c_11 = (a2 + a2 > TINY) ? (1 - exp(-(a2 + a2) * M)) / (a2 + a2) : M;
        double c_01 = (a1 + a2 > TINY) ? (1 - exp(-(a1 + a2) * M)) / (a1 + a2) : M;
        double c_10 = c_01;

        // These are the eigen values
        double lambda_0 = (c_00 + c_11) / 2. + sqrt(::pow(c_00 + c_11, 2) / 4. - (c_00 * c_11 - c_01 * c_10));
        double lambda_1 = (c_00 + c_11) / 2. - sqrt(::pow(c_00 + c_11, 2) / 4. - (c_00 * c_11 - c_01 * c_10));

        // These are the eigen vectors
        double u_00 = 1 / sqrt(1 + ::pow( (c_00 - lambda_0) / c_01 , 2) );
        double u_01 = (c_00 - lambda_0) / c_01 / sqrt(1 + ::pow( (c_00 - lambda_0) / c_01 , 2) );
        double u_11 = u_00;
        double u_10 = -u_01;

        // This is the fundamental solution
        double k_00 = u_00 / sqrt(lambda_0);
        double k_01 = u_10 / sqrt(lambda_0);
        double k_10 = u_01 / sqrt(lambda_1);
        double k_11 = u_11 / sqrt(lambda_1);

        // This is the rotation matrix
        double r_00 = cos(angle);
        double r_01 = sin(angle);
        double r_10 = -sin(angle);
        double r_11 = cos(angle);

        // This is the general solution
        (*rK)[0][0] = r_00 * k_00 + r_01 * k_10;
        (*rK)[0][1] = r_00 * k_01 + r_01 * k_11;
        (*rK)[1][0] = r_10 * k_00 + r_11 * k_10;
        (*rK)[1][1] = r_10 * k_01 + r_11 * k_11;

        return rK;
    }
    catch(exception& e) { throw ModelException(e, "VnfmT::setrKmatrix"); }
}
        
// find the angle in the PC3 representation that corresponds to given Fix3 paramters. 
double VnfmT::findAngle(double angle, void* data) // (I) angle in PC3 parameters
{
    DoubleArray &fix3param(*((DoubleArray*)data));
    try 
    {
        double beta1        = fix3param[0];
        double beta2        = fix3param[1];
        double alpha1       = fix3param[2];
        double alpha2       = fix3param[3];
        double rho          = fix3param[4];

        CDoubleMatrixSP rK = setrKmatrix(beta1, beta2, angle);

        double result = alpha1*alpha1 * (*rK)[0][1]*(*rK)[1][1] + alpha2*alpha2 * (*rK)[0][0]*(*rK)[1][0] 
                        - rho*alpha1*alpha2 * ((*rK)[0][0]*(*rK)[1][1] + (*rK)[1][0]*(*rK)[0][1]);

        return result;
    }
    catch(exception& e) { throw ModelException(e, "VnfmT::findAngle"); }
}

/*************************************************************************************************/


/***** VNFM PLUS *********************************************************************************/
// Excel interface.
void VnfmTPlus::load(CClassSP& clazz) {
    clazz->setPublic();                         // make visible to EAS/spreadsheet
    clazz->setDescription("My description");
    REGISTER(VnfmTPlus, clazz);
    SUPERCLASS(VnfmT);                           // base class
    FIELD(weightMatrix,"Weight matrix: swap vols");
    FIELD(tenorPairs,"Tenor pairs (in years)");
    FIELD(marketCorrelation,"Market correlation");
    FIELD(calibrateCorrelation,"Boolean vector");
    FIELD(calibrator,"Quasi-Newton or Levenberg-Marquardt"); FIELD_MAKE_OPTIONAL(calibrator);
    FIELD(humptime,"Location (in years)"); FIELD_MAKE_OPTIONAL(humptime);
    EMPTY_SHELL_METHOD(defaultConstructor);

    Addin::registerObjectMethod(
            "CORRELATION_T", Addin::UTILITIES,
            "Help message",
            false, Addin::expandSimple,
            &VnfmTPlus::getCorrelation);    

    Addin::registerObjectMethod(
            "CALIBRATE_PC3_T", Addin::UTILITIES,
            "Help message",
            false, Addin::expandSimple,
            &VnfmTPlus::getPC3Calibration); 
}

CClassConstSP const VnfmTPlus::TYPE = CClass::registerClassLoadMethod("VnfmTPlus", typeid(VnfmTPlus), load);

// constructor
VnfmTPlus::VnfmTPlus(CDoubleMatrix &param, DayCountConventionSP dcc, int freq, CIntSP ref_tenor, 
                   IRVolSP ir, YieldCurveSP yc, CDoubleMatrix &weights, CDoubleMatrix &tp,
                   CDoubleArray &mc, CBoolArray &cc, CDoubleSP hp, CClassConstSP const &type) :
VnfmT(param,dcc,freq,ref_tenor,ir,yc,TYPE),
weightMatrix(weights),
tenorPairs(tp),
marketCorrelation(mc),
calibrateCorrelation(cc),
humptime(hp)
{
}

// calculates root-mean-square-error between model and market volatilities
DoubleArray VnfmTPlus::rmse_swap_vol(const CDoubleArray&  x)
{
    try
    {       
        DoubleArray rmse(nbFunc, 1.0e+10);

        double mr1 = x[0]; double mr2 = x[1];
        // set pC3 parameters
        bool valid_corr = true;
        if(nbVar == 2 || nbVar == 3 )
        {
            // we are doing the PC3 calibration
			for(int i=0; i<sv_data.NbSwaptionExpiries; i++)
			{
				pC3param[0][i]=x[0];
				pC3param[1][i]=x[1];
				if(nbVar ==3)
				{
					pC3param[4][i]=x[2];
				}
			}

            // set factor weights
            set_factorWeights();
        }
        else
        {
            throw ModelException("invalid model");
        }
        
        int row_size = weightMatrix.numRows(); 
        int col_size = weightMatrix.numCols(); 
        double model_squared_vol; double market_squared_vol; int cc = 0; rmse[0]=0.;
        for(int i = 0; i < row_size; i++)
            for(int j = 0; j < col_size; j++)
            {
                market_squared_vol = ::pow(sv_data.VolMatrix[i][j] , 2); 
                model_squared_vol  = ::pow(swap_vol(sv_data.SwaptionExpiries[i],sv_data.SwapTenors[j]) , 2);

                if(weightMatrix[j][i]>TINY)
                    if(nbFunc==1)
					{
                        rmse[0] += ::pow(model_squared_vol - market_squared_vol, 2) * weightMatrix[j][i];
					}
                    else
                    {
                        rmse[cc] = ( model_squared_vol - market_squared_vol ) * sqrt(weightMatrix[j][i]);
                        cc++;
                    }
            }       

        // make sure constraints are satisfied
        double a1       = pC3param[0][0];
        double a2       = pC3param[1][0];

		CDoubleMatrixSP fix3param = pC3toFix3(pC3param);
		double c1=pC3param[2][0]; double c2=pC3param[3][0]; double rho=fabs((*fix3param)[4][0]);
		for(int i=0; i<sv_data.NbSwaptionExpiries; i++)
		{
			c1       = min(c1,pC3param[2][i]);
			c2       = min(c2,pC3param[3][i]);
			rho      = max(rho,fabs((*fix3param)[4][i]));
		}

        // values below might be changed if needed
		bool valid_rho	= rho <= 0.95;
		bool valid_mr	= fabs(a2-a1) >= 0.01;
		bool valid_fact	= (c2 > TINY) && (c1 > TINY);
		bool valid_hump = true;
		if(humptime.get()!= NULL) 
		{ 
			valid_hump = fabs((*hump(pC3param))[0] - humptime->doubleValue()) < TINY; 
		}
        bool c_ok = valid_fact && valid_rho && valid_mr && valid_corr && valid_hump;

        if(!c_ok)
        {
            // force back to good region
            for(cc = 0; cc < nbFunc; cc++)
                rmse[cc] = 1.0e+10;
        }

        return (rmse);

    }
    catch(exception& e) { throw ModelException(e, "VnfmT::rmse_swap_vol_PC3"); }
}

// Calculates the model correlation between two swap rates. This is a public function.
double VnfmTPlus::correlation(int expiry, int tenor1, int tenor2)
{
	try{

    double b1_1 = b_factor(expiry, tenor1, 0);
    double b2_1 = b_factor(expiry, tenor1, 1);
    double b1_2 = b_factor(expiry, tenor2, 0);
    double b2_2 = b_factor(expiry, tenor2, 1);

    CDoubleMatrixSP   fix3param = pC3toFix3(pC3param);
	// mean reversion parameters are not allowed to be time dependent
    double mr1  = (*fix3param)[0][0];  double mr2  = (*fix3param)[1][0]; 

	double cov = 0.; double var1 = 0.; double var2 = 0.;
    double j11, j12, j22;

	double t1 = 0; double t2 = sv_data.SwaptionExpiries[0]/12.; double t3 = expiry/12.; 
    for(int i=0; sv_data.SwaptionExpiries[i]<=expiry && i<sv_data.NbSwaptionExpiries; i++)
    {
		double w1 = (*fix3param)[2][i]; double w2 = (*fix3param)[3][i]; double rho = (*fix3param)[4][i];

        j11 = ((fabs(2 * mr1) > TINY) ? (exp(-2*mr1*(t3-t2)) - exp(-2*mr1*(t3-t1)) )/(2*mr1) : t2-t1 )*w1*w1;
        j22 = ((fabs(2 * mr2) > TINY) ? (exp(-2*mr2*(t3-t2)) - exp(-2*mr2*(t3-t1)) )/(2*mr2) : t2-t1 )*w2*w2;
        j12 = ((fabs(mr1+mr2) > TINY) ? (exp(-(mr1+mr2)*(t3-t2)) - exp(-(mr1+mr2)*(t3-t1)) ) 
                                                                                / (mr1+mr2) : t2-t1 )*w1*w2*rho;
		double sv2	= ::pow((*spotvol)[i],2);

        cov += sv2 * (b1_1 * b1_2 * j11 + (b1_1 * b2_2 + b1_2 * b2_1) * j12 + b2_1 * b2_2 * j22);
        var1+= sv2 * (b1_1 * b1_1 * j11 + 2 * b1_1 * b2_1 * j12 + b2_1 * b2_1 * j22);
        var2+= sv2 * (b1_2 * b1_2 * j11 + 2 * b1_2 * b2_2 * j12 + b2_2 * b2_2 * j22);

        t1 = t2;
        t2 = sv_data.SwaptionExpiries[i+1]/12.;
    }

    return (cov/sqrt(var1*var2));   

    }
    catch(exception& e) { throw ModelException(e, "VnfmTPlus::correlation"); }
}

// see technical document for definition of this function. 
double VnfmTPlus::g_function(int start, int end, int expiry, int tenor1, int tenor2, int index)
{
	try
	{
		double b1_1 = b_factor(expiry, tenor1, 0);
		double b2_1 = b_factor(expiry, tenor1, 1);
		double b1_2 = b_factor(expiry, tenor2, 0);
		double b2_2 = b_factor(expiry, tenor2, 1);

		double j11 = j_function(start, end, expiry, 0, 0, index);
		double j12 = j_function(start, end, expiry, 0, 1, index);
		double j22 = j_function(start, end, expiry, 1, 1, index);

		double out = (b1_1 * b1_2 * j11 + (b1_1 * b2_2 + b1_2 * b2_1) * j12 + b2_1 * b2_2 * j22);

		return ( out / expiry * 12. );  
    }
    catch(exception& e) { throw ModelException(e, "VnfmTPlus::g_function"); }
}

// see technical document for definition of this function. 
double VnfmTPlus::j_function(int start, int end, int expiry, int n, int m, int index)
{
	try
	{
		// mean reversion parameters are not allowe to be time dependent
		double mr1  = pC3param[n][0];  double mr2  = pC3param[m][0]; 
		double j=0.;
		
		if((reference_tenor.get()==NULL)&&(index==0))
		{
			// flat parameters
			j = ( (fabs(mr1+mr2) > TINY) ? (1 - exp(-(mr1+mr2) * expiry / 12)) / (mr1+mr2) : expiry );
		}
		else
		{
			// curve parameters
			int offset = 0;
			while( start > sv_data.SwaptionExpiries[offset] ) { offset++; }

			double t1 = start/12.; double t2 = sv_data.SwaptionExpiries[offset]/12.; double t3 = expiry/12.; 
			for(int i=offset; sv_data.SwaptionExpiries[i]<=end && i<sv_data.NbSwaptionExpiries; i++)
			{
				double vol2= ::pow((pC3param)[2+index][i],2) * ::pow((*spotvol)[i],2);
				j += vol2 * ( (fabs(mr1+mr2) > TINY) ? 
							(exp(-(mr1+mr2)*(t3-t2)) - exp(-(mr1+mr2)*(t3-t1)) ) / (mr1+mr2) : t2-t1 );

				t1 = t2; t2 = sv_data.SwaptionExpiries[i+1]/12.;
			}
		}

		return ( j * (*rKmatrix)[index][n] * (*rKmatrix)[index][m] );   

    }
    catch(exception& e) { throw ModelException(e, "VnfmTPlus::j_function"); }
}

// calculates the weights of each principal factor. 
void VnfmTPlus::set_factorWeights()
{
	try
	{
		// given mean reversion values and the rotation angles, set the RK matrix
		update_rKmatrix();

		double q=0;
		// no spotvol bootstrapping
		if(reference_tenor.get()==NULL)
		{
			int start = 0;
			int end   = 0;
			int j =0;
			for(int i=0; i<sv_data.NbSwaptionExpiries; i++)
			{
				pC3param[2][i]=1;
				pC3param[3][i]=1;
				if(calibrateCorrelation[i])
				{
					int expiry = sv_data.SwaptionExpiries[i];
	
					double g1_01 = g_function(0,expiry,expiry,12*(int)tenorPairs[0][i],12*(int)tenorPairs[1][i],0);
					double g1_00 = g_function(0,expiry,expiry,12*(int)tenorPairs[0][i],12*(int)tenorPairs[0][i],0);
					double g1_11 = g_function(0,expiry,expiry,12*(int)tenorPairs[1][i],12*(int)tenorPairs[1][i],0);

					double g2_01 = g_function(0,end,expiry,12*(int)tenorPairs[0][i],12*(int)tenorPairs[1][i],1);
					double g2_00 = g_function(0,end,expiry,12*(int)tenorPairs[0][i],12*(int)tenorPairs[0][i],1);
					double g2_11 = g_function(0,end,expiry,12*(int)tenorPairs[1][i],12*(int)tenorPairs[1][i],1);

					double G2_01 = g_function(end,expiry,expiry,12*(int)tenorPairs[0][i],12*(int)tenorPairs[1][i],1);
					double G2_00 = g_function(end,expiry,expiry,12*(int)tenorPairs[0][i],12*(int)tenorPairs[0][i],1);
					double G2_11 = g_function(end,expiry,expiry,12*(int)tenorPairs[1][i],12*(int)tenorPairs[1][i],1);

					double x1 = ::pow(g1_01 + g2_01,2);
					double x2 = ::pow(G2_01,2);
					double x3 = (g1_01 + g2_01)*G2_01;
	                
					double cor = ::pow( marketCorrelation[i] , 2 );

					double y1 = cor * (g1_00 + g2_00) * (g1_11 + g2_11);
					double y2 = cor * G2_00 * G2_11;
					double y3 = cor * (G2_00 * (g1_11 + g2_11) + G2_11 * (g1_00 + g2_00));
				

					// define coeeficients in 2nd order equation a*q^2+b*q+c=0, where q=c2^2/c1^2
					double a = x2-y2; double b = 2.*x3-y3; double c = x1-y1; 
					for(int z=j;z<=i;z++)
					{
						if(b*b-4.*a*c>0)
						{
							double r1 = ( -b + sqrt( b*b - 4.*a*c )) /2./a;
							double r2 = ( -b - sqrt( b*b - 4.*a*c )) /2./a;
							
							double root =0.;
							if(r1*r2<0)
								root = max(r1,r2);
							else
								root = min(r1,r2);
							
							root = max(root,0.0);
							(pC3param)[3][z] = sqrt(root);
												
						}
						else
						{
							(pC3param)[3][z] = 0.0;//max(-b/2./a , 0.);
						}
					}

					// update end point for bootstrapping
					end = sv_data.SwaptionExpiries[i];
					j = i+1;
				}
			}

			// flat extrapolation
			for(int i=j; i<sv_data.NbSwaptionExpiries; i++)
			{
				(pC3param)[3][i]=(pC3param)[3][j-1];
			}

			// solve for c1 by ordinary least square
			double xx=0; double xy=0.;
			for(int i = 0; i < sv_data.NbSwaptionExpiries; i++)
			{
				q = ::pow((pC3param)[3][i],2);
				int expiry = sv_data.SwaptionExpiries[i];
				for(int j = 0; j < sv_data.NbSwapTenors; j++)
				{
					double w = weightMatrix[j][i];
					double v = ::pow(sv_data.VolMatrix[i][j],2);

					xx += w * ::pow( g_function(0,expiry,expiry,sv_data.SwapTenors[j],sv_data.SwapTenors[j],0)
							+ g_function(0,expiry,expiry,sv_data.SwapTenors[j],sv_data.SwapTenors[j],1) , 2 );
					xy += w * v * ( g_function(0,expiry,expiry,sv_data.SwapTenors[j],sv_data.SwapTenors[j],0)
							+ g_function(0,expiry,expiry,sv_data.SwapTenors[j],sv_data.SwapTenors[j],1) );
				}
			}

			//set weights
			for(int i = 0; i < sv_data.NbSwaptionExpiries; i++)
			{
				pC3param[2][i] = sqrt(xy/xx);
				pC3param[3][i] = pC3param[3][i] * sqrt(xy/xx);
			}
		}
		else
		// do spot vol bootstrapping
		{
			// set start point for bootstrapping
			spotStart=0;
			for(int i = 0; i < sv_data.NbSwaptionExpiries; i++)
			{	
				if(calibrateCorrelation[i])
				{
					// set end point for bootstrapping
					spotEnd =i;
					bool not_ready = true; double high =0.1; int counter = 0;
					while(not_ready)
					{
						counter++;
						double start = find_q(TINY,this);
						double end = find_q(high,this);
						not_ready = ((start * end >= 0) && (counter < 1000));
						if(not_ready) { high += 0.1; }
					}
					// for some parameter values, we basically cannot find the value of q that makes
					// model correlation and market correlation equal. 
					// below is a hacker fix, and should be improved for the future
					if(high<10)
					{
						try 
						{
							q = zbrentUseful(
									find_q,             // (I) The function to find the root of
									this,               // (I) class input
									TINY,                   // (I) Low value for x
									high,               // (I) High value for x
									1.e-10);            // (I) Tolerance
						} 
						catch (exception&) 
						{
							throw ModelException(   // e, // hide "zbrent" message which confuses users
													"Brent algorithm failed");
						}               

					}
					// when setting q (i.e. c2) to zero, we violate the constraints in the optimizer
					// this will force us back to the valid region.
					else
					{
						q=0.0;
					}

					// set weights
					for(int z = spotStart; z < sv_data.NbSwaptionExpiries; z++)
					{
						pC3param[2][z] = 1.;
						pC3param[3][z] = sqrt(q);
					}
					
					// update start point for bootstrapping 
					spotStart = i+1;
				}
			}
		}
    }
    catch(exception& e) { throw ModelException(e, "VnfmTPlus::set_factorWeights"); }
}

// find c2 when we have a spot volatility 
double VnfmTPlus::find_q(double q, void* data) 
{
    VnfmTPlus *vnfmPlus = (VnfmTPlus*) data;
	for(int j=vnfmPlus->spotStart; j<vnfmPlus->sv_data.NbSwaptionExpiries; j++)
	{
		vnfmPlus->pC3param[2][j] = 1.;
		vnfmPlus->pC3param[3][j] = sqrt(q);
	}
    vnfmPlus->set_spotvol(vnfmPlus->reference_tenor->intValue());
    try 
    {
        int i=vnfmPlus->spotEnd;
        double result = vnfmPlus->correlation(vnfmPlus->sv_data.SwaptionExpiries[i],
                        12*vnfmPlus->tenorPairs[0][i],12*vnfmPlus->tenorPairs[1][i])
						- vnfmPlus->marketCorrelation[i];
		return result;
    }
    catch(exception& e) { throw ModelException(e, "VnfmTPlus::find_q"); }
}

// update the R*K matrix
void VnfmTPlus::update_rKmatrix()
{
	// mean reversion parameters and angle are not allowed to be time dependent
	double angle=0.;
	if(humptime.get() != NULL) 
	{
		bool not_ready = true; double low =-0.7854; double high =low+0.01; int counter = 0;
		while(not_ready)
		{
			counter++;
			double start = getAnglefromHumpTime(low,this);
			double end = getAnglefromHumpTime(high,this);
			not_ready = ((start * end >= 0) && (counter < 1000));
			if(not_ready) { high += 0.01; }
		}
		// for some parameter values, we basically cannot find the value of q that makes
		// model correlation and market correlation equal. 
		// below is a hacker fix, and should be improved for the future
		if(!not_ready)
		{
			try 
			{
				angle = zbrentUseful(
						getAnglefromHumpTime,	// (I) The function to find the root of
						this,					// (I) class input
						low,					// (I) Low value for x
						high,					// (I) High value for x
						1.e-10);				// (I) Tolerance
			} 
			catch (exception&) 
			{
				throw ModelException(   // e, // hide "zbrent" message which confuses users
										"Brent algorithm failed");
			}               
		}

		for(int i = 0; i < sv_data.NbSwaptionExpiries; i++)
		{
	        pC3param[4][i]=angle;
		}
	}

	rKmatrix = setrKmatrix(pC3param[0][0], pC3param[1][0], pC3param[4][0]);
}

// find the angle in the PC3 representation that corresponds to given Fix3 paramters. 
double VnfmTPlus::getAnglefromHumpTime(double angle, void* data) // (I) angle in PC3 parameters
{
    VnfmTPlus *vnfmPlus = (VnfmTPlus*) data;
    try 
    {
		for(int i = 0; i < vnfmPlus->sv_data.NbSwaptionExpiries; i++)
		{
	        vnfmPlus->pC3param[4][i]=angle;
		}

        double result = (*hump(vnfmPlus->pC3param))[0] - vnfmPlus->humptime->doubleValue();

        return result;
    }
    catch(exception& e) { throw ModelException(e, "VnfmT::getAnglefromHumpTime"); }
}

// Populate the model swaption vol grid. This is a public function.
CDoubleArraySP VnfmTPlus::getCorrelation()
{
	try
	{
		// local variables for function call
		ExpiryArray selectedTenors, selectedExpiries; 
		IrConverter::to_SWAPVOL_DATA(sv_data, selectedTenors, selectedExpiries, irvol.get());  

		CDoubleArraySP corr(new CDoubleArray(sv_data.NbSwaptionExpiries));

		for(int i = 0; i < sv_data.NbSwaptionExpiries; i++)
			(*corr)[i] = correlation(sv_data.SwaptionExpiries[i],12*tenorPairs[0][i],12*tenorPairs[1][i]);

		return corr;
    }
    catch(exception& e) { throw ModelException(e, "VnfmTPlus::getCorrelation"); }
}

// Calibrate the PC3 model
CDoubleMatrixSP VnfmTPlus::getPC3Calibration()
{
	(humptime.get() == NULL) ? nbVar = 3 : nbVar = 2; 
    return (getCalibration("PC3"));
}

// Calibrate the model
CDoubleMatrixSP VnfmTPlus::getCalibration(string model)
{
    try
    {
		// local variables for function call
		ExpiryArray selectedTenors, selectedExpiries; 
		IrConverter::to_SWAPVOL_DATA(sv_data, selectedTenors, selectedExpiries, irvol.get());  

        CDoubleMatrixSP param(new CDoubleMatrix(pC3param.numCols(),pC3param.numRows()));

        // initialise optimizer//
        int row_size = weightMatrix.numRows(); 
        int col_size = weightMatrix.numCols(); 
        if(calibrator == "Quasi-Newton")
        {
            nbFunc = 1;
        }
        else if(calibrator == "Levenberg-Marquardt")
        {
            int nb = 0;
            for(int i = 0; i < row_size; i++)
            {
                for(int j = 0; j < col_size; j++)
                    if(weightMatrix[j][i]>TINY) { nb++; }   
            }

            nbFunc = nb;
        }
        else 
        {
            throw ModelException("optimizer not identified");
        }

        DoubleArray x(nbVar); DoubleArray guess(nbVar); RangeArray ranges(nbVar);
        if(model == "PC3")
        {
            guess[0] = pC3param[0][0]; guess[1] = pC3param[1][0]; 
			if(humptime.get() == NULL) { guess[2] = pC3param[4][0]; }

            ranges[0] = RangeSP(new Range (OpenBoundary(0.0),    OpenBoundary(3.0)));
            ranges[1] = RangeSP(new Range (OpenBoundary(0.0),    OpenBoundary(3.0)));
			if(humptime.get() == NULL) 
			{
				ranges[2] = RangeSP(new Range (OpenBoundary(-0.7854), OpenBoundary(1.57)));   
			}	
		}
        else
        {
            throw ModelException("invalid model");
        }

        ObjFuncT objFunc(ranges,nbVar,nbFunc,*this,model);
        if(calibrator == "Quasi-Newton")
        {
            QuasiNewton().minimize(objFunc,guess,x);
        }
        else if(calibrator == "Levenberg-Marquardt")
        {
            LevenbergMarquardt().minimize(objFunc,guess,x);
        }
        else 
        {
            throw ModelException("optimizer not identified"); 
        }

		double angle;
		if(humptime.get() != NULL) 
		{				
			bool not_ready = true; double low =-0.7854; double high =low + 0.01; int counter = 0;
			while(not_ready)
			{
				counter++;
				double start = getAnglefromHumpTime(low,this);
				double end = getAnglefromHumpTime(high,this);
				not_ready = ((start * end >= 0) && (counter < 1000));
				if(not_ready) { high += 0.01; }
			}
			// for some parameter values, we basically cannot find the value of q that makes
			// model correlation and market correlation equal. 
			// below is a hacker fix, and should be improved for the future
			if(!not_ready)
			{
				try 
				{
					angle = zbrentUseful(
							getAnglefromHumpTime,	// (I) The function to find the root of
							this,					// (I) class input
							low,					// (I) Low value for x
							high,					// (I) High value for x
							1.e-10);				// (I) Tolerance
				} 
				catch (exception&) 
				{
					throw ModelException(   // e, // hide "zbrent" message which confuses users
											"Brent algorithm failed");
				}               
			}

		}
		for(int i = 0; i < sv_data.NbSwaptionExpiries; i++)
		{
	        pC3param[0][i]=x[0]; pC3param[1][i]=x[1]; 
			humptime.get() != NULL ? pC3param[4][i] = angle : pC3param[4][i] = x[2]; 
		}
        set_factorWeights();

        // assign values
        for(int i=0; i<pC3param.numCols(); i++)
		{
	        for(int j=0; j<pC3param.numRows(); j++)
			{
		        ((*param)[i][j])=pC3param[i][j];
			}
		}

        return param;
    }
    catch(exception& e) { throw ModelException(e, "VnfmTPlus::getCalibration"); }
}

/*************************************************************************************************/


/***** VNFM PC MOVE ******************************************************************************/
// Excel interface.

void VnfmTPCMove::load(CClassSP& clazz) {
    clazz->setPublic();                         // make visible to EAS/spreadsheet
    clazz->setDescription("My description");
    REGISTER(VnfmTPCMove, clazz);
    SUPERCLASS(VnfmTPlus);                           // base class
    FIELD(pC_move,"Perturbation of sensitivity"); FIELD_MAKE_OPTIONAL(pC_move);
    FIELD(pC_order,"Order of orthogonality");
    FIELD(pC_bump,"Bump size to calculate sensitivities"); FIELD_MAKE_OPTIONAL(pC_bump);
    EMPTY_SHELL_METHOD(defaultConstructor);

    Addin::registerObjectMethod(
            "PC_MOVE_T", Addin::UTILITIES,
            "Help message",
            false, Addin::expandSimple,
            &VnfmTPCMove::getPCMove); 

    Addin::registerObjectMethod(
            "PC_SENSITIVITY_T", Addin::UTILITIES,
            "Help message",
            false, Addin::expandSimple,
            &VnfmTPCMove::getPCSensitivity);  
}

CClassConstSP const VnfmTPCMove::TYPE = CClass::registerClassLoadMethod("VnfmTPCMove", typeid(VnfmTPCMove), load);

// constructor
VnfmTPCMove::VnfmTPCMove(CDoubleMatrix &param, DayCountConventionSP dcc, int freq, CIntSP ref_tenor, 
                   IRVolSP ir, YieldCurveSP yc, CDoubleMatrix &weights, CDoubleMatrix &tp,
                   CDoubleArray &mc, CBoolArray &cc, CDoubleSP hp, CDoubleArraySP &pert, 
				   IntArray &ortho, CDoubleArraySP &bs) :
VnfmTPlus(param,dcc,freq,ref_tenor,ir,yc,weights,tp,mc,cc,hp,TYPE),
pC_move(pert),
pC_order(ortho),
pC_bump(bs)
{
}

CDoubleArraySP VnfmTPCMove::getPCMove()
{
	try
	{
		// local variables for function call
		ExpiryArray selectedTenors, selectedExpiries; 
		IrConverter::to_SWAPVOL_DATA(sv_data, selectedTenors, selectedExpiries, irvol.get());  

		// base case
		CDoubleMatrixSP pC = pComponents(pC3param);

		// various perturbed cases
		double a1 = pC3param[0][0]; double a2 = pC3param[1][0]; double angle = pC3param[4][0]; 
		CDoubleMatrix rK = (*setrKmatrix(a1, a2, angle));   

		pC3param[0][0] = a1 + (*pC_bump)[0]; pC3param[1][0] = a2; pC3param[4][0] = angle;
		set_factorWeights();
		CDoubleMatrixSP pC_1 = pComponents(pC3param);   

		pC3param[0][0] = a1; pC3param[1][0] = a2 + (*pC_bump)[1]; pC3param[4][0] = angle;
		set_factorWeights();
		CDoubleMatrixSP pC_2 = pComponents(pC3param);   

		pC3param[0][0] = a1; pC3param[1][0] = a2; pC3param[4][0] = angle + (*pC_bump)[2];
		set_factorWeights();
		CDoubleMatrixSP pC_3 = pComponents(pC3param);   

		// reset
		pC3param[0][0] = a1; pC3param[1][0] = a2; pC3param[4][0] = angle;


		// calculate the length of the sensitivities
		CDoubleMatrix LPC1(3,3);
		CDoubleMatrix LPC2(3,3);
		DoubleArray deltaPC(6);
		for(int i=0; i<NUM_TO_M; i++)
		{
			// calculate i'th sensitivity and order it.
			deltaPC[pC_order[0]-1]      = ((*pC_1)[1][i]-(*pC)[1][i])/(*pC_bump)[0];
			deltaPC[pC_order[1]-1]      = ((*pC_2)[1][i]-(*pC)[1][i])/(*pC_bump)[1];
			deltaPC[pC_order[2]-1]      = ((*pC_3)[1][i]-(*pC)[1][i])/(*pC_bump)[2];
			deltaPC[3 + pC_order[0]-1]  = ((*pC_1)[2][i]-(*pC)[2][i])/(*pC_bump)[0];
			deltaPC[3 + pC_order[1]-1]  = ((*pC_2)[2][i]-(*pC)[2][i])/(*pC_bump)[1];
			deltaPC[3 + pC_order[2]-1]  = ((*pC_3)[2][i]-(*pC)[2][i])/(*pC_bump)[2];

			// calculate the norms
			LPC1[0][0] += deltaPC[0]*deltaPC[0];
			LPC1[0][1] += deltaPC[0]*deltaPC[1];
			LPC1[0][2] += deltaPC[0]*deltaPC[2];
			LPC1[1][1] += deltaPC[1]*deltaPC[1];
			LPC1[1][2] += deltaPC[1]*deltaPC[2];
			LPC1[2][2] += deltaPC[2]*deltaPC[2];

			LPC2[0][0] += deltaPC[3+0]*deltaPC[3+0];
			LPC2[0][1] += deltaPC[3+0]*deltaPC[3+1];
			LPC2[0][2] += deltaPC[3+0]*deltaPC[3+2];
			LPC2[1][1] += deltaPC[3+1]*deltaPC[3+1];
			LPC2[1][2] += deltaPC[3+1]*deltaPC[3+2];
			LPC2[2][2] += deltaPC[3+2]*deltaPC[3+2];
		}

		// useful variables for Gram-Schmidt 
		double x_1 = ( LPC1[0][1] / sqrt(LPC1[0][0]*LPC1[1][1]) );
		double y_1 = ( LPC1[0][2] / sqrt(LPC1[0][0]*LPC1[2][2]) );
		double z_1 = 1./sqrt(1.-x_1*x_1)*LPC1[1][2]/sqrt(LPC1[1][1]*LPC1[2][2]) 
						-x_1/sqrt(1.-x_1*x_1)*LPC1[0][2]/sqrt(LPC1[0][0]*LPC1[2][2]);
		double x_2 = ( LPC2[0][1] / sqrt(LPC2[0][0]*LPC2[1][1]) );
		double y_2 = ( LPC2[0][2] / sqrt(LPC2[0][0]*LPC2[2][2]) );
		double z_2 = 1./sqrt(1.-x_2*x_2)*LPC2[1][2]/sqrt(LPC2[1][1]*LPC2[2][2]) 
						-x_2/sqrt(1.-x_2*x_2)*LPC2[0][2]/sqrt(LPC2[0][0]*LPC2[2][2]);

		// based on G-S we derive the transition matrix from base moves to PC3 moves
		CDoubleMatrix DPC(6,3);
		DPC[0][0] = 1. / sqrt(LPC1[0][0]);
		DPC[0][1] = (-x_1 / sqrt(1. - x_1*x_1)) / sqrt(LPC1[0][0]);
		DPC[1][1] = (1. / sqrt(1. - x_1*x_1)) / sqrt(LPC1[1][1]);
		if(humptime.get() == NULL) 
		{				
			DPC[0][2] = (x_1*z_1 / sqrt(1. - x_1*x_1) / sqrt(1. - y_1*y_1 - z_1*z_1) 
						- y_1 / sqrt(1. - y_1*y_1 - z_1*z_1)) / sqrt(LPC1[0][0]);
			DPC[1][2] = (-z_1 / sqrt(1. - x_1*x_1) / sqrt(1. - y_1*y_1 - z_1*z_1)) / sqrt(LPC1[1][1]);
			DPC[2][2] = (1. / sqrt(1. - y_1*y_1 - z_1*z_1)) / sqrt(LPC1[2][2]);
		}
		else
		{
			DPC[0][2] = 0.0;
			DPC[1][2] = 0.0;
			DPC[2][2] = 0.0;
		}
		DPC[3][0] = 1. / sqrt(LPC2[0][0]);
		DPC[3][1] = (-x_2 / sqrt(1. - x_2*x_2)) / sqrt(LPC2[0][0]);
		DPC[4][1] = (1. / sqrt(1. - x_2*x_2)) / sqrt(LPC2[1][1]);
		if(humptime.get() == NULL) 
		{				
			DPC[3][2] = (x_2*z_2 / sqrt(1. - x_2*x_2) / sqrt(1. - y_2*y_2 - z_2*z_2) 
						- y_2 / sqrt(1. - y_2*y_2 - z_2*z_2)) / sqrt(LPC2[0][0]);
			DPC[4][2] = (-z_2 / sqrt(1. - x_2*x_2) / sqrt(1. - y_2*y_2 - z_2*z_2)) / sqrt(LPC2[1][1]);
			DPC[5][2] = (1. / sqrt(1. - y_2*y_2 - z_2*z_2)) / sqrt(LPC2[2][2]); 
		}
		else
		{
			DPC[3][2] = 0.0;
			DPC[4][2] = 0.0;
			DPC[5][2] = 0.0;
		}

		// we focus on the 1st principal component
		int dummy = LPC1[0][0] > LPC2[0][0] ? 0 : 1; dummy = 0;


		CDoubleArraySP out(new DoubleArray(3));
		if(pC_move.get()!=NULL)
		{
			for(int i=0;i<3;i++)
			{
				(*out)[i] =0;
				for(int j=0;j<3;j++)
					(*out)[i] += DPC[pC_order[i]-1+3*dummy][j] * (*pC_move)[j];
			}
		}

		(*out)[0] += a1;
		(*out)[1] += a2;
		(*out)[2] += angle;

		return (out);
    }
    catch(exception& e) { throw ModelException(e, "VnfmTPlus::getPCMove"); }
}

CDoubleMatrixSP VnfmTPCMove::getPCSensitivity()
{
	try
	{
		// local variables for function call
		ExpiryArray selectedTenors, selectedExpiries; 
		IrConverter::to_SWAPVOL_DATA(sv_data, selectedTenors, selectedExpiries, irvol.get());  

		// base case
		CDoubleMatrixSP pC = pComponents(pC3param);

		// various perturbed cases
		double a1 = pC3param[0][0]; double a2 = pC3param[1][0]; double angle = pC3param[4][0]; 
		CDoubleMatrix rK = (*setrKmatrix(a1, a2, angle));   

		pC3param[0][0] = a1 + (*pC_bump)[0]; pC3param[1][0] = a2; pC3param[4][0] = angle;
		set_factorWeights();
		CDoubleMatrixSP pC_1 = pComponents(pC3param);   

		pC3param[0][0] = a1; pC3param[1][0] = a2 + (*pC_bump)[1]; pC3param[4][0] = angle;
		set_factorWeights();
		CDoubleMatrixSP pC_2 = pComponents(pC3param);   

		pC3param[0][0] = a1; pC3param[1][0] = a2; pC3param[4][0] = angle + (*pC_bump)[2];
		set_factorWeights();
		CDoubleMatrixSP pC_3 = pComponents(pC3param);   

		// reset
		pC3param[0][0] = a1; pC3param[1][0] = a2; pC3param[4][0] = angle;


		// calculate the length of the sensitivities
		CDoubleMatrix LPC1(3,3);
		CDoubleMatrix LPC2(3,3);

		CDoubleMatrix deltaPC(7,NUM_TO_M);
		for(int i=0; i<NUM_TO_M; i++)
		{
			// calculate i'th sensitivity and order it.
			deltaPC[0][i]               = (*pC)[0][i];
			deltaPC[pC_order[0]][i]     = ((*pC_1)[1][i]-(*pC)[1][i])/(*pC_bump)[0];
			deltaPC[pC_order[1]][i]     = ((*pC_2)[1][i]-(*pC)[1][i])/(*pC_bump)[1];
			deltaPC[pC_order[2]][i]     = ((*pC_3)[1][i]-(*pC)[1][i])/(*pC_bump)[2];
			deltaPC[3 + pC_order[0]][i] = ((*pC_1)[2][i]-(*pC)[2][i])/(*pC_bump)[0];
			deltaPC[3 + pC_order[1]][i] = ((*pC_2)[2][i]-(*pC)[2][i])/(*pC_bump)[1];
			deltaPC[3 + pC_order[2]][i] = ((*pC_3)[2][i]-(*pC)[2][i])/(*pC_bump)[2];

			// calculate the norms
			LPC1[0][0] += deltaPC[1][i]*deltaPC[1][i];
			LPC1[0][1] += deltaPC[1][i]*deltaPC[2][i];
			LPC1[0][2] += deltaPC[1][i]*deltaPC[3][i];
			LPC1[1][1] += deltaPC[2][i]*deltaPC[2][i];
			LPC1[1][2] += deltaPC[2][i]*deltaPC[3][i];
			LPC1[2][2] += deltaPC[3][i]*deltaPC[3][i];

			LPC2[0][0] += deltaPC[3+1][i]*deltaPC[3+1][i];
			LPC2[0][1] += deltaPC[3+1][i]*deltaPC[3+2][i];
			LPC2[0][2] += deltaPC[3+1][i]*deltaPC[3+3][i];
			LPC2[1][1] += deltaPC[3+2][i]*deltaPC[3+2][i];
			LPC2[1][2] += deltaPC[3+2][i]*deltaPC[3+3][i];
			LPC2[2][2] += deltaPC[3+3][i]*deltaPC[3+3][i];
		}

		// useful variables for Gram-Schmidt 
		double x_1 = ( LPC1[0][1] / sqrt(LPC1[0][0]*LPC1[1][1]) );
		double y_1 = ( LPC1[0][2] / sqrt(LPC1[0][0]*LPC1[2][2]) );
		double z_1 = 1./sqrt(1.-x_1*x_1)*LPC1[1][2]/sqrt(LPC1[1][1]*LPC1[2][2]) 
						-x_1/sqrt(1.-x_1*x_1)*LPC1[0][2]/sqrt(LPC1[0][0]*LPC1[2][2]);
		double x_2 = ( LPC2[0][1] / sqrt(LPC2[0][0]*LPC2[1][1]) );
		double y_2 = ( LPC2[0][2] / sqrt(LPC2[0][0]*LPC2[2][2]) );
		double z_2 = 1./sqrt(1.-x_2*x_2)*LPC2[1][2]/sqrt(LPC2[1][1]*LPC2[2][2]) 
						-x_2/sqrt(1.-x_2*x_2)*LPC2[0][2]/sqrt(LPC2[0][0]*LPC2[2][2]);

		// create the ortogonal vectors
		CDoubleMatrixSP DPC(new CDoubleMatrix(7,NUM_TO_M));
		for(int i=0;i<NUM_TO_M;i++)
		{
			(*DPC)[0][i] = deltaPC[0][i];
			(*DPC)[1][i] = deltaPC[1][i] / sqrt(LPC1[0][0]);
			(*DPC)[2][i] = ( deltaPC[2][i] / sqrt(LPC1[1][1]) - (*DPC)[1][i]*x_1 ) / sqrt(1.-x_1*x_1);
			if(humptime.get() == NULL) 
			{				
				(*DPC)[3][i] = ( deltaPC[3][i] / sqrt(LPC1[2][2]) - (*DPC)[1][i]*y_1 - (*DPC)[2][i]*z_1 ) 
							/ sqrt(1.-y_1*y_1-z_1*z_1);
			}
			else
			{
				(*DPC)[3][i] = 0.0;
			}
			(*DPC)[4][i] = deltaPC[4][i] / sqrt(LPC2[0][0]);
			(*DPC)[5][i] = ( deltaPC[5][i] / sqrt(LPC2[1][1]) - (*DPC)[4][i]*x_2 ) / sqrt(1.-x_2*x_2);
			if(humptime.get() == NULL) 
			{							
				(*DPC)[6][i] = ( deltaPC[6][i] / sqrt(LPC2[2][2]) - (*DPC)[4][i]*y_2 - (*DPC)[5][i]*z_2 ) 
							/ sqrt(1.-y_2*y_2-z_2*z_2);
			}
			else
			{
				(*DPC)[6][i] = 0.0;
			}
		}

		return (DPC);
    }
    catch(exception& e) { throw ModelException(e, "VnfmTPlus::getPCSensitivity"); }
}

/*************************************************************************************************/

// Finally we need this function to force the compiler to register our class
bool VnfmTestTimeLoad(void) {
    return (VnfmT::TYPE != 0);
}

DRLIB_END_NAMESPACE

