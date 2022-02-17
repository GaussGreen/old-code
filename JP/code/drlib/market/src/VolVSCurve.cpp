//GAD, 15/02/2006

#include "edginc/config.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/Format.hpp"
#include "edginc/Maths.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Complex.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include "edginc/VolVSCurve.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/Nrfns.hpp"

//length for which the FVSCPC is built
static const double XMAX = 5.0;
//number of observations per year (hard coded for now)
static const int PPY = 252;

DRLIB_BEGIN_NAMESPACE

/********************/
/*** class FVSCPC ***/

//piecewise constant forward variance swap curve

//empty constructor
FVSCPC::FVSCPC():CObject(TYPE){
}

//wrap around empty constructor
IObject* FVSCPC::defaultCtor(){
    return new FVSCPC();
}

//constructor
FVSCPC::FVSCPC(const DoubleArray&  inTenors,
               const DoubleArray&  inValues):CObject(TYPE){
    static const string method = "FVSCPC::FVSCPC";
    if( inTenors.size() != inValues.size() ){
        throw ModelException("tenor and value arrays don't match", method);
    }
    else{
        tenors.resize(inTenors.size());
		values.resize(inValues.size());
		for(int i=0; i<inTenors.size(); i++){
            tenors[i]=inTenors[i];
			values[i]=inValues[i];
		}
    }
}

//destructor
FVSCPC::~FVSCPC(){
}


//find the left tenor (piecewise constant curve is assumed to be right continuous)
int FVSCPC::tenorLocation(double tenor) const{
    static const string method = "FVSCPC::tenorLocation";
    try{
        //result
        int res;

        if( tenor < tenors.front() ){
            throw ModelException("input passed was negative", method);
        }
        else{
            if( tenor >= tenors.back() ){
                res = tenors.size() - 1;
            }
            else{
                int compt = 0;
                while( tenor >= tenors[compt] ){
                    compt++;
                }
                res = compt - 1;
            }
        }
        return res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}
/*
int FVSCPC::tenorLocation(double tenor) const{
    static const string method = "FVSCPC::tenorLocation";
    try{
        //result
        unsigned long res;
        unsigned long n = tenors.size();
        double * newtenors = new double[n];
        for (int i = 0; i < tenors.size(); i++) {
            newtenors[i] = tenors[i];
        }
        locate(newtenors,n,tenor,&res);
        return (int) res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}
*/
//yield using piecewise constant interpolation
double FVSCPC::yield(double tenor) const{
    static const string method = "FVSCPC::yield";
    try{
        //result
        double res;
        
        if( tenor < tenors.front() ){
            throw ModelException("input passed was negative", method);
        }
        else{
            int position = tenorLocation(tenor);
            res = values[position];
        }
        return res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

//find the floor fvscpcFloor on range [0,xMax]
smartPtr<FVSCPC> FVSCPC::floor(double   limit,
                               double   xMax) const{
    static const string method = "FVSCPC::floor";
    try{
        //time subdivision
        DoubleArray tenorsFloor(1);
	    tenorsFloor[0] = 0.0;
        while( Maths::isPositive( limit - xMax - tenorsFloor.back() - 0.0000000001 ) ){
            double delay1 = tenors[ tenorLocation( tenorsFloor.back() ) + 1 ] - tenorsFloor.back();
		    double delay2 = tenors[ tenorLocation( tenorsFloor.back() + xMax ) + 1 ] - tenorsFloor.back() - xMax;

            if( Maths::isPositive(delay1-delay2) ){
                tenorsFloor.push_back( (tenorsFloor.back() + delay2) );
            }
            else{
                tenorsFloor.push_back( (tenorsFloor.back() + delay1) );
            }
        }

    	//associated values
        DoubleArray valuesFloor(tenorsFloor.size());
	    int nUp, nDown;
        for(int k=0; k< valuesFloor.size(); k++){
		    nDown = tenorLocation(tenorsFloor[k]);
            nUp = tenorLocation(tenorsFloor[k] + xMax);
            valuesFloor[k] = values[nDown];
            for(int j=nDown; j<nUp+1; j++){
                if( valuesFloor[k] > values[j] ){
                    valuesFloor[k] = values[j];
                }
		    }
	    }
    
        //declaration of curvePtr
        smartPtr<FVSCPC> curvePtr( new FVSCPC(tenorsFloor, valuesFloor) ) ;
        return curvePtr;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

//TYPE used for registration
CClassConstSP const FVSCPC::TYPE = CClass::registerClassLoadMethod("FVSCPC", typeid(FVSCPC), load);
//template<> CClassConstSP const FVSCPCArray::TYPE = CClass::registerClassLoadMethod("FVSCPCArray", typeid(FVSCPCArray), load);

//load method
void FVSCPC::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FVSCPC, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultCtor);

    //makes inline first field
	FIELD(tenors, "tenors");
    //makes inline second field
    FIELD(values, "values");
}

/*** end of class FVSCPC ***/
/***************************/


/********************/
/*** class VSCAff ***/

//piecewise affine variance swap curve

//empty constructor
VSCAff::VSCAff():CObject(TYPE){
}

//wrap around empty constructor
IObject* VSCAff::defaultCtor(){
    return new VSCAff();
}

//constructor
VSCAff::VSCAff(const DoubleArray&   inTenors,
               const DoubleArray&   inValues):CObject(TYPE){
    static const string method = "VSCAff::VSCAff";
    if( inTenors.size() != inValues.size() ){
        throw ModelException("tenor and value arrays don't match", method);
    }
    else{
        tenors.resize(inTenors.size());
		values.resize(inValues.size());
		for(int i=0; i<inTenors.size(); i++){
            tenors[i]=inTenors[i];
			values[i]=inValues[i];
		}
    }
}

//destructor
VSCAff::~VSCAff(){
}

//find the left tenor
int VSCAff::tenorLocation(double tenor) const{
    static const string method = "VSCAff::tenorLocation";
    try{
        //result
        int res;

        if( tenor < tenors.front() ){
            throw ModelException("input passed was negative", method);
        }
        else{
            if( tenor >= tenors.back() ){
                res = tenors.size() - 1;
            }
            else{
                int compt = 0;
                while( tenor >= tenors[compt] ){
                    compt++;
                }
                res = compt-1;
            }           
        }
        return res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}
/*
int VSCAff::tenorLocation(double tenor) const{
    static const string method = "VSCAff::tenorLocation";
    try{
        //result
        unsigned long res;
        unsigned long n = tenors.size();
        double * newtenors = new double[n];
        for (int i = 0; i < tenors.size(); i++) {
            newtenors[i] = tenors[i];
        }
        locate(newtenors,n,tenor,&res);
        return (int) res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}
*/    
//yield using linear interpolation
double VSCAff::yield(double tenor) const{
    static const string method = "VSCAff::yield";
    try{
        //result
        double res;

        if( tenor < tenors.front() ){
            throw ModelException("input passed was negative", method);
        }
        else{
            //last derivative is assumed to be the same as the previous one
            if(tenor >= tenors.back() ){
                int size = tenors.size();
                double slope = (values[size-1] - values[size-2]) / (tenors[size-1] - tenors[size-2]); 
                res = values.back() + slope * (tenor - tenors.back());
            }
            else{
                int position = tenorLocation(tenor);
                double slope = (values[position+1] - values[position]) / (tenors[position+1] - tenors[position]);
                res = values[position] + slope * (tenor - tenors[position]);
            }
        }
        return res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

//TYPE used for registration
CClassConstSP const VSCAff::TYPE = CClass::registerClassLoadMethod("VSCAff", typeid(VSCAff), load);
//template<> CClassConstSP const VSCAffArray::TYPE = CClass::registerClassLoadMethod("VSCAffArray", typeid(VSCAffArray), load);

//load method
void VSCAff::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VSCAff, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultCtor);
    
    //makes inline first field
	FIELD(tenors, "tenors");
    //makes inline second field
    FIELD(values, "values");
}

/*** end of class VSCAff ***/
/***************************/


/*****************************************/
/*** class VolVSCurve::VSCurveVolParam ***/

//volVSCurve volParam

//constructor
VolVSCurve::VSCurveVolParam::VSCurveVolParam(): CVolParam(TYPE){
}	

/*** end of class VolVSCurve::VSCurveVolParam ***/
/************************************************/


/******************************************/
/*** class CalcAlphaBeta (helper class) ***/

//components alpha and beta entering the Laplace transform of volVSCurve

class CalcAlphaBeta{
private:
    //beta
    static Complex calcBetaUniDim(double            tau, //tau = T-t
		    			          double	        mRR, //mean reverse rate
			    		          double	        vVol, //vol of vol
                                  double            corr, //correlation
					              const Complex&	uGam1, //used in the computation of gamma (spot)
                                  const Complex&    uGam2, //used in the computation of gamma (real var)
                                  const Complex&    uT){ //used as a terminal condition
        static const string method = "CalcAlphaBeta::calcBetaUniDim";
        try{
            //result
            Complex res;

			//case vVol = 0 
            if( Maths::isZero(vVol) ){
				//coefficients
				throw ModelException(method, "zero vol of vol not supported");
	        }
            //vVol != 0
            else{
                //gamma
				double mRRSQ = mRR * mRR;
				double vVolSQ = vVol * vVol;
				double corrSQ = corr * corr;
				
                double coeff1 = corrSQ - 1.0;
				double coeff2 = 1.0 - 2.0 * corr * mRR / vVol;

                Complex gamPt1 = coeff1 * uGam1 * uGam1;
                Complex gamPt2 = coeff2 * uGam1;
                Complex gamPt3 = - 2.0 * uGam2;
                Complex gamTot = gamPt1 + gamPt2 + gamPt3;

				Complex gammaSQR = mRRSQ + vVolSQ * gamTot;
		        Complex gamma = sqrt(gammaSQR);
		
                //case gamma = 0
                if( Maths::isZero(gamma) ){
                    //coefficients
                    double a = mRR / vVolSQ;
	        	    Complex b = uT - a;
                    Complex c = 0.5 * vVolSQ * tau;
                    Complex d = 1.0 - c;
                    Complex e = b / d;
                    
                    res = a + e;
                }
                //gamma != 0
		        else{
			        //coefficients
	        	    Complex a = ( mRR - gamma ) / vVolSQ;
                    Complex b = uT - a;
                    Complex c = 0.5 * vVolSQ / gamma;
                    
                    Complex e = - gamma * tau;
                    Complex expUp = exp(e);
                    Complex expDown = 1.0 - exp(e);

                    Complex Up = b * expUp;
                    Complex Down = 1.0 - c * b * expDown;
                    
    			    res = a + (Up / Down); 
                }
            }
            return res;
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

	//beta2 Extended Int used in the computation of alpha from a to b (vVolSQ * beta2 - corrSQ * un+1)
    //beta2 Extended Int Pt1
	static Complex calcBeta2IntegrUniDimPt1(double          tau1, //tau1 = T-a 
	        							    double          tau2, //tau2 = T-b
                                            double          mRR, // mean reverse rate
									        double          vVol, //vol of vol
                                            double          corr, //correlation
									        const Complex&  uGam1, //used in the computation of gamma (spot)
                                            const Complex&  uGam2, //used in the computation of gamma (real var)
                                            const Complex&  uT){ //used as a terminal condition
        static const string method = "CalcAlphaBeta::calcBeta2IntegrUniDimPt1";
        try{
            //result
            Complex res;

			//case vVol = 0
            if( Maths::isZero(vVol) ){
				throw ModelException(method, "zero vol of vol not supported"); 
            }
            //vVol != 0
            else{
                //gamma
				double mRRSQ = mRR * mRR;
				double vVolSQ = vVol * vVol;
				double corrSQ = corr * corr;
				
                double coeff1 = corrSQ - 1.0;
				double coeff2 = 1.0 - 2.0 * corr * mRR / vVol;

                Complex gamPt1 = coeff1 * uGam1 * uGam1;
                Complex gamPt2 = coeff2 * uGam1;
                Complex gamPt3 = - 2.0 * uGam2;
                Complex gamTot = gamPt1 + gamPt2 + gamPt3;

				Complex gammaSQR = mRRSQ + vVolSQ * gamTot;
		        Complex gamma = sqrt(gammaSQR);
		
                //case gamma = 0
                if( Maths::isZero(gamma) ){
                    //coefficients
                    double resPt1 = mRRSQ / vVolSQ;
                    Complex resPt2 = corrSQ * uGam1 * uGam1;
                    Complex resTot = resPt1 - resPt2;

                    double tau = tau1 - tau2;
                    
					res = resTot * tau;
                }
                //gamma != 0
                else{
                    //coefficients
                    Complex a = ( mRR - gamma ) / vVolSQ;
                    Complex resPt1 = vVolSQ * a * a;
                    Complex resPt2 = corrSQ * uGam1 * uGam1;
                    Complex resTot = resPt1 - resPt2;

                    double tau = tau1 - tau2;
                    
					res = resTot * tau;
                }
            }
            return res;
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
	}

    //beta2 Extended Int Pt2
    static Complex calcBeta2IntegrUniDimPt2(double              tau1, //tau1 = T-a
                                            double              tau2, //tau2 = T-b
				    			   	        double              mRR, //mean reverse rate
					    			        double              vVol, //vol of vol
			    					        double              corr, //correlation
							    	        const Complex&      uGam1, //used in the computation of gamma (spot)
                                            const Complex&      uGam2, //used in the computation of gamma (real var)
                                            const Complex&      uT){ //used as a terminal condition
        static const string method = "CalcAlphaBeta::calcBeta2IntegrUniDimPt2";
        try{
            //result
            Complex res;

            //case vVol = 0
			if( Maths::isZero(vVol) ){
                throw ModelException(method, "zero vol of vol not suppported");
            }
            //vVol != 0
            else{
                //gamma
				double mRRSQ = mRR * mRR;
				double vVolSQ = vVol * vVol;
				double corrSQ = corr * corr;
				
                double coeff1 = corrSQ - 1.0;
				double coeff2 = 1.0 - 2.0 * corr * mRR / vVol;

                Complex gamPt1 = coeff1 * uGam1 * uGam1;
                Complex gamPt2 = coeff2 * uGam1;
                Complex gamPt3 = - 2.0 * uGam2;
                Complex gamTot = gamPt1 + gamPt2 + gamPt3;

				Complex gammaSQR = mRRSQ + vVolSQ * gamTot;
		        Complex gamma = sqrt(gammaSQR);
            
				//case gamma = 0
		        if( Maths::isZero(gamma) ){
		            //coefficients
					double a = mRR / vVolSQ;
                    Complex b = uT - a;
                    Complex c = 0.5 * vVolSQ * b;
                     
					Complex g1 = 1.0 - c * tau1;
                    Complex g2 = 1.0 - c * tau2;
			
                    Complex resPt1 = 1.0 / g1;
                    Complex resPt2 = 1.0 / g2;
                    Complex resCoeff = 2.0 * b;

					res = resCoeff * ( resPt1 - resPt2 );                    
	            }
                //gamma != 0
                else{
                    //coefficients
					Complex a = ( mRR - gamma ) / vVolSQ;
                    Complex b = uT - a;
                    
                    Complex resCoeff1 = 4.0 * gamma / vVolSQ;
                    Complex resCoeff2 = 0.5 * vVolSQ / gamma;
                    Complex resCoeff3 = resCoeff2 * b;
                    Complex resCoeff4 = 1.0 - resCoeff3;
                    Complex resCoeff5 = resCoeff1 * resCoeff4;

                    Complex g1 = - gamma * tau1; 
					Complex g2 = - gamma * tau2;
					Complex c1 = resCoeff4 + resCoeff3 * exp(g1);
					Complex c2 = resCoeff4 + resCoeff3 * exp(g2);
                    
                    Complex resPt1 = 1.0 / c1;
                    Complex resPt2 = 1.0 / c2;
                    
					res = resCoeff5 * ( resPt2 - resPt1 );
		        }
	        }
            return res;
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

    //beta2 Extended Int Pt3
    static Complex calcBeta2IntegrUniDimPt3(double              tau1, //tau1 = T-a
	    							        double              tau2, //tau2 = T-b
			    				            double              mRR, //mean reverse rate
                                            double              vVol, //vol of vol
                                            double              corr, //correlation
						    		        const Complex&      uGam1, //used for the computation of gamma (spot)
                                            const Complex&      uGam2, //used for the computation of gamma (real var)
                                            const Complex&      uT){ //used as a terminal condition
        static const string method = "CalcAlphaBeta::calcBeta2IntegrUniDimPt3";
        try{
            //result
            Complex res;

            //case vVol = 0
			if( Maths::isZero(vVol) ){
                throw ModelException(method, "zero vol of vol not supported");
            }
            //vVol != 0
            else{                
                //Jackel method for regularization of log complex ratio
              
                //gamma
				double mRRSQ = mRR * mRR;
				double vVolSQ = vVol * vVol;
				double corrSQ = corr * corr;
				
                double coeff1 = corrSQ - 1.0;
				double coeff2 = 1.0 - 2.0 * corr * mRR / vVol;

                Complex gamPt1 = coeff1 * uGam1 * uGam1;
                Complex gamPt2 = coeff2 * uGam1;
                Complex gamPt3 = - 2.0 * uGam2;
                Complex gamTot = gamPt1 + gamPt2 + gamPt3;

				Complex gammaSQR = mRRSQ + vVolSQ * gamTot;
		        Complex gamma = sqrt(gammaSQR);
		
				//case gamma = 0
                if( Maths::isZero(gamma) ){
                    //tau1 = tau2
                    if( Maths::isZero(tau2 - tau1) ){
                        res = 0.0;
                    }
                    else{
                        //c in the model;
                        Complex cCoeff = vVolSQ * uT - mRR;
                                            
                        //denominator
                        Complex cDenom = 0.5 * cCoeff * tau1 - 1.0;
                        //test for infinite values
                        double kiDown;
                        if( !Maths::isZero( cDenom.imag() ) && Maths::isZero( cDenom.real() )  ){
                            kiDown = Maths::sign(cDenom.imag()) * Maths::PI * 0.5; 
                        }
                        else{
                            kiDown = atan( cDenom.imag() / cDenom.real() );
                        }
                        double rDownSQ = Complex::absSquare(cDenom);
                        
                        //numerator
                        Complex cNum = 0.5 * cCoeff * tau2 - 1.0;
                        //test for infinite values
                        double kiUp;
                        if( !Maths::isZero( cNum.imag() ) && Maths::isZero( cNum.real() )  ){
                            kiUp = Maths::sign(cNum.imag()) * Maths::PI * 0.5;    
                        }
                        else{
                            kiUp = atan( cNum.imag() / cNum.real() );
                        }
                        double rUpSQ = Complex::absSquare(cNum);

                        //result
                        double rCoeff = 0.5 * log(rUpSQ / rDownSQ);
                        double argCoeff = kiUp - kiDown;
                        double resCoeff = 4.0 * mRR / vVolSQ;
                        
                        res = resCoeff * Complex(rCoeff, argCoeff); //corresponds with limit of the general formula
                        // res = resCoeff * log(tau2/tau1);
                      
                    }
                }
                //gamma != 0
                else{
                    //tau1 = tau2
                    if( Maths::isZero(tau2 - tau1) ){
                        res = 0.0;
                    }
                    //tau1 != tau2
                    else{
                        //c in the model is c1/c2 * exp()   ABSOLUTELY NECESSARY TO WRITE IT LIKE THIS!!!!!!
    	        	    Complex c1 = (mRR - gamma) - vVolSQ * uT;
	    	            Complex c2 = (mRR + gamma) - vVolSQ * uT;
                        //test for infinite values
                        double tc1;
                        if( !Maths::isZero( c1.imag() ) && Maths::isZero( c1.real() )  ){
                            tc1 = Maths::sign(c1.imag()) * Maths::PI * 0.5;
                        }
                        else{
                            tc1 = atan( c1.imag()/c1.real() );
                        }
                        //test for infinite values
                        double tc2;
                        if( !Maths::isZero( c2.imag() ) && Maths::isZero( c2.real() )  ){
                            tc2 = Maths::sign(c2.imag()) * Maths::PI * 0.5;
                        }
                        else{
                            tc2 = atan( c2.imag()/c2.real() );
                        }
                        double tc = tc1 - tc2 - gamma.imag() * tau1;

                        //denominator
                        //cDown = c1 / c2 * exp(-gamma*tau1) - 1.0; 
                        double m = floor( (tc + Maths::PI) / (2.0 * Maths::PI) );
                        Complex cDown = (c1 / c2) * exp(- gamma * tau1) - 1.0;
                        //test for infinite values
                        double kiDown;
                        if( !Maths::isZero( cDown.imag() ) && Maths::isZero( cDown.real() )  ){
                            kiDown = Maths::sign(cDown.imag()) * Maths::PI * 0.5;
                        }
                        else{
                            kiDown = atan( cDown.imag() / cDown.real() );
                        }
                        double rDownSQ = Complex::absSquare(cDown);

                        //numerator
                        //cUp = c1 / c2 * exp(-gamma*tau2) - 1.0;
                        double bd = gamma.imag();
                        double tau = tau1 - tau2;
                        double n = floor( (tc + bd * tau + Maths::PI) / (2.0 * Maths::PI) );
                        Complex cUp = (c1 / c2) * exp(- gamma * tau2) - 1.0;
                        //test for infinite values
                        double kiUp;
                        if( !Maths::isZero( cUp.imag() ) && Maths::isZero( cUp.real() )  ){
                            kiUp = Maths::sign(cUp.imag()) * Maths::PI * 0.5;
                        }
                        else{
                            kiUp = atan( cUp.imag() / cUp.real() );
                        }
                        double rUpSQ = Complex::absSquare(cUp);
        
                        //result
                        double resCoeff = 4.0 * mRR / vVolSQ;
                        double rCoeff = 0.5 * log(rUpSQ / rDownSQ);
                        // double argCoeff = kiUp - kiDown + 2.0 * Maths::PI * (n - m);
                        double argCoeff = kiUp - kiDown;

                        res = resCoeff * Complex(rCoeff, argCoeff);
                    }
                }
            }
            return res;   
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

    static Complex calcBeta2IntegrUniDimPt4(double              tau, //tau1 = T-t
			    				            double              mRR, //mean reversion rate
                                            double              vVol, //vol of vol
                                            double              corr, //correlation
						    		        const Complex&      uGam1, //used for the computation of gamma (spot)
                                            const Complex&      uGam2, //used for the computation of gamma (real var)
                                            const Complex&      uT){ //used as a terminal condition
        static const string method = "CalcAlphaBeta::calcBeta2IntegrUniDimPt4";
        // int_t^T log((a exp(gamma * tau) -1) / (a-1)) du
        try{
            Complex res;

            if( Maths::isZero(vVol) ){
                throw ModelException(method, "zero vol of vol not supported");
            }
            //vVol != 0
            else{                
              
                //gamma
				double mRRSQ = mRR * mRR;
				double vVolSQ = vVol * vVol;
				double corrSQ = corr * corr;
				
                double coeff1 = corrSQ - 1.0;
				double coeff2 = 1.0 - 2.0 * corr * mRR / vVol;

                Complex gamPt1 = coeff1 * uGam1 * uGam1;
                Complex gamPt2 = coeff2 * uGam1;
                Complex gamPt3 = - 2.0 * uGam2;
                Complex gamTot = gamPt1 + gamPt2 + gamPt3;

				Complex gammaSQR = mRRSQ + vVolSQ * gamTot;
		        Complex gamma = sqrt(gammaSQR);

                if (Maths::isZero(gamma)){
                    Complex phi = 0.5 * (uT * vVolSQ - mRR);
                    Complex simplecoeff = - 1.0 / (phi * tau - 1.0);
                    double simpsq = Complex::absSquare(simplecoeff);
                    //test for infinite values
                    double argsimp;
                    if( !Maths::isZero( simplecoeff.imag() ) && Maths::isZero( simplecoeff.real() )  ){
                        argsimp = Maths::sign(simplecoeff.imag()) * Maths::PI * 0.5;
                    }
                    else{
                        argsimp = atan(simplecoeff.imag() / simplecoeff.real());
                    }
                    double rsimp = 0.5 * log(simpsq);
                    Complex simpleLog = Complex(rsimp,argsimp);
                    res = 1.0 / phi * (simpleLog - 1.0);
                }
                else{
                    Complex c1 = (mRR - gamma) - vVolSQ * uT;
	    	        Complex c2 = (mRR + gamma) - vVolSQ * uT;
                    Complex c3 = exp(-gamma * tau);
                    Complex a = c1 / c2 * c3;
                    if (! Maths::isZero(a)) {
                        Complex coeff1num = a * exp(gamma * tau) - 1.0;
                        Complex coeff1den = a - 1.0;
                        Complex Dcoeff1 = coeff1num / coeff1den;
                        double coeff1sq = Complex::absSquare(Dcoeff1);
                        //test for infinite values
                        double argcoeff1;
                        if( !Maths::isZero( Dcoeff1.imag() ) && Maths::isZero( Dcoeff1.real() )  ){
                            argcoeff1 = Maths::sign(Dcoeff1.imag()) * Maths::PI * 0.5;
                        }
                        else{
                            argcoeff1 = atan(Dcoeff1.imag() / Dcoeff1.real());
                        }
                        double rcoeff1 = 0.5 * log(coeff1sq);
                        Complex logcoeff1 = Complex(rcoeff1,argcoeff1);
                        Complex Dcoeff2 = a * exp(gamma * tau);
                        double coeff2sq = Complex::absSquare(Dcoeff2);
                        //test for infinite values
                        double argcoeff2;
                        if( !Maths::isZero( Dcoeff2.imag() ) && Maths::isZero( Dcoeff2.real() )  ){
                            argcoeff2 = Maths::sign(Dcoeff2.imag()) * Maths::PI * 0.5;
                        }
                        else{
                            argcoeff2 = atan(Dcoeff2.imag() / Dcoeff2.real());
                        }
                        double rcoeff2 = 0.5 * log(coeff2sq);
                        Complex logcoeff2 = Complex(rcoeff2,argcoeff2);
                        Complex logcoeff = logcoeff1 * logcoeff2 / gamma;

                        Complex dilog1 = CalcAlphaBeta::diLogarithm(a);
                        Complex dilogcoeff = a * exp(gamma * tau);
                        Complex dilog2 = CalcAlphaBeta::diLogarithm(dilogcoeff);
                        res = logcoeff - 1.0 / gamma * (dilog1 - dilog2);
                    } else{
                        res = 0.0;
                    }

                }


            }

            return res;
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

    static Complex calcBeta2IntegrUniDimPt5(double              tau, //tau1 = T-t
			    				            double              mRR, //mean reversion rate
                                            double              vVol, //vol of vol
                                            double              corr, //correlation
						    		        const Complex&      uGam1, //used for the computation of gamma (spot)
                                            const Complex&      uGam2, //used for the computation of gamma (real var)
                                            const Complex&      uT){ //used as a terminal condition
        static const string method = "CalcAlphaBeta::calcBeta2IntegrUniDimPt5";
        try{
            Complex res;

            if( Maths::isZero(vVol) ){
                throw ModelException(method, "zero vol of vol not supported");
            }
            //vVol != 0
            else{                
              
                //gamma
				double mRRSQ = mRR * mRR;
				double vVolSQ = vVol * vVol;
				double corrSQ = corr * corr;
				
                double coeff1 = corrSQ - 1.0;
				double coeff2 = 1.0 - 2.0 * corr * mRR / vVol;

                Complex gamPt1 = coeff1 * uGam1 * uGam1;
                Complex gamPt2 = coeff2 * uGam1;
                Complex gamPt3 = - 2.0 * uGam2;
                Complex gamTot = gamPt1 + gamPt2 + gamPt3;

				Complex gammaSQR = mRRSQ + vVolSQ * gamTot;
		        Complex gamma = sqrt(gammaSQR);

                if (Maths::isZero(gamma)){
                    Complex coeff1 = 0.5 * vVolSQ * (uT - mRR / vVolSQ);
                    Complex coeff2 = 1.0 / (1.0 - coeff1 * tau);
                    double coeff2sq = Complex::absSquare(coeff2);
                    //test for infinite values
                    double argcoeff;
                    if( !Maths::isZero( coeff2.imag() ) && Maths::isZero( coeff2.real() )  ){
                        argcoeff = Maths::sign(coeff2.imag()) * Maths::PI * 0.5;
                    }
                    else{
                        argcoeff = atan(coeff2.imag() / coeff2.real());
                    }
                    double rcoeff = 0.5 * log(coeff2sq);
                    Complex coeff3 = Complex(rcoeff,argcoeff);
                    res = coeff3 / coeff1;
                }
                else{
                    Complex coeffnum = vVolSQ * uT - (mRR - gamma);
                    Complex coeffden = 2.0 * gamma;
                    Complex coeff = coeffnum / coeffden;
                    if (Maths::isZero(coeff - 1.0)) {
                        Complex decay = exp(gamma * tau);
                        res = 1.0 / gamma * (decay - 1.0);
                    }
                    else {
                        Complex decay = exp(-gamma * tau);
                        Complex dcoeff1 = 1.0 - coeff + coeff * decay;
                        double coeffsq = Complex::absSquare(dcoeff1);
                        //test for infinite values
                        double argcoeff;
                        if( !Maths::isZero( dcoeff1.imag() ) && Maths::isZero( dcoeff1.real() )  ){
                            argcoeff = Maths::sign(dcoeff1.imag()) * Maths::PI * 0.5;
                        }
                        else{
                            argcoeff = atan(dcoeff1.imag() / dcoeff1.real());
                        }
                        double rcoeff = 0.5 * log(coeffsq);
                        Complex logcoeff = Complex(rcoeff,argcoeff);
                        Complex den = gamma * (1.0 - coeff);
                        res = tau / (1.0 - coeff) + logcoeff / den;
                    }

                }


            }

            return res;
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }
                            


    static Complex calcBeta2IntegrUniDimPt6(double              tau, //tau1 = T-t
			    				            double              mRR, //mean reversion rate
                                            double              vVol, //vol of vol
                                            double              corr, //correlation
						    		        const Complex&      uGam1, //used for the computation of gamma (spot)
                                            const Complex&      uGam2, //used for the computation of gamma (real var)
                                            const Complex&      uT){ //used as a terminal condition
        static const string method = "CalcAlphaBeta::calcBeta2IntegrUniDimPt6";
        // computing \int_t^T F(u) du where F(u) = \int_t^u a^2 \beta^2(s) ds
        try{
            Complex res;

            if( Maths::isZero(vVol) ){
                throw ModelException(method, "zero vol of vol not supported");
            }
            //vVol != 0
            else{                
              
                //gamma
				double mRRSQ = mRR * mRR;
				double vVolSQ = vVol * vVol;
				double corrSQ = corr * corr;
				
                double coeff1 = corrSQ - 1.0;
				double coeff2 = 1.0 - 2.0 * corr * mRR / vVol;

                Complex gamPt1 = coeff1 * uGam1 * uGam1;
                Complex gamPt2 = coeff2 * uGam1;
                Complex gamPt3 = - 2.0 * uGam2;
                Complex gamTot = gamPt1 + gamPt2 + gamPt3;

				Complex gammaSQR = mRRSQ + vVolSQ * gamTot;
		        Complex gamma = sqrt(gammaSQR);

                if (Maths::isZero(gamma)) {
                    Complex coeff = 0.5 * vVolSQ * (uT - mRR / vVolSQ);
                    Complex coeff8 = 2.0 * (mRR / vVolSQ - uT);
                    Complex term1 = - coeff8 * tau / (1.0 - coeff * tau);
                    Complex coeff9 = CalcAlphaBeta::calcBeta2IntegrUniDimPt5(tau,mRR,vVol,corr,uGam1,uGam2,uT);
                    Complex term2 = coeff8 * coeff9;
                    Complex coeff10 = 4.0 * mRR / vVolSQ;
                    Complex coeff11 = CalcAlphaBeta::calcBeta2IntegrUniDimPt4(tau,mRR,vVol,corr,uGam1,uGam2,uT);
                    Complex term3 = coeff10 * coeff11;
                    Complex coeff12 = (mRR - gamma) / vVol;
                    Complex term4 = 0.5 * Maths::square(coeff12 * tau) ;
                    res = term1 + term2 + term3 + term4;  
                } else{
                    Complex coeff = (vVolSQ * uT - (mRR - gamma)) / (2.0 * gamma);
                    Complex decay = exp(-gamma * tau);
                    Complex coeff4 = 1.0 / (1.0 - coeff + coeff * decay);
                    Complex term1 = - 4.0 * gamma / vVolSQ * tau * (1.0-coeff) * coeff4;
                    Complex coeff8 = 4.0 * gamma / vVolSQ * (1.0 - coeff);
                    Complex coeff9 = CalcAlphaBeta::calcBeta2IntegrUniDimPt5(tau,mRR,vVol,corr,uGam1,uGam2,uT);
                    Complex term2 = coeff8 * coeff9;
                    Complex coeff10 = 4.0 * mRR / vVolSQ;
                    Complex coeff11 = CalcAlphaBeta::calcBeta2IntegrUniDimPt4(tau,mRR,vVol,corr,uGam1,uGam2,uT);
                    Complex term3 = coeff10 * coeff11;
                    Complex coeff12 = (mRR - gamma) / vVol;
                    Complex term4 = 0.5 * Maths::square(coeff12 * tau) ;
                    res = term1 + term2 + term3 + term4;  
                }
            }
            return res;
        }
                                                                                                                  
        catch(exception& e){
            throw ModelException(e, method);
        }
      }
                                            




    //beta2 Extended Int
    static Complex calcBeta2IntegrUniDim(double             tau1, //tau1 = T-a
	            	    				 double             tau2, //tau2 = T-b
		                   		         double             mRR, //mean reverse rate
				    			         double             vVol, //vol of vol
		    			                 double             corr, //correlation
						    	         const Complex&     uGam1, //used for the computation of gamma (spot)
                                         const Complex&     uGam2, //used for the computation of gamma (real var)
                                         const Complex&     uT){ //used as a terminal condition
        static const string method = "CalcAlphaBeta::calcBeta2IntegrUniDim";
        try{
            //result
            Complex res;
            
            //computation of beta2 Extended Int
            Complex a = CalcAlphaBeta::calcBeta2IntegrUniDimPt1(tau1, tau2, mRR, vVol, corr, uGam1, uGam2, uT);
            Complex b = CalcAlphaBeta::calcBeta2IntegrUniDimPt2(tau1, tau2, mRR, vVol, corr, uGam1, uGam2, uT);
            Complex c = CalcAlphaBeta::calcBeta2IntegrUniDimPt3(tau1, tau2, mRR, vVol, corr, uGam1, uGam2, uT);
            res = a + b + c;
        

            return res;
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }


   static Complex diLogarithm(const Complex& z){
       static const string method = "CalcAlphaBeta::diLogarithm";
       try{
           Complex res = 0.0;
           // series expansion inside a disc of radius 1
           if (! Maths::isZero(Complex::absSquare(z)-1.0)) {
           // dilog(1) = 0.0
            if (Complex::absSquare(z) < 1.0) {
                for (int k = 1; k < 10; k++) {
                    Complex arg1 = pow(-1.0,k);
                    Complex argz = z - 1.0;
                    Complex arg2 = pow(argz,k);
                    Complex arg3 = k*k;
                    res += arg1 * arg2 / arg3;                
                }
            }
            else{
            // functional relationship between dilog(x) and dilog(1/x)
                Complex t = 1.0 / z;
                Complex scaledRes = 0.0;
                for (int j =1; j<10; j++){
                    Complex arg1 = pow(-1.0,j);
                    Complex argz = t - 1.0;
                    Complex arg2 = pow(argz,j);
                    Complex arg3 = j*j;
                    scaledRes += arg1 * arg2 / arg3; 
                }
                double coeffsq = Complex::absSquare(t);
                
                //test for infinite values
                double argCoeff;
                if( !Maths::isZero( t.imag() ) && Maths::isZero( t.real() )  ){
                    argCoeff = Maths::sign(t.imag()) * Maths::PI * 0.5;
                }
                else{
                    argCoeff = atan(t.imag() / t.real());
                }
                double rCoeff = 0.5 * log(coeffsq);
                Complex coeff = Complex(rCoeff,argCoeff);
                res = - 0.5 * Maths::square(coeff) - scaledRes;
            }
           }
           return res;

       }
       catch(exception& e){
           throw ModelException(e,method);
       }
   }


public:
    //Alpha and Beta for Laplace transform
    static void calcJointLapAlphaBeta(const DateTime&       t, //date of evaluation = today
                                      const DateTime&	    T, //maturity
				                      const ComplexArray&   u, //initial conditions
							          Complex&	            alpha, //AJD
							          ComplexArray&         beta, //AJD
                                      const VolVSCurve*     volvscur){ //used to pass parameters
        static const string method = "CalcAlphaBeta::calcJointLapAlphaBeta";
        try {
            //tenors x are passed directly from volVSCurve-> because of MaturityPeriod copy constructor declared as private
		    
            //size = number of factors
            int size = volvscur->n;

            double k_J = volvscur->crashMRR;
            double lambda_V = volvscur->crashRate;
            double mu_V = volvscur->crashSize;

            double lambda_S = volvscur->spotcrashRate;
            double k_S = volvscur->spotcrashSize;
            double delta_S = volvscur->spotcrashUncertainty;

            double gamma_S = log(1.0 + k_S) - 0.5 * Maths::square(delta_S);
            double chi_S = log(1.0 + k_S) + 0.5 * Maths::square(delta_S);

            //delays between dates expressed in years
            double tempt = (volvscur->timeMetric)->yearFrac(volvscur->baseDate, t);
            double tempT = (volvscur->timeMetric)->yearFrac(volvscur->baseDate, T);
		    double tau = tempT - tempt; //in years

            FVSCPCSP fvscpcFloor = volvscur->fvscpcFloor;
            

            if( (int)(u.size()) != size+3)
                throw ModelException("size of u does not correspond to the number of factors used");
            else{
                //beta coefficients
    	        beta.resize(size+3);

                //constant beta coefficients
                beta[size] = u[size]; //spot
                beta[size+1] = u[size+1]; // quadratic variation   

                //compute beta coeff for Poisson state variable
                Complex phi = 0.5 * u[size] - u[size+1] - 0.5 * Maths::square(u[size]);
                if ( ! Maths::isZero(k_J)) {
                    Complex Ccoeff = - phi / k_J;
                    Complex Acoeff = (u[size+2] - Ccoeff); 
                    beta[size+2] = Ccoeff + Acoeff * exp(- k_J * tau);                 
                } else {
                    beta[size+2] = u[size+2] - phi * tau;
                }
                
                //beta functions
                for(int i=0; i<size; i++){
                    beta[i] = calcBetaUniDim(tau,
                                             volvscur->meanReversRate[i],
                                             volvscur->volVol[i],
                                             volvscur->correlation[i],
                                             u[size], //corresponding to spot beta
                                             u[size+1], //corresponding to real var beta
                                             u[i]); //used as a terminal condition)
                }

                //Dummy coefficient used for integrated variance
                double VS_0T = (volvscur->vscaff)->yield( tempT ); 
                double VS_0t = (volvscur->vscaff)->yield( tempt ); 
		    
		        //alpha coefficient Pt1
                Complex alphaCoeff = 0.5 * u[size] * (u[size] - 1.0) + u[size+1];
		        alpha = alphaCoeff * ( VS_0T - VS_0t );

                // adding jump components (factor-independent)
                double int1 = 0.5 * (Maths::square(tempT) - Maths::square(tempt));
                double sum1 = 0.0;
                Complex alpha_J;
                // \int_t^T \beta_{n+3}(s) ds
                Complex LastBetaInt;
                if (! Maths::isZero(k_J)) {
                     Complex Ccoeff = - phi / k_J;
                     Complex Acoeff = (u[size+2]-Ccoeff);
                     LastBetaInt = Ccoeff * tau + Acoeff / k_J * (1.0 - exp(-k_J * tau));
                } else {
                    Complex Ccoeff = u[size+2] - phi * tempT;
                    Complex Acoeff = phi;
                    LastBetaInt = Ccoeff * tau + int1 * Acoeff;
                }
                Complex coeffRatio = (1.0 - mu_V * u[size+2]) / (1.0 - mu_V * u[size+2] + mu_V * phi * tau);
                double coeffsq = Complex::absSquare(coeffRatio);
                //test for infinite values
                double argcoeff;
                if( !Maths::isZero( coeffRatio.imag() ) && Maths::isZero( coeffRatio.real() )  ){
                    argcoeff = Maths::sign(coeffRatio.imag()) * Maths::PI * 0.5;
                }
                else{
                    argcoeff = atan(coeffRatio.imag() / coeffRatio.real());
                }
                double rcoeff = 0.5 * log(coeffsq);
                Complex logcoeff = Complex(rcoeff,argcoeff);
                if (! Maths::isZero(k_J)) {
                    Complex Ccoeff = - phi / k_J;
                    Complex jumpLogCoeff = 1.0 - mu_V * Ccoeff + mu_V * (Ccoeff - u[size+2]) * exp(-k_J * tau);
                    double jumpLogSq = Complex::absSquare(jumpLogCoeff);
                    //test for infinite values
                    double jumpLogArg;
                    if( !Maths::isZero( jumpLogCoeff.imag() ) && Maths::isZero( jumpLogCoeff.real() )  ){
                        jumpLogArg = Maths::sign( jumpLogCoeff.imag() ) * Maths::PI * 0.5;
                    }
                    else{
                        jumpLogArg = atan(jumpLogCoeff.imag() / jumpLogCoeff.real());
                    }
                    double rjumpLog = 0.5 * log(jumpLogSq);
                    Complex jumpLog = Complex(rjumpLog,jumpLogArg);
                    Complex alpha_J = - lambda_V * tau + lambda_V * tau / (1.0 - mu_V * Ccoeff) + lambda_V / (k_J * (1.0 - mu_V * Ccoeff)) * jumpLog;
                } else {
                    if (! Maths::isZero(mu_V) && ! Maths::isZero(phi)){
                        alpha_J = - lambda_V * tau - lambda_V / (mu_V * phi) * logcoeff;                 
                    } else{
                        alpha_J = (lambda_V * mu_V * u[size+2]) / (1.0 - mu_V * u[size+2]) * tau;
                    }
                }

                // computing alpha contribution from jumps in the spot
                Complex alpha_JS;
                Complex spotjump1 = exp(gamma_S * u[size] + 0.5 * Maths::square(delta_S * u[size]));
                Complex spotjump2 = 1.0 / sqrt(1.0 - 2.0 * u[size+1] * Maths::square(delta_S)) * exp(u[size+1] * Maths::square(gamma_S) / (1.0 - 2.0 * u[size+1] * Maths::square(delta_S)));
                alpha_JS = lambda_S * tau * (spotjump1 + spotjump2 - 2.0);
                for (int j=0; j<size; j++) {
                    double corr = volvscur->correlation[j];
                    double lamb = volvscur->lambda[j];
                    sum1 += lamb * (Maths::square(corr));
                }

                alpha += int1 * lambda_V * mu_V * (0.5 * sum1 *Maths::square(u[size])) - 0.5 * lambda_S * chi_S * sum1 * Maths::square(u[size]) * tau - lambda_V * mu_V * LastBetaInt + alpha_JS + alpha_J;


                //gives the lower int and higher int for the PC interpolation of the initial FVSC
		        int nLow = fvscpcFloor->tenorLocation( tempt ); //computation from t 
		        int nUp = fvscpcFloor->tenorLocation( tempT );  //to T 
		    
                if (Maths::isZero(nLow-nUp)){
                    //range of integration used throughout the computation
                    double ta, tb;

                    for(int k=0; k<size; k++){
                        //k-th factor parameters
                        double corr = volvscur->correlation[k]; //correlation
                        double mRR = volvscur->meanReversRate[k]; //mean reverse rate
                        double vVol = volvscur->volVol[k]; //vol of vol
                        double lamb = volvscur->lambda[k]; //lambda
                    
                        //range of integration              
                        ta = tempT - tempt; 
                        tb = 0.0; 
                        
                        //alpha
                        alpha+= 0.5 * fvscpcFloor->values[nLow] * lamb * calcBeta2IntegrUniDim(ta,
                                                                                           tb,
                                                                                           mRR,
                                                                                           vVol,
                                                                                           corr,
                                                                                           u[size], //corresponding to spot
                                                                                           u[size+1], //corresponding to real var
                                                                                           u[k]); //used as a terminal condition    
                        
                    // adding jump components
                        // T * F(T)
                        Complex jump1 = tempT * (calcBeta2IntegrUniDim(tau,0.0,mRR,vVol,corr,u[size],u[size+1],u[k]) + tau * Maths::square(corr * u[size]));
                        // \int_t^T a^2 u \beta_j^2(u) du   =  T * F(T) - \int_t^T F(u) du
                        Complex jump2 = jump1 - (calcBeta2IntegrUniDimPt6(tau,mRR,vVol,corr,u[size],u[size+1],u[k])); 
                        // final contribution from jumps in the spot
                        Complex jump3 = (-0.5) * lambda_S * chi_S * lamb * (calcBeta2IntegrUniDim(tau,0.0,mRR,vVol,corr,u[size],u[size+1],u[k]));

                        alpha += (-0.5) * lambda_V * mu_V * lamb * jump2 + jump3;
                    
                    }
                }
                else{
                    //range of integration used throughout the computation
                    double ta, tb;

                    for(int k=0; k<size; k++){
                        //k-th factor parameters
                        double corr = volvscur->correlation[k]; //correlation
                        double mRR = volvscur->meanReversRate[k]; //mean reverse rate
                        double vVol = volvscur->volVol[k]; //vol of vol
                        double lamb = volvscur->lambda[k]; //lambda
                    
                        //range of integration              
                        ta = tempT - tempt; 
                        tb = tempT - fvscpcFloor->tenors[nLow+1]; 
                        
                        //alpha
                        alpha+= 0.5 * fvscpcFloor->values[nLow] * lamb * calcBeta2IntegrUniDim(ta,
                                                                                                   tb,
                                                                                                   mRR,
                                                                                                   vVol,
                                                                                                   corr,
                                                                                                   u[size], //corresponding to spot
                                                                                                   u[size+1], //corresponding to real var
                                                                                                   u[k]); //used as a terminal condition
                        int compt = nLow + 1;
                        while(compt<nUp){
                            //range of integration                
                            ta = tempT - fvscpcFloor->tenors[compt]; 
                            tb = tempT - fvscpcFloor->tenors[compt+1]; 
                            
                                //alpha
                                alpha += 0.5 * fvscpcFloor->values[compt] * lamb * calcBeta2IntegrUniDim(ta,
                                                                                                 tb,
                                                                                                 mRR,
                                                                                                 vVol,
                                                                                                 corr,
                                                                                                 u[size], //corresponding to spot
                                                                                                 u[size+1], //corresponding to real var
                                                                                                 u[k]); //used as a terminal condition
                                compt++;
                        }
                        //range of integration
                            ta = tempT - fvscpcFloor->tenors[nUp]; 
                            tb = 0.0; 

                            //alpha
                            alpha += 0.5 * fvscpcFloor->values[nUp] * lamb * calcBeta2IntegrUniDim(ta,
                                                                                           tb,
                                                                                               mRR,
                                                                                               vVol,
                                                                                           corr,
                                                                                           u[size], //corresponding to spot
                                                                                           u[size+1], //corresponding to real var
                                                                                           u[k]); //used as a terminal condition
                        // adding jump components
                        // T * F(T)
                        Complex jump1 = tempT * (calcBeta2IntegrUniDim(tau,0.0,mRR,vVol,corr,u[size],u[size+1],u[k]) + tau * Maths::square(corr * u[size]));
                        // \int_t^T a^2 u \beta_j^2(u) du   =  T * F(T) - \int_t^T F(u) du
                        Complex jump2 = jump1 - (calcBeta2IntegrUniDimPt6(tau,mRR,vVol,corr,u[size],u[size+1],u[k])); 
                        // final contribution from jumps in the spot
                        Complex jump3 = (-0.5) * lambda_S * chi_S * lamb * (calcBeta2IntegrUniDim(tau,0.0,mRR,vVol,corr,u[size],u[size+1],u[k]));

                        alpha += (-0.5) * lambda_V * mu_V * lamb * jump2 + jump3;
                    
                    }
                }
            }
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }   
};

/*** end of class CalcAlphaBeta (helper class) ***/
/*************************************************/


/*************************/
/*** class VolVSCurve  ***/

//MODIFIED GAD 15-Oct-2006

//default values for parameters
const double VolVSCurve::DefaultVal::correlation = -0.5;
const double VolVSCurve::DefaultVal::volVol = 0.1;
const double VolVSCurve::DefaultVal::meanReversRate = 0.5;
const double VolVSCurve::DefaultVal::lambda = 1.0;
const double VolVSCurve::DefaultVal::crashMRR = 0.5;
const double VolVSCurve::DefaultVal::crashRate = 0.0;
const double VolVSCurve::DefaultVal::crashSize = 0.0;
const double VolVSCurve::DefaultVal::spotcrashRate = 0.0;
const double VolVSCurve::DefaultVal::spotcrashSize = 0.0;
const double VolVSCurve::DefaultVal::spotcrashUncertainty = 0.0;

//range for parameters
const Range VolVSCurve::RangeDef::correlation = Range( OpenBoundary(-1.0), OpenBoundary(1.0) );
const Range VolVSCurve::RangeDef::volVol = Range( OpenBoundary(0.0), Infinity(Infinity::Plus) ); //modified
const Range VolVSCurve::RangeDef::meanReversRate = Range( OpenBoundary(0.0), Infinity(Infinity::Plus) ); //modified
const Range VolVSCurve::RangeDef::lambda = Range( ClosedBoundary(0.0), ClosedBoundary(1.0) ); //modified
const Range VolVSCurve::RangeDef::crashMRR = Range( ClosedBoundary(0.0), Infinity(Infinity::Plus) );
const Range VolVSCurve::RangeDef::crashRate = Range( ClosedBoundary(0.0), Infinity(Infinity::Plus) );
const Range VolVSCurve::RangeDef::crashSize = Range( ClosedBoundary(0.0), Infinity(Infinity::Plus) );
const Range VolVSCurve::RangeDef::spotcrashRate = Range( ClosedBoundary(0.0), Infinity(Infinity::Plus) );
const Range VolVSCurve::RangeDef::spotcrashSize = Range( Infinity(Infinity::Minus), ClosedBoundary(0.0) ); //modified
const Range VolVSCurve::RangeDef::spotcrashUncertainty = Range( ClosedBoundary(0.0), Infinity(Infinity::Plus) );

//END MODIFIED GAD

//registration
CClassConstSP const VolVSCurve::TYPE = CClass::registerClassLoadMethod("VolVSCurve", typeid(VolVSCurve), load);
DEFINE_TEMPLATE_TYPE(VolVSCurveArray);
CClassConstSP const VolVSCurve::VSCurveVolParam::TYPE = CClass::registerClassLoadMethod("VolVSCurve::VSCurveVolParam", typeid(VSCurveVolParam), load);


//MODIFIED GAD 15-Oct-2006

//validation
//now empty; simple range tests moved to update()
void VolVSCurve::validatePop2Object(){
    //no test on var jump size, additional test could be added to guarantee that
    //diffusion state variables Qi can become negative and take the FVSC on the downside
}

//END MODIFIED GAD

//update for calibration
void VolVSCurve::update(){
    
    try{
        // range test, i.e. code from old validate
        
        //MODIFIED GAD 15-Oct-2006

        //double volVolMin = volVol[0], volVolMax = volVol[0];
        double lambdaMin = lambda[0], lambdaMax = lambda[0], lambdaSum = lambda[0];
        //double correlationMin = correlation[0], correlationMax = correlation[0];

        for(int i=1; i<n; i++){
            //volVolMin = min(volVolMin, volVol[i]);
            //volVolMax = max(volVolMax, volVol[i]);

            //lambdaMin = min(lambdaMin, lambda[i]);
            //lambdaMax = max(lambdaMax, lambda[i]);
            lambdaSum += lambda[i];

            //correlationMin = min(correlationMin, correlation[i]);
            //correlationMax = max(correlationMax, correlation[i]);    
        }

        Calibrator::IAdjustable::checkRange(this);
        if( (n<1) || (n>4) )
            throw ModelException("VolVSCurve::update","number of factors must be in 1 - 4");
        if( !Maths::isZero(lambdaSum - 1.0) )
            throw ModelException("VolVSCurve::update","weights must sum up to 1");        
        
        //not necessary if checked from calibrator range
        /*
        
        if( Maths::isPositive( - volVolMin) || Maths::isZero(volVolMin) ) //modified
		    throw ModelException("VolVSCurve::update","vol of vol must be strictly positive");
        if( Maths::isPositive(lambdaMax - 1.0) || Maths::isPositive(-lambdaMin) ) //modified
            throw ModelException("VolVSCurve::update","weight outside correct range");
        if(  Maths::isPositive(-1.0 + correlationMax) || Maths::isPositive(-correlationMin - 1.0) ) //modified
            throw ModelException("VolVSCurve::update","correlation outside correct range");
        if( Maths::isPositive(-crashSize) ) //modified
		    throw ModelException("VolVSCurve::update","crashSize must be positive");
        if(Maths::isPositive(-crashRate)) //modified
		    throw ModelException("VolVSCurve::update","crashRate must be positive");
        if( Maths::isPositive(-crashMRR) ) //modified
		    throw ModelException("VolVSCurve::update","crashMRR must be positive");
        if( Maths::isPositive(-spotcrashRate) ) //modified
		    throw ModelException("VolVSCurve::update","spotcrashRate must be positive");
        if( Maths::isPositive(-spotcrashUncertainty) ) //modified
		    throw ModelException("VolVSCurve::update","spotcrashUncertainty must be strictly positive");
        if( Maths::isPositive(spotcrashSize) ) //added
		    throw ModelException("VolVSCurve::update","spotcrashSize must be negative");
        */

        //END MODIFIED GAD

        // code from old getMarket
        //construction of transient fields
    
        //tenorCurve expressed in DateTime
        DateTimeArray tenorCurveBis(tenorCurve.size());
        for(int i=0; i<tenorCurve.size(); i++){
            ExpirySP thisExpiry = tenorCurve[i];
            MaturityPeriodSP period = MaturityPeriodSP::dynamicCast(thisExpiry);
            if(period.get()){
                //Convert to EOD
                thisExpiry = MaturityTimePeriodSP(new MaturityTimePeriod(period->toString(), DateTime::END_OF_DAY_TIME)); 
            } 
            tenorCurveBis[i] = thisExpiry->toDate(baseDate);
        }
        
        // Validate input tenors`
        DateTime::ensureIncreasing(tenorCurveBis, "Variance Swap maturities", true);
        if(baseDate.getDate() >= tenorCurveBis[0].getDate() ) {
            throw ModelException( "Base date (excluding time) " + baseDate.toString() + 
                " has to be before first Variance swap maturity " +
                tenorCurveBis[0].toString() );
        }

        //tenorCurve expressed in double and used to build VSCAff and FVSCPC
        DoubleArray tenorCurveBisBis(tenorCurve.size());
        for(int j=0; j<tenorCurve.size(); j++){
            tenorCurveBisBis[j] = timeMetric->yearFrac(baseDate, tenorCurveBis[j]);    
        }

        //trading time vs PPY
        const HolidayWrapper assetHols = timeMetric->getHolidays();

        //internal computation of VSCurve
        if( Maths::isZero(tenorCurve.size()) ){
            //VSCAFF Construction (first value is added / left end of the curve)
            DoubleArray vscaffTenorTemp(tenorCurve.size()+1), vscaffValueTemp(tenorCurve.size()+1); 
            //vscaff construction
            vscaffTenorTemp[0] = 0.0;
            vscaffValueTemp[0] = 0.0; //TO BE CHANGED
            for(int l=1; l<tenorCurve.size()+1; l++){
                vscaffTenorTemp[l] = tenorCurveBisBis[l-1];
                vscaffValueTemp[l] = 0.15 * 0.15 * tenorCurveBisBis[l-1]; //TO BE CHANGED
            }
            //VS curve
            vscaff = VSCAffSP( new VSCAff(vscaffTenorTemp, vscaffValueTemp) );
        }
        //VSCurve is input
        else{
            //reference 0D EOD point
            MaturityTimePeriod ref0("0D", DateTime::END_OF_DAY_TIME);
            //into date
            DateTime ref0Date = (ref0).toDate(baseDate);
            //into double
            double ref0Double = timeMetric->yearFrac(baseDate, ref0Date); 

            //fwd starting VS
            if(isFwdStartingVS){
                //not a real fwd starting VS, baseDate = 0D EOD
                if(Maths::isZero(ref0Double)){
                    //VSCAFF Construction (first value is added / left end of the curve)
                    DoubleArray vscaffTenorTemp(tenorCurve.size()+1), vscaffValueTemp(tenorCurve.size()+1);
                    //trading time scaling factors, i.e. expN / PPY * 1/T , for each input tenor T
                    DoubleArray tradingTimeScaling(tenorCurve.size());
                    
                    //zero point
                    vscaffTenorTemp[0] = 0.0;
                    vscaffValueTemp[0] = 0.0;

                    //other points
                     for(int m=1; m<tenorCurve.size()+1; m++){
                         int expectedN = assetHols->businessDaysDiff(baseDate,tenorCurveBis[m-1]);
                         tradingTimeScaling[m-1] = (double)expectedN / (double) PPY;
                         vscaffTenorTemp[m] = tenorCurveBisBis[m-1];
                         vscaffValueTemp[m] = valueCurve[m-1] * valueCurve[m-1] * tradingTimeScaling[m-1];
                     }
                     //VS curve
                     vscaff = VSCAffSP( new VSCAff(vscaffTenorTemp, vscaffValueTemp) );
                }
                //real fwd starting VS
                else{
                    //VSCAFF Construction (first value is added / left end of the curve)
                    DoubleArray vscaffTenorTemp(tenorCurve.size()+2), vscaffValueTemp(tenorCurve.size()+2);
                    //trading time scaling factors, i.e. expN / PPY * 1/T , for each input tenor T
                    DoubleArray tradingTimeScaling(tenorCurve.size());

                    //zero point
                    vscaffTenorTemp[0] = 0.0;
                    vscaffValueTemp[0] = 0.0;
                
                    //first point (1D VS)
                    vscaffTenorTemp[1] = ref0Double;
                    vscaffValueTemp[1] = valueCurve[0] * valueCurve[0] * ref0Double;
                    //other points
                    for(int m=2; m<tenorCurve.size()+2; m++){
                        int expectedN = assetHols->businessDaysDiff(baseDate,tenorCurveBis[m-2]);
                        tradingTimeScaling[m-2] = (double)expectedN / (double) PPY;
                        vscaffTenorTemp[m] = tenorCurveBisBis[m-2];
                        vscaffValueTemp[m] = valueCurve[m-2] * valueCurve[m-2] * tradingTimeScaling[m-2] + valueCurve[0] * valueCurve[0] * ref0Double;
                    }
                    //VS curve
                    vscaff = VSCAffSP( new VSCAff(vscaffTenorTemp, vscaffValueTemp) );
                }
            }
            //starting VS
            else{
                //VSCAFF Construction (first value is added / left end of the curve)
                DoubleArray vscaffTenorTemp(tenorCurve.size()+1), vscaffValueTemp(tenorCurve.size()+1);
                //trading time scaling factors, i.e. expN / PPY * 1/T , for each input tenor T
                DoubleArray tradingTimeScaling(tenorCurve.size());

                //zero point
                vscaffTenorTemp[0] = 0.0;
                vscaffValueTemp[0] = 0.0;

                for(int m=1; m<tenorCurve.size()+1; m++){
                    int expectedN = assetHols->businessDaysDiff(baseDate,tenorCurveBis[m-1]);
                    tradingTimeScaling[m-1] = (double)expectedN / (double) PPY;
                    vscaffTenorTemp[m] = tenorCurveBisBis[m-1];
                    vscaffValueTemp[m] = valueCurve[m-1] * valueCurve[m-1] * tradingTimeScaling[m-1];
                }
                //VS curve
                vscaff = VSCAffSP( new VSCAff(vscaffTenorTemp, vscaffValueTemp) );
            }
        }

        //FVSCPC Construction (last value is added / right end of the curve)
        int sizeVSC = vscaff->tenors.size();
        DoubleArray fvscpcTenorTemp(sizeVSC), fvscpcValueTemp(sizeVSC);
        for(int lBis=0; lBis < sizeVSC-1; lBis++){
            fvscpcTenorTemp[lBis] = vscaff->tenors[lBis];
            fvscpcValueTemp[lBis] = ( (vscaff->values[lBis+1]) - (vscaff->values[lBis]) ) / ( (vscaff->tenors[lBis+1]) - (vscaff->tenors[lBis]) );
        }
        //for the last value, we use the same value so as to have a continuous derivative
        fvscpcTenorTemp[sizeVSC-1] = vscaff->tenors[sizeVSC-1];
        fvscpcValueTemp[sizeVSC-1] = fvscpcValueTemp[sizeVSC-2];

        //definition
        fvscpc = FVSCPCSP( new FVSCPC(fvscpcTenorTemp, fvscpcValueTemp) );

        double limit = timeMetric->yearFrac( baseDate, (tenorCurve.back())->toDate(baseDate) );
        fvscpcFloor = fvscpc->floor(limit,XMAX);      
    }
    catch(exception& e){
        throw ModelException(e, "VolVSCurve::update");
    }
}

bool VolVSCurve::sensShift(Theta* shift){
    static const string method = "VolVSCurve::sensShift";
    try{
        // Call parent
        return VolBaseParam::sensShift(shift);

        // Build DateTimeArray from baseDate, tenors
        // Delete input (tenors , values) if baseDate.getDate() <= tenors.getDate()
        int nSize = tenorCurve.size();
        for(int i = 0; i<nSize; i++){
            ExpirySP thisExpiry = tenorCurve.front();
            DateTime thisExpiryDate = thisExpiry->toDate(baseDate);
            bool isGreater = baseDate.getDate() <= thisExpiryDate.getDate();
            if( isGreater ){
                tenorCurve.erase(tenorCurve.begin());
                valueCurve.erase(valueCurve.begin());
            }
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

//implied vol in VolParam
void VolVSCurve::VSCurveVolParam::ComputeImpVol(const CVolBase*          vol,
                                                const CLatticeDouble&    strikes,
                                                const DateTimeArray&     maturities,
                                                CLatticeDouble&          impV) const{
    // turn the vol into what we must have
    const VolVSCurve* myVol = static_cast<const VolVSCurve*>(vol);
    // then just pass through the parameterised vol
    myVol->ComputeImpVol(strikes, maturities, impV);
}

//implementation in VolVSCurve
void VolVSCurve::ComputeImpVol(const CLatticeDouble&      strikes,
                               const DateTimeArray&       maturities,
                               CLatticeDouble&            impV) const{
    static const string routine("VolVSCurve::ComputeImpVol");
    throw ModelException(routine, "Not supported");
    if ((maturities.size() != strikes.size()) || (maturities.size() != impV.size())){
        throw ModelException(routine, "Size mismatch between strikes ("+ 
                             Format::toString(strikes.size()) +
                             "), maturities ("+ 
                             Format::toString(maturities.size())+
                             ") and impV ("+ 
                             Format::toString(impV.size())+ ")");
    }
    for (int iMat = 0; iMat < maturities.size(); iMat++){
        if (strikes[iMat].size() != impV[iMat].size()){
            throw ModelException(routine, "Size mismatch between strikes"
                                 " & maturities for Mat " +
                                 maturities[iMat].toString() +
                                 " (n "+ Format::toString(iMat) + ")");
        }
        for (int iStrike = 0; iStrike < strikes[iMat].size(); iStrike ++){
            impV[iMat][iStrike] = 0.0;
        }
    }
}

//Matrix Vol Surface in VolPAram
VolSurface* VolVSCurve::VSCurveVolParam::spotVolSurfaceFromStrikes(const CVolBase*       vol,
														           const CDoubleArray&   strikes) const{
	// turn the vol into what we must have
    const VolVSCurve* myVol = static_cast<const VolVSCurve*>(vol);
    // then just pass through the parameterised vol
    return myVol->spotVolSurfaceFromStrikes(strikes);
}

//implementation in VolVSCurve
VolSurface* VolVSCurve::spotVolSurfaceFromStrikes(const CDoubleArray& strikes) const{
    static const string routine("VolVSCurve::spotVolSurfaceFromStrikes");
    try{
        throw ModelException(routine, "Not supported");

        const VolSurface* backbone = getBackboneSurface();
        const DateTimeArray& dates = backbone->getDates();
        CDoubleMatrix matrix(strikes.size(),
                             dates.size());
        for (int iStrike = 0; iStrike < strikes.size(); iStrike++){
            for (int iMat = 0; iMat < dates.size(); iMat++){
                matrix[iStrike][iMat] = 0.0;
            }
        }
        /** for performance need constructor that takes in
            cached values (to do) */
        VolSurface* volSurf = new VolSurface(getName(),
                                             timeMetric.get(),
                                             strikes,
                                             matrix,
                                             backbone->getExpiries().get(),
                                             baseDate);
        return volSurf;
    }
    catch (exception& e){
        throw ModelException(e, routine);
    }
}

//load function
void VolVSCurve::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VolVSCurve, clazz);
    SUPERCLASS(VolBaseParam);
    IMPLEMENTS(Calibrator::IAdjustable);
    EMPTY_SHELL_METHOD(defaultCtor);

	FIELD(n, "n");
    FIELD(meanReversRate, "meanReversRate");
	FIELD(volVol, "volVol");
	FIELD(correlation, "correlation");
	FIELD(lambda, "lambda");
    FIELD(crashMRR,"crashMRR");
    FIELD_MAKE_OPTIONAL(crashMRR);
    FIELD(crashRate,"crashRate");
    FIELD_MAKE_OPTIONAL(crashRate);
    FIELD(crashSize,"crashSize");
    FIELD_MAKE_OPTIONAL(crashSize);
    FIELD(spotcrashRate,"spotcrashRate");
    FIELD_MAKE_OPTIONAL(spotcrashRate);
    FIELD(spotcrashSize,"spotcrashSize");
    FIELD_MAKE_OPTIONAL(spotcrashSize);
    FIELD(spotcrashUncertainty,"spotcrashUncertainty");
    FIELD_MAKE_OPTIONAL(spotcrashUncertainty);
    FIELD(tenorCurve, "tenorCurve");
    FIELD(isFwdStartingVS,"isFwdStartingVS");
    FIELD_MAKE_OPTIONAL(isFwdStartingVS);

    //optional (if not rebuild internally)
    FIELD(valueCurve, "valueCurve");
    FIELD_MAKE_OPTIONAL(valueCurve);
 
    //transient, declared and used without appearing on the spreadsheet 
    FIELD(baseDate, "");
    FIELD_MAKE_TRANSIENT(baseDate);
    
    FIELD(fvscpc,"");
    FIELD_MAKE_TRANSIENT(fvscpc);
    
    FIELD(vscaff,"");
    FIELD_MAKE_TRANSIENT(vscaff);

    FIELD(fvscpcFloor,"");
    FIELD_MAKE_TRANSIENT(fvscpcFloor);


    
	//add our fields and their ranges to central list
    Calibrator::IAdjustable::registerField(
        clazz, "meanReversRate",
        new Range(VolVSCurve::RangeDef::meanReversRate));

    Calibrator::IAdjustable::registerField(
        clazz, "volVol", 
        new Range(VolVSCurve::RangeDef::volVol));
    
	Calibrator::IAdjustable::registerField(
        clazz, "correlation",
        new Range(VolVSCurve::RangeDef::correlation));

    Calibrator::IAdjustable::registerField(
        clazz, "lambda",
        new Range(VolVSCurve::RangeDef::lambda));

    Calibrator::IAdjustable::registerField(
        clazz, "crashMRR",
        new Range(VolVSCurve::RangeDef::crashMRR));

    Calibrator::IAdjustable::registerField(
        clazz, "crashRate",
        new Range(VolVSCurve::RangeDef::crashRate));

    Calibrator::IAdjustable::registerField(
        clazz, "crashSize",
        new Range(VolVSCurve::RangeDef::crashSize));

    Calibrator::IAdjustable::registerField(
        clazz, "spotcrashRate",
        new Range(VolVSCurve::RangeDef::spotcrashRate));

    Calibrator::IAdjustable::registerField(
        clazz, "spotcrashSize",
        new Range(VolVSCurve::RangeDef::spotcrashSize));

    Calibrator::IAdjustable::registerField(
        clazz, "spotcrashUncertainty",
        new Range(VolVSCurve::RangeDef::spotcrashUncertainty));
}



VolVSCurve::VolVSCurve() :
VolBaseParam(TYPE),
crashMRR(VolVSCurve::DefaultVal::crashMRR),
crashRate(VolVSCurve::DefaultVal::crashRate),
crashSize(VolVSCurve::DefaultVal::crashSize),
spotcrashRate(VolVSCurve::DefaultVal::spotcrashRate),
spotcrashSize(VolVSCurve::DefaultVal::spotcrashSize),
spotcrashUncertainty(VolVSCurve::DefaultVal::spotcrashUncertainty),
isFwdStartingVS(true){}



/*
METHOD IMPLEMENTED IN VOLSVJ

// Build the parameterised vol and cache any values
void VolVSCurve::buildCache(){
    const VolSurface* backbone = getBackboneSurface();
    baseDate = backbone->getBaseDate();
}

*/

//METHOD IMPLEMENTED IN SVCJ

// populate from market cache
void VolVSCurve::getMarket(const IModel* model, const MarketData* market){
    try{
        // calls parent method first
        // populates parent
        VolBaseParam::getMarket(model, market);
    
        //const VolSurface* backbone = getBackboneSurface();
        //returns a constant reference to surface to be used for the backbone
        //VolSurfaceSP backbone = VolSurfaceSP(VolSurfaceSP::dynamicCast(market->GetData(getName(), VolSurface::TYPE)));
        //backbone->getMarket(model,market);  // put holiday and so on into volsurface.
        
        //populates itself
        baseDate = market->GetReferenceDate();       
        
        update();
    }
    catch(exception& e){
        throw ModelException(e, "VolVSCurve::getMarket curve construction");
    }
}

/** Called after adjustments have been made to fields (eg calibrator) */
void VolVSCurve::fieldsUpdated(const CFieldArray& fields){
    update();
}

//Laplace transform to be computed

/* Started Log Return */
Complex VolVSCurve::scalelessCumulant(const StFourierProcessLogRtn& process,
                                      const StFourierProductLogRtn&	product, 
                                      const Complex&		        z, 
                                      const DateTime&      			matDate) const{
    static const string method = "VolVSCurve::scalelessCumulant StFourierProcessLogRtn";
    Complex res;
    try{
        //used in CalcAlpha
        Complex alpha;
        ComplexArray beta(n+3);
        
        //terminal conditions
        ComplexArray u(n+3);
        
        //variance curve state variables 
        for(int i=0; i<n; i++){
			//case vVol[i] = 0
			if( Maths::isZero(volVol[i]) ){
				throw ModelException(method, "zero vol of vol not supported");
			}
			//vVol[i] != 0
			else{
				u[i] = z * correlation[i] / volVol[i];
			}
        }
        u[n] =  z; //spot 
        u[n+1] = 0.0; //real var
        u[n+2] = 0.0; //jumps

        CalcAlphaBeta::calcJointLapAlphaBeta(baseDate,
                                             matDate,
                                             u,
							                 alpha,
                                             beta,
                                             this);
        res = alpha;
        return res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Fwd starting Log Return  */
Complex VolVSCurve::scalelessCumulant(const FwdStFourierProcessLogRtn&  process,
                                      const FwdStFourierProductLogRtn&  product, 
                                      const Complex&                    z, 
                                      const DateTime&					matDate) const{
    static const string method = "VolVSCurve::scalelessCumulant FwdStFourierProcessLogRtn";
    Complex res;
    try{
        /* Calculate alpha and beta from start date till maturity */
        
        //used in CalcAlpha
        Complex alpha_tT;
        ComplexArray beta_tT(n+3);

        //Terminal conditions
        ComplexArray u_tT(n+3);
        for(int i=0; i<n; i++){
			//case vVol[i] = 0
			if( Maths::isZero(volVol[i]) ){
				throw ModelException(method, "zero vol of vol not supported");			
			}
			//vVol[i] != 0
			else{
				u_tT[i] = z * correlation[i] / volVol[i];			
			}
        }
        u_tT[n] = z; //spot
        u_tT[n+1] = 0.0; //real var
        u_tT[n+2] = 0.0; //jumps

        CalcAlphaBeta::calcJointLapAlphaBeta(product.getStartDate(),
                                             matDate,
                                             u_tT,
                                             alpha_tT,
                                             beta_tT,
                                             this);
        
        /* Then, from today till start date */
        
        //used in CalcAlpha
        Complex alpha_0t;
        ComplexArray beta_0t(n+3);

        //Terminal conditions
        ComplexArray u_0t(n+3);
        for(int j=0; j<n; j++){
			//case vVol[j] = 0
			if( Maths::isZero(volVol[j]) ){
				throw ModelException(method, "zero vol of vol not supported");
			}
			//vVol[j] != 0
			else{
				u_0t[j] = beta_tT[j] - z * correlation[j] / volVol[j];
			}
        }
        u_0t[n] = beta_tT[n] - z; //spot
        u_0t[n+1] = beta_tT[n+1]; //real var
        u_0t[n+2] = beta_tT[n+2]; //jumps

        CalcAlphaBeta::calcJointLapAlphaBeta(baseDate,
                                             product.getStartDate(),
                                             u_0t,
                                             alpha_0t,
                                             beta_0t,
                                             this);

		res = alpha_tT + alpha_0t;
        return res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Started integrated variance */
Complex VolVSCurve::cumulant(const StFourierProcessIntVar&  process,
                             const StFourierProductIntVar&  product, 
                             const Complex&				    z, 
                             const DateTime&				matDate) const{
    static const string method = "VolVSCurve::cumulant StFourierProcessIntVar";
    Complex res;
    try{
        //used in CalcAlpha
        Complex alpha;
        ComplexArray beta (n+3);

        //Terminal conditions
		ComplexArray u(n+3);
        for(int j=0; j<n; j++){
            u[j] = 0.0;
        }
        u[n] = 0.0; //spot
        u[n+1] = z; //real var, no maturity adjustment
        u[n+2] = 0.0;

        CalcAlphaBeta::calcJointLapAlphaBeta(baseDate,
                                             matDate,
                                             u,
                                             alpha,
                                             beta,
                                             this);

		res = alpha;
        return res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Forward Starting integrated variance */
Complex VolVSCurve::cumulant(const FwdStFourierProcessIntVar&	process,
                             const FwdStFourierProductIntVar&	product, 
                             const Complex& 					z, 
                             const DateTime&					matDate) const{        
    static const string method = "VolVSCurve::cumulant FwdStFourierProcessIntVar";
    Complex res;
    try{
        /* Calculate alpha and beta from start date till maturity */
        
        //used in CalcAlphaBeta 
        Complex alpha_tT;
        ComplexArray beta_tT(n+3);

        //Terminal conditions
        ComplexArray u_tT(n+3);
        for(int j=0; j<n; j++){
            u_tT[j] = 0.0;
        }
        u_tT[n] = 0.0; //spot
        u_tT[n+1] = z; //real var, no maturity adjustment
        u_tT[n+2] = 0.0; //jumps

		CalcAlphaBeta::calcJointLapAlphaBeta(product.getStartDate(),
                                             matDate,
                                             u_tT,
                                             alpha_tT,
                                             beta_tT,
                                             this);

        /* Then, from today till start date */

        //used in CalcAlphaBeta
		Complex alpha_0t;
        ComplexArray beta_0t(n+3);
        
		//Terminal conditions
        ComplexArray u_0t(n+3);
        for(int k=0; k<n; k++){
            u_0t[k] = beta_tT[k];
        }
        u_0t[n] = beta_tT[n]; //spot
        u_0t[n+1] = beta_tT[n+1] - z; //real var, no maturity adjustment
        u_0t[n+2] = beta_tT[n+2]; //jumps
        
		CalcAlphaBeta::calcJointLapAlphaBeta(baseDate,
                                             product.getStartDate(),
                                             u_0t,
                                             alpha_0t,
                                             beta_0t,
                                             this);

        res = alpha_tT + alpha_0t;
        return res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Started Quadratic Variation  */
Complex VolVSCurve::cumulant(const StFourierProcessQuadVar&    process,
                             const StFourierProductQuadVar&    product, 
                             const Complex&	                   z, 
                             const DateTime&				   matDate) const{
    static const string method = "VolVSCurve::cumulant StFourierProcessQuadVar";
    Complex res;
    try{
        //used in CalcAlpha
        Complex alpha;
        ComplexArray beta(n+3);

        //Terminal conditions
		ComplexArray u(n+3);
        for(int j=0; j<n; j++){
            u[j] = 0.0;
        }
        u[n] = 0.0; //spot
        u[n+1] = z; //real var, no maturity adjustment
        u[n+2] = 0.0;

        CalcAlphaBeta::calcJointLapAlphaBeta(baseDate,
                                             matDate,
                                             u,
                                             alpha,
                                             beta,
                                             this);
		res = alpha;
        return(res);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

//expectation
double VolVSCurve::expectation(const StFourierProcessQuadVar&   process,
                               const StFourierProductQuadVar&   product, 
                               const DateTime&                  matDate) const{    
    static const string method = "VolVSCurve::cumulant";
    try{
        //result
        double res;
        
        double matDateDouble = timeMetric->yearFrac(baseDate, matDate);
        if( Maths::isZero(matDateDouble) ){
            throw ModelException("matDate can't be zero", method);
        }
        else{
            res = vscaff->yield(matDateDouble) / (matDateDouble);
            double lambda_S = spotcrashRate;
            double k_S = spotcrashSize;
            double delta_S = spotcrashUncertainty;

            double chi_S = log(1.0 + k_S) + 0.5 * Maths::square(delta_S);
            res += lambda_S * chi_S;
            return res;
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Forward Starting Quadratic Variation */
Complex VolVSCurve::cumulant(const FwdStFourierProcessQuadVar&	process,
                             const FwdStFourierProductQuadVar&	product, 
                             const Complex&						z, 
                             const DateTime&					matDate) const {
    static const string method = "VolVSCurve::cumulant FwdStFourierProcessQuadVar";
    Complex res;
    try{
        /* Calculate alpha and beta from start date till maturity */
        
        //used in CalcAlphaBeta 
        Complex alpha_tT;
        ComplexArray beta_tT(n+3);

        //Terminal conditions
        ComplexArray u_tT(n+3);
        for(int j=0; j<n; j++){
            u_tT[j] = 0.0;
        }
        u_tT[n] = 0.0; //spot
        u_tT[n+1] = z; //real var, no maturity adjustment
        u_tT[n+2] = 0.0; //jumps

		CalcAlphaBeta::calcJointLapAlphaBeta(product.getStartDate(),
                                             matDate,
                                             u_tT,
                                             alpha_tT,
                                             beta_tT,
                                             this);

        /* Then, from today till start date */

        //used in CalcAlphaBeta
		Complex alpha_0t;
        ComplexArray beta_0t(n+3);
        
		//Terminal conditions
        ComplexArray u_0t(n+3);
        for(int k=0; k<n; k++){
            u_0t[k] = beta_tT[k];
        }
        u_0t[n] = beta_tT[n]; //spot
        u_0t[n+1] = beta_tT[n+1] - z; //real var, no maturity adjustment
        u_0t[n+2] = beta_tT[n+2];
        
		CalcAlphaBeta::calcJointLapAlphaBeta(baseDate,
                                             product.getStartDate(),
                                             u_0t,
                                             alpha_0t,
                                             beta_0t,
                                             this);

        res = alpha_tT + alpha_0t;
        return res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

//expectation
double VolVSCurve::expectation(const FwdStFourierProcessQuadVar& process,
                               const FwdStFourierProductQuadVar& product, 
                               const DateTime& matDate) const{    
    static const string method = "VolVSCurve::cumulant FwdStFourierProcessQuadVar";
    try{
        //result
        double res;

        double startDateDouble = timeMetric->yearFrac(baseDate, product.getStartDate());
        double matDateDouble = timeMetric->yearFrac(baseDate, matDate);
        if( Maths::isZero(matDateDouble) || Maths::isPositive(startDateDouble - matDateDouble) ){
            throw ModelException("matDate can't be zero", method);
        }
        else{
            double resMat = vscaff->yield(matDateDouble);
            double resStart = vscaff->yield(startDateDouble);
            double resCoeff = matDateDouble - startDateDouble;
            res = 1 / (resCoeff) * (resMat - resStart);
            double lambda_S = spotcrashRate;
            double k_S = spotcrashSize;
            double delta_S = spotcrashUncertainty;

            double chi_S = log(1.0 + k_S) + 0.5 * Maths::square(delta_S);
            res += lambda_S * chi_S;
            return res;
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}


// Laplace transform for expected quad var
// Pricing VIX futures with VSCurve (ARNAUD)
Complex VolVSCurve::cumulant(const FwdStFourierProcessExpQuadVar&	process,
                             const FwdStFourierProductExpQuadVar&	product, 
                             const Complex&						z, 
                             const DateTime&					matDate) const {
    static const string method = "VolVSCurve::cumulant FwdStFourierProcessExpQuadVar";
    try{
        //Time difference from today till start date and from start till maturity
        double t = timeMetric->yearFrac(baseDate, product.getStartDate());
        double tau = timeMetric->yearFrac(product.getStartDate(), matDate);

        //variance swap
        double VS_0T = vscaff->yield( timeMetric->yearFrac(baseDate, matDate) );
        double VS_0t = vscaff->yield( timeMetric->yearFrac(baseDate, product.getStartDate()) );
        /* Calculate alpha component from start date till maturity */
        Complex alpha_tT;
		alpha_tT = z * (VS_0T - VS_0t); 
		//Terminal conditions
        ComplexArray u_tT(n+3);
        for(int j=0; j<n; j++){
            u_tT[j] =  z * (1 - exp(-meanReversRate[j]*tau)) / meanReversRate[j];
        }
        u_tT[n] = Complex(0.0, 0.0); //spot
        u_tT[n+1] = Complex(0.0, 0.0); //real var
        if (Maths::isZero(crashMRR)) {
            throw ModelException(method, "zero crashMRR not supported for VSW calculation");
        } else {
            u_tT[n+2] = z * (1 - exp(-crashMRR*tau)) / crashMRR; // jumps
        }
        /* Then, from today till start date */
		Complex alpha_t;
        ComplexArray beta_t(n+3);
        CalcAlphaBeta::calcJointLapAlphaBeta(baseDate,
                                             product.getStartDate(),
                                             u_tT,
                                             alpha_t,
                                             beta_t,
                                             this);
    
        return (alpha_t + alpha_tT);

    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}



//get params
DoubleArray VolVSCurve::getVSCurveParam(const EnumParam param) const{
    DoubleArray resCorr(n), resVol(n), resMrr(n), resTen(n);
    DoubleArray rescrashRate(1), rescrashSize(1), rescrashMRR(1), resspotcrashRate(1), resspotcrashSize(1), resspotcrashUncertainty(1);
    for(int j=0; j<n; j++){
            resCorr[j] = correlation[j];
            resVol[j] = volVol[j];
            resMrr[j] = meanReversRate[j];
        }

    DoubleArray resLamb(n);
    double lambSum = 0.0;
    for(int k=0; k<n-1; k++){
        resLamb[k] = lambda[k];
        lambSum += lambda[k];
    }
    resLamb[n-1] = 1.0-lambSum;
    rescrashMRR[0] = crashMRR;
    rescrashRate[0] = crashRate;
    rescrashSize[0] = crashSize;
    resspotcrashRate[0] = spotcrashRate;
    resspotcrashSize[0] = spotcrashSize;
    resspotcrashUncertainty[0] = spotcrashUncertainty;
    
    switch (param)
    {
    case CORRELATION:
        return(resCorr);

    case VOL_VOL:
        return(resVol);

    case MEAN_REVERS_RATE:
        return(resMrr);
        
    case LAMBDA:
        return(resLamb);

    case CRASH_MRR:
        return(rescrashMRR);

    case CRASH_RATE:
        return(rescrashRate);

    case CRASH_SIZE:
        return(rescrashSize);

    case SPOT_CRASH_RATE:
        return(resspotcrashRate);

    case SPOT_CRASH_SIZE:
        return(resspotcrashSize);

    case SPOT_CRASH_UNCERTAINTY:
        return(resspotcrashUncertainty);

        
    default:
        throw ModelException("VolVSCurve::getVSCurveParam", "unkown param request"); // cannot reach here
        DoubleArray Zerror(1);
        Zerror[0]=0.0;
        return(Zerror);
    }
}

/** Needed for IAdjustable interface. Returns market data name for vol */
string VolVSCurve::getName() const{
    return CVolBaseParamSurface::getName();
}

DRLIB_END_NAMESPACE