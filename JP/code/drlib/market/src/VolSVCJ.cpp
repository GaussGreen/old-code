//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolSVCJ.cpp
//
//   Description : CommonStockVolDiff + CommonStockVolJump
//
//   Date        : 20 Apr 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolSV.hpp"
#include "edginc/VolSVJ.hpp"
#include "edginc/VolSVJJ.hpp"
#include "edginc/VolSVCJ.hpp"
#include "edginc/mathlib.hpp"

#include "edginc/VolProcessed.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/VolRequestTime.hpp"
#include "edginc/FlatVol.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"


DRLIB_BEGIN_NAMESPACE

/****** GAD 19/02/2006 ******/

//Order of the expansion
static const int NInfExpan = 5;

/*** helper functions for SVCJ quadratic variation alpha ***/
class CommonAlphaHelper{
private:
    //integral term for Beta inverse
    static Complex IntegralBetaInv(int            l, //order of the integral (l from 0 to ...)
                                   double         mRR, //mean reversion rate
                                   double         vVol, //vol of vol
                                   double         mu, //jump size mean
                                   double         delta, //jump size variance
                                   double         lambda, //jump intensity
                                   double         rhoJ, //correlation rhoJ
                                   double         tau, //tau=T-t
                                   const Complex& u1, //input parameter corresponding to quad var
                                   const Complex& u2){ //input parameter corresponding to spot variance
        static const string method = "CommonAlphaHelper::IntegralBetaInv";
        try{
            //result
            Complex res;

            //coefficients
            double mRRSQ = mRR * mRR;
            double vVolSQ = vVol * vVol;
            //gamma
            Complex gammaSQ = mRRSQ - 2.0 * vVolSQ * u1;
            Complex gamma = sqrt(gammaSQ);

            //case gamma = 0
            if( Maths::isZero(gamma) ){
                //case vVol = 0
                if( Maths::isZero(vVol) ){
                    Complex cComp = 1.0 - mu * u2;
                    res = pow(cComp, l+1) * tau;
                }
                else{
                    //coefficients in beta
                    double coeff0 = mRR / vVolSQ;
                    Complex coeff1 = u2 - coeff0;
                    Complex coeff2 = 0.5 * vVolSQ * coeff1;
                    double coeff3 = vVolSQ - mu * mRR;

                    //coefficients a, b, c, d
                    //betaInv = (a+b*tau) / (c+d*tau);
                    double a = 1.0;
                    Complex b = - coeff2;
                    Complex c = 1.0 - mu * u2;
                    Complex d = -0.5 * coeff3 * coeff1;

                    //tests for limiting cases
                    bool bZero = Maths::isZero(b);  
                    bool dZero = Maths::isZero(d);
                    bool bdZero = bZero && dZero;

                    //both b and d are zero
                    if(bdZero){
                        //case c = 0
                        if( Maths::isZero(c) ){
                            throw ModelException(method, "c can't be zero");
                        }
                        //c != 0
                        else{
                            Complex resCoeff = a / c; 
                            res = pow(resCoeff, l+1) * tau;          
                        }
                    }
                    //b or d is zero
                    else{
                        //b != 0, d = 0
                        if(dZero){
                            //case c = 0
                            if( Maths::isZero(c) ){
                                throw ModelException(method, "c can't be zero");
                            }
                            //c != 0
                            else{
                                Complex resCoeff = a / c;
                                Complex bTil = b / a;
                                Complex resPow = 1.0 + bTil * tau;
                                res = pow(resCoeff, l+1) * (pow(resPow,l+2) - 1.0) / (double(l+2) * bTil);
                            }
                        }
                        //b, d != 0
                        else{
                            //coefficients aTil, cTil
                            Complex aTil = a / b;
                            Complex cTil = c / d;

                            //tests for limiting cases
                            bool cZero = Maths::isZero(c);
                            
                            //c = 0
                            if(cZero){
                                throw ModelException(method, "c can't be zero");
                            }
                            //c != 0
                            else{
                                //coefficients
                                Complex resCoeff = b / d;
                                Complex acTil = aTil - cTil;
                        
                                //l = 0 (integrals Ii from 0 to 1)
                                Complex resTemp = tau; //I0

                                Complex resLog = 1.0 + tau / cTil;
                                double rLog = Complex::absSquare(resLog);
                                double argLog = Complex::argument(resLog);
                                Complex logResult = Complex( 0.5 * log(rLog), argLog);
        
                                resTemp += (double)(l+1) * acTil * logResult; //I1
    
                                //for l>0 (integrals Ii from 2 to...)
                                double binom, binomP;
                                Complex cTilTau;
                                for(int compt=2; compt<l+2; compt++){
                                    binom = combineBinom(l+1, compt);
                                    binomP = binom / ( (double)compt - 1.0 );
                                    cTilTau = cTil + tau;
                                    resTemp += binomP * pow(acTil, compt) * ( 1.0 / pow(cTil, compt-1) - 1.0 / pow(cTilTau, compt-1) ); //Icompt
                                }
                                res = pow(resCoeff, l+1) * resTemp;
                            }
                        }
                    }
                }    
            }
            //case gamma != 0
            else{
                //case vVol = 0
                if( Maths::isZero(vVol) ){
                    //case mu = 0
                    if(Maths::isZero(mu)){
                        res = tau;
                    }
                    //mu != 0
                    else{
                        Complex cTil = - mu * u2;
                        Complex resCoeff = 1.0 / pow(cTil, l+1);
                        res = resCoeff * Integrali(l+1,
                                                   cTil,
                                                   tau,
                                                   mRR);
                    }
                }
                //vVol != 0
                else{
                    //coefficients in beta
                    Complex coeff0 = mRR - gamma;
                    Complex coeff1 = coeff0 / vVolSQ;
                    Complex coeff2 = u2 - coeff1;
                    Complex coeff3 = 0.5 * (vVolSQ / gamma) * coeff2;
                    Complex coeff4 = ( vVolSQ - mu * (mRR - gamma) ) / vVolSQ;
                    Complex coeff5 = 0.5 * ( vVolSQ - mu * (mRR + gamma) ) / gamma;

                    //coefficients a, b, c, d
                    //betaInv = (a + b * exp(- gamma*tau) ) / (c + d * exp(- gamma*tau) )
                    Complex a = 1.0 - coeff3;
                    Complex b = coeff3;
                    Complex c = coeff4 * a;
                    Complex d = coeff5 * coeff2;

                    //tests for limiting cases
                    bool bZero = Maths::isZero(b);  
                    bool dZero = Maths::isZero(d);
                    bool bdZero = bZero && dZero;
                
                    //both b and d are zero
                    if(bdZero){
                        //case c = 0
                        if( Maths::isZero(c) ){
                            throw ModelException(method, "c can't be zero");
                        }
                        //c != 0
                        else{
                            Complex resCoeff = a / c;
                            res = pow(resCoeff, l+1) * tau;
                        }                    
                    }
                    //b or d is zero
                    else{
                        //b != 0, d = 0
                        if(dZero){
                            //case c = 0
                            if( Maths::isZero(c) ){
                                throw ModelException(method, "c can't be zero");
                            }
                            //c != 0
                            else{
                                //coefficients
                                Complex resCoeff = a / c;
                        
                                //l = 0
                                Complex resTemp = tau; //I0
                                Complex bTil = b / a;
                                resTemp += (double)(l+1) * bTil * ( 1.0 - exp(- gamma * tau) ) / gamma; //I1

                                //for l>0 (integrals Ii from 2 to...)
                                double binom, binomP;
                                Complex resExp;
                                for(int compt=2; compt<l+2; compt++){
                                    binom = combineBinom(l+1, compt);
                                    binomP = binom / (double)compt;
                                    resExp = ( 1.0 - exp(-(double)compt * gamma * tau) ) / gamma;
                                    resTemp += binomP * pow(bTil, compt) * resExp; //Icompt
                                }
                                res = pow(resCoeff, l+1) * resTemp;
                            }
                        }
                        // b, d != 0
                        else{
                            //tests for limiting cases
                            bool cZero = Maths::isZero(c);

                            //c = 0
                            if(cZero){
                                //coefficients
                                Complex resCoeff = b / d;
                            
                                //l = 0
                                Complex resTemp = tau; // I0
                                Complex aTil = a / b;
                                resTemp += (double)(l+1) * aTil * ( exp(gamma * tau) - 1.0 ) / gamma; //I1
                            
                                //for l>0 (integrals Ii from 2 to...)
                                double binom, binomP;
                                Complex resExp;
                                for(int compt=2; compt<l+2; compt++){
                                    binom = combineBinom(l+1, compt);
                                    binomP = binom / (double)compt;
                                    resExp = (exp((double)compt * gamma * tau) - 1.0) / gamma;
                                    resTemp += binomP * pow(aTil, compt) * resExp; //Icompt
                                }
                                res = pow(resCoeff, l+1) * resTemp;
                            }
                            //c != 0
                            else{
                                //coefficients
                                Complex resCoeff = b / d;

                                //coefficients aTil, cTil
                                Complex aTil = a / b;
                                Complex cTil = c / d;
                                Complex acTil = aTil - cTil;

                                Complex resTemp = tau;
                                double resFactTemp;
                                for(int compt = 1; compt < l+2; compt++ ){
                                    resFactTemp = combineBinom(l+1, compt);
                                    resTemp += resFactTemp * pow(acTil, compt) * Integrali(compt,
                                                                                           cTil,
                                                                                           tau,
                                                                                           gamma);
                                }   
                                res = pow(resCoeff, l+1) * resTemp;   
                            }
                        }
                    } 
                }
            }
            return res;
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

    //partial exponential sum (in l in paper)
    static Complex ExpSum(int               l, //start integer
                          const Complex&    x){ //point x
        static const string method = "CommonAlphaHelper::ExpSum";
        try{
            //result
            Complex res = 0.0;

            //coefficients
            Complex resCoeff = 2.0 * x;
            switch(l){
            case 0:
                res = exp(x);
                break;

            case 1:
                res = pow(resCoeff,1) * exp(x);
                break;
            case 2:
                res = pow(resCoeff,1) * (2.0*x + 1.0) * exp(x);
                break;

            case 3:
                res = pow(resCoeff,2) * (2.0*x + 3.0) * exp(x);
                break;
            case 4:
                res = pow(resCoeff,2) * (4.0*pow(x,2) + 12.0*x + 3.0) * exp(x);
                break;

            case 5:
                res = pow(resCoeff,3) * (4.0*pow(x,2) + 20.0*x + 15.0) * exp(x);
                break;
            case 6:
                res = pow(resCoeff,3) * (8.0*pow(x,3) + 60.0*pow(x,2) + 90.0*x + 15.0) * exp(x);
                break;
            
            case 7:
                res = pow(resCoeff,4) * (8.0*pow(x,3) + 84.0*pow(x,2) + 210.0*x + 105.0) * exp(x);
                break;            
            case 8:
                res = pow(resCoeff,4) * (16.0*pow(x,4) + 224.0*pow(x,3) + 840.0*pow(x,2) + 840.0*x + 105.0) * exp(x);
                break;

            case 9:
                res = pow(resCoeff,5) * (16.0*pow(x,4) + 288.0*pow(x,3) + 1512.0*pow(x,2) + 2520.0*x + 945.0) * exp(x);
                break;
            case 10:
                res = pow(resCoeff,5) * (32.0*pow(x,5) + 720.0*pow(x,4) + 5040.0*pow(x,3) + 12600.0*pow(x,2) + 9450.0*x + 945.0) * exp(x);
                break;

            case 11:
                res = pow(resCoeff,6) * (32.0*pow(x,5) + 880.0*pow(x,4) + 7920.0*pow(x,3) + 27720.0*pow(x,2) + 34650.0*x + 10395.0) * exp(x);
                break;
            case 12:
                res = pow(resCoeff,6) * (64.0*pow(x,6) + 2112.0*pow(x,5) + 23760.0*pow(x,4) + 110880.0*pow(x,3) + 207900.0*pow(x,2) + 124740.0*x + 10395.0) * exp(x);
                break;

            case 13:
                res = pow(resCoeff,7) * (64.0*pow(x,6) + 2496.0*pow(x,5) + 34320.0*pow(x,4) + 205920.0*pow(x,3) + 540540.0*pow(x,2) + 540540.0*x + 135135.0) * exp(x);
                break;
            case 14:
                res = pow(resCoeff,7) * (128.0*pow(x,7) + 5824.0*pow(x,6) + 96096.0*pow(x,5) + 720720.0*pow(x,4) + 2522520.0*pow(x,3) + 3783780.0*pow(x,2) + 1891890.0*x + 135135.0) * exp(x);
                break;

            case 15:
                res = pow(resCoeff,8) * (128.0*pow(x,7) + 6720.0*pow(x,6) + 131040.0*pow(x,5) + 1201200.0*pow(x,4) + 5405400.0*pow(x,3) + 11351340.0*pow(x,2) + 9459450.0*x + 2027025.0) * exp(x);
                break;
            case 16:
                res = pow(resCoeff,8) * (256.0*pow(x,8) + 15360.0*pow(x,7) + 349440.0*pow(x,6) + 3843840.0*pow(x,5) + 21621600.0*pow(x,4) + 60540480.0*pow(x,3) + 75675600.0*pow(x,2) + 32432400.0*x + 2027025.0) * exp(x);
                break;

            case 17: pow(resCoeff,9) * (256.0*pow(x,8) + 17408.0*pow(x,7) + 456960.0*pow(x,6) + 5940480.0*pow(x,5) + 40840800.0*pow(x,4) + 147026880.0*pow(x,3) + 257297040.0*pow(x,2) + 183783600.0*x + 34459425.0) * exp(x);
                break;
            case 18: pow(resCoeff,9) * (512.0*pow(x,9) + 39168.0*pow(x,8) + 1175040.0*pow(x,7) + 17821440.0*pow(x,6) + 147026880.0*pow(x,5) + 661620960.0*pow(x,4) + 1543782240.0*pow(x,3) + 1654052400.0*pow(x,2) + 620269650.0*x + 34459425.0) * exp(x);
                break;

            case 19: pow(resCoeff,10) * (512.0*pow(x,9) + 43776.0*pow(x,8) + 1488384.0*pow(x,7) + 26046720.0*pow(x,6) + 253955520.0*pow(x,5) + 1396755360.0*pow(x,4) + 4190266080.0*pow(x,3) + 6285399120.0*pow(x,2) + 3928374450.0*x + 654729075.0) * exp(x);
                break;
            case 20: pow(resCoeff,10) * (1024.0*pow(x,10) + 97280.0*pow(x,9) + 3720960.0*pow(x,8) + 74419200.0*pow(x,7) + 846518400.0*pow(x,6) + 5587021440.0*pow(x,5) + 20951330400.0*pow(x,4) + 41902660800.0*pow(x,3) + 39283744500.0*pow(x,2) + 13094581500.0*x + 654729075.0) * exp(x);
                break;

            case 21: pow(resCoeff,11) * (1024.0*pow(x,10) + 107520.0*pow(x,9) + 4596480.0*pow(x,8) + 104186880.0*pow(x,7) + 1367452800.0*pow(x,6) + 10666131840.0*pow(x,5) + 48886437600.0*pow(x,4) + 125707982400.0*pow(x,3) + 164991726900.0*pow(x,2) + 91662070500.0*x + 13749310575.0) * exp(x);
                break;
            case 22: pow(resCoeff,11) * (2048.0*pow(x,11) + 236544.0*pow(x,10) + 11235840.0*pow(x,9) + 286513920.0*pow(x,8) + 4297708800.0*pow(x,7) + 39109150080.0*pow(x,6) + 215100325440.0*pow(x,5) + 691393903200.0*pow(x,4) + 1209939330600.0*pow(x,3) + 1008282775500.0*pow(x,2) + 302484832650.0*x + 13749310575.0) * exp(x);
                break;

            case 23: pow(resCoeff,12) * (2048.0*pow(x,11) + 259072.0*pow(x,10) + 13601280.0*pow(x,9) + 387636480.0*pow(x,8) + 6589820160.0*pow(x,7) + 69193111680.0*pow(x,6) + 449755225920.0*pow(x,5) + 1766895530400.0*pow(x,4) + 3975514943400.0*pow(x,3) + 4638100767300.0*pow(x,2) + 2319050383650.0*x + 316234143225.0) * exp(x);
                break;
            case 24: pow(resCoeff,12) * (4096.0*pow(x,12) + 565248.0*pow(x,11) + 32643072.0*pow(x,10) + 1033697280.0*pow(x,9) + 19769460480.0*pow(x,8) + 237233525760.0*pow(x,7) + 1799020903680.0*pow(x,6) + 8481098545920.0*pow(x,5) + 23853089660400.0*pow(x,4) + 37104806138400.0*pow(x,3) + 27828604603800.0*pow(x,2) + 7589619437400.0*x + 316234143225.0) * exp(x);
                break;

            default:
                throw ModelException(method, "coefficient not implemented yet");
            }
            return res;
        }    
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

    //factorial function
    static double fact(int n){
        //result
        int res = 1;
        if(n == 0){
            //nothing
        }
        else{
            for(int i=1; i<n+1; i++){
                res *= i; 
            }
        }
        return (double)res;
    }

    //Anp function : n! / p!
    static double combineBinom(int n,
                               int p){
        static const string method = "CommonAlphaHelper::combineBinom";
        if(p > n){
            throw ModelException(method, "n must be higher than p");
        }
        else{
            try{
                //result
                double res = 1.0;
                
                if(p == 0){
                    //nothing
                }
                else{
                    for(int i=1; i<(p+1); i++){
                        res *= (double)(n+1-i) / (double)i;
                    }
                }
                return res;
            }
            catch(exception& e){
                throw ModelException(e, method);
            }
        }    
    }

    //integral Ii in the general case (a, b, c, d != 0)
    static Complex Integrali(int                i, //order
                             const Complex&     cTil, //cTilde
                             double             tau, //tau
                             const Complex&     gamma){ //gamma != 0
        static const string method = "CommonAlphaHelper::Integrali";
        try{
            //result
            Complex res = 0.0;

            //I0
            if(i == 0){
                res = tau;
            }
            //Ii for i from 1 to...
            else{
                //coefficients
                Complex resCoeff = pow(cTil, i);
                
                //1 tau term
                Complex resTemp = tau;             
                
                //2 logComplex using Jackel method
                //complex c
                Complex cComp = - (1.0 / cTil);
                double tc = Complex::argument(cComp);
                //phase change
                double bd = - gamma.imag();

                //numerator
                Complex cNum = cComp * exp(- gamma * tau) - 1.0;
                double rNum = Complex::absSquare(cNum);
                double kiNum = Complex::argument(cNum);
                double nNum = floor( (tc + Maths::PI + bd * tau) / (2.0 * Maths::PI) );

                //denominator
                Complex cDenom = cComp - 1.0;
                double rDenom = Complex::absSquare(cDenom);
                double kiDenom = Complex::argument(cDenom);
                double nDenom = floor( (tc + Maths::PI) / (2.0 * Maths::PI) );

                double rCoeff = 0.5 * log(rNum / rDenom);
                //double argCoeff = kiNum - kiDenom + 2.0 * Maths::PI * (nNum - nDenom);
                double argCoeff = kiNum - kiDenom;
                
                Complex logComp = Complex( rCoeff, argCoeff);       
                resTemp += logComp / gamma;

                //3 sum of rational fractions (i>=2)
                Complex resCoeff1, resCoeff2, resPt1, resPt2;
                for(int compt=1; compt<i; compt++ ){
                    resCoeff1 = cTil / ( cTil + 1.0 );
                    resCoeff2 = cTil / ( cTil + exp(-gamma*tau) );
                    resPt1 = pow(resCoeff1, compt);
                    resPt2 = pow(resCoeff2, compt);
                    resTemp += ( resPt1 - resPt2 ) / ( (double)compt * gamma );
                }
                
                //final
                res = resTemp / resCoeff;    
            }
            return res;
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

public:
    //total sum to order N (in k in paper)
    static Complex AlphaSumN(int                N, //order of the integral 
                             double             mRR, //mean reversion rate
                             double             vVol, //vol of vol
                             double             mu, //jump size
                             double             delta, //variance
                             double             lambda, //jump intensity
                             double             rhoJ, //correlation rhoJ
                             double             gammaBar, //gammaBAR in the computation
                             double             tau, //tau=T-t
                             const Complex&     u1, //input parameter corresponding to quad var
                             const Complex&     u2){ //input parameter corresponding to spot variance
        static const string method = "CommonAlphaHelper::AlphaSumN";
        try{
            //result
            Complex res = 0.0;
            
            //case gammaBar = 0
            if( Maths::isZero(gammaBar) ){
                //coefficient
                double rhoJSQ = rhoJ * rhoJ;
                double muSQ = mu * mu;
                double deltaSQ = delta * delta;
                
                //U
                Complex UUp = u1 * rhoJSQ * muSQ;
                Complex UDown = 1.0 - 2.0 * u1 * deltaSQ;
                Complex U = UUp / UDown;

                Complex resTemp;
                //odd or even number of terms
                //for(int l=0; l<(2*N+1); l++){
                for(int l=0; l<(2*N+2); l++){
                    resTemp = combineBinom(2*l, l) * pow(U, l) * IntegralBetaInv(2*l,
                                                                                 mRR,
                                                                                 vVol,
                                                                                 mu,
                                                                                 delta,
                                                                                 lambda,
                                                                                 rhoJ,
                                                                                 tau,
                                                                                 u1,
                                                                                 u2);
                    res += resTemp;
                }
            }
            //gammaBar != 0
            else{
                //coefficient
                double gammaBarSQ = gammaBar * gammaBar;
                double deltaSQ = delta * delta; 
    
                //U
                Complex UUp = u1 * gammaBarSQ;
                Complex UDown = 1.0 - 2.0 * u1 * deltaSQ;
                Complex U = UUp / UDown;

                //GAM
                double GAM = rhoJ * mu / gammaBar;
            
                Complex resTemp;
                //Odd / Even number of terms
                //for(int l=0; l<(2*N+1); l++){
                for(int l=0; l<(2*N+2); l++){
                    resTemp = pow(GAM,l) * ExpSum(l,U) * IntegralBetaInv(l,
                                                                         mRR,
                                                                         vVol,
                                                                         mu,
                                                                         delta,
                                                                         lambda,
                                                                         rhoJ,
                                                                         tau,
                                                                         u1,
                                                                         u2);                   
                    res += resTemp;
                }
            }
            return res;
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

    //total sum to order N using Pade equivalent of the 2N Taylor expansion (in k in paper)
    static Complex AlphaSumNPade(int                N, //order of the integral 
                                 double             mRR, //mean reversion rate
                                 double             vVol, //vol of vol
                                 double             mu, //jump size
                                 double             delta, //variance
                                 double             lambda, //jump intensity
                                 double             rhoJ, //correlation rhoJ
                                 double             gammaBar, //gammaBAR in the computation
                                 double             tau, //tau=T-t
                                 const Complex&     u1, //input parameter corresponding to quad var
                                 const Complex&     u2){ //input parameter corresponding to spot variance
        static const string method = "CommonAlphaHelper::AlphaSumNPade";
        try{
            //result
            Complex res = 0.0;

            //case gammaBar = 0
            if( Maths::isZero(gammaBar) ){
                //coefficients
                double rhoJSQ = rhoJ * rhoJ;
                double muSQ = mu * mu;
                double deltaSQ = delta * delta;
                
                //U
                Complex UUp = u1 * rhoJSQ * muSQ;
                Complex UDown = 1.0 - 2.0 * u1 * deltaSQ;
                Complex U = UUp / UDown;

                //u1 = 0, no inversion can be performed in Pade
                if( Maths::isZero(u1) ){
                    res = IntegralBetaInv(0,
                                          mRR,
                                          vVol,
                                          mu,
                                          delta,
                                          lambda,
                                          rhoJ,
                                          tau,
                                          u1,
                                          u2);
                }
                //u1 != 0, Pade can be used
                else{
                    //array (matrix) a (a1 -> a2N-1)
                    vector<d_complex> a(N*N);
                    for(int i=0; i<N; i++){
                        for(int j = 0; j<N; j++){
                            int l = N + (i+1) - (j+1);
                            Complex coeff = combineBinom(2*l, l) * pow(U,l) * IntegralBetaInv(2*l,
                                                                                              mRR,
                                                                                              vVol,
                                                                                              mu,
                                                                                              delta,
                                                                                              lambda,
                                                                                              rhoJ,
                                                                                              tau,
                                                                                              u1,
                                                                                              u2);                    
                            //definition of term ( (i*N)+j )
                            a [(i*N) + j].re = coeff.real();
                            a [(i*N) + j].im = coeff.imag();
                        }
                    }

                    //array b (aN+1 -> a2N)
                    vector<d_complex> b(N);
                    for(int m = 0; m<N; m++){
                        int n = N + (m+1);
                        Complex coeff = combineBinom(2*n, n) * pow(U,n) * IntegralBetaInv(2*n,
                                                                                          mRR,
                                                                                          vVol,
                                                                                          mu,
                                                                                          delta,
                                                                                          lambda,
                                                                                          rhoJ,
                                                                                          tau,
                                                                                          u1,
                                                                                          u2);
                        b[m].re = coeff.real();
                        b[m].im = coeff.imag();
                    }

                    //denominator B
                    vector<d_complex> B_imsl(N);
                    imsl_z_lin_sol_gen(N, &a[0], &b[0], IMSL_RETURN_USER, &B_imsl[0], 0);
                    //transformation into a Complex usable in qLib
                    ComplexArray B(N+1);
                    B[0] = 1.0;
                    for(int p=0; p<N; p++){
                        B[p+1] = Complex(B_imsl[p].re, B_imsl[p].im);
                    }

                    //numerator A
                    ComplexArray A(N+1);
                    for(int q = 0; q < N+1; q++){
                        A[q] = 0.0;
                        for(int r = 0; r<q+1; r++){
                            int l = q - r;
                            A[q] += B[r] * combineBinom(2*l, l) * pow(U,l) * IntegralBetaInv(l,
                                                                                             mRR,
                                                                                             vVol,
                                                                                             mu,
                                                                                             delta,
                                                                                             lambda,
                                                                                             rhoJ,
                                                                                             tau,
                                                                                             u1,
                                                                                             u2);
                        }
                    }

                    //construction of Pade approximant at order N / N
                    Complex resTempNum = 0.0;
                    Complex resTempDenom = 0.0;
                    for(int compt = 0; compt<(N+1); compt++){
                        resTempNum += pow(U,compt) * A[compt];                 
                        resTempDenom += pow(U,compt) * B[compt];
                    }
                    res = resTempNum / resTempDenom;            
                }
            }
            //gammaBar != 0
            else{
                //coefficients
                double gammaBarSQ = gammaBar * gammaBar;
                double deltaSQ = delta * delta; 

                //U
                Complex UUp = u1 * gammaBarSQ;
                Complex UDown = 1.0 - 2.0 * u1 * deltaSQ;
                Complex U = UUp / UDown;
                
                //GAM
                double GAM = rhoJ * mu / gammaBar;

                //u1 = 0, no inversion can be performed in Pade
                if( Maths::isZero(u1) ){
                    res = IntegralBetaInv(0,
                                          mRR,
                                          vVol,
                                          mu,
                                          delta,
                                          lambda,
                                          rhoJ,
                                          tau,
                                          u1,
                                          u2);
                }
                //u1 != 0, Pade can be used
                else{
                    //array (matrix) a (a1 -> a2N-1)
                    vector<d_complex> a(N*N);
                    for(int i=0; i<N; i++){
                        for(int j = 0; j<N; j++){
                            int l = N + (i+1) - (j+1);
                            Complex coeff= ExpSum(l,U) * IntegralBetaInv(l,
                                                                         mRR,
                                                                         vVol,
                                                                         mu,
                                                                         delta,
                                                                         lambda,
                                                                         rhoJ,
                                                                         tau,
                                                                         u1,
                                                                         u2);                    
                            //definition of term ( (i*N)+j )
                            a [(i*N) + j].re = coeff.real();
                            a [(i*N) + j].im = coeff.imag();
                        }
                    }

                    //array b (aN+1 -> a2N)
                    vector<d_complex> b(N);
                    for(int m = 0; m<N; m++){
                        int n = N + (m+1);
                        Complex coeff = - ExpSum(n,U) * IntegralBetaInv(n,
                                                                        mRR,
                                                                        vVol,
                                                                        mu,
                                                                        delta,
                                                                        lambda,
                                                                        rhoJ,
                                                                        tau,
                                                                        u1,
                                                                        u2);
                        b[m].re = coeff.real();
                        b[m].im = coeff.imag();
                    }

                    //denominator B
                    vector<d_complex> B_imsl(N);
                    imsl_z_lin_sol_gen(N, &a[0], &b[0], IMSL_RETURN_USER, &B_imsl[0], 0);
                    //transformation into a Complex usable in qLib
                    ComplexArray B(N+1);
                    B[0] = 1.0;
                    for(int p=0; p<N; p++){
                        B[p+1] = Complex(B_imsl[p].re, B_imsl[p].im);
                    }

                    //numerator A
                    ComplexArray A(N+1);
                    for(int q = 0; q < N+1; q++){
                        A[q] = 0.0;
                        for(int r = 0; r<q+1; r++){
                            int l = q - r;
                            A[q] += B[r] * ExpSum(l,U) * IntegralBetaInv(l,
                                                                         mRR,
                                                                         vVol,
                                                                         mu,
                                                                         delta,
                                                                         lambda,
                                                                         rhoJ,
                                                                         tau,
                                                                         u1,
                                                                         u2);
                        }
                    }

                    //construction of Pade approximant at order N / N
                    Complex resTempNum = 0.0;
                    Complex resTempDenom = 0.0;
                    for(int compt = 0; compt<(N+1); compt++){
                        resTempNum += pow(GAM,compt) * A[compt];                 
                        resTempDenom += pow(GAM,compt) * B[compt];
                    }
                    res = resTempNum / resTempDenom;
                }
            }
            return res;
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }
};

/****** end of GAD 19/02/2006 ******/

/************************ helper functions ***************/
class CalcAlpha{
public:

    /** simplified for VolSVCJ. Refer to comments for same method in Heston.cpp */
    static void calcAlphaBeta(double         tau,    // tau = T - t (in years)
                              const Complex& u1,
                              const Complex& u2,
                              const Complex& u3,
                              Complex&       alpha,
                              Complex&       beta,
                              const VolSVCJ* svcj){

        double dCorr = svcj->correlation;

        double mr = svcj->meanReversRate;
        double vVol = svcj->volVol;

        bool zeroVolVolZeromeanReversRate = Maths::isZero(vVol) && Maths::isZero(mr);
        bool zeroVolVolNonZeromeanReversRate  = Maths::isZero(vVol) && !Maths::isZero(mr);
        double mrSqMVol = mr * Maths::square(svcj->meanVol);
        double sqVolVol = Maths::square(vVol);

        Complex a = -2.0 * u3 + u1 * (1.0 - u1);
        if (zeroVolVolZeromeanReversRate){
            /* Calculate beta */
            beta = u2 - a * (0.5 * tau);
            /* Calculate alpha */
            alpha = 0.0;
        }
        else if (zeroVolVolNonZeromeanReversRate){
            /* Calculate beta */
            double amortFact = exp(-mr * tau);
            beta = u2 * amortFact;
            amortFact = (1.0 - amortFact) / mr;
            beta -=  a * (0.5 * amortFact);

            /* Calculate alpha */
            alpha = a * (0.5 * (amortFact - tau) / mr);
            alpha = (alpha + u2 * amortFact) * mrSqMVol;
        }
        else{
            Complex b = dCorr * vVol * u1 - mr;
            Complex gamma = sqrt(Maths::square(b) + a * sqVolVol);
            if (Maths::isZero(gamma)){
                Complex beta_minus = - b / sqVolVol;
                Complex u2_less_beta_minus = u2 - beta_minus;

                /* Calculate numerator of beta - beta_minus */
                Complex num = u2_less_beta_minus;

                /* Calculate denominator of beta - beta_minus */
                Complex den = 1.0 - u2_less_beta_minus * (tau * (0.5 * sqVolVol));

                /* Calculate beta */
                beta = beta_minus + num / den;

                /* Calculate stochastic vol contribution to alpha */
                alpha = (beta_minus * tau - log(den) * (2.0 / sqVolVol)) * mrSqMVol;
            }
            else{
                Complex beta_minus = - (b + gamma) / sqVolVol;
                Complex u2_less_beta_minus = u2 - beta_minus;

                /* Calculate numerator of beta - beta_minus */
                Complex amortFact = exp(- gamma * tau);
                Complex num = u2_less_beta_minus * amortFact;

                /* Calculate denominator of beta - beta_minus */
                amortFact = (1.0 - amortFact) / gamma;
                Complex den = 1.0 - u2_less_beta_minus * amortFact * (0.5 * sqVolVol);

                /* Calculate beta */
                beta = beta_minus + num / den;

                /* Calculate stochastic vol contribution to alpha */
                alpha = (beta_minus * tau - log(den) * (2.0 / sqVolVol)) * mrSqMVol;
            }
        }
    }

    /** simplified for VolSVCJ. Refer to comments for same method in Heston.cpp */
    static Complex calcAlpha(double  tau,    // tau = T - t (in years)
                           const Complex& u1,
                           const Complex& u2,
                           const Complex& u3,
                           const VolSVCJ* svcj){

        double dCorr = svcj->correlation;
        double jCorr = svcj->stockVolCrashSizeCorrelation;
        double jRate = svcj->commonCrashRate;
        double jVolMean = svcj->commonVolCrashSizeMean;
        double jStkMean = svcj->commonStockCrashSizeMean;
        double jStkWidth = svcj->commonStockCrashSizeUncertainty;
        double mr = svcj->meanReversRate;
        double vVol = svcj->volVol;

        /* If zero rate, alpha = 0.0 */
        if (!Maths::isPositive(jRate)){
            return 0.0;
        }
        /* Compute growth adjustement. NB: chi is close to zero */
        double chi = 1.0 - 1.0 / (1.0 - jCorr * jVolMean);
        double growthAdj = jRate * (jStkMean * (1.0 - chi) - chi) * tau;
        /* Compute stock crash contribution */
        double stockCrashGamma = log(1.0 + jStkMean);
        Complex a = -2.0 * u3 + u1 * (1.0 - u1);
        Complex alpha = exp(u1 * stockCrashGamma - 0.5 * Maths::square(jStkWidth) * a);
        /* Compute vol crash contribution */
        Complex d;
        if (Maths::isZero(jVolMean)){
            d = jRate * tau;
        }
        else{
            Complex c = 1.0 - jCorr * jVolMean * u1;
            if (Maths::isZero(vVol) && Maths::isZero(mr)){
                if (Maths::isZero(a)){
                    d = jRate * tau / (c - jVolMean * u2);
                }
                else{
                    Complex temp = 0.5 * a * jVolMean;
                    d = jRate / temp * log(1.0 + temp * tau / (c - jVolMean * u2));
                }
            }
            else if (Maths::isZero(vVol) && !Maths::isZero(mr)){
                double amortFact = (1.0 - exp(-mr * tau)) / mr;
                d = log(1.0 + jVolMean * (mr * u2 + 0.5 * a) / (c - jVolMean * u2) * amortFact);
                d = jRate / (mr * c + 0.5 * a * jVolMean) * (mr * tau + d);
            }
            else{
                Complex b = dCorr * vVol * u1 - mr;
                double sqVolVol = Maths::square(vVol);
                Complex gamma = sqrt(Maths::square(b) + a * sqVolVol);
                if (Maths::isZero(gamma)){
                    Complex beta_minus = - b / sqVolVol;
                    Complex e = c - jVolMean * beta_minus;
                    Complex f = gamma * jVolMean - 0.5 * sqVolVol * e;
                    d = jRate / e * (tau + jVolMean / f * log(1.0 + f * (u2 - beta_minus) 
                                                                      / (c - jVolMean * u2)
                                                                      * tau));
                }
                else{
                    Complex beta_minus = - (b + gamma) / sqVolVol;
                    Complex amortFact = (1.0 - exp(-gamma * tau)) / gamma;
                    Complex e = c - jVolMean * beta_minus;
                    Complex f = gamma * jVolMean - 0.5 * sqVolVol * e;
                    d = jRate / e * (tau + jVolMean / f * log(1.0 + f * (u2 - beta_minus) 
                                                                     / (c - jVolMean * u2)
                                                                      * amortFact));
                }
            }
        }
        return (alpha * d - jRate * tau - growthAdj * u1);
    }

	/**	Calculates the Laplace Transform of the stationnary distribution of the variance
		and plug the log of it into the cumulant function.	
		This version stick to what MC is doing. It indeed assumes a Gamma distribution for the 
		initial Var and the mean of the distribution is the square of the SVCJ parameter initialVol.
		The jumpPart is also not included since it is not taken into account in the MC.
		Actually the distribution given to the intial Var in the MC is a mere Gamma and not the real 
		long term stationary distribution.
		The exact long term distribution is given in the function calcAlphaStationaryVolTerm below. */
	static Complex calcAlphaStationary (const Complex& u2, const VolSVCJ* svcj) { 
	
		double mr = svcj->meanReversRate;
        double vVol = svcj->volVol;
		double sqIVol = Maths::square(svcj->initialVol);

		Complex hestonPart = 0.0;

		/* In case no vol of vol the term beta*initialVar of standard Heston is needed */ 
		if (Maths::isZero(vVol)) {
			hestonPart = u2 * sqIVol;
		}
		/* If Mean Reversion is zero then we should not have a stationary distribution term */
		else if (!Maths::isZero(mr)){
			double root = Maths::square(vVol) / (2.0 * mr);
			hestonPart = log(1.0 - u2 * root);
			hestonPart *= - sqIVol / root;
		}
		return hestonPart;
	}

	/****** GAD 19/02/2006 ******/

    //computes betaH (see AJD paper by Regis)
    static Complex calcAlphaBetaCommon_Beta(double         tau, // tau = T - t (in years)
                                            const Complex& u1, //quad var
                                            const Complex& u2, //spot variance
                                            const VolSVCJ* svcj){
        static const string method = "CalcAlpha::calcAlphaBetaCommon_Beta";
        try{
            //result
            Complex res = 0.0;

            double mRR = svcj->meanReversRate;
            double vVol = svcj->volVol;
            double mRRSQ = mRR * mRR;
            double vVolSQ = vVol * vVol;

            //gamma
            Complex gammaSQ = mRRSQ - 2.0 * vVolSQ * u1;
            Complex gamma = sqrt(gammaSQ);
        
            //case gamma = 0
            if( Maths::isZero(gamma) ){
                //case vVol = 0
                if( Maths::isZero(vVol) ){
                    throw ModelException(method, "vVol can't be zero");
                }
                //vVol != 0
                else{
                    //coefficients
                    double coeff1 = mRR / vVolSQ;
                    Complex coeff2 = u2 - coeff1;
                    Complex coeff3 = 0.5 * vVolSQ * coeff2 * tau;

                    res = coeff1 + coeff2 / (1.0 - coeff3);
                }
            }
            //gamma != 0
            else{
                //case vVol = 0
                if( Maths::isZero(vVol) ){
                    res = u2 * exp(- mRR * tau);
                }
                //vVol != 0
                else{
                    //coefficients
                    Complex coeff1 = (mRR - gamma) / vVolSQ;
                    Complex coeff2 = u2 - coeff1;
                    Complex coeff3 = - gamma * tau;
                    Complex coeff4 = 0.5 * vVolSQ / gamma;
                    Complex coeff5 = exp(coeff3);
                    Complex coeff6 = 1.0 - coeff5;
        
                    res = coeff1 + (coeff2 * coeff5) / (1.0 - coeff4 * coeff2 * coeff6);
                }
            }
            return res;
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

    //computes alphaH (see Regis paper)
    static Complex calcAlphaBetaCommon_Alpha(double         tau, // tau = T - t (in years)
                                             const Complex& u1, //quad var
                                             const Complex& u2, //spot variance
                                             const VolSVCJ* svcj){
        static const string method = "CalcAlpha::calcAlphaBetaCommon_Alpha";
        
        //Jackel method
        try{
            //result
            Complex res = 0.0;

            double mRR = svcj->meanReversRate;
            double vVol = svcj->volVol;
            double mean = svcj->meanVol;
            double mRRSQ = mRR * mRR;
            double vVolSQ = vVol * vVol;
            double meanSQ = mean * mean;

            //gamma
            Complex gammaSQ = mRRSQ - 2.0 * vVolSQ * u1;
            Complex gamma = sqrt(gammaSQ);
        
            //case gamma = 0
            if( Maths::isZero(gamma) ){
                //case vVol = 0
                if( Maths::isZero(vVol) ){
                    throw ModelException(method, "vVol can't be zero");
                }
                //vVol != 0
                else{
                    //coefficients
                    double coeff0 = (vVolSQ / 2.0) * tau;
                    Complex coeff1 = u2 - (mRR / vVolSQ);
                    Complex coeffLog = 1.0 - coeff0 * coeff1;
                
                    double resPt1 = (mRR / vVolSQ) * tau;
                    Complex resPt2 = (2.0 / vVolSQ) * log(coeffLog);
                    double coeffRes = mRR * meanSQ;
                    
                    res = coeffRes * (resPt1 + resPt2);
                }
            }
            //gamma != 0
            else{
                //case vVol = 0
                if( Maths::isZero(vVol) ){
                    res = 0.0;
                }
                //vVol != 0
                else{
                    //coefficients c
				    Complex coeff0 = mRR - gamma;
                    Complex coeff1 = mRR + gamma;
				    Complex coeff2 = vVolSQ * u2;
				    Complex coeffUp = coeff2 - coeff0 ;
				    Complex coeffDown = coeff2 - coeff1;
				
    				//logComplex using Jackel method
                    //c = 0
                    if( Maths::isZero(coeffUp) ){
                        res = 0.0;
                    }
                    else{
                        //limiting case, no log complex
                        if( Maths::isZero(coeffDown) ){
                            res = - gamma * tau;
                        }
                        else{
                            //complex c
                            Complex cComp = coeffUp / coeffDown;
                            //arg
                            double tc = Complex::argument(cComp);                 
                            //phase change
                            double bd = - gamma.imag();

                            //numerator
                            Complex gamCoeff = - gamma * tau;
                            Complex expCoeff = exp(gamCoeff);
                            Complex cNum = cComp * expCoeff - 1.0;
                            double rNum = Complex::absSquare(cNum);
                            double kiNum = Complex::argument(cNum);
                            double nNum = floor( (tc + Maths::PI + bd * tau) / (2.0 * Maths::PI) );

                            //denominator
                            Complex cDenom = cComp - 1.0;
                            double rDenom = Complex::absSquare(cDenom);
                            double kiDenom = Complex::argument(cDenom);
                            double nDenom = floor( (tc + Maths::PI) / (2.0 * Maths::PI) );
    
                            //double argCoeff = kiNum - kiDenom + 2.0 * Maths::PI * (nNum - nDenom);
                            double argCoeff = kiNum - kiDenom;
                            double rCoeff = 0.5 * log(rNum / rDenom);
                            Complex logComp = Complex(rCoeff, argCoeff);
                            Complex logCoeff = - 2.0 * mRR * meanSQ / vVolSQ;       
		    		        Complex resPt1 = logComp * logCoeff; 
                    
    				        Complex resPt2 = (coeff0 / vVolSQ) * mRR * meanSQ * tau;
				
	    			        res = resPt1 + resPt2;
                        }   
                    }
                }
            }
            return res;
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

    //computes alphaCDPS-I (see Regis paper)
    //infinite expansion formula
    static Complex calcAlphaBetaCommon_AlphaJumpNew(double tau, //tau = T - t (in years)
                                                    const Complex& u1, //quad var
                                                    const Complex& u2, //spot variance
                                                    const VolSVCJ* svcj){
        static const string method = "CalcAlpha::CalcAlphaBetaCommon_AlphaJumpNew";
        try{
            //result
            Complex res = 0.0;
            
            double dCorr = svcj->correlation;
            double jCorr = svcj->stockVolCrashSizeCorrelation;
            double jRate = svcj->commonCrashRate;
            double jVolMean = svcj->commonVolCrashSizeMean;
            double jStkMean = svcj->commonStockCrashSizeMean;
            double jStkWidth = svcj->commonStockCrashSizeUncertainty;
            double mr = svcj->meanReversRate;
            double vVol = svcj->volVol;

            //coefficients
            double delta = jStkWidth; //delta in Regis paper
            double mu = jVolMean; //mu in Regis paper
            double deltaSQ = delta * delta;
            double gammaBar = log(1.0 + jStkMean) - 0.5 * deltaSQ;
            double rhoJ = jCorr;
            double lambda = jRate;

            //resPt1
            double resPt1 = - tau;
            
            //resPt2
            Complex SQCoeff = 1.0 - 2.0 * u1 * deltaSQ;
            Complex resPt2Coeff = 1.0 / sqrt(SQCoeff);
            Complex resPt2 = resPt2Coeff * CommonAlphaHelper::AlphaSumN(NInfExpan,
                                                                        mr,
                                                                        vVol,
                                                                        mu,
                                                                        delta,
                                                                        lambda,
                                                                        rhoJ,
                                                                        gammaBar,
                                                                        tau,
                                                                        u1,
                                                                        u2);
            res = lambda * (resPt1 + resPt2);
            return res;
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

    //computes alphaCDPS-I using Pade approximant(see Regis paper)
    //infinite expansion formulae transformed into its Pade equivalent
    static Complex calcAlphaBetaCommon_AlphaJumpNewPade(double tau, //tau = T - t (in years)
                                                        const Complex& u1, //quad var
                                                        const Complex& u2, //spot variance
                                                        const VolSVCJ* svcj){
        static const string method = "CalcAlpha::CalcAlphaBetaCommon_AlphaJumpNew";
        try{
            //result
            Complex res = 0.0;
            
            double dCorr = svcj->correlation;
            double jCorr = svcj->stockVolCrashSizeCorrelation;
            double jRate = svcj->commonCrashRate;
            double jVolMean = svcj->commonVolCrashSizeMean;
            double jStkMean = svcj->commonStockCrashSizeMean;
            double jStkWidth = svcj->commonStockCrashSizeUncertainty;
            double mr = svcj->meanReversRate;
            double vVol = svcj->volVol;

            //coefficients
            double delta = jStkWidth; //delta in Regis paper
            double mu = jVolMean; //mu in Regis paper
            double deltaSQ = delta * delta;
            double gammaBar = log(1.0 + jStkMean) - 0.5 * deltaSQ;
            double rhoJ = jCorr;
            double lambda = jRate;

            //resPt1
            double resPt1 = - tau;
            
            //resPt2
            Complex SQCoeff = 1.0 - 2.0 * u1 * deltaSQ;
            Complex resPt2Coeff = 1.0 / sqrt(SQCoeff);
            Complex resPt2 = resPt2Coeff * CommonAlphaHelper::AlphaSumNPade(NInfExpan,
                                                                            mr,
                                                                            vVol,
                                                                            mu,
                                                                            delta,
                                                                            lambda,
                                                                            rhoJ,
                                                                            gammaBar,
                                                                            tau,
                                                                            u1,
                                                                            u2);
            res = lambda * (resPt1 + resPt2);
            return res;
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

    /****** end of GAD 19/02/2006 ******/
};


/************************ end helper functions ***************/

VolSVCJ::VolSVCJ(const string& name,
                 double initialVol,
                 double correlation,
                 double volVol,
                 double meanVol,
                 double meanReversRate,
                 double commonCrashRate,
                 double commonStockCrashSizeMean,
                 double commonStockCrashSizeUncertainty,
                 double commonVolCrashSizeMean,
                 double stockVolCrashSizeCorrelation,
                 bool randInitVol):
VolBaseParam(TYPE, name), initialVol(initialVol), correlation(correlation), volVol(volVol), meanVol(meanVol),
meanReversRate(meanReversRate), commonCrashRate(commonCrashRate), 
commonStockCrashSizeMean(commonStockCrashSizeMean), commonStockCrashSizeUncertainty(commonStockCrashSizeUncertainty),
commonVolCrashSizeMean(commonVolCrashSizeMean), stockVolCrashSizeCorrelation(stockVolCrashSizeCorrelation),
randInitVol(randInitVol) {
    validatePop2Object();
}

CVolProcessed* VolSVCJ::getProcessedVol(const CVolRequest* volRequest,
                                       const CAsset*      /*asset*/) const{
    if (VolRequestRaw::TYPE->isInstance(volRequest) || 
        VolRequestTime::TYPE->isInstance(volRequest)){
        // it's ours or can just use this
        return const_cast<VolSVCJ*>(this);
    }
    else if (ATMVolRequest::TYPE->isInstance(volRequest) ||
             LinearStrikeVolRequest::TYPE->isInstance(volRequest)) {
        // FWD_AT_MAT request for protected assets will ask for this
        // Looking just for something that doesn't break
        HolidaySP noHols(Holiday::noHolidays());
        TimeMetricSP tm(new TimeMetric(1.0, noHols.get()));
        //double flatFXVol = !compVol.empty()? compVol[0] : spotVol[0];
        double flatFXVol = initialVol;
        FlatVolSP flatVol(new FlatVol(this->getName(), 
                                      baseDate, // baseDate
                                      tm.get(),
                                      flatFXVol));
        return flatVol->getProcessedVol(volRequest, 0);
    }
    throw ModelException("VolSVCJ:getProcessedVol", 
                         "Request of type "+
                         volRequest->getClass()->getName()+
                         " not supported for " +this->getName());
}

VolSVSP VolSVCJ::convert(VolSV* p) const {
    VolSVSP volSV;
    if(!randInitVol) {
        if( Maths::isZero(commonCrashRate) || 
            (Maths::isZero(commonStockCrashSizeMean) && Maths::isZero(commonStockCrashSizeUncertainty) &&
             Maths::isZero(commonVolCrashSizeMean))) {
                // Convert if no jumps
                volSV = VolSVSP(new VolSV(getName(), initialVol, meanVol, meanReversRate, volVol, correlation));
        }
    }
    
    return volSV;
}
    

VolSVJSP VolSVCJ::convert(VolSVJ* p) const {
    VolSVJSP volSVJ;
    if(!randInitVol) {
        if(Maths::isZero(commonCrashRate) || Maths::isZero(commonVolCrashSizeMean)) {
            // If no jumps at all or if no jumps in vol
            volSVJ = VolSVJSP(new VolSVJ(
                getName(), initialVol, correlation, volVol, meanVol, meanReversRate,
                commonCrashRate, commonStockCrashSizeMean, commonStockCrashSizeMean, 0.0));
        }
    }
    return volSVJ;
}


VolSVJJSP VolSVCJ::convert(VolSVJJ* p) const {
    VolSVJJSP volSVJJ;
    if(!randInitVol) {
        volSVJJ = VolSVJJSP(new VolSVJJ(
            getName(), initialVol, meanVol, meanReversRate, volVol, correlation, 0.0, 0.0, 0.0, 0.0, 0.0, 
            commonCrashRate, commonStockCrashSizeMean, commonStockCrashSizeMean, commonVolCrashSizeMean, 
            stockVolCrashSizeCorrelation));
    }

    return volSVJJ;
}


void VolSVCJ::SVCJVolParam::ComputeImpVol(const CVolBase*          vol,
                                          const CLatticeDouble&    strikes,
                                          const DateTimeArray&     maturities,
                                          CLatticeDouble&          impV) const{
    // turn the vol into what we must have
    const VolSVCJ* myVol = static_cast<const VolSVCJ *>(vol);
    // then just pass through the parameterised vol
    myVol->ComputeImpVol(strikes, maturities, impV);
}

/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolSVCJ::SVCJVolParam::spotVolSurfaceFromStrikes(
    const CVolBase*       vol,
    const CDoubleArray&   strikes) const{
    // turn the vol into what we must have
    const VolSVCJ* myVol = 
        static_cast<const VolSVCJ *>(vol);
    // then just pass through the parameterised vol
    return myVol->spotVolSurfaceFromStrikes(strikes);
}

void VolSVCJ::validatePop2Object(){
    try {
        Calibrator::IAdjustable::checkRange(this);
        if (randInitVol && meanReversRate < 1.0e-10)
            throw ModelException("VolSVCJ::validatePop2Object",
                        "meanReversionRate too small for random initial variance.");
        // build diffusion + stockCrash + volCrash + commonCrash
    }
    catch(exception& e){
        throw ModelException(e, "VolSVCJ::validatePop2Object");
    }
}

void VolSVCJ::ComputeImpVol(const CLatticeDouble&      strikes,
                            const DateTimeArray&       maturities,
                            CLatticeDouble&            impV) const {
    static const string routine("VolSVCJ::ComputeImpVol");
    throw ModelException(routine, "Not supported");
    if ((maturities.size() != strikes.size()) ||
        (maturities.size() != impV.size())) {
        throw ModelException(routine, "Size mismatch between strikes ("+ 
                             Format::toString(strikes.size()) +
                             "), maturities ("+ 
                             Format::toString(maturities.size())+
                             ") and impV ("+ 
                             Format::toString(impV.size())+ ")");
    }
    
    for (int iMat = 0; iMat < maturities.size(); iMat++) {
        if (strikes[iMat].size() != impV[iMat].size()){
            throw ModelException(routine, "Size mismatch between strikes"
                                 " & maturities for Mat " +
                                 maturities[iMat].toString() +
                                 " (n "+ Format::toString(iMat) + ")");
        }
        for (int iStrike = 0; iStrike < strikes[iMat].size(); iStrike ++) {
            impV[iMat][iStrike] = 0.0;
        }
    }
}

/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
VolSurface* VolSVCJ::spotVolSurfaceFromStrikes(
    const CDoubleArray&   strikes) const{
    static const string routine("VolSVCJ::spotVolSurfaceFromStrikes");
    try{
        throw ModelException(routine, "Not supported");

        const VolSurface* backbone = getBackboneSurface();
        const DateTimeArray& dates = backbone->getDates();
        CDoubleMatrix matrix(strikes.size(),
                             dates.size());
        for (int iStrike = 0; iStrike < strikes.size(); iStrike++) {
            for (int iMat = 0; iMat < dates.size(); iMat++) {
                matrix[iStrike][iMat] = 0.0;
            }
        }
        /** for performance need constructor that takes in
            cached values (to do) */
        VolSurface* volSurf = 
            new VolSurface(getName(),
                           timeMetric.get(),
                           strikes,
                           matrix,
                           backbone->getExpiries().get(),
                           baseDate);
        return volSurf;
    } catch (exception& e){
        throw ModelException(e, routine);
    }
} 

/** Invoked when Class is 'loaded' */
void VolSVCJ::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VolSVCJ, clazz);
    SUPERCLASS(VolBaseParam);
    IMPLEMENTS(Calibrator::IAdjustable);
    IMPLEMENTS(IDynamicsParameter);
    IMPLEMENTS(MarketDataConvert::IConvert<VolSV>);
    IMPLEMENTS(MarketDataConvert::IConvert<VolSVJ>);
    IMPLEMENTS(MarketDataConvert::IConvert<VolSVJJ>);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(initialVol, "initialVol");
    FIELD(meanVol, "meanVol");
    FIELD(meanReversRate, "meanReversRate");
    FIELD(volVol, "volVol");
    FIELD(correlation, "correlation");
    FIELD(commonCrashRate, "commonCrashRate");
    FIELD(commonStockCrashSizeMean, "commonStockCrashSizeMean");
    FIELD(commonStockCrashSizeUncertainty, "commonStockCrashSizeUncertainty");
    FIELD(commonVolCrashSizeMean, "commonVolCrashSizeMean");
    FIELD(stockVolCrashSizeCorrelation, "stockVolCrashSizeCorrelation");
    FIELD(randInitVol, "if true, use random initial variance (gamma distribution)");
    FIELD_MAKE_OPTIONAL(randInitVol);

    // transient        
    FIELD(baseDate, "");
    FIELD_MAKE_TRANSIENT(baseDate);
    // add our fields and their ranges to central list
    Calibrator::IAdjustable::registerField(
        clazz, "initialVol",
        new Range(Heston::RangeDef::initialVol));
    Calibrator::IAdjustable::registerField(
        clazz, "meanVol", 
        new Range(Heston::RangeDef::meanVol));
    Calibrator::IAdjustable::registerField(
        clazz, "meanReversRate",
        new Range(Heston::RangeDef::meanReversRate));
    Calibrator::IAdjustable::registerField(
        clazz, "volVol", 
        new Range(Heston::RangeDef::volVol));
    Calibrator::IAdjustable::registerField(
        clazz, "correlation", 
        new Range(Heston::RangeDef::correlation));

    Calibrator::IAdjustable::registerField(
        clazz, "commonCrashRate",
        new Range(Heston::CommonCrash::RangeDef::crashRate));
    Calibrator::IAdjustable::registerField(
        clazz, "commonStockCrashSizeMean", 
        new Range(Heston::CommonCrash::RangeDef::stockCrashSizeMean));
    Calibrator::IAdjustable::registerField(
        clazz, "commonStockCrashSizeUncertainty", 
        new Range(Heston::CommonCrash::RangeDef::stockCrashSizeUncertainty));
    Calibrator::IAdjustable::registerField(
        clazz, "commonVolCrashSizeMean", 
        new Range(Heston::CommonCrash::RangeDef::volCrashSizeMean));
    Calibrator::IAdjustable::registerField(
        clazz, "stockVolCrashSizeCorrelation", 
        new Range(Heston::CommonCrash::RangeDef::crashCorrelation));
}

VolSVCJ::VolSVCJ(CClassConstSP clazz) :
VolBaseParam(clazz),
initialVol(Heston::DefaultVal::initialVol),
correlation(-0.9999),
volVol(Heston::DefaultVal::volVol),
meanVol(Heston::DefaultVal::meanVol),
meanReversRate(Heston::DefaultVal::meanReversRate),
commonCrashRate(Heston::CommonCrash::DefaultVal::crashRate),
commonStockCrashSizeMean(Heston::CommonCrash::DefaultVal::stockCrashSizeMean),
commonStockCrashSizeUncertainty(Heston::CommonCrash::DefaultVal::stockCrashSizeUncertainty),
commonVolCrashSizeMean(Heston::CommonCrash::DefaultVal::volCrashSizeMean),
stockVolCrashSizeCorrelation(Heston::CommonCrash::DefaultVal::crashCorrelation),
randInitVol(true){}

/** populate from market cache */
void VolSVCJ::getMarket(const IModel* model, const MarketData* market) {

    // call parent method first
    VolBaseParam::getMarket(model, market);
    
    // const VolSurface* backbone = getBackboneSurface();
    /** returns a constant reference to surface to be used for the backbone */
    VolSurfaceSP backbone = VolSurfaceSP(VolSurfaceSP::dynamicCast(
                market->GetData(getName(),VolSurface::TYPE)));
    backbone->getMarket(model,market);  // put holiday and so on into volsurface.
    baseDate = backbone->getBaseDate();
}

/* Started Log Return */
Complex VolSVCJ::scalelessCumulant(const StFourierProcessLogRtn& process,
                                   const StFourierProductLogRtn& product, 
                                   const Complex& z, 
                                   const DateTime& matDate) const{
    static const string method = "VolSVCJ::scalelessCumulant";
    try{
        double tau = timeMetric->yearFrac(baseDate,
                                          matDate);

        Complex alpha, beta;
        calcJointLapAlphaBeta(tau,
                              z,
                              0.0,
                              0.0,
                              alpha,
                              beta);
        
		/* Add invariant component to get random initial variance (gamma distribution) */
		if (randInitVol) {
			alpha += CalcAlpha::calcAlphaStationary(beta, this);
			beta = 0.0;
		}
        
        return (alpha + beta * Maths::square(initialVol));
        
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Fwd starting Log Return  */
Complex VolSVCJ::scalelessCumulant(const FwdStFourierProcessLogRtn& process,
                                   const FwdStFourierProductLogRtn& product, 
                                   const Complex& z, 
                                   const DateTime& matDate) const{
    static const string method = "VolSVCJ::scalelessCumulant";
    try{
        /* Time fractions from today till start of option and then from
           start of option till maturity */
        double t = timeMetric->yearFrac(baseDate,
                                        product.getStartDate());
        double tau = timeMetric->yearFrac(product.getStartDate(),
                                          matDate);

        /* Calculate alpha, beta component from start date till maturity */
        Complex alpha_tT, beta_tT;
        calcJointLapAlphaBeta(tau,
                              z,
                              0.0,
                              0.0,
                              alpha_tT,
                              beta_tT);

        /* Then, from today till start date */
        Complex alpha_t, beta_t;
        calcJointLapAlphaBeta(t,
                              0.0,
                              beta_tT,
                              0.0,
                              alpha_t,
                              beta_t);
        
		/* Add invariant component to get random initial variance (gamma distribution) */
		if (randInitVol) {
			alpha_t += CalcAlpha::calcAlphaStationary(beta_t, this);
			beta_t = 0.0;
		}

        return (alpha_tT + alpha_t + beta_t * Maths::square(initialVol));

    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/* Started integrated variance */
Complex VolSVCJ::cumulant(const StFourierProcessIntVar& process,
                          const StFourierProductIntVar& product, 
                          const Complex z, 
                          const DateTime& matDate) const {    
    static const string method = "VolSVCJ::cumulant";
    try{
        double tau = timeMetric->yearFrac(baseDate,
                                          matDate);

        Complex alpha, beta;
        calcJointLapAlphaBeta(tau,
                              0.0,
                              0.0,
                              z,
                              alpha,
                              beta);
       
		/* Add invariant component to get random initial variance (gamma distribution) */
		if (randInitVol) {
			alpha += CalcAlpha::calcAlphaStationary(beta, this);
			beta = 0.0;
		}

        return (alpha + beta * Maths::square(initialVol));

    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}
    
/* Forward Starting integrated variance */
Complex VolSVCJ::cumulant(const FwdStFourierProcessIntVar& process,
                          const FwdStFourierProductIntVar& product, 
                          const Complex z, 
                          const DateTime& matDate) const {        
    static const string method = "VolSVCJ::cumulant";
    try{
        /* Time fractions from today till start of option and then from
           start of option till maturity */
        double t = timeMetric->yearFrac(baseDate,
                                        product.getStartDate());
        double tau = timeMetric->yearFrac(product.getStartDate(),
                                          matDate);

        /* Calculate alpha, beta component from start date till maturity */
        Complex alpha_tT, beta_tT;          
		calcJointLapAlphaBeta(tau,
                              0.0,
                              0.0,
                              z,
                              alpha_tT,
                              beta_tT);

        /* Then, from today till start date */
        Complex alpha_t, beta_t;
        calcJointLapAlphaBeta(t,
                              0.0,
                              beta_tT,
                              0.0,
                              alpha_t,
                              beta_t);
        
        /* Add invariant component to get random initial variance (gamma distribution) */
		if (randInitVol) {
			alpha_t += CalcAlpha::calcAlphaStationary(beta_t, this);
			beta_t = 0.0;
		}

        return (alpha_tT + alpha_t + beta_t * Maths::square(initialVol));

    }
    catch(exception& e){
        throw ModelException(e, method);
    }
} 

/****** GAD 19/02/2006 ******/

/* Started quadratic variation */
Complex VolSVCJ::cumulant(const StFourierProcessQuadVar& process,
                          const StFourierProductQuadVar& product, 
                          const Complex&                 z, 
                          const DateTime&                matDate) const {    
    static const string method = "VolSVCJ::cumulant";
    try{
        Complex res;
        double tau = timeMetric->yearFrac(baseDate,
                                          matDate);

        Complex alpha, beta;
        calcJointLapAlphaBetaCommon(tau,
                                    z, //quad var
                                    0.0, //spot variance
                                    alpha,
                                    beta);
        
		/* Add invariant component to get random initial variance (gamma distribution) */
		if (randInitVol) {
			alpha += CalcAlpha::calcAlphaStationary(beta, this);
			beta = 0.0;
		}
        res = alpha + beta * Maths::square(initialVol);
        return res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

//expectation
double VolSVCJ::expectation(const StFourierProcessQuadVar& process,
                            const StFourierProductQuadVar& product, 
                            const DateTime& matDate) const{    
    static const string method = "VolSVCJ::cumulant";
    try{
        //result
        double res;

        //coefficients
        double nu0 = initialVol * initialVol;
        double theta = meanVol * meanVol;
        double ka = meanReversRate;
        double lamb = commonCrashRate;
        double k = commonStockCrashSizeMean;
        double delta = commonStockCrashSizeUncertainty;
        double mu = commonVolCrashSizeMean;
        double rho = stockVolCrashSizeCorrelation;
        
        double t = timeMetric->yearFrac(baseDate, matDate);
        
        double rhoSQ = rho * rho;
        double muSQ = mu * mu;
        double deltaSQ = delta * delta;

        double gamma = log(1.0 + k) - 0.5 * deltaSQ;
        
        double aPt1 = gamma * gamma;
        double aPt2 = 2.0 * rho * mu * gamma;
        double aPt3 = 2.0 * rhoSQ * muSQ + deltaSQ;
        double a = aPt1 + aPt2 + aPt3;

        double resPt1 = theta + lamb * mu / ka;
        double resPt2 = nu0 - resPt1;
        double kaT = ka * t;
        double resPt3 = ( 1.0 - exp(-kaT) ) / kaT;
        double resPt4 = lamb * a;

        res = resPt1 + resPt2 * resPt3 + resPt4;

        return res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}
    
/* Forward Starting quadratic variation */
Complex VolSVCJ::cumulant(const FwdStFourierProcessQuadVar& process,
                          const FwdStFourierProductQuadVar& product, 
                          const Complex&                    z, 
                          const DateTime&                   matDate) const {        
    static const string method = "VolSVCJ::cumulant";
    try{
        Complex res;
        /* Time fractions from today till start of option and then from
           start of option till maturity */
        double t = timeMetric->yearFrac(baseDate,
                                        product.getStartDate());
        double tau = timeMetric->yearFrac(product.getStartDate(),
                                          matDate);

        /* Calculate alpha, beta component from start date till maturity */
        Complex alpha_tT, beta_tT;
        calcJointLapAlphaBetaCommon(tau,
                                    z, //quad var
                                    0.0, //spot variance
                                    alpha_tT,
                                    beta_tT);

        /* Then, from today till start date */
        Complex alpha_0t, beta_0t;
        calcJointLapAlphaBetaCommon(t,
                                    0.0, //quad var
                                    beta_tT, //spot variance
                                    alpha_0t,
                                    beta_0t);

        /* Add invariant component to get random initial variance (gamma distribution) */
		if (randInitVol) {
			alpha_0t += CalcAlpha::calcAlphaStationary(beta_0t, this);
			beta_0t = 0.0;
		}
        res = alpha_tT + alpha_0t + beta_0t * Maths::square(initialVol);
        return(res);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

//expectation
double VolSVCJ::expectation(const FwdStFourierProcessQuadVar& process,
                            const FwdStFourierProductQuadVar& product, 
                            const DateTime& matDate) const{    
    static const string method = "VolVSCurve::cumulant FwdStFourierProcessQuadVar";
    try{
        //result
        double res;

        //coefficients
        double nu0 = initialVol * initialVol;
        double theta = meanVol * meanVol;
        double ka = meanReversRate;
        double lamb = commonCrashRate;
        double k = commonStockCrashSizeMean;
        double delta = commonStockCrashSizeUncertainty;
        double mu = commonVolCrashSizeMean;
        double rho = stockVolCrashSizeCorrelation;
        double T = timeMetric->yearFrac(baseDate, matDate);
        double t = timeMetric->yearFrac(baseDate, product.getStartDate());

        double rhoSQ = rho * rho;
        double muSQ = mu * mu;
        double deltaSQ = delta * delta;

        double gamma = log(1.0 + k) - 0.5 * deltaSQ;
        
        double aPt1 = gamma * gamma;
        double aPt2 = 2.0 * rho * mu * gamma;
        double aPt3 = 2.0 * rhoSQ * muSQ + deltaSQ;
        double a = aPt1 + aPt2 + aPt3;

        double res_tPt1 = theta + lamb * mu / ka;
        double res_tPt2 = nu0 - res_tPt1;
        double res_tPt3 = ( 1.0 - exp(- ka * t) ) / (ka * t);
        double res_tPt4 = lamb * a;
        double resCoeff_t = t/(T-t);

        double res_TPt1 = theta + lamb * mu / ka;
        double res_TPt2 = nu0 - res_TPt1;
        double res_TPt3 = ( 1.0 - exp(- ka * T) ) / (ka * T);
        double res_TPt4 = lamb * a;
        double resCoeff_T = T/(T-t);

        double res_t = res_tPt1 + res_tPt2 * res_tPt3 + res_tPt4;
        double res_T = res_TPt1 + res_TPt2 * res_TPt3 + res_TPt4;

        res = resCoeff_T * res_T - resCoeff_t * res_t;

        return res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/*** end of GAD 19/02/2006 ***/

/** Calculates the components alpha and beta that appear in 
    the exponent of the time-t joint Bilateral Laplace 
    of X_T = (Y_T, V_T, I_T) where Y_T = ln(S_T / F(0, T)) is the dimension-less
    log spot at time T, V_T is the instantaneous variance at time T
    and I_T is the instantaneous variance from 0 to time T
    The Bilateral Laplace transform is of the form
        exp(alpha(tau, u1, u2, u3) 
            + u1 * Y_t 
            + beta(tau, u1, u2, u3) * V_t 
            + u3 * I_t)
    where tau = T - t and the complex numbers u1, u2, u3 are the frequencies 
    wrt Y_T, V_T and I_T, respectively. */
void VolSVCJ::calcJointLapAlphaBeta(double         tau,    // tau = T - t (in years)
                                    const Complex& u1,
                                    const Complex& u2,
                                    const Complex& u3,
                                    Complex&       alpha,
                                    Complex&       beta) const{
    static const string method = "VolSVCJ::calcJointLapAlphaBeta";

    try{    // in case Complex operations fail for some reason
        /* Calculate Heston contribution */
		CalcAlpha::calcAlphaBeta(tau,
								 u1,
								 u2,
								 u3,
								 alpha,
								 beta,
								 this);
		
        /* Calculate common jump contribution to alpha */
        alpha += CalcAlpha::calcAlpha(tau,
										u1,
                                        u2,
                                        u3,
                                        this);

    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

/****** GAD 19/02/2006 ******/

//common alpha component for quadratic variation
void VolSVCJ::calcJointLapAlphaBetaCommon(double         tau, // tau = T - t (in years)
                                          const Complex& u1, // quad var
                                          const Complex& u2, //spot variance
                                          Complex&       alpha,
                                          Complex&       beta) const{
    static const string method = "VolSVCJ::calcJointLapAlphaBetaCommon";
    try{       
        //betaH
        beta = CalcAlpha::calcAlphaBetaCommon_Beta(tau, // tau = T - t
                                                   u1, // quad var
                                                   u2, // spot variance
                                                   this);
        //alphaH
        Complex alphaDiff = CalcAlpha::calcAlphaBetaCommon_Alpha(tau, // tau = T - t
                                                                 u1, // quad var
                                                                 u2, // spot variance
                                                                 this);
/*
        //alphaDPS-I
        Complex alphaJump = CalcAlpha::calcAlphaBetaCommon_AlphaJumpNew(tau, //tau = T - t
                                                                        u1, // quad var
                                                                        u2, // spot variance
                                                                        this);
*/
        //alphaDPS-I Pade
        Complex alphaJump = CalcAlpha::calcAlphaBetaCommon_AlphaJumpNewPade(tau, //tau = T-t
                                                                            u1, // quad var
                                                                            u2, //spot variance
                                                                            this);
        //alpha
        alpha = alphaDiff + alphaJump;
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

/****** end of GAD 19/02/2006 ******/

double VolSVCJ::getSVCJParam(const EnumParam param) const
{
    switch (param)
    {
    case INITIAL_VOL:
        return initialVol;
    case CORRELATION:
        return correlation;
    case VOL_VOL:
        return volVol;
    case MEAN_VOL:
        return meanVol;
    case MEAN_REVERS_RATE:
        return meanReversRate;
    case COMMON_CRASH_RATE:
        return commonCrashRate;
    case COMMON_STOCK_CRASH_SIZE_MEAN:
        return commonStockCrashSizeMean;
    case COMMON_STOCK_CRASH_SIZE_UNCERTAINTY:
        return commonStockCrashSizeUncertainty;
    case COMMON_VOL_CRASH_SIZE_MEAN:
        return commonVolCrashSizeMean;
    case STOCK_VOL_CRASH_SIZE_CORRELATION:
        return stockVolCrashSizeCorrelation;
    default:
        throw ModelException("VolSVCJ::getSVCJParam", "unkown param request"); // cannot reach here
        return 0.0;
    }
}

//  -----------------  MC support   ---------------------------------------------------------------

// *********** MC support for simulating arrival times *********************/

/** compute jump contributions for var and log spots
jumpsAdded[0] -> tradYrs, [1] -> instVars, [2] -> integratedVars, [3] -> logSpots */
void VolSVCJ::addJumps(list<double>& tradYrs,
                        const double* jumpTimes,
                        int           numJumps,
                        BoolArray     jumpParticipation,
                        const double* jumpDeviatesVar,
                        const double* jumpDeviatesSpot,
                        list<STEP_TYPE>& stepTypes,
                        list<double>& instVars,
                        list<double>& integratedVars,
                        list<double>& logSpots,
                        vector<double>& varJumps,
                        vector<STEP_TYPE_ITER>& jumpTypesAdded,
                        vector<vector<LIST_ITER> >&jumpsAdded,
                        double          randomInitVar) const
{
    ASSERT(jumpsAdded.size() == 4);

    LIST_ITER t0 = tradYrs.begin();
    // start from step 1
    LIST_ITER t = ++tradYrs.begin();
    LIST_ITER v = ++instVars.begin();
    LIST_ITER w = ++integratedVars.begin();
    LIST_ITER s = ++logSpots.begin();
    
    STEP_TYPE_ITER type = ++stepTypes.begin();

    // init
    double startVar = initialVol*initialVol;
    if (randInitVol){
        startVar = randomInitVar*volVol*volVol/meanReversRate/2.0;
    }

    *instVars.begin() = startVar;
    *integratedVars.begin() = 0.0;
    *logSpots.begin() = 0.0;

    int nbSteps = tradYrs.size();

    double logMean = log(1.0 + commonStockCrashSizeMean) 
                    - 0.5 * Maths::square(commonStockCrashSizeUncertainty);

    //  insert jumpTimes
    for (int i = 0, j=0; i < nbSteps-1; i++, t0=t++, v++, w++, type++){

        // init
        *v = 0.0;
        *s = 0.0;

        if (j<numJumps && *t0 < jumpTimes[j] && jumpTimes[j] <= *t)
        {
            if (Maths::isZero(jumpTimes[j] - *t)){
                // do we ever land in here?
                t = tradYrs.insert(t, jumpTimes[j] - 1.0/365.0/24.0/60.0); // move by 1 min or less
            }
            else
                t = tradYrs.insert(t, jumpTimes[j]);

            double uniform = Maths::max(N1(jumpDeviatesVar[j]), 0.00000001);
            double vJump = 0.0;
            double sJump = 0.0;
            
            if (jumpParticipation[j] == true)
            {
                vJump = -commonVolCrashSizeMean*log(uniform);
                sJump = logMean + stockVolCrashSizeCorrelation*vJump
                    + commonStockCrashSizeUncertainty * jumpDeviatesSpot[j];
            }

            varJumps.push_back(vJump);
            jumpTypesAdded.push_back(type = stepTypes.insert(type, JUMP_DATE));
            jumpsAdded[0].push_back(t);
            jumpsAdded[1].push_back(v = instVars.insert(v, vJump));
            jumpsAdded[2].push_back(w = instVars.insert(w, 0.0));
            jumpsAdded[3].push_back(s = logSpots.insert(s, sJump));

            nbSteps++;
            j++;
        }
        if (*type == SAMPLE_DATE || *type == JUMP_DATE){
            s++;
        }
    }
}

/** Given an array of gaussian deviates, populate the arrays 'instVars' 
    and 'integratedVars' with simulated values using an Euler
    discretization scheme */
void VolSVCJ::generateVarPathsEuler(LIST_ITER          tradYears,
                                    int                nbSteps,
                                    const double*      deviates,
                                    LIST_ITER          instVars,
                                    LIST_ITER          integratedVars) const
{
    LIST_ITER t0, v0, w0;

    // double meanVar = Maths::square(meanVol); now 
    for (int iStep = 0; iStep < nbSteps - 1; iStep++) {
        t0 = tradYears++;
        v0 = instVars ++;
        w0 = integratedVars ++;
        double dt = *tradYears - *t0;
        double meanVar = meanVol*meanVol*dt;
        // drift contribution
        //double drift = meanReversRate * (computeMeanVar(*t0) - *v0) * dt;
        double drift = meanReversRate * (meanVar - *v0 * dt);
        // diffusion contribution
        double instVarPlus = Maths::max(0.0, *v0);
        double diffusion = volVol * sqrt(instVarPlus * dt) * deviates[iStep];
        // compute inst variance
        *instVars += *v0 + drift + diffusion;
        // compute integrated var using an Euler scheme too
        *integratedVars = *w0 + instVarPlus * dt;
    }
}

/** Given an array of gaussian deviates, populate the arrays 'instVars' 
    and 'integratedVars' with simulated values using a variable transform
    scheme together with an Euler discretization scheme */
void VolSVCJ::generateVarPathsTransEuler(LIST_ITER          tradYears,
                                         int                nbSteps,
                                         const double*      deviates,
                                         LIST_ITER          instVars,
                                         LIST_ITER          integratedVars) const
{
    LIST_ITER t0, v0, w0;

    for (int iStep = 0; iStep < nbSteps - 1; iStep++) {
        t0 = tradYears++;
        v0 = instVars ++;
        w0 = integratedVars ++;
        double dt = *tradYears - *t0;
        double meanVar = meanVol*meanVol*dt;
        // integrated var is function of last step values only
        *integratedVars = *w0 + *v0 * dt;
        // Euler approx for the vol
        double vol = sqrt(*v0);
        ASSERT(!Maths::isZero(vol));
        double alpha = meanReversRate * meanVar /* computeMeanVar(*t0)*dt */ 
                       - Maths::square(volVol)*dt / 4.0;
        double drift = alpha / vol - meanReversRate * vol * dt;
        double diffusion = volVol * sqrt(dt) * deviates[iStep];
        vol += 0.5*(drift + diffusion);
        *instVars += Maths::square(vol);
    }
}

/** Given two arrays of (independent) gaussian deviates, populate the 
    arrays 'instVars', 'integratedVars' and 'exJumpLogSpots' with simulated 
    values using an Euler discretization scheme for the spot */
void VolSVCJ::generateSpotPathsEuler(LIST_ITER          tradYears,
                                     int                nbSteps,
                                     STEP_TYPE_ITER     stepTypes,
                                     const double*      spotDeviates,
                                     const double*      varDeviates,
                                     LIST_ITER          logSpots,
                                     LIST_ITER          instVars,
                                     LIST_ITER          integratedVars) const
{
    LIST_ITER t0, w0;

    double beta1 = correlation;
    double beta2 = sqrt(1.0 - beta1 * beta1);
    double logSpot = 0.0;
    double x = stockVolCrashSizeCorrelation * commonVolCrashSizeMean;
    double jumpDrift = - commonCrashRate*((commonStockCrashSizeMean+ x)/(1.0-x));
    for (int iStep = 0; iStep < nbSteps-1; iStep++, instVars++) {
        t0 = tradYears++;
        w0 = integratedVars ++;
        stepTypes++;

        double dt = *tradYears - *t0;
        double drift = -0.5 * (*integratedVars - *w0);
        double instVarPlus = Maths::max(0.0, *instVars);
        double diffusion = sqrt(instVarPlus * dt) 
                           * (beta1 * varDeviates[iStep] 
                           + beta2 * spotDeviates[iStep]);
        logSpot += drift + diffusion + jumpDrift*dt;
        if (*stepTypes == SAMPLE_DATE || *stepTypes == JUMP_DATE){
            logSpots++;
            logSpot += *logSpots; // add jump contribution if any
            *logSpots = logSpot;
        }
    }
}

/** Given two arrays of (independent) gaussian deviates, populate the 
    arrays 'instVars', 'integratedVars' and 'exJumpLogSpots' with simulated 
    values. The log spots are simulated exactly */
void VolSVCJ::generateSpotPathsExact(LIST_ITER          tradYears,
                                     int                nbSteps,
                                     STEP_TYPE_ITER     stepTypes,
                                     const double*      spotDeviates,
                                     vector<double>&    varJumps,
                                     LIST_ITER          logSpots,
                                     LIST_ITER          instVars,
                                     LIST_ITER          integratedVars) const
{
    LIST_ITER s0;
    LIST_ITER t0 = tradYears;
    LIST_ITER v0 = instVars; 
    LIST_ITER w0 = integratedVars;

    double mean, sqrtVar;

    double beta1 = correlation;
    double beta2 = sqrt(1.0 - beta1 * beta1);
    double x = stockVolCrashSizeCorrelation * commonVolCrashSizeMean;
    double jumpDrift = - commonCrashRate*((commonStockCrashSizeMean+ x)/(1.0-x));
    for (int iStep = 0, iSpot=0, iJump=0; iStep < nbSteps-1; iStep++) {

        tradYears++; stepTypes++; instVars++; integratedVars++;

        if (*stepTypes == SAMPLE_DATE || *stepTypes == JUMP_DATE){
            s0 = logSpots++;
            double dt = *tradYears - *t0;
            double integratedVarDiff = *integratedVars - *w0;
            // special treatment for the deterministic case
            if(Maths::isZero(volVol)){
                mean = -0.5 * integratedVarDiff;
                sqrtVar = sqrt(integratedVarDiff);
            }
            // general case
            else{
                double instVarDiff = *instVars - *v0; 
                if (*stepTypes == JUMP_DATE){
                    instVarDiff -= varJumps[iJump];
                    iJump++;
                }
                double meanVar = meanVol*meanVol*dt;
                double meanAdjTerm = (instVarDiff 
                                      - meanReversRate * (meanVar - integratedVarDiff))
                                     / volVol;
                mean = -0.5 * integratedVarDiff + beta1 * meanAdjTerm;
                sqrtVar = beta2 * sqrt(integratedVarDiff);
            }
            *logSpots += *s0 + mean + sqrtVar * spotDeviates[iSpot] + jumpDrift*dt;
            iSpot++;
            t0 = tradYears;
            v0 = instVars; 
            w0 = integratedVars;
        }
    }
}

DEFINE_TEMPLATE_TYPE(VolSVCJArray);

CClassConstSP const VolSVCJ::TYPE =
CClass::registerClassLoadMethod("VolSVCJ", typeid(VolSVCJ), load);

CClassConstSP const VolSVCJ::SVCJVolParam::TYPE =
CClass::registerClassLoadMethod("VolSVCJ::SVCJVolParam", typeid(SVCJVolParam), load);



/************************/
/*** SVCJ + LV scheme ***/
/************************/

/** Needed for IAdjustable interface. Returns market data name for vol */
string VolSVCJ::getName() const{
    return CVolBaseParamSurface::getName();
}

/* external symbol to allow class to be forced to be linked in */
bool VolSVCJLinkIn(){
    return (VolSVCJ::TYPE != 0);
}


void VolSVCJLV::validatePop2Object(){
    static const string method = "VolSVCJLV::validatePop2Object";
    try {
	    VolSVCJ::validatePop2Object();

        // check that LV A parameters are provided
        if (!locVolParamsA.get())
        {
            throw ModelException(method,
                                 "Local Volatility A parameters were not found.");
        }
        // check that LV B parameters are provided
        if (!locVolParamsB.get())
        {
            throw ModelException(method,
                                 "Local Volatility B parameters were not found.");
        }
        // check that LV C parameters are provided
        if (!locVolParamsC.get())
        {
            throw ModelException(method,
                                 "Local Volatility C parameters were not found.");
        }

        int locVolParamsSize = locVolExpiries->size(); // nber of expiries
        if (locVolParamsSize <= 0)
        {
            throw ModelException(method,
                                 "the number of Local Volatility expiaries should be > 0");
        }

        if (locVolParamsSize != locVolParamsA->size())
        {
            throw ModelException(method,
                                 "the number of Local Volatility A parameters is different from "
                                 "the number of expiries");
        }
        if (locVolParamsSize != locVolParamsB->size())
        {
            throw ModelException(method,
                                 "the number of Local Volatility B parameters is different from "
                                 "the number of expiries");
        }
        if (locVolParamsSize != locVolParamsC->size())
        {
            throw ModelException(method,
                                 "the number of Local Volatility C parameters is different from "
                                 "the number of expiries");
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

/** return a copy of LV expiries */
ExpiryArraySP VolSVCJLV::getLocVolExpiries(const IObject* obj) {
    const VolSVCJLV* vol = dynamic_cast<const VolSVCJLV*>(obj);
    if (obj) {
        if (!vol->locVolExpiries){
            throw ModelException("VolSVCJLV::getLocVolExpiries",
                                 "Internal error: null lv expiries.");
        }
        return ExpiryArraySP(copy(vol->locVolExpiries.get()));
    }
    // should never happen
    throw ModelException("VolSVCJLV::getLocVolExpiries",
                         "Internal error: obj is not of the desired type.");
}

/** compute local volatility function for SVCJ + LV scheme */
double VolSVCJLV::calcLocVol(const double logSpot, const double paramA, const double paramB, const double paramC) const
{
    /* SRM functional form */
    /*
    double locVol = (1.0 + (paramA * tanh(paramC * logSpot)) + (paramB * (1.0 - (1.0 / cosh(paramC * logSpot)))));
    return Maths::max(locVol, 0.0);*/

    /* polynomial functional form */
    double locVol = ((paramA + paramB * logSpot) * logSpot) + paramC;
    double maxLocVol = 4.0 * paramC;
    locVol = Maths::min(Maths::max(locVol, 0.000001), maxLocVol);
    
    if (Maths::isPositive(paramA))
    {
        double disc = paramB * paramB - 4.0 * paramA * paramC; // polynomial discriminant
        if (Maths::isPositive(disc))
        {
            // locVol is always positive except between the polynomial roots
            // two roots on the same side than x = 0
            double minAbsRoot; // root with lowest absolute value
            if (paramB >= 0)
            {
                minAbsRoot = (-paramB + sqrt(disc)) / (2 * paramA);
                locVol = (logSpot <= minAbsRoot) ? 0.0 : locVol;
            }
            else
            {
                minAbsRoot = (-paramB - sqrt(disc)) / (2 * paramA);
                locVol = (logSpot >= minAbsRoot) ? 0.0 : locVol;
            }
        }
    }
    
    if (Maths::isNegative(paramA))
    {
        double disc = paramB * paramB - 4.0 * paramA * (paramC - maxLocVol); // polynomial discriminant
        if (Maths::isPositive(disc))
        {
            // (locVol - maxLocVol) is always negative except between the polynomial roots
            // two roots on the same side than x = 0
            double minAbsRoot; // root with lowest absolute value
            if (paramB >= 0)
            {
                minAbsRoot = (-paramB + sqrt(disc)) / (2 * paramA);
                locVol = (logSpot <= minAbsRoot) ? 0.0 : locVol;
            }
            else
            {
                minAbsRoot = (-paramB - sqrt(disc)) / (2 * paramA);
                (logSpot >= minAbsRoot) ? 0.0 : locVol;
            }
        }
    }
    
    return locVol;
}

/** extend local volatility parameters arrays
    in order to get array with size equal to the number of steps in MC */
void VolSVCJLV::extendLocVolParams(LIST_ITER tradYears,
                                   int nbSteps,
                                   DoubleArraySP locVolParamsAExt,
                                   DoubleArraySP locVolParamsBExt,
                                   DoubleArraySP locVolParamsCExt) const
{
    // trading time fractions for local volatility adjustment in SVCJ + LV scheme
    DoubleArraySP locVolYearFracs = DoubleArraySP(new DoubleArray(locVolExpiries->size()));
    int iDate;
    for (iDate = 0 ; iDate < locVolExpiries->size() ; iDate++)
    {
        // date corresponding to the local volatility expiry
        DateTime thisDate = (*locVolExpiries)[iDate]->toDate(baseDate);
        // trading time fraction from base date to date corresponding to the local volatility expiry
        double timeFrac = timeMetric->yearFrac(baseDate,
                                               thisDate);
        (*locVolYearFracs)[iDate] = timeFrac;
    }

    LIST_ITER t0; 
    int currentIndex = 0; // index of local volatility parameters currently used
    double currentYearFrac = 0.0; // year fraction from base date to current date

    // first step treated separately
    (*locVolParamsAExt)[0] = (*locVolParamsA)[currentIndex];
    (*locVolParamsBExt)[0] = (*locVolParamsB)[currentIndex];
    (*locVolParamsCExt)[0] = (*locVolParamsC)[currentIndex];

    // other steps
    int iStep = 0;
    for (iStep = 0; iStep < nbSteps - 1; iStep++) {
        t0 = tradYears++;
        currentYearFrac += (*tradYears - *t0);

        // if the current year fraction is greater than the current LV yearFrac
        // increment currentIndex by 1
        if ((currentYearFrac > (*locVolYearFracs)[currentIndex]) &&
            (currentIndex < (locVolExpiries->size() - 1)))
        {
            if (currentYearFrac > (*locVolYearFracs)[currentIndex+1])
            {
                throw ModelException("VolSVCJLV::extendLocVolParams",
                                     "LV parameters for maturity Format::toString(currentIndex-1) not used.");
            }
            currentIndex++;
        }

        (*locVolParamsAExt)[iStep+1] = (*locVolParamsA)[currentIndex];
        (*locVolParamsBExt)[iStep+1] = (*locVolParamsB)[currentIndex];
        (*locVolParamsCExt)[iStep+1] = (*locVolParamsC)[currentIndex];
    }
}

/** Given two arrays of (independent) gaussian deviates, populate the 
    arrays 'instVars', 'integratedVars' and 'exJumpLogSpots' with simulated 
    values using an Euler discretization scheme for the spot */
void VolSVCJLV::generateSpotPathsEuler(LIST_ITER          tradYears,
                                       int                nbSteps,
                                       STEP_TYPE_ITER     stepTypes,
                                       const double*      spotDeviates,
                                       const double*      varDeviates,
                                       LIST_ITER          logSpots,
                                       LIST_ITER          instVars,
                                       LIST_ITER          integratedVars) const
{
    // extend local volatility parameters
    DoubleArraySP locVolParamsAExt = DoubleArraySP(new DoubleArray(nbSteps));
    DoubleArraySP locVolParamsBExt = DoubleArraySP(new DoubleArray(nbSteps));
    DoubleArraySP locVolParamsCExt = DoubleArraySP(new DoubleArray(nbSteps));
    extendLocVolParams(tradYears, nbSteps, locVolParamsAExt, locVolParamsBExt, locVolParamsCExt);
    LIST_ITER t0, w0;

    double beta1 = correlation;
    double beta2 = sqrt(1.0 - beta1 * beta1);
    double logSpot = 0.0;
    double x = stockVolCrashSizeCorrelation * commonVolCrashSizeMean;
    double jumpDrift = - commonCrashRate*((commonStockCrashSizeMean+ x)/(1.0-x));
    for (int iStep = 0; iStep < nbSteps-1; iStep++, instVars++) {
        t0 = tradYears++;
        w0 = integratedVars ++;
        stepTypes++;

        // SVCJ + LV scheme
        double locVol =
            calcLocVol(logSpot, (*locVolParamsAExt)[iStep], (*locVolParamsBExt)[iStep], (*locVolParamsCExt)[iStep]);
        double dt = *tradYears - *t0;
        double drift = -0.5 * (*integratedVars - *w0) * locVol * locVol;
        double instVarPlus = Maths::max(0.0, *instVars);
        double diffusion = sqrt(instVarPlus * locVol * locVol * dt)
                           * (beta1 * varDeviates[iStep] 
                           + beta2 * spotDeviates[iStep]);
        logSpot += drift + diffusion + jumpDrift*dt;
        if (*stepTypes == SAMPLE_DATE || *stepTypes == JUMP_DATE){
            logSpots++;
            logSpot += *logSpots; // add jump contribution if any
            *logSpots = logSpot;
        }
    }
}

/** Given two arrays of (independent) gaussian deviates, populate the 
    arrays 'instVars', 'integratedVars' and 'exJumpLogSpots' with simulated 
    values. The log spots are simulated exactly */
void VolSVCJLV::generateSpotPathsExact(LIST_ITER          tradYears,
                                       int                nbSteps,
                                       STEP_TYPE_ITER     stepTypes,
                                       const double*      spotDeviates,
                                       vector<double>&    varJumps,
                                       LIST_ITER          logSpots,
                                       LIST_ITER          instVars,
                                       LIST_ITER          integratedVars) const
{
    // extend local volatility parameters
    DoubleArraySP locVolParamsAExt = DoubleArraySP(new DoubleArray(nbSteps));
    DoubleArraySP locVolParamsBExt = DoubleArraySP(new DoubleArray(nbSteps));
    DoubleArraySP locVolParamsCExt = DoubleArraySP(new DoubleArray(nbSteps));
    extendLocVolParams(tradYears, nbSteps, locVolParamsAExt, locVolParamsBExt, locVolParamsCExt);
    LIST_ITER s0;
    LIST_ITER t0 = tradYears;
    LIST_ITER v0 = instVars; 
    LIST_ITER w0 = integratedVars;

    double mean, sqrtVar;

    double beta1 = correlation;
    double beta2 = sqrt(1.0 - beta1 * beta1);
    double x = stockVolCrashSizeCorrelation * commonVolCrashSizeMean;
    double jumpDrift = - commonCrashRate*((commonStockCrashSizeMean+ x)/(1.0-x));
    for (int iStep = 0, iSpot=0, iJump=0; iStep < nbSteps-1; iStep++) {

        tradYears++; stepTypes++; instVars++; integratedVars++;

        if (*stepTypes == SAMPLE_DATE || *stepTypes == JUMP_DATE){
            s0 = logSpots++;
            double dt = *tradYears - *t0;
            double integratedVarDiff = *integratedVars - *w0;
            // special treatment for the deterministic case
            if(Maths::isZero(volVol)){
                mean = -0.5 * integratedVarDiff;
                sqrtVar = sqrt(integratedVarDiff);
            }
            // general case
            else{
                double instVarDiff = *instVars - *v0; 
                if (*stepTypes == JUMP_DATE){
                    instVarDiff -= varJumps[iJump];
                    iJump++;
                }
                double meanVar = meanVol*meanVol*dt;
                double meanAdjTerm = (instVarDiff 
                                      - meanReversRate * (meanVar - integratedVarDiff))
                                     / volVol;
                mean = -0.5 * integratedVarDiff + beta1 * meanAdjTerm;
                sqrtVar = beta2 * sqrt(integratedVarDiff);
            }
            *logSpots += *s0 + mean + sqrtVar * spotDeviates[iSpot] + jumpDrift*dt;
            iSpot++;
            t0 = tradYears;
            v0 = instVars; 
            w0 = integratedVars;
        }
    }
}

VolSVCJLV::VolSVCJLV() :
    VolSVCJ(TYPE)
{}

/** Invoked when Class is 'loaded' */
void VolSVCJLV::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VolSVCJLV, clazz);
    SUPERCLASS(VolSVCJ);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(locVolExpiries, "expiries for local volatility adjustment in SVCJ + LV scheme");
    FIELD(locVolParamsA, "parameters A for local volatility adjustment in SVCJ + LV scheme");
    FIELD(locVolParamsB, "parameters B for local volatility adjustment in SVCJ + LV scheme");
    FIELD(locVolParamsC, "parameters C for local volatility adjustment in SVCJ + LV scheme");

    // add our fields and their ranges to central list
    Calibrator::IAdjustable::registerBootstrappableField(
            clazz, "locVolParamsA",
            new Range(Infinity(Infinity::Minus), Infinity(Infinity::Plus)),
            getLocVolExpiries);
    Calibrator::IAdjustable::registerBootstrappableField(
            clazz, "locVolParamsB",
            new Range(Infinity(Infinity::Minus), Infinity(Infinity::Plus)),
            getLocVolExpiries);
    Calibrator::IAdjustable::registerBootstrappableField(
            clazz, "locVolParamsC",
            new Range(ClosedBoundary(0), Infinity(Infinity::Plus)),
            getLocVolExpiries);
}

CClassConstSP const VolSVCJLV::TYPE =
CClass::registerClassLoadMethod("VolSVCJLV", typeid(VolSVCJLV), load);

template<>
const CClassConstSP VolSVCJConvert::TYPE = 
CClass::registerInterfaceLoadMethod("MarketDataConvert::IConvert<VolSVCJ>", typeid(VolSVCJConvert), load);

/* external symbol to allow class to be forced to be linked in */
bool VolSVCJLVLinkIn(){
    return (VolSVCJLV::TYPE != 0) && (VolSVCJConvert::TYPE != 0);
}

DRLIB_END_NAMESPACE
