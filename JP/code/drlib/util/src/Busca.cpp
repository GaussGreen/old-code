//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Busca.cpp
//
//   Description : Implementation of Busca equation given VHT local vol
//
//   Date        : 20 December 2004
//
//
//   $Log:$
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/PDESolverMoL.hpp"
#include "edginc/Format.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Function.hpp"
#include "edginc/Integrator.hpp"
#include "edginc/Addin.hpp"


DRLIB_BEGIN_NAMESPACE

#define xxxSECONDWAY

class Busca: public PDESolverMoL::PDE {    
public:
    static CClassConstSP const TYPE;

    class ILocalVol:virtual public IObject {
    public:
        static CClassConstSP const TYPE;
        ILocalVol():IObject(){}
        ILocalVol(smartConstPtr<ILocalVol> sigma);
        virtual double getValue(double    x,
                                double    t) const = 0;
        virtual double getFirstDer(double    x,
                                   double    t) const = 0;
        virtual double getSecondDer(double    x,
                                    double    t) const = 0;
        virtual double getThirdDer(double    x,
                                   double    t) const = 0;
        virtual double getValueAtSpacePositiveInfinity(double t) const = 0;
        virtual double getValueAtSpaceNegativeInfinity(double t) const = 0;
    private:
        static void load(CClassSP& clazz){
            clazz->setPublic(); // make visible to EAS/spreadsheet
            REGISTER_INTERFACE(ILocalVol, clazz);
            EXTENDS(IObject);
        }
    };

    typedef smartConstPtr<ILocalVol> ILocalVolConstSP;
    typedef smartPtr<ILocalVol> ILocalVolSP;


    class BuscaIntegrand:public Function1DDouble {
    public:
        static CClassConstSP const TYPE;

        // constructor
        BuscaIntegrand(int                  type,
                       double               x,
                       double               t,
                       ILocalVolConstSP     sigma):
        Function1DDouble(Range(ClosedBoundary(0.0), ClosedBoundary(1.0))), 
        integrationType(type),
        t(t),
        x(x),
        sigma(sigma){}

        //function
        virtual double operator()(double  s) const {
            double value = 0;

            switch(integrationType)
            {
                case 0:
                    value = 1.0/sigma->getValue(s*x,0.0);
                    break;

                case 1:
                    value = Maths::square(sigma->getValueAtSpacePositiveInfinity(s*t));
                    break;

                case -1:
                    value = Maths::square(sigma->getValueAtSpaceNegativeInfinity(s*t));
                    break;

                case 2:
                    value = 1.0/sigma->getValue(s*x,t);
                    break;
            }

            return value;
        }

    protected:
        int                 integrationType;   //0 for intial, 1 or positive infinity, -1 for negative infinity,
        double              x;
        double              t;
        ILocalVolConstSP    sigma;
    };
        
    /*  Given the value of the k-th function u_k and of its 1st order
        space derivatives ux_k at the node (x,t), evaluate the quantity 
        P,Q and Q defined by 
                    Sum_j P_i_j*du_j/dt + Q_i = dR_i/dx
        where i,j = 0,..., npdes-1 */
    void func(double                x,
              double                t,
              const CDoubleArray&   u,
              const CDoubleArray&   ux,
              CDoubleMatrix&        P,
              CDoubleArray&         Q,
              CDoubleArray&         R) const {
        static const string method = "Busca::func";
        int n = u.size();
        try{
            if (n != 1){
                throw ModelException(method,
                                    "Busca PDE must be dimension 1; not "
                                    + Format::toString(n));
            }
        }
        catch (exception& e){
            throw ModelException(e,method);
        }

        double localVolSqr = Maths::square(sigma->getValue(x,t));

#ifdef SECONDWAY
        P[0][0] = 2*t/localVolSqr;
        Q[0] = u[0]/localVolSqr - Maths::square(1-x*ux[0]/u[0])/u[0] + Maths::square(0.5*t*u[0]*ux[0])/u[0]*weight2;
        R[0] = t*ux[0]*weight;
#else
        P[0][0] = 2*t*u[0]/localVolSqr;
        Q[0] = Maths::square(u[0])/localVolSqr - Maths::square(1-x*ux[0]/u[0]) + t*Maths::square(ux[0])*weight + Maths::square(0.5*t*u[0]*ux[0])*weight2;
        R[0] = t*u[0]*ux[0]*weight;
#endif
    }

   
    /** Evaluate the coefficients of the lower boundary condition at time t. 
        The boundary condition has to be of the form
                        beta * R = gamma         */
    void lowerBound(double              t,
                    const CDoubleArray& u,
                    const CDoubleArray& ux,
                    CDoubleArray&       beta,
                    CDoubleArray&       gamma) const {
        // set finite boundary value at infinity, otherwise set R zero.
        if (boundOption == 0){
            beta[0] = 0.0;
            gamma[0] = u[0] - solutionAtSpaceNegativeInfinity(t);
        } else {
            beta[0] = 1.0;
            gamma[0] = 0.0;
        }
    }


    /** Evaluate the coefficients of the upper boundary condition at time t. 
        The boundary condition has to be of the form
                        beta * R = gamma         */
    void upperBound(double              t,
                    const CDoubleArray& u,
                    const CDoubleArray& ux,
                    CDoubleArray&       beta,
                    CDoubleArray&       gamma) const {
        // set finite boundary value at infinity, otherwise set R zero.
        if (boundOption == 0){
            beta[0] = 0.0;
            gamma[0] = u[0] - solutionAtSpacePositiveInfinity(t);
        } else {
            beta[0] = 1.0;
            gamma[0] = 0.0;
        }
    }
        
    /* Given x, returns the initial values u_k */
    void init(double        x,
              DoubleArray&  u) const {

        u[0] = solutionAtTimeZero(x);

        // to be consistent with boundary condition
        if (boundOption == 0){
            if (Maths::isZero(x - lowerBoundary)){
                u[0] = solutionAtSpaceNegativeInfinity(0.0);
            }
            if (Maths::isZero(x - upperBoundary)){
                u[0] = solutionAtSpacePositiveInfinity(0.0);
            }
        }
    }

    /* Asympototics of Busca equation*/
    //solution at t=0
    virtual double solutionAtTimeZero(double x) const{

        BuscaIntegrand  myIntegrand(0,x,0.0,sigma);
        string          integrationMethod = "ClosedRomberg1D";
        double          absPrecision = 0.0001;
        double          relPrecision = 0.0001;

        // create an integrator by name and pass on precision in terms of integrator
        Integrator1DSP myIntegrator = Integrator1D::createIntegrator(integrationMethod, 
                                                                    absPrecision,
                                                                    relPrecision);

        return 1.0/myIntegrator->integrate(myIntegrand);
    }

    //solution at x=+infinity
    virtual double solutionAtSpacePositiveInfinity(double t) const {

        BuscaIntegrand  myIntegrand(1,1000,t,sigma);
        string          integrationMethod = "ClosedRomberg1D";
        double          absPrecision = 0.0001;
        double          relPrecision = 0.0001;

        // create an integrator by name and pass on precision in terms of integrator
        Integrator1DSP myIntegrator = Integrator1D::createIntegrator(integrationMethod, 
                                                                    absPrecision,
                                                                    relPrecision);

        return sqrt(myIntegrator->integrate(myIntegrand));
    }

    //solution at x=-infinity
    virtual double solutionAtSpaceNegativeInfinity(double t) const {
        
        BuscaIntegrand  myIntegrand(-1,1000.0,t,sigma);
        string          integrationMethod = "ClosedRomberg1D";
        double          absPrecision = 0.0001;
        double          relPrecision = 0.0001;

        // create an integrator by name and pass on precision in terms of integrator
        Integrator1DSP myIntegrator = Integrator1D::createIntegrator(integrationMethod, 
                                                                    absPrecision,
                                                                    relPrecision);

        return sqrt(myIntegrator->integrate(myIntegrand));
    }

    // similar to Busca initial solution, but taking local vol at time t
    virtual double getIntegratedVolAtTime(double x,
                                          double t) const{

        BuscaIntegrand  myIntegrand(2,x,t,sigma);
        string          integrationMethod = "ClosedRomberg1D";
        double          absPrecision = 0.0001;
        double          relPrecision = 0.0001;

        // create an integrator by name and pass on precision in terms of integrator
        Integrator1DSP myIntegrator = Integrator1D::createIntegrator(integrationMethod, 
                                                                    absPrecision,
                                                                    relPrecision);

        return 1.0/myIntegrator->integrate(myIntegrand);
    }

    //analytical approximations
    double localExpansionApprox(double   x,
                                double   t);

    double globalIterationApprox(double   x,
                                 double   t);

    CDoubleArray globalIterationApprox(double           x,
                                       CDoubleArray     t,
                                       int              method);

    double GatheralFirstApprox(double   x,
                               double   t);

    double SolutionOfZeroOrderBuscaWithSeparableVHT(double x, 
                                                    double t);

    // approximation: set finite lower and upper boundaries. 
    void setLowerBound(double x) {
        lowerBoundary = x;
    }

    void setUpperBound(double x) {
        upperBoundary = x;
    }

    void setBoundOption(int i) {
        boundOption = i;
    }

private:
    Busca():PDE(TYPE),boundOption(1), weight(1.0), weight2(1.0){}

    static IObject* defaultBusca(){
        return new Busca(); 
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Busca, clazz);
        SUPERCLASS(PDESolverMoL::PDE);
        EMPTY_SHELL_METHOD(defaultBusca);
        FIELD(sigma, "local vol");
        FIELD_INLINE(lowerBoundary, "lower boundary");
        FIELD_INLINE(upperBoundary, "upper boundary");
        FIELD_INLINE(boundOption, "how to specify boundary condition");
        FIELD_MAKE_OPTIONAL(boundOption);
        FIELD_INLINE(weight, "weight of first order term of time on the RHS of Busca equation");
        FIELD_MAKE_OPTIONAL(weight);
        FIELD_INLINE(weight2, "weight of second order term of time on the RHS of Busca equation");
        FIELD_MAKE_OPTIONAL(weight2);
    }

protected:
    ILocalVolConstSP    sigma;
    double              lowerBoundary;  // replace infinity boundary by finite ones
    double              upperBoundary;
    int                 boundOption;    //0 specifying u; otherwise specifying R
    double              weight;         //weight of first order term (in time) on the RHS of busca equation
    double              weight2;        //weight of second order term (in time) on the RHS of busca equation
};

typedef smartConstPtr<Busca> BuscaConstSP;
typedef smartPtr<Busca> BuscaSP;

CClassConstSP const Busca::TYPE = CClass::registerClassLoadMethod(
                                    "Busca", typeid(Busca), Busca::load);

CClassConstSP const Busca::ILocalVol::TYPE = CClass::registerInterfaceLoadMethod("Busca::ILocalVol",
                                                                                 typeid(Busca::ILocalVol), 
                                                                                 Busca::ILocalVol::load);

// * for class loading (avoid having header file) */
bool BuscaLoad() {
    return (Busca::TYPE != 0);
}


class VHT: public CObject,
           virtual public Busca::ILocalVol{
public:
    friend class SeparableVHT;
    static CClassConstSP const TYPE;

    double getValue(double    x,
                    double    t) const {
        double  C = getParamAtTime(paramC,t),
                D = getParamAtTime(paramD,t),
                sigma0 = getParamAtTime(atmVol,t),
                vol = 1.0 + paramA * tanh(C*x) + paramB * (1.0-1.0/cosh(D*x));
        return sigma0 * sqrt(vol);
    }

    double getFirstDer(double    x,
                       double    t) const {
        double  C = getParamAtTime(paramC,t),
                D = getParamAtTime(paramD,t),
                sigma0 = getParamAtTime(atmVol,t),
                value = getValue(x,t),
                sechcx = 1.0/cosh(C*x),
                term = paramA*C*sechcx*sechcx + paramB*D/cosh(D*x)*tanh(D*x);
        return 0.5*sigma0*sigma0/value*term;
    }

    double getSecondDer(double    x,
                        double    t) const {
        double  C = getParamAtTime(paramC,t),
                D = getParamAtTime(paramD,t),
                sigma0 = getParamAtTime(atmVol,t),
                value = getValue(x,t),
                sechcx = 1.0/cosh(C*x),
                sechdx = 1.0/cosh(D*x),
                tanhcx = tanh(C*x),
                tanhdx = tanh(D*x),
                term1 = 2.0*paramA*C*C*sechcx*sechcx*tanhcx+paramB*D*D*sechdx*(tanhdx*tanhdx-sechdx*sechdx),
                term2 = paramA*C*sechcx*sechcx+paramB*D*sechdx*tanhdx;
        return -(0.5*term1+0.25*sigma0*sigma0/value/value*term2*term2)*sigma0*sigma0/value;
    }

    double getThirdDer(double    x,
                       double    t) const {
        double  C = getParamAtTime(paramC,t),
                D = getParamAtTime(paramD,t),
                sigma0 = getParamAtTime(atmVol,t),
                value = getValue(x,t),
                sechcx = 1.0/cosh(C*x),
                sechdx = 1.0/cosh(D*x),
                tanhcx = tanh(C*x),
                tanhdx = tanh(D*x),
                term1,term2,term3, term4;

        term1 = paramA*C*sechcx*sechcx+paramB*D*sechdx*tanhdx;
        term2 = term1*(-2*paramA*C*C*sechcx*sechcx*tanhcx+paramB*D*D*sechdx*(sechdx*sechdx-tanhdx*tanhdx));
        term3 = 2.0*paramA*C*C*C*sechcx*sechcx*(-sechcx*sechcx+2.0*tanhcx*tanhcx)
                    +paramB*D*D*D*sechdx*tanhdx*(-5.0*sechdx*sechdx+tanhdx*tanhdx);
        term4 = sigma0*sigma0/value;

        term1 *= term4;
        term1 = 0.375*term1*term1*term1/value/value;

        term2 *= -0.75*term4*term4/value;
        term3 *= 0.5*term4;
        return term1+term2+term3;
    }

    double getValueAtSpacePositiveInfinity(double t) const {
        double  vol     = 1.0 + paramA + paramB,
                sigma0  = getParamAtTime(atmVol,t);
        return sigma0 * sqrt(vol);
    }

    double getValueAtSpaceNegativeInfinity(double t) const {
        double  vol     = 1.0 - paramA + paramB,
                sigma0  = getParamAtTime(atmVol,t);
        return sigma0 * sqrt(vol);
    }

    int getTenorIndexAtTime(double  t) const {
        static const string method = "VHT::getTenorIndexAtTime";

        if (Maths::isNegative(t)){
            throw ModelException(method, 
                             "time can not be negative.");
        }

        int     i = 0,
                n = tenor.size();

        while (Maths::isPositive(t-tenor[i])){
            i++;
            if (i==n) {
                return i-1;
            }
        }

        return i;
    }

    double getParamAtTime(CDoubleArray param,
                          double t) const {
        return param[getTenorIndexAtTime(t)];
    }

    CDoubleArray getTenor() const {
        return tenor;
    }

    void validatePop2Object(){
        static const string method = "VHT::validatePop2Object";
        
        int sizeAtmVol = atmVol.size(),
            sizeC = paramC.size(),
            sizeD = paramD.size(),
            size = tenor.size();

        Maths::checkPositive(size, "size of tenor");

        if (sizeAtmVol!=size){
            throw ModelException(method, 
                             "size of atmVol "+ Format::toString(sizeAtmVol) + 
                             "is not equal to size of tenor "+ Format::toString(size));
        }

        if (sizeC!=size){
            throw ModelException(method, 
                             "size of paramC "+ Format::toString(sizeC) + 
                             "is not equal to size of tenor "+ Format::toString(size));
        }

        if (sizeD!=size){
            throw ModelException(method, 
                             "size of paramD "+ Format::toString(sizeD) + 
                             "is not equal to size of tenor "+ Format::toString(size));
        }

        if (!Maths::isPositive(tenor[0])){
            throw ModelException(method, 
                             "tenor should start with positive number.");
        }

        for (int i=1; i<size; i++){
            if (!Maths::isPositive(tenor[i]-tenor[i-1])){
                throw ModelException(method, 
                                 "tenor should be in an ascending order.");
            }
        }
    }

    double GatheralFirstApprox(double x,
                            double t) const {
        if (Maths::isZero(t)) {
            return getValue(x,0.0);
        }

        const double largeNumber = 100.0;
        const double smallNumber = 1.0e-10;

        int     i, finalIndex = getTenorIndexAtTime(t)+1;
        double  tdiff, sum = 0.0;

        DoubleArray tempTenor(finalIndex+1), xbkpts(finalIndex+1);
        for (i=0;i<finalIndex-1;i++){
            tempTenor[i+1] = tenor[i];
        }
        tempTenor[0] = 0.0;
        tempTenor[finalIndex] = t;

        for (i=0;i<finalIndex+1;i++){
            xbkpts[i] = tempTenor[i]*x/t;
        }
        
        // atm term 
        for (i=0; i<finalIndex; i++){
            tdiff = tempTenor[i+1] - tempTenor[i];
            sum += atmVol[i]*atmVol[i]*tdiff;
        }
        if ((x>smallNumber)||(x<-smallNumber)){
            // skew and convexity terms
            double cx, logCosh1, logCosh2, dx, atanExp1, atanExp2;
            for (i=0; i<finalIndex; i++){  
                logCosh1 = logCosh(paramC[i]*xbkpts[i+1]);
                logCosh2 = logCosh(paramC[i]*xbkpts[i]);
                cx = paramC[i]*(xbkpts[i+1]-xbkpts[i]);

                atanExp1 = atanExp(paramD[i]*xbkpts[i+1]);
                atanExp2 = atanExp(paramD[i]*xbkpts[i]);
                dx = paramD[i]*(xbkpts[i+1]-xbkpts[i]);

                tdiff = tempTenor[i+1] - tempTenor[i];
                sum += tdiff*atmVol[i]*atmVol[i]*(paramA*(logCosh1-logCosh2)/cx + paramB*(1.0-2.0*(atanExp1-atanExp2)/dx));
            }
        }
        return sqrt(sum/t);
    }

    VHT(CDoubleArray    atmVol,
        double          paramA,
        double          paramB,
        CDoubleArray    paramC,
        CDoubleArray    paramD,
        CDoubleArray    tenor):
        CObject(TYPE),
        atmVol(atmVol),
        paramA(paramA),
        paramB(paramB),
        paramC(paramC),
        paramD(paramD),
        tenor(tenor){}
                       
protected:
    /** VHT local volatility specified by parameters atmVol, A, B, C and D,
        where atmVol, C and D are step functions of t. */
    CDoubleArray    atmVol;
    double          paramA;
    double          paramB;
    CDoubleArray    paramC;
    CDoubleArray    paramD;
    CDoubleArray    tenor;     

private:
    VHT():CObject(TYPE){}
    static IObject* defaultVHT(){
        return new VHT();
    }

    /**Invoke when Class is 'load' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(VHT,clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultVHT);
        IMPLEMENTS(Busca::ILocalVol);
        FIELD_INLINE(atmVol, "at the money vols");
        FIELD_INLINE(paramA, "paramA of VHT local vol");
        FIELD_INLINE(paramB, "paramB of VHT local vol");
        FIELD_INLINE(paramC, "paramC of VHT local vol");
        FIELD_INLINE(paramD, "paramD of VHT local vol");
        FIELD_INLINE(tenor, "tenor of VHT local vol");
    }

    double logCosh(double x) const {
        double val = 0.0, largeNumber = 100.0;
        
        if (x>largeNumber){
            val = x - log(2.0);
        } else {
            if (x<-largeNumber) {
                val = -x - log(2.0);
            } else {
                val = log(cosh(x));   
            }
        }
        return val;
    }

    double atanExp(double x) const {
        double val = 0.0, largeNumber = 100.0;
        if(x>largeNumber){
            val = 0.5*Maths::PI;
        } else {
            val = atan(exp(x));
        }
        return val;
    }

};

typedef smartConstPtr<VHT> VHTConstSP;
typedef smartPtr<VHT> VHTSP;

CClassConstSP const VHT::TYPE = CClass::registerClassLoadMethod(
    "VHT", typeid(VHT), VHT::load);

// * for class loading (avoid having header file) */
bool VHTLoad() {
    return (VHT::TYPE != 0);
}


class SeparableVHT: public CObject,
                    virtual public Busca::ILocalVol{
public:
    static CClassConstSP const TYPE;

     void validatePop2Object(){
        static const string method = "SeparableVHT::validatePop2Object";
        
        int sizeAtmVol = atmVol.size(),
            size = tenor.size();

        Maths::checkPositive(size, "size of tenor");

        if (sizeAtmVol!=size){
            throw ModelException(method, 
                             "size of atmVol "+ Format::toString(sizeAtmVol) + 
                             "is not equal to size of tenor "+ Format::toString(size));
        }

        if (!Maths::isPositive(tenor[0])){
            throw ModelException(method, 
                             "tenor should start with positive number.");
        }

        for (int i=1; i<size; i++){
            if (!Maths::isPositive(tenor[i]-tenor[i-1])){
                throw ModelException(method, 
                                 "tenor should be in an ascending order.");
            }
        }

        // set transient VHTSP
        int             n = tenor.size();
        DoubleArray     arrayC(n,paramC), arrayD(n,paramD);
        vhtPtr = VHTSP(new VHT(atmVol, paramA, paramB, arrayC, arrayD, tenor));
    }

    double getValue(double    x,
                    double    t) const {
        return vhtPtr->getValue(x,t);
    }
    
    double getFirstDer(double    x,
                       double    t) const {
        return vhtPtr->getFirstDer(x,t);
    }

    double getSecondDer(double    x,
                        double    t) const {
        return vhtPtr->getSecondDer(x,t);
    }

    double getThirdDer(double    x,
                        double    t) const {
        return vhtPtr->getThirdDer(x,t);
    }

    double getValueAtSpacePositiveInfinity(double t) const {
        return vhtPtr->getValueAtSpacePositiveInfinity(t);
    }

    double getValueAtSpaceNegativeInfinity(double t) const {
        return vhtPtr->getValueAtSpaceNegativeInfinity(t);
    }

    double GatheralFirstApprox(double x,
                                double t) const {
        return vhtPtr->GatheralFirstApprox(x,t);
    }

    double getAtmVol(double t) const {
        double vol = vhtPtr->getParamAtTime(atmVol,t);
        return vol;
    }

    VHTSP getVHTSP() const {
        return vhtPtr;
    }

    void setAtmVol(double vol) {
        for (int i=0; i<atmVol.size(); i++){
            atmVol[i] = vol;
        }
    }

    void setAtmVol(CDoubleArray vol) {
        atmVol = vol;
    }

    double integrateTimeComponent(double t){

        //time part
        int     finalIndex = vhtPtr->getTenorIndexAtTime(t) + 1;

        DoubleArray tempTenor(finalIndex+1);
        int i;
        for (i=0;i<finalIndex-1;i++){
            tempTenor[i+1] = tenor[i];
        }
        tempTenor[0] = 0.0;
        tempTenor[finalIndex] = t;

        double sum = 0.0;
        for (i=0;i<finalIndex;i++){
            sum += atmVol[i]*atmVol[i]*(tempTenor[i+1]-tempTenor[i]);
        }

        return sqrt(sum/t);
    }

protected:
    /** VHT local volatility specified by parameters atmVol, A, B, C and D,
        where atmVol is step functions of t. */
    CDoubleArray    atmVol;
    double          paramA;
    double          paramB;
    double          paramC;
    double          paramD;
    CDoubleArray    tenor; 
    
    //transient
    VHTSP           vhtPtr;

private:
    SeparableVHT():CObject(TYPE){}
    static IObject* defaultSeparableVHT(){
        return new SeparableVHT();
    }

    /**Invoke when Class is 'load' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(SeparableVHT,clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultSeparableVHT);
        IMPLEMENTS(Busca::ILocalVol);
        FIELD_INLINE(atmVol, "at the money vols");
        FIELD_INLINE(paramA, "paramA of VHT local vol");
        FIELD_INLINE(paramB, "paramB of VHT local vol");
        FIELD_INLINE(paramC, "paramC of VHT local vol");
        FIELD_INLINE(paramD, "paramD of VHT local vol");
        FIELD_INLINE(tenor, "tenor of VHT local vol");
        FIELD(vhtPtr,"VHT pointer");
        FIELD_MAKE_TRANSIENT(vhtPtr);
    }        
};

typedef smartConstPtr<SeparableVHT> SeparableVHTConstSP;
typedef smartPtr<SeparableVHT> SeparableVHTSP;

CClassConstSP const SeparableVHT::TYPE = CClass::registerClassLoadMethod(
    "SeparableVHT", typeid(SeparableVHT), SeparableVHT::load);

// * for class loading (avoid having header file) */
bool SeparableVHTLoad() {
    return (SeparableVHT::TYPE != 0);
}

//method based on locally expansing total variance
double Busca::localExpansionApprox(double   x,
                                   double   t){
    static const string method = "Busca::localExpansionApprox";
    try{
        if (!(VHT::TYPE->isInstance(sigma)||SeparableVHT::TYPE->isInstance(sigma))){
            throw ModelException(method, 
                    "analytical approximation is not available for local Vol "
                    + sigma->getClass()->getName());
        }
        
        double              delta = 0.005;
        double              epsilon = 0.05;        // for approximation around x = 0
        double              temp, xtemp, tdiff;

        const VHT*      localVol = NULL;
        if (VHT::TYPE->isInstance(sigma)){
            localVol = dynamic_cast<const VHT*>(sigma.get());
        } 
        if (SeparableVHT::TYPE->isInstance(sigma)){
             localVol = (dynamic_cast<const SeparableVHT*>(sigma.get()))->getVHTSP().get();
        }
        CDoubleArray        tenor = localVol->getTenor();

        int                 n = 3;
        int                 tIndex = localVol->getTenorIndexAtTime(t);
        int                 nx = 2*tIndex+1;
        int                 i;
        CDoubleArray        xarray(nx);
        CDoubleMatrix       prev(nx,n);
        CDoubleMatrix       curr(nx,n);
       
        double  sigma0 = localVol->getValue(0.0,0.0),
                firstDer = localVol->getFirstDer(0.0,0.0),
                secondDer = localVol->getSecondDer(0.0,0.0),
                thirdDer = localVol->getThirdDer(0.0,0.0),
                zeroOrder = sigma0*(secondDer*sigma0-0.5*firstDer*firstDer)/12.0,
                firstOrder = sigma0*sigma0*thirdDer/24.0+0.5*firstDer/sigma0*zeroOrder;

        for (int i=0; i<nx; i++){
            xarray[i] = x+(i-tIndex)*delta;
            curr[i][0] = solutionAtTimeZero(xarray[i]);

            xtemp = xarray[i];
            temp = curr[i][0];
            if (Maths::isPositive(xtemp+epsilon)&&Maths::isNegative(xtemp-epsilon)){
                curr[i][1] = zeroOrder+firstOrder*xarray[i];
            } else {
                curr[i][1] = temp*temp*temp/xtemp/xtemp*log(sqrt(localVol->getValue(xtemp,0.0)*sigma0)/temp);
            }

            curr[i][2] = 0.0;
        }

        if(t<tenor[1]){
            tdiff = t-tenor[0];
            return curr[tIndex][0]+curr[tIndex][1]*tdiff+curr[tIndex][2]*tdiff*tdiff;
        } else {
            for (int i=0; i<nx; i++){
/*                double temp1, temp2;
                xtemp = xarray[i];
                tdiff = tenor[1]-tenor[0];
                temp1 = curr[i][0]+curr[i][1]*tdiff+curr[i][2]*tdiff*tdiff;
                temp2 = temp1*temp1+2.0*temp1*(curr[i][1]+2.0*curr[i][2]*tdiff)*tdiff;
                curr[i][0] = temp1*temp1*tdiff;
                temp = localVol->getValue(xarray[i],tenor[1])/localVol->getValue(xarray[i],tenor[0]);
                curr[i][1] = temp*temp*temp2;
*/
                curr[i][2] = 2.0*curr[i][0]*curr[i][1];
                curr[i][1] = curr[i][0]*curr[i][0];
                curr[i][0] = 0.0;
            }
        }


        for (int k=1; k<=tIndex; k++){
            double ratio0, ratio1, secondDer;
            prev = curr;

            tdiff = tenor[k]-tenor[k-1];

            for (i=k-1;i<nx-k+1;i++){
                curr[i][0] = prev[i][0] + prev[i][1]*tdiff + prev[i][2]*tdiff*tdiff;
                temp = localVol->getValue(xarray[i],tenor[k])/localVol->getValue(xarray[i],tenor[k-1]);
                curr[i][1] = temp*temp*(prev[i][1]+2.0*prev[i][2]);
            }

            for (i=k;i<nx-k;i++){
                ratio0 = (curr[k+1][0]-curr[k-1][0])/2./delta/curr[k][0];
                ratio1 = (curr[k+1][1]-curr[k-1][1])/2./delta/curr[k][1];
                secondDer = (curr[k+1][1]+curr[k-1][1]-2.*curr[k][1])/delta/delta;
                temp = localVol->getValue(xarray[i],tenor[k]);

                curr[i][2] = -xarray[i]*(1-0.5*xarray[i]*ratio0)*curr[i][1]/curr[i][0]*(ratio1-ratio0);
                curr[i][2] += 0.5*secondDer-0.25*curr[i][1]*ratio0*(2.*ratio1-ratio0);
                curr[i][2] *= 0.5*temp*temp;
            }
        }
        tdiff = t - tenor[tIndex];
        return sqrt((curr[tIndex][0]+curr[tIndex][1]*tdiff+curr[tIndex][2]*tdiff*tdiff)/t);
    }
    catch (exception& e){
        throw ModelException(e,method);
    }
}

CDoubleArray Busca::globalIterationApprox(double            x,
                                          CDoubleArray      t,
                                          int               methodIdx){
    static const string method = "Busca::globalIterationApprox";
    try{
        if (!(VHT::TYPE->isInstance(sigma)||SeparableVHT::TYPE->isInstance(sigma))){
            throw ModelException(method, 
                    "analytical approximation is not available for local Vol "
                    + sigma->getClass()->getName());
        }

        if (methodIdx>2||methodIdx<0){
            throw ModelException(method, 
                    "methodIdx must be 0, 1, or 2");
        }

        const VHT*      localVol = NULL;
        if (VHT::TYPE->isInstance(sigma)){
            localVol = dynamic_cast<const VHT*>(sigma.get());
        } 
        if (SeparableVHT::TYPE->isInstance(sigma)){
             localVol = (dynamic_cast<const SeparableVHT*>(sigma.get()))->getVHTSP().get();
        }

        int             tSize = t.size();
        int             tIndex = localVol->getTenorIndexAtTime(t[tSize-1]);
        CDoubleArray    u(tSize), totalVar(tIndex+1);
        CDoubleArray    tenor = localVol->getTenor(); 

        //for convenience, define temporary time arrays starting with zero
        CDoubleArray    tmpTenor(tIndex+2);
        int             i,k;
        for (k=1; k<tIndex+2; k++){
            tmpTenor[k] = tenor[k-1];
        }
        
        //vo
        double vol, tdiff;
        for (k=0; k<tIndex; k++){
            vol = getIntegratedVolAtTime(x,tenor[k]);
            tdiff = tmpTenor[k+1]-tmpTenor[k];
            totalVar[k+1] = totalVar[k] + vol*vol*tdiff;           
        }
        for (int i=0; i<tSize; i++){
            k = localVol->getTenorIndexAtTime(t[i]);
            vol = getIntegratedVolAtTime(x,tenor[k]);
            tdiff = t[i]-tmpTenor[k];
            u[i] += totalVar[k] + vol*vol*tdiff; 
        }

        //v1
        if (methodIdx>0){
            CDoubleArray    V1(tIndex+1), asum(tIndex+1), bsum(tIndex+1);
            double  aterm, bterm, a, b, c, tdiffpre;
            k = 0;
            vol = getIntegratedVolAtTime(x,tenor[k]);
            aterm = vol*vol*vol/localVol->getValue(x,tenor[k]);
            bterm = vol*vol;
            for (k=1; k<tIndex; k++){
                tdiffpre = tmpTenor[k]-tmpTenor[k-1];
                asum[k] = asum[k-1] + aterm*tdiffpre;
                bsum[k] = bsum[k-1] + bterm*tdiffpre;

                vol = getIntegratedVolAtTime(x,tenor[k]);
                aterm = vol*vol*vol/localVol->getValue(x,tenor[k]);
                bterm = vol*vol;

                a = aterm/asum[k];
                b = bterm/bsum[k];
                c = asum[k]/bsum[k];
    
                tdiff = tmpTenor[k+1]-tmpTenor[k];
                V1[k+1] = V1[k] + vol*vol*((b-a)*(b-a)/a/a*tdiff/(1.0+b*tdiff)+2.0*(b-a)/a/b*log(1.0+b*tdiff));           
            }
            for (i=0; i<tSize; i++){
                k = localVol->getTenorIndexAtTime(t[i]);
                if (k>0){
                    tdiffpre = tmpTenor[k]-tmpTenor[k-1];
                    double asumtmp = asum[k] + aterm*tdiffpre;
                    double bsumtmp = bsum[k] + bterm*tdiffpre;

                    vol = getIntegratedVolAtTime(x,tenor[k]);
                    aterm = vol*vol*vol/localVol->getValue(x,tenor[k]);
                    bterm = vol*vol;

                    a = aterm/asumtmp;
                    b = bterm/bsumtmp;
                    c = asumtmp/bsumtmp;

                    tdiff = t[i]-tenor[k-1];
                    u[i] += V1[k] + vol*vol*((b-a)*(b-a)/a/a*tdiff/(1.0+b*tdiff)+2.0*(b-a)/a/b*log(1.0+b*tdiff));    
                }
            }
        }
        
        //v2
        if (methodIdx>1){
            CDoubleArray    V2(tIndex+1,0.0);
            double  temp, local;
            double  firstDer, secondDer, thirdDer, zeroOrder, firstOrder;
            double  delta = 0.05;
            for (k=0; k<tIndex; k++){
                vol = getIntegratedVolAtTime(x,tenor[k]);
                local = localVol->getValue(x,tenor[k]);
                firstDer = localVol->getFirstDer(x,tenor[k]);
                tdiff = tmpTenor[k+1]*tmpTenor[k+1]-tmpTenor[k]*tmpTenor[k];
                if ((x<delta)&&(x>-delta)){
                    firstDer = localVol->getFirstDer(0.0,tenor[k]),
                    secondDer = localVol->getSecondDer(0.0,tenor[k]),
                    thirdDer = localVol->getThirdDer(0.0,tenor[k]),
                    zeroOrder = (secondDer-0.5*firstDer*firstDer/local)/3.0,
                    firstOrder = 0.25*thirdDer-0.5*firstDer/local*zeroOrder;
                    temp = zeroOrder + firstOrder * x;
                    V2[k+1] = V2[k] + 0.5*vol*vol*vol*temp*tdiff;
                } else {
                    V2[k+1] = V2[k] + 0.5*vol*vol*vol/x/x*(x*firstDer-2.0*local+2.0*vol)*tdiff; 
                }
            }
            for (i=0; i<tSize; i++){
                k = localVol->getTenorIndexAtTime(t[i]);
                vol = getIntegratedVolAtTime(x,tenor[k]);
                local = localVol->getValue(x,tenor[k]);
                firstDer = localVol->getFirstDer(x,tenor[k]);
                tdiff = t[i]*t[i]-tmpTenor[k]*tmpTenor[k];
                if ((x<delta)&&(x>-delta)){
                    firstDer = localVol->getFirstDer(0.0,tenor[k]),
                    secondDer = localVol->getSecondDer(0.0,tenor[k]),
                    thirdDer = localVol->getThirdDer(0.0,tenor[k]),
                    zeroOrder = (secondDer-0.5*firstDer*firstDer/local)/3.0,
                    firstOrder = 0.25*thirdDer-0.5*firstDer/local*zeroOrder;
                    temp = zeroOrder + firstOrder * x;
                    u[i] += V2[k] + 0.5*vol*vol*vol*temp*tdiff;
                } else {
                    u[i] += V2[k] + 0.5*vol*vol*vol/x/x*(x*firstDer-2.0*local+2.0*vol)*tdiff; 
                }
            }
        }

        for (i=0; i<tSize; i++){
            u[i] = sqrt(u[i]/t[i]);
        }
        return u;
    }
    catch (exception& e){
        throw ModelException(e,method);
    }
}

//method based on Gatheral approximation
double Busca::GatheralFirstApprox(double   x,
                                  double   t){
    static const string method = "Busca::GatheralFirstApprox";
    try{
        if (!(VHT::TYPE->isInstance(sigma)||SeparableVHT::TYPE->isInstance(sigma))){
            throw ModelException(method, 
                    "analytical approximation is not available for local Vol "
                    + sigma->getClass()->getName());
        }

        const VHT*      localVol = NULL;
        if (VHT::TYPE->isInstance(sigma)){
            localVol = dynamic_cast<const VHT*>(sigma.get());
        } 
        if (SeparableVHT::TYPE->isInstance(sigma)){
             localVol = (dynamic_cast<const SeparableVHT*>(sigma.get()))->getVHTSP().get();
        }

        return localVol->GatheralFirstApprox(x,t); 
    }
    catch (exception& e){
        throw ModelException(e,method);
    }
}

//axact solution of first order Busca with separable VHT
double Busca::SolutionOfZeroOrderBuscaWithSeparableVHT(double x, 
                                                       double t){
    static const string method = "Busca::SolutionOfZeroOrderBuscaWithSeparableVHT";
    try{
        if (!(SeparableVHT::TYPE->isInstance(sigma.get()))){
            throw ModelException(method,
                "analytical solution is not available for local vol "
                + sigma->getClass()->getName());
        }

        const SeparableVHT*  localVol = dynamic_cast<const SeparableVHT*>(sigma.get());

        SeparableVHT         sigma = static_cast<SeparableVHT>(*localVol);
        
        //spatial component
        double sol = getIntegratedVolAtTime(x,t);
        sol /= sigma.getAtmVol(t);

        //time component
        sol *= sigma.integrateTimeComponent(t);

        return sol;
    }      
    catch (exception& e){
        throw ModelException(e,method);
    }
}


class BuscaApproxAddin:public CObject{
public:
    static CClassConstSP const TYPE;

    /**the 'addin function' */ 
    static IObjectSP BuscaApprox(BuscaApproxAddin* params){
   
        int             numCols = params->xbkpts.size();
        int             numRows = params->tend.size();
        DoubleMatrix    ubkpts(numCols, numRows);

        for (int i=0; i<numCols; i++){    

            if(params->method == 0){
                for (int j=0; j<numRows; j++){
                    ubkpts[i][j] = params->pde->localExpansionApprox(params->xbkpts[i],params->tend[j]);
                }
            } else {
                CDoubleArray u = params->pde->globalIterationApprox(params->xbkpts[i],params->tend,params->methodIdx);
                for (int j=0; j<numRows; j++){
                    ubkpts[i][j] = u[j];
                }
            }
        }

        return IObjectSP(copy(&(ubkpts))); 
    }

private:
    CDoubleArray    tend;
    CDoubleArray    xbkpts;       
    BuscaSP         pde;  
    int             method;     // 0 for local expansion, otherwise for global iteration
    int             methodIdx;  //for global iteration method, i for ith order approximation, i = 0,1,2

    /** for reflection */ 
    BuscaApproxAddin():CObject(TYPE), method(1), methodIdx(2) {}

    static IObject* defaultBuscaApproxAddin(){
        return new BuscaApproxAddin();
    } 

    /**Invoke when Class is 'load' */ 
    static void load(CClassSP& clazz){
        REGISTER(BuscaApproxAddin,clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultBuscaApproxAddin);
        FIELD_INLINE(tend, "ending time.");
        FIELD_INLINE(xbkpts, "spatial points.");
        FIELD(pde, "Busca PDE");
        FIELD_INLINE(method, "method choice");
        FIELD_MAKE_OPTIONAL(method);
        FIELD_INLINE(methodIdx, "method index");
        FIELD_MAKE_OPTIONAL(methodIdx);

        Addin::registerClassObjectMethod("BUSCA_APPROX",
                            Addin::UTILITIES,
                            "Analytical approx of a Busca PDE",
                            TYPE,
                            false,
                            Addin::expandMulti, // returnHandle ?
                            (Addin::ObjMethod*)BuscaApprox);

    }
};

CClassConstSP const BuscaApproxAddin::TYPE = CClass::registerClassLoadMethod(
    "BuscaApproxAddin", typeid(BuscaApproxAddin), BuscaApproxAddin::load);

bool BuscaApproxAddinLinkIn() {
    return BuscaApproxAddin::TYPE != 0;
}

class GatheralApproxAddin:public CObject{
public:
    static CClassConstSP const TYPE;

    /**the 'addin function' */ 
    static IObjectSP GatheralApprox(GatheralApproxAddin* params){
   
        int             numCols = params->xbkpts.size();
        int             numRows = params->tend.size();
        DoubleMatrix    ubkpts(numCols, numRows);

        for (int i=0; i<numCols; i++){    

            for (int j=0; j<numRows; j++){
               
                ubkpts[i][j] = params->pde->GatheralFirstApprox(params->xbkpts[i],params->tend[j]);
            }
        }
        
        return IObjectSP(copy(&(ubkpts)));
    }

private:
    CDoubleArray    tend;
    CDoubleArray    xbkpts;       
    BuscaSP         pde;  

    /** for reflection */ 
    GatheralApproxAddin():CObject(TYPE) {}

    static IObject* defaultGatheralApproxAddin(){
        return new GatheralApproxAddin();
    }

    /**Invoke when Class is 'load' */ 
    static void load(CClassSP& clazz){
        REGISTER(GatheralApproxAddin,clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGatheralApproxAddin);
        FIELD_INLINE(tend, "ending time.");
        FIELD_INLINE(xbkpts, "spatial points.");
        FIELD(pde, "Busca PDE");
       
        Addin::registerClassObjectMethod("GATHERAL_APPROX",
                            Addin::UTILITIES,
                            "Gatheral approx of a Busca PDE",
                            TYPE,
                            false,
                            Addin::expandMulti, // returnHandle ?
                            (Addin::ObjMethod*)GatheralApprox);

    }
};

CClassConstSP const GatheralApproxAddin::TYPE = CClass::registerClassLoadMethod(
    "GatheralApproxAddin", typeid(GatheralApproxAddin), GatheralApproxAddin::load);

bool GatheralApproxAddinLinkIn() {
    return GatheralApproxAddin::TYPE != 0;
}


class SeparableVHTAddin:public CObject{
public:
    static CClassConstSP const TYPE;

    /**the 'addin function' */ 
    static IObjectSP SolutionOfZeroOrderBusca(SeparableVHTAddin* params){
   
        int             numCols = params->xbkpts.size();
        int             numRows = params->tend.size();
        DoubleMatrix    ubkpts(numCols, numRows);

        for (int i=0; i<numCols; i++){    

            for (int j=0; j<numRows; j++){
               
                ubkpts[i][j] = params->pde->SolutionOfZeroOrderBuscaWithSeparableVHT(params->xbkpts[i],params->tend[j]);
            }
        }
        
        return IObjectSP(copy(&(ubkpts)));
    }

private:
    CDoubleArray    tend;
    CDoubleArray    xbkpts;       
    BuscaSP         pde;  

    /** for reflection */ 
    SeparableVHTAddin():CObject(TYPE) {}

    static IObject* defaultSeparableVHTAddin(){
        return new SeparableVHTAddin();
    }

    /**Invoke when Class is 'load' */ 
    static void load(CClassSP& clazz){
        REGISTER(SeparableVHTAddin,clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultSeparableVHTAddin);
        FIELD_INLINE(tend, "ending time.");
        FIELD_INLINE(xbkpts, "spatial points.");
        FIELD(pde, "Busca PDE");
       
        Addin::registerClassObjectMethod("ZERO_ORDER_BUSCA",
                            Addin::UTILITIES,
                            "Solution of zero order Busca with separable VHT",
                            TYPE,
                            false,
                            Addin::expandMulti, // returnHandle ?
                            (Addin::ObjMethod*)SolutionOfZeroOrderBusca);

    }
};


CClassConstSP const SeparableVHTAddin::TYPE = CClass::registerClassLoadMethod(
    "SeparableVHTAddin", typeid(SeparableVHTAddin), SeparableVHTAddin::load);

bool SeparableVHTAddinLinkIn() {
    return SeparableVHTAddin::TYPE != 0;
}

    
class BuscaAsymptoticsAddin:public CObject{
public:
    static CClassConstSP const TYPE;

    /**the 'addin function' */ 
    static IObjectSP BuscaAsymptotics(BuscaAsymptoticsAddin* params){
   
        int             i, num = params->bkpts.size();
        CDoubleArray    u(num);

        switch(params->type){
        case 0:
            for (i=0; i<num; i++){      
                u[i] = params->pde->solutionAtTimeZero(params->bkpts[i]);
            }
            break;
        case 1:
            for (i=0; i<num; i++){      
                u[i] = params->pde->solutionAtSpacePositiveInfinity(params->bkpts[i]);
            }
            break;
        case -1:
            for (i=0; i<num; i++){      
                u[i] = params->pde->solutionAtSpaceNegativeInfinity(params->bkpts[i]);
            }
            break;
        }

        return IObjectSP(copy(&(u)));
    }

private:
    int             type;  // 0 for solution at time zero, 1 for solution at positive infinity and -1 for solution at negative infinity
    CDoubleArray    bkpts; // maturity or moneyness array      
    BuscaSP         pde;  

    /** for reflection */ 
    BuscaAsymptoticsAddin():CObject(TYPE) {}

    static IObject* defaultBuscaAsymptoticsAddin(){
        return new BuscaAsymptoticsAddin();
    }

    /**Invoke when Class is 'load' */ 
    static void load(CClassSP& clazz){
        REGISTER(BuscaAsymptoticsAddin,clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultBuscaAsymptoticsAddin);
        FIELD_INLINE(type, "asymptotics type");
        FIELD_INLINE(bkpts, "maturity or moneyness array");
        FIELD(pde, "Busca PDE");
       
        Addin::registerClassObjectMethod("BUSCA_ASYMPTOTICS",
                            Addin::UTILITIES,
                            "Asymptotics of a Busca PDE",
                            TYPE,
                            false,
                            Addin::expandMulti, // returnHandle ?
                            (Addin::ObjMethod*)BuscaAsymptotics);

    }
};

CClassConstSP const BuscaAsymptoticsAddin::TYPE = CClass::registerClassLoadMethod(
    "BuscaAsymptoticsAddin", typeid(BuscaAsymptoticsAddin), BuscaAsymptoticsAddin::load);

bool BuscaAsymptoticsAddinLinkIn() {
    return BuscaAsymptoticsAddin::TYPE != 0;
}
  
DRLIB_END_NAMESPACE
