//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : NumericalIntegrationLN.hpp
//
//   Description : Numerical Integration Algorithm 
//
//   Author      : Andrew J Swain
//
//   Date        : 5 April 2002
//
//
//----------------------------------------------------------------------------

#ifndef NUMERICALINTEGRATIONLN_HPP
#define NUMERICALINTEGRATIONLN_HPP
#include "edginc/ModelLN.hpp"
#include "edginc/PDFCalculator.hpp"
#include "edginc/Theta.hpp"

DRLIB_BEGIN_NAMESPACE 

/** Numerical Integration Algorithm */
class PRODUCTS_DLL NumericalIntegrationLN: public CModelLN, virtual public Theta::RestorableShift {
public:
    static CClassConstSP const TYPE;
    
    virtual void validatePop2Object();
    
    /** Override default createMDF in order to set the right MDF */
    virtual MarketDataFetcherSP createMDF() const;
    
    /** Does special things for Theta-type tweaks */
    virtual bool sensShift(Theta* shift);
    virtual void sensRestore(Theta* shift);

    /** the class that the product must be able to create */
    class PRODUCTS_DLL IProduct{
    public:
        /** invoke the pricing */
        virtual void price(NumericalIntegrationLN* model,
                           Control*                control, 
                           CResults*               results) = 0;

        /** following methods drive the numerical integration */
        /** gives the distribution of the fwd */
        virtual PDFCalculator* pdfCalculator() const = 0;

        /** given a spot level, return the payoff */
        virtual double payoff(double spot) const = 0;

        /** centre point of distribution */
        virtual double centre(const DateTime& date) const = 0;

        /** variance at T for given K */
        virtual double variance(double strike, const DateTime& date) const = 0;

        /** when are we getting the distrubution for */
        virtual DateTime time() const = 0;

        virtual ~IProduct() {};
    };

    /** interface that the instrument must implement */
    class PRODUCTS_DLL IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        friend class NumericalIntegrationLNHelper;
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(NumericalIntegrationLN* model) const = 0;
    };

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       Control*      control, 
                       CResults*     results);
 
    /** integrate payoff using Simpson's rule */
    virtual double integrate(NumericalIntegrationLN::IProduct* product) const;

    /** Constructor takes type of vol to use */
    NumericalIntegrationLN(const string& volType);
    
    /** are we doing a time shift? */
    bool doingTimeShift() const;
private:
    friend class NumericalIntegrationLNHelper;
    NumericalIntegrationLN();
    NumericalIntegrationLN(const NumericalIntegrationLN &rhs);
    NumericalIntegrationLN& operator=(const NumericalIntegrationLN& rhs);

    double lowBound(IProduct* product) const;
    double highBound(IProduct* product) const;

    int    steps;
    double hiCutOff;
    double loCutOff;

    // transient
    bool isTimeShift;    
};

DRLIB_END_NAMESPACE
#endif
