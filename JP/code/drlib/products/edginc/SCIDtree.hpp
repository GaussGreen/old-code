/**
 * @file SCIDtree.hpp
 */

#ifndef QR_SCIDtree_HPP
#define QR_SCIDtree_HPP

#include "edginc/Model.hpp"
#include "edginc/Results.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/CDSIndexParSpreads.hpp" //required for wrapper array somehow
#include "edginc/SCIDtreeParameters.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(SCIDtree)
FORWARD_DECLARE(CInstrument)
FORWARD_DECLARE(Control)

/** sCID model contains all parameters used by sCID in all types of Monte Carlo simulations. 
 *  The model diffuses spreads, simulates defaults, and computes a curve of expected losses
 FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
 */
  
class PRODUCTS_DLL SCIDtree: public CModel {
public:
    static CClassConstSP const TYPE;
    friend class SCIDtreeHelper;

    /** the class that the product must be able to create */
    class PRODUCTS_DLL IProduct{
    public:
        virtual void price(SCIDtree* model,
                           Control*        control, 
                           CResults*       results) const = 0;
        virtual ~IProduct() {};
    };

    /** interface that the instrument must implement */
    class PRODUCTS_DLL IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        friend class SCIDtreeHelper;
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(SCIDtree* model) const = 0;
    };

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);

 
    void getMarket(const MarketData*  market,
                           IInstrumentCollectionSP instruments);
    /** get the number of paths for FAST simulation **/
    IntArray    getNrFastPaths(){return nbPathsFast;}
    SCIDtree();
    SCIDtree(CClassConstSP clazz);
    virtual ~SCIDtree();
	string getMCAlgorithm() { return MCmethod; };
	int getSeed() { return seed; };
	double getCFtimeSteps() { return timeSteps; }
	MaturityPeriodSP getFreqFastMC() { return freqFastMC; };
	MaturityPeriodSP getFreqFullMC() { return freqFullMC; };
	int getNbPathsFastNoJump() { return nbPathsFast[0]; };
	int getNbPathsFastAtLeastOneJump() { return nbPathsFast[1]; };
	int getNbPathsFull() { return nbPathsFull; };
	int getConvolutionNoJump() { return ConvolutionMethod[0]; };
	int getConvolutionAtLeastOneJump() { return ConvolutionMethod[1]; };
    SCIDtreeParametersSP getSCIDtreeParameters(){return SCIDtreeParam.getSP();};
	double getTest(int i) { return test[i]; }

private:
	string							MCmethod;  
	int								seed;
	MaturityPeriodSP				freqFastMC, freqFullMC;
	double							timeSteps;
	IntArray						nbPathsFast, ConvolutionMethod;
	int								nbPathsFull;
	SCIDtreeParametersWrapper		SCIDtreeParam;
	DoubleArray						test;
};

DRLIB_END_NAMESPACE
#endif
