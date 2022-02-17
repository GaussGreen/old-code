//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PDFCalculator.cpp
//
//   Description : 
//
//   Author      : Andrew J Swain
//
//   Date        : 15 March 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PDFCalculator.hpp"
#include "edginc/Addin.hpp"
#include "edginc/PDFRequest.hpp"
#include "edginc/Model.hpp"
#include "edginc/Asset.hpp"

DRLIB_BEGIN_NAMESPACE

IPDFCalculator::IPDFCalculator(){}
IPDFCalculator::~IPDFCalculator(){}

void IPDFCalculator::load(CClassSP& clazz){
    REGISTER_INTERFACE(IPDFCalculator, clazz);
    EXTENDS(IObject);
}    

CClassConstSP const IPDFCalculator::TYPE = CClass::registerInterfaceLoadMethod(
    "IPDFCalculator", typeid(IPDFCalculator), load);

PDFCalculator::~PDFCalculator(){}

class PDFCalculatorHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PDFCalculator, clazz);
        SUPERCLASS(CObject);
    }
};

PDFCalculator::PDFCalculator(const CClassConstSP& clazz):CObject(clazz){}

CClassConstSP const PDFCalculator::TYPE = CClass::registerClassLoadMethod(
    "PDFCalculator", typeid(PDFCalculator), PDFCalculatorHelper::load);


class PVolPDFAddin: public CObject{
    static CClassConstSP const TYPE;

    IModelSP      model;
    CMarketDataSP market;
    CAssetSP      asset;
    PDFRequestSP  request;
    DoubleArray   strikes;
    DateTime      maturity;   
    bool          failProbs;

    static IObjectSP probs(PVolPDFAddin* params){
        static const string routine = "PVolPDFAddin::probs";
        try {
            params->asset->getMarket(params->model.get(),params->market.get());

            PDFCalculatorSP pdf(params->asset->pdfCalculator(params->request.get()));

            DoubleArray output(params->strikes.size());

            try {
                pdf->probabilities(params->strikes, params->maturity, output);
            } catch(PDFCalculatorException e) {
                if(params->failProbs) {
                    throw ModelException(e);
                }
            }
                              
            return IObjectSP(output.clone());

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    static IObjectSP local(PVolPDFAddin* params){
        static const string routine = "PVolPDFAddin::local";
        try {
            params->asset->getMarket(params->model.get(),params->market.get());

            PDFCalculatorSP pdf(params->asset->pdfCalculator(params->request.get()));

            DoubleArray output(params->strikes.size());

            try {
                pdf->localDensity(params->strikes, params->maturity, output);
            } catch(PDFCalculatorException e) {
                if(params->failProbs) {
                    throw ModelException(e);
                }
            }
                              
            return IObjectSP(output.clone());

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    static IObjectSP integrated(PVolPDFAddin* params){
        static const string routine = "PVolPDFAddin::integrated";
        try {
            params->asset->getMarket(params->model.get(),params->market.get());

            PDFCalculatorSP pdf(params->asset->pdfCalculator(params->request.get()));

            DoubleArray output(params->strikes.size());

            try {
                pdf->integratedDensity(params->strikes, params->maturity, output);
            } catch(PDFCalculatorException e) {
                if(params->failProbs) {
                    throw ModelException(e);
                }
            }
                              
            return IObjectSP(output.clone());

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    PVolPDFAddin():  CObject(TYPE), failProbs(false) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PVolPDFAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultPVolPDFAddin);
        FIELD(model, "model");
        FIELD(market, "market");
        FIELD(asset, "asset");
        FIELD(request, "request");
        FIELD(strikes, "strikes");
        FIELD(maturity, "maturity");
        FIELD(failProbs, "Fail probabilities");
        FIELD_MAKE_OPTIONAL(failProbs);

        Addin::registerClassObjectMethod("PDF_PROBABILITY",
                                         Addin::RISK,
                                         "PDF probability distribution",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)probs);

        Addin::registerClassObjectMethod("PDF_LOCAL_DENSITY",
                                         Addin::RISK,
                                         "PDF local density",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)local);

        Addin::registerClassObjectMethod("PDF_INTEGRATED_DENSITY",
                                         Addin::RISK,
                                         "PDF integrated density",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)integrated);
    }

    static IObject* defaultPVolPDFAddin(){
        return new PVolPDFAddin();
    }   
};

CClassConstSP const PVolPDFAddin::TYPE = CClass::registerClassLoadMethod(
    "PVolPDFAddin", typeid(PVolPDFAddin), load);

DRLIB_END_NAMESPACE
