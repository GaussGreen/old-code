//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ResultsSet.cpp
//
//   Description : Return type of a composite instrument
//
//   Author      : Andrew J Swain
//
//   Date        : 24 May 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ResultsSet.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

/** Does not clone supplied parameters */
ResultsSet::ResultsSet(Results*       composite,
                       CResultsArray* components) :
    CObject(TYPE), composite(composite), components(components) {
    // empty 
}


// for reflection
ResultsSet::ResultsSet() : 
    CObject(TYPE), composite(0), components(0) {
    // empty
}

class ResultsSetHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ResultsSet, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultResultsSet);
        FIELD(composite, "composite");
        FIELD(components, "components");       
    }

    static IObject* defaultResultsSet(){
        return new ResultsSet();
    }
};

CClassConstSP const ResultsSet::TYPE = 
CClass::registerClassLoadMethod("ResultsSet", typeid(ResultsSet), ResultsSetHelper::load);

/** Addin for getting results back */
class ResultsSetAddin: public CObject{
    static CClassConstSP const TYPE;

    /** addin takes two parameters - the name of the result to get
        and the result object */
    ResultsSetSP        results; // holds the results
    string              id1;
    string              id2;
    string              id3;
    int                 index;   // which component

    /** the 'addin function' - builds array of correct type */
    static IObjectSP getResult(ResultsSetAddin* params){
        static const string routine = "ResultsSetAddin::getResult";
        try{
            int size = 3;
            ResultsSP component;

            if (params->id3.empty()){
                size--;
                if (params->id2.empty()){
                    size--;
                }
            }
            if (params->id1.empty()){
                // perhaps do something else here
                throw ModelException(routine, "The first string id is empty");
            }
            if (params->index > 0) {
                CResultsArraySP components = params->results->components;
                component = (*components)[params->index-1];
            }
            
            if (size == 1 && params->id1 == Results::VALUE){
                double value;
                if (params->index == 0) {
                    value = params->results->composite->retrievePrice();
                }
                else {
                    value = component->retrievePrice();
                }

                return IObjectSP(CDouble::create(value));
            } 
            if (size < 2){
                // to do: return all values for given packet name
                throw ModelException(routine, "Supplying just the packet name"
                                     " is not supported yet");
            }
            OutputNameSP name = OutputNameSP(size == 3?
                                             new OutputName(params->id2,
                                                            params->id3):
                                             new OutputName(params->id2));
            IObjectConstSP result;

            if (params->index == 0) {
                result = 
                    params->results->composite->retrieveGreek(params->id1, name);
            }
            else {
                result = component->retrieveGreek(params->id1, name);
            }
            // return a copy of the results to the addin
            return IObjectSP(result.clone());
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    ResultsSetAddin():  CObject(TYPE), index(0) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ResultsSetAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultResultsSetAddin);
        FIELD(results, "The results from pricing");
        FIELD(id1, "Packet name or VALUE");
        FIELD(id2, "First identifier");
        FIELD_MAKE_OPTIONAL(id2);
        FIELD(id3, "Second identifier");
        FIELD_MAKE_OPTIONAL(id3);
        FIELD(index, "Component index (starts at 1)");
        FIELD_MAKE_OPTIONAL(index);
        Addin::registerInstanceObjectMethod("GET_COMPOSITE_RESULT",
                                            Addin::RISK,
                                            "Returns a specific result",
                                            TYPE,
                                            false,
                                            Addin::expandMulti,
                                            (Addin::ObjMethod*)getResult);
    }

    static IObject* defaultResultsSetAddin(){
        return new ResultsSetAddin();
    }
    
};

CClassConstSP const ResultsSetAddin::TYPE = CClass::registerClassLoadMethod(
    "ResultsSetAddin", typeid(ResultsSetAddin), load);




DRLIB_END_NAMESPACE
