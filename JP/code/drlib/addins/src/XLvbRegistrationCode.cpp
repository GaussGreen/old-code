//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLvbGetRegistrationCode.cpp
//
//   Description : autmatically generated registration code which will add the
//                 a new type to the Visual Basic type table
//
//   Author      : André Segger
//
//   Date        : 21 Sep 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Modifier.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE

/* non static variable to force linkage of this file. Outside of class to
   avoid necessity of header file */
bool XLvbRegistrationCodeRegistered = true;


// addin support
class GetVBRegistrationCode: public CObject{
    static CClassConstSP const TYPE;

    string    typeKey;

    /** create a market data cache */
    static IObjectSP getCode(GetVBRegistrationCode* params){
        static const string routine("GetVBRegistrationCode::getCode");

        CClassConstSP clazz = CClass::forName(params->typeKey);
        if (clazz->isInterface()){
            return CStringSP(CString::create(params->typeKey + "Error: is an inter"
                                             "face. It has no components"));
        }
        if (!Modifier::isPublic(clazz->getModifiers())){
            return CStringSP(CString::create(params->typeKey + "Error: is not a"
                                             "public type"));
        }

        CStringArraySP    codeFragment(new StringArray(0));

        CFieldArray fields(Addin::getDataClassFields(clazz));
        int numRows = fields.size();

        codeFragment->push_back("Public Sub RegisterInstrument()");
        codeFragment->push_back("   Dim typeDescriptor(0 To " + Format::toString(numRows) + ") As Variant");
        codeFragment->push_back("   typeDescriptor(0) = Array(\"" + params->typeKey + "\", " + Format::toString(numRows) + ")" );

        for (unsigned int i = 0; i < fields.size(); i++) {

            codeFragment->push_back("   ' " + fields[i]->getDescription() + 
                                    (fields[i]->isOptional()? " (Optional)": " (Mandatory)"));
            if (fields[i]->isOptional()){
                codeFragment->push_back("   typeDescriptor(" + 
                                        Format::toString(int(i+1)) + 
                                        ") = Array(\"" + fields[i]->getName() + "\", \"" + 
                                        fields[i]->getType()->getName() + "\", \"OPTIONAL\", EXPAND_NONE, " +
                                        "\"<RANGE_NAME>\", \"<RANGE_NAME>\")");
            }
            else{
                codeFragment->push_back("   typeDescriptor(" + 
                                        Format::toString(int(i+1)) + 
                                        ") = Array(\"" + fields[i]->getName() + "\", \"" + 
                                        fields[i]->getType()->getName() + "\", \"<RANGE_NAME>\", EXPAND_NONE )");
            }
        }

        codeFragment->push_back("   Call EQD_AddToTypeTable(typeDescriptor)");
        codeFragment->push_back("End Sub");

        return codeFragment;
    }
 
    GetVBRegistrationCode(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(GetVBRegistrationCode, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetVBRegistrationCode);
        FIELD(typeKey,          "Object type key");
        // register addin for adding objects to cache
        Addin::registerClassObjectMethod(
            "GET_VB_TYPETABLE_CODE",
            Addin::MARKET,
            "Creates vb code to add an object type to the VB type table",
            TYPE,
            false,
            Addin::expandSimple,
            (Addin::ObjMethod*)getCode);
    }

    static IObject* defaultGetVBRegistrationCode(){
        return new GetVBRegistrationCode();
    }
};

CClassConstSP const GetVBRegistrationCode::TYPE= CClass::registerClassLoadMethod(
    "GetVBRegistrationCode", typeid(GetVBRegistrationCode), GetVBRegistrationCode::load);

DRLIB_END_NAMESPACE

