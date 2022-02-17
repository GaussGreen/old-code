//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLHelp.cpp
//
//   Description : Help type addin functions
//
//   Author      : Mark A Robson
//
//   Date        : 28 Feb 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XLAddin.hpp"
#include "edginc/Version.hpp"
#include "edginc/Modifier.hpp"
#include <fstream>
#include <algorithm>

DRLIB_BEGIN_NAMESPACE

class XLHelp: public CObject{
public:
    static CClassConstSP const TYPE;

    // single addin parameter (optional). If present detail addin function else
    // return list of all available addin functions
    string   addinName;

    /** return either a list of all addins or a description about a specific
        one */
    static IObjectSP addinHelp(XLHelp *params){
        if (params->addinName.empty()){
            CStringArraySP names(Addin::names());
            // now add prefix
            for (int i = 0; i < names->size(); i++){
                (*names)[i] = XLAddin::ADDIN_PREFIX + (*names)[i];
            }
            return names;
        }
        const Addin* addin = 
            Addin::lookUp(XLAddin::stripAddinPrefix(params->addinName));
        /* return something like :
           1. <parameter name>  <parameter type> <optional or not> <info>
           2.      "                  "               "              "
           and so on */
        CFieldArray fields(addin->getDataClassFields());
        ObjectArraySP output(new ObjectArray(5)); // 5 columns
        int offset = addin->hasHandleNameParam()? 1: 0;
        int numRows = fields.size() + offset;
        // construct our arrays
        CIntArraySP paramsIdx(new IntArray(numRows));
        CStringArraySP paramName(new StringArray(numRows));
        CStringArraySP paramType(new StringArray(numRows));
        CStringArraySP paramOptional(new StringArray(numRows));
        CStringArraySP paramInfo(new StringArray(numRows));
        // save in the output
        (*output)[0] = paramsIdx;
        (*output)[1] = paramName;
        (*output)[2] = paramType;
        (*output)[3] = paramOptional;
        (*output)[4] = paramInfo;
        if (addin->hasHandleNameParam()){
            (*paramsIdx)[0] = 1;
            (*paramName)[0] = "userHandleName";
            (*paramType)[0] = CString::TYPE->getName();
            (*paramOptional)[0] = "Mandatory";
            (*paramInfo)[0] = "Prefix for new handle";
        }
        for (unsigned int i = 0; i < fields.size(); i++){
            (*paramsIdx)[i+offset] = i + 1 + offset;
            (*paramName)[i+offset] = fields[i]->getName();
            (*paramType)[i+offset] = fields[i]->getType()->getName();
            (*paramOptional)[i+offset] = fields[i]->isOptional()? 
                "Optional": "Mandatory";
            (*paramInfo)[i+offset] = fields[i]->getDescription();
        }
        return output;
    }

    /** Provide detailed help about a specific addin */
    static string addinInfo(XLHelp *params){
        if (params->addinName.empty()){
            return "Please supply an addin name";
        }
        const Addin* addin = 
            Addin::lookUp(XLAddin::stripAddinPrefix(params->addinName));
        return addin->getDescription();
    }

    /** What category does an addin belong to */
    static string addinCategory(XLHelp *params){
        if (params->addinName.empty()){
            return "Please supply an addin name";
        }
        const Addin* addin = 
            Addin::lookUp(XLAddin::stripAddinPrefix(params->addinName));
        return addin->getCategory();
    }

    XLHelp(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLHelp, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLHelp);
        FIELD(addinName, "Name of addin function");
        FIELD_MAKE_OPTIONAL(addinName);
        Addin::registerClassObjectMethod("HELP",
                                         Addin::UTILITIES,
                                         "Provides help for addin functions",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)addinHelp);

        Addin::registerClassStringMethod("ADDIN_INFO",
                                         Addin::UTILITIES,
                                         "Returns a description of an addin",
                                         TYPE,
                                         (Addin::StringMethod*)addinInfo);

        Addin::registerClassStringMethod("ADDIN_CATEGORY",
                                         Addin::UTILITIES,
                                         "What category does an addin belong to",
                                         TYPE,
                                         (Addin::StringMethod*)addinCategory);
    }

    static IObject* defaultXLHelp(){
        return new XLHelp();
    }
};

CClassConstSP const XLHelp::TYPE = CClass::registerClassLoadMethod(
    "XLHelp", typeid(XLHelp), load);

/** class for dealing with addin parameters that take no parameters */
class XLNoParams: public CObject{
public:
    static CClassConstSP const TYPE;

    //// list all types that aren't private
    static IObjectSP listTypes(XLNoParams* params){
        const CClassVec& allClasses = CClass::allClasses();
        CStringArraySP types(new CStringArray());
        types->reserve(allClasses.size());
        for (unsigned int i = 0; i < allClasses.size(); i++){
            // hide private classes from spreadsheet
            if (!(Modifier::isPrivate(allClasses[i]->getModifiers()))){
                types->push_back(allClasses[i]->getName());
            }
        }
        // now sort them
        sort(types->begin(), types->end());
        return types;
    }

    static string version(XLNoParams* params){
        return CVersion::DRLibVersion();
    }

    XLNoParams(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLNoParams, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLNoParams);
        Addin::registerClassObjectMethod("LIST_TYPES",
                                         Addin::UTILITIES,
                                         "Lists all known types",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)listTypes);
        Addin::registerClassStringMethod("VERSION",
                                         Addin::UTILITIES,
                                         "Returns the version of the library",
                                         TYPE,
                                         (Addin::StringMethod*)version);
     }

    static IObject* defaultXLNoParams(){
        return new XLNoParams();
    }
};

CClassConstSP const XLNoParams::TYPE = CClass::registerClassLoadMethod(
    "XLNoParams", typeid(XLNoParams), load);


class TypeHelp: public CObject{
public:
    static CClassConstSP const TYPE;
    // single addin parameter 
    string   typeName;

    /** return description about the fields making up the type */
    static IObjectSP typeHelp(TypeHelp *params){
        CClassConstSP clazz = CClass::forName(params->typeName);
        if (clazz->isInterface()){
            return CStringSP(CString::create(params->typeName + " is an inter"
                                             "face. It has no components"));
        }
        if (!Modifier::isPublic(clazz->getModifiers())){
            return CStringSP(CString::create(params->typeName + " is not a "
                                             "public type"));
        }

        /* return something like :
           <parameter name>  <parameter type> <optional or not> <info>
                 "                  "               "              "
           and so on */
        CFieldArray fields(Addin::getDataClassFields(clazz));
        ObjectArraySP output(new ObjectArray(4)); // 4 columns
        int numRows = fields.size();
        // construct our arrays
        CStringArraySP paramName(new StringArray(numRows));
        CStringArraySP paramType(new StringArray(numRows));
        CStringArraySP paramOptional(new StringArray(numRows));
        CStringArraySP paramInfo(new StringArray(numRows));
        // save in the output
        (*output)[0] = paramName;
        (*output)[1] = paramType;
        (*output)[2] = paramOptional;
        (*output)[3] = paramInfo;
        for (unsigned int i = 0; i < fields.size(); i++){
            (*paramName)[i] = fields[i]->getName();
            (*paramType)[i] = fields[i]->getType()->getName();
            (*paramOptional)[i] = fields[i]->isOptional()? 
                "Optional": "Mandatory";
            (*paramInfo)[i] = fields[i]->getDescription();
        }
        return output;
    }

    /** return string representation of modifiers to supplied class */
    static string modifiersHelp(TypeHelp *params){
        CClassConstSP clazz = CClass::forName(params->typeName);
        return Modifier::toString(clazz->getModifiers());
    }

    /** Returns an array of strings which are the types which can used
        to build an object of the supplied type via the data
        dictionary route */
    static IObjectSP typeConstructors(TypeHelp *params){
        CClassConstSP clazz = CClass::forName(params->typeName);
        CClassVec classes(clazz->CClass::listAllConstructorClasses());
        CStringArraySP types(new CStringArray(classes.size()));
        for (unsigned int i = 0; i < classes.size(); i++){
            (*types)[i] = classes[i]->getName();
        }
        sort(types->begin(), types->end());
        return types;
    }

    /** Returns a pair of array of strings which are the possible
        values for the enum of the specified type together with their comment */
    IObjectSP enumValues(){
        CClassConstSP clazz = CClass::forName(typeName);
        const vector<CClass::EnumValue>* enumValues = clazz->getEnumValues();
        if (!enumValues){
            throw ModelException("TypeHelp::enumValues",
                                 typeName+" is not an enum");
        }
        ObjectArraySP output(new ObjectArray(2));
        StringArraySP values(new StringArray(enumValues->size()));
        StringArraySP comments(new StringArray(enumValues->size()));
        (*output)[0] = values;
        (*output)[1] = comments;
        for (unsigned int i = 0; i < enumValues->size(); i++){
            (*values)[i] = (*enumValues)[i].valueAsString;
            (*comments)[i] = (*enumValues)[i].comment;
        }
        return output;
    }

    static IObjectSP driTypeHelp(TypeHelp *params){
        CClassConstSP clazz = CClass::forName(params->typeName);
        return clazz->getDRIType();
    }


    TypeHelp(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(TypeHelp, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultTypeHelp);
        FIELD(typeName, "Type name");
        Addin::registerClassObjectMethod("TYPE_COMPONENTS",
                                         Addin::UTILITIES,
                                         "Provides details of the components"
                                         " that makes up a type",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)typeHelp);
        
        Addin::registerClassStringMethod("TYPE_PROPERTIES",
                                         Addin::UTILITIES,
                                         "Provides details of the type's"
                                         " properties",
                                         TYPE,
                                         (Addin::StringMethod*)modifiersHelp);
        Addin::registerClassObjectMethod("TYPE_CONSTRUCTORS",
                                         Addin::UTILITIES,
                                         "Lists types that can be used to"
                                         " build object of supplied type",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)typeConstructors);
        
        Addin::registerClassObjectMethod("TYPE_INFO",
                                         Addin::UTILITIES,
                                         "Gives information in DRI style",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)driTypeHelp);
        Addin::registerObjectMethod("ENUM_VALUES",
                                    Addin::UTILITIES,
                                    "Lists possible values for enums of the "
                                    "specified type",
                                    false,
                                    Addin::expandMulti,
                                    &TypeHelp::enumValues);
                                    
    }
    
    static IObject* defaultTypeHelp(){
        return new TypeHelp();
    }
};

CClassConstSP const TypeHelp::TYPE = CClass::registerClassLoadMethod(
    "TypeHelp", typeid(TypeHelp), load);

/** class for addin function which returns the default value for fields (if
    they are optional) */
class XLDefaultValues: public CObject{
public:
    static CClassConstSP const TYPE;

    string   className;
    string   fieldName;

    static IObjectSP defaultValues(XLDefaultValues *params){
        // 1. look up class
        CClassConstSP clazz = CClass::forName(params->className);
        // 2. look up field if it exists
        CClassConstSP c = clazz;
        CFieldConstSP field(0);
        do{
            if (c->hasDeclaredField(params->fieldName)){
                field = c->getDeclaredField(params->fieldName);
            }
            // recurse over parents
        } while (!field && (c = c->getSuperClass()));
        if (!field){
            throw ModelException("Class "+params->className+" does not have "
                                 "a field called "+params->fieldName);
        }
        if (Modifier::isAbstract(clazz->getModifiers())){
            throw ModelException("Class "+params->className+" is abstract. "
                                 "Abstract classes are not supported. Choose "
                                 "a concrete instance instead");
        }

        IObjectSP returnVal;
        // if the field is mandatory then there is no default value. However
        // some people like to have "suggested" values - since we're at the
        // spreadsheet level we'll be generous and allow access to non
        // optional fields. A better way would be to have some marker interface
        // that classes implement if they want to provide suggested values.
        // Implementing this interface => non default values are populated
        // NB restricted to atomic values 
        // Perhaps better to have a different addin eg suggested value
        CClassConstSP type = field->getType();
        if (!field->isOptional() && 
            !CDouble::TYPE->isAssignableFrom(type) &&
            !CInt::TYPE->isAssignableFrom(type) && 
            !CBool::TYPE->isAssignableFrom(type) &&
            !CString::TYPE->isAssignableFrom(type)){
            returnVal = IObjectSP(CString::create("No default"));
        } else {
            // 3. create empty instance of class
            IObjectSP newInstance(clazz->newInstance());
            IObjectSP defaultVal(field->get(newInstance));
            returnVal = IObjectSP(defaultVal.clone()); /* best to clone to 
                                                          avoid memory issues*/
        }
        return returnVal;
    }


    XLDefaultValues(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLDefaultValues, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLDefaultValues);
        FIELD(className, "Name of Class");
        FIELD(fieldName, "Name of Field");
        Addin::registerClassObjectMethod("DEFAULT_VALUE",
                                         Addin::UTILITIES,
                                         "Gives default values for parameters",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)defaultValues);

    }

    static IObject* defaultXLDefaultValues(){
        return new XLDefaultValues();
    }
};

CClassConstSP const XLDefaultValues::TYPE = CClass::registerClassLoadMethod(
    "XLDefaultValues", typeid(XLDefaultValues), load);


class XLDoc: public CObject{
public:
    static CClassConstSP const TYPE;

    static const string HTML_START;
    static const string HEADER_START;
    static const string HEADER_END;
    static const string BODY_START;
    static const string BODY_END;
    static const string HTML_END;
    static const string MAIN_TITLE_START;
    static const string NAME_START;
    static const string NAME_END;
    static const string TITLE_START;
    static const string TITLE_END;
    static const string MAIN_DESC_START;
    static const string MAIN_DESC_END;
    static const string TABLE_START;
    static const string TABLE_END;
    static const string CONTENTS_ROW_START;
    static const string TABLE_ROW_START;
    static const string TABLE_ROW_END;
    static const string OUT_TABLE_ROW_END;
    static const string COL1_START;
    static const string COL1_END;
    static const string COL2_START;
    static const string COL3_START;
    static const string COL2OR3_END;
    static const string CONTENTS_COL1_START;
    static const string CONTENTS_COL1_MID;
    static const string CONTENTS_COL1_END;
    static const string CONTENTS_COL2_START;

    string   filename;

    /** write out documentation for all addins to file name supplied */
    static bool addinDoc(XLDoc *params){

        // open doc to write to
        ofstream doc(params->filename.c_str());
        if (!doc) {
            throw ModelException("XLDoc::addinDoc",
                                 "Couldn't open " + params->filename);
        }

        // get list of all addins
        CStringArraySP names(Addin::names());

        int i;

        doc << HTML_START << endl;
        doc << MAIN_TITLE_START << "Quantitative Research QLib Add-in" << "<br>"; 
        doc << "(Version " << CVersion::DRLibVersion() << ")" << TITLE_END << endl;
        doc << "<br>" << endl;

        // table of contents - TO DO: links to addin pages
        doc << TABLE_START << endl;
        doc << CONTENTS_ROW_START << "Function Name" << OUT_TABLE_ROW_END;
        for (i = 0; i < names->size(); i++){
            const Addin* addin = 
                Addin::lookUp(XLAddin::stripAddinPrefix((*names)[i]));
            doc << "<tr>" << CONTENTS_COL1_START << (*names)[i] << CONTENTS_COL1_MID << (*names)[i] << CONTENTS_COL1_END;
            doc << CONTENTS_COL2_START << addin->getDescription() << COL2OR3_END << "</tr>";
        }
        doc << TABLE_END << endl;

        // write each addin out in turn
        for (i = 0; i < names->size(); i++){
            doc << HEADER_START << (*names)[i] << HEADER_END << endl;
            doc << BODY_START << endl;
            doc << NAME_START << (*names)[i] << NAME_END << endl;
            doc << TITLE_START << (*names)[i] << TITLE_END << endl;

            const Addin* addin = 
                Addin::lookUp(XLAddin::stripAddinPrefix((*names)[i]));

            doc << MAIN_DESC_START << addin->getDescription() << MAIN_DESC_END << endl;

            CFieldArray fields(addin->getDataClassFields());
            ObjectArraySP output(new ObjectArray(5)); // 5 columns
            int offset = addin->hasHandleNameParam()? 1: 0;
            int numRows = fields.size() + offset;
            // construct our arrays
            CIntArraySP paramsIdx(new IntArray(numRows));
            CStringArraySP paramName(new StringArray(numRows));
            CStringArraySP paramType(new StringArray(numRows));
            CStringArraySP paramOptional(new StringArray(numRows));
            CStringArraySP paramInfo(new StringArray(numRows));
            if (addin->hasHandleNameParam()){
                (*paramsIdx)[0] = 1;
                (*paramName)[0] = "Handle Name";
                (*paramType)[0] = CString::TYPE->getName();
                (*paramOptional)[0] = "Mandatory";
                (*paramInfo)[0] = "Prefix for new handle";
            }
            int j;
            for (j = 0; j < (int)fields.size(); j++){
                (*paramsIdx)[j+offset] = j + 1 + offset;
                (*paramName)[j+offset] = fields[j]->getName();
                (*paramType)[j+offset] = fields[j]->getType()->getName();
                (*paramOptional)[j+offset] = fields[j]->isOptional()? 
                    "Optional": "Mandatory";
                (*paramInfo)[j+offset] = fields[j]->getDescription();
            }

            // write out inputs
            if (fields.size() > 0) {
                doc << TABLE_START << endl;
                doc << TABLE_ROW_START << "Input" << TABLE_ROW_END;
                for (j = 0; j < numRows; j++){
                    doc << "<tr>" << COL1_START << j+1 << COL1_END << COL2_START;
                    doc << (*paramName)[j] << "</b> [" << (*paramType)[j] << "]";
                    doc << "<br>" << (*paramOptional)[j];
                    doc << COL2OR3_END << COL3_START;
                    doc << (*paramInfo)[j];
                    doc << COL2OR3_END << "</tr>";
                }
                doc << TABLE_END << endl;
            }
            
            doc << "<br>" << endl;

            // now output
            doc << TABLE_START << endl;
            doc << TABLE_ROW_START << "Returns" << OUT_TABLE_ROW_END;
            doc << "<tr>" << COL1_START << "1" << COL1_END << COL2_START;
            doc << "[" << addin->getReturnType()->getName() << "]";
            doc << COL2OR3_END << COL3_START;
            doc << addin->getDescription() << COL2OR3_END << "</tr>";
            doc << TABLE_END << endl;
                
            doc << BODY_END << endl;
            doc << HEADER_END << endl;
        }

        doc << HTML_END << endl;

        return true;
    }

    XLDoc(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLDoc, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLDoc);
        FIELD(filename, "Where to write document");

        Addin::registerClassBoolMethod("DOC",
                                       Addin::UTILITIES,
                                       "Generates HTML addin documentation",
                                       TYPE,
                                       (Addin::BoolMethod*)addinDoc);
    }

    static IObject* defaultXLDoc(){
        return new XLDoc();
    }
};

const string XLDoc::HTML_START = "<HTML>";
const string XLDoc::HEADER_START = "<HEAD><TITLE>";
const string XLDoc::HEADER_END = "</TITLE></HEAD>";
const string XLDoc::BODY_START = "<BODY bgcolor=\"#FFFFFF\" link=\"#800000\" vlink=\"#800080\">";
const string XLDoc::BODY_END = "</BODY>";
const string XLDoc::HTML_END = "</HTML>";
const string XLDoc::MAIN_TITLE_START = "<p align=\"center\"><font size=\"7\" face=\"Tahoma\"><b>";
const string XLDoc::NAME_START = "<a name=\"";
const string XLDoc::NAME_END = "\"></a>";
const string XLDoc::TITLE_START = "<p align=\"center\"><font size=\"5\" face=\"Tahoma\"><b>";
const string XLDoc::TITLE_END = "</b></font></p>";
const string XLDoc::MAIN_DESC_START = "<p><font size=\"3\" face=\"Tahoma\" color=\"#0000A0\"><I>";
const string XLDoc::MAIN_DESC_END = "</I></b></font></p>";;
const string XLDoc::TABLE_START = "<table border=\"1\" cellpadding=\"7\" cellspacing=\"1\" width=\"100%\" bordercolor=\"#000000\">";
const string XLDoc::TABLE_END = "</table>";
const string XLDoc::CONTENTS_ROW_START = "<tr><td valign=\"top\" width=\"50%\" bgcolor=\"#3161FF\" height=\"18\">\n"
                                         "<font color=\"#FFFFFF\" size=\"3\" face=\"Tahoma\"><b>";
const string XLDoc::TABLE_ROW_START = "<tr><td valign=\"top\" width=\"6%\" height=\"18\">&nbsp;</td>\n"
                                    "<td valign=\"top\" width=\"63%\" bgcolor=\"#3161FF\" height=\"18\">\n"
                                    "<font color=\"#FFFFFF\" size=\"3\" face=\"Tahoma\"><b>";
const string XLDoc::TABLE_ROW_END   = " Parameters</b></font></td>\n"
                                    "<td valign=\"top\" width=\"31%\" bgcolor=\"#3161FF\" height=\"18\">\n"
                                    "<font color=\"#FFFFFF\" size=\"3\" face=\"Tahoma\"><b>Description</b></font></td></tr>\n";
const string XLDoc::OUT_TABLE_ROW_END = "</b></font></td>\n"
                                        "<td valign=\"top\" width=\"31%\" bgcolor=\"#3161FF\" height=\"18\">\n"
                                        "<font color=\"#FFFFFF\" size=\"3\" face=\"Tahoma\"><b>Description</b></font></td></tr>\n";

const string XLDoc::COL1_START = "<td align=\"center\" valign=\"top\" width=\"6%\""
                                " bgcolor=\"#3161FF\"><font color=\"#FFFFFF\" size=\"2\" face=\"Tahoma\"><strong>";
const string XLDoc::COL1_END = "</strong></font></td>";
const string XLDoc::COL2_START = "<td valign=\"top\" width=\"63%\"><font size=\"2\" face=\"Tahoma\"><b>";
const string XLDoc::COL3_START = "<td valign=\"top\" width=\"31%\"><font size=\"2\" face=\"Tahoma\">";
const string XLDoc::COL2OR3_END = "</font></td>";
const string XLDoc::CONTENTS_COL1_START = "<td valign=\"top\" width=\"50%\"><font size=\"2\" face=\"Tahoma\"><b><a href=\"#";
const string XLDoc::CONTENTS_COL1_MID = "\">";
const string XLDoc::CONTENTS_COL1_END = "</a></font></td>";
const string XLDoc::CONTENTS_COL2_START = "<td valign=\"top\" width=\"50%\"><font size=\"2\" face=\"Tahoma\">";


CClassConstSP const XLDoc::TYPE = CClass::registerClassLoadMethod(
    "XLDoc", typeid(XLDoc), load);


// symbol (referenced by Addin.cpp) to ensure file gets linked in
bool XLHelpExists = true;


DRLIB_END_NAMESPACE
