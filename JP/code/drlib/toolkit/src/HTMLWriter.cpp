//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : HTMLWriter.hpp
//
//   Description : The HTMLWriter writes class documentation in HTML format
//
//   Author      : Andrew J Swain
//
//   Date        : 20 December 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/HTMLWriter.hpp"
#include "edginc/Modifier.hpp"
#include "edginc/Format.hpp"
#include "edginc/Addin.hpp"
#include "edginc/ClientRunnable.hpp"
#include <fstream>

DRLIB_BEGIN_NAMESPACE

const string HTMLWriter::HTML_START = "<HTML>";
const string HTMLWriter::HEADER_START = "<HEAD><TITLE>";
const string HTMLWriter::HEADER_END = "</TITLE></HEAD>";
const string HTMLWriter::BODY_START = "<BODY bgcolor=\"#FFFFFF\" link=\"#800000\" vlink=\"#800080\">";
const string HTMLWriter::BODY_END = "</BODY>";
const string HTMLWriter::HTML_END = "</HTML>";
const string HTMLWriter::MAIN_TITLE_START = "<p align=\"center\"><font size=\"7\" face=\"Tahoma\"><b>";
const string HTMLWriter::NAME_START = "<a name=\"";
const string HTMLWriter::NAME_END = "\"></a>";
const string HTMLWriter::TITLE_START = "<p align=\"center\"><font size=\"5\" face=\"Tahoma\"><b>";
const string HTMLWriter::TITLE_END = "</b></font></p>";
const string HTMLWriter::MAIN_DESC_START = "<p><font size=\"3\" face=\"Tahoma\" color=\"#0000A0\">";
const string HTMLWriter::MAIN_DESC_END = "</b></font></p>";;
const string HTMLWriter::TABLE_START = "<table border=\"1\" cellpadding=\"7\" cellspacing=\"1\" width=\"100%\" bordercolor=\"#000000\">";
const string HTMLWriter::TABLE_END = "</table>";
const string HTMLWriter::CONTENTS_ROW_START = "<tr><td valign=\"top\" width=\"50%\" bgcolor=\"#3161FF\" height=\"18\">\n"
                                         "<font color=\"#FFFFFF\" size=\"3\" face=\"Tahoma\"><b>";
const string HTMLWriter::TABLE_ROW_START = "<tr><td valign=\"top\" width=\"6%\" height=\"18\">&nbsp;</td>\n"
                                    "<td valign=\"top\" width=\"63%\" bgcolor=\"#3161FF\" height=\"18\">\n"
                                    "<font color=\"#FFFFFF\" size=\"3\" face=\"Tahoma\"><b>";
const string HTMLWriter::TABLE_ROW_END   = " Parameters</b></font></td>\n"
                                    "<td valign=\"top\" width=\"31%\" bgcolor=\"#3161FF\" height=\"18\">\n"
                                    "<font color=\"#FFFFFF\" size=\"3\" face=\"Tahoma\"><b>Description</b></font></td></tr>\n";
const string HTMLWriter::OUT_TABLE_ROW_END = "</b></font></td>\n"
                                        "<td valign=\"top\" width=\"31%\" bgcolor=\"#3161FF\" height=\"18\">\n"
                                        "<font color=\"#FFFFFF\" size=\"3\" face=\"Tahoma\"><b>Description</b></font></td></tr>\n";

const string HTMLWriter::COL1_START = "<td align=\"center\" valign=\"top\" width=\"6%\""
                                " bgcolor=\"#3161FF\"><font color=\"#FFFFFF\" size=\"2\" face=\"Tahoma\"><strong>";
const string HTMLWriter::COL1_END = "</strong></font></td>";
const string HTMLWriter::COL2_START = "<td valign=\"top\" width=\"63%\"><font size=\"2\" face=\"Tahoma\"><b>";
const string HTMLWriter::COL3_START = "<td valign=\"top\" width=\"31%\"><font size=\"2\" face=\"Tahoma\">";
const string HTMLWriter::COL2OR3_END = "</font></td>";
const string HTMLWriter::CONTENTS_COL1_START = "<td valign=\"top\" width=\"50%\"><font size=\"2\" face=\"Tahoma\"><b><a href=\"#";
const string HTMLWriter::CONTENTS_COL1_MID = "\">";
const string HTMLWriter::CONTENTS_COL1_END = "</a></font></td>";
const string HTMLWriter::CONTENTS_COL2_START = "<td valign=\"top\" width=\"50%\"><font size=\"2\" face=\"Tahoma\">";

HTMLWriter::~HTMLWriter(){}

/** start a comment */
void HTMLWriter::commentStart() {
    write("<!--\n");
}

/** end a comment */
void HTMLWriter::commentEnd() {
    write("-->\n");
}

/** write free-form data */
void HTMLWriter::write(const string& data) {
    *stream << data;
}
 
HTMLWriter* HTMLWriter::htmlFile(const string& filename) {
    HTMLWriter* writer = new HTMLWriter;
    writer->stream = new ofstream(filename.c_str());
    return writer;
}

HTMLWriter::HTMLWriter() : stream(0) {}

static CClassVec allInterfaces(const CClass* clazz) {
    CClassVec interfaces;
    do {
        CClassVec ifaces = clazz->getInterfaces();
        for (int i = 0; i < (int)ifaces.size(); i++) {
            interfaces.push_back(ifaces[i]);
        }
    } while ((clazz = clazz->getSuperClass()));
    return interfaces;
}

void HTMLWriter::classStart(const CClass* clazz) {
    write(HTML_START);
    write(HEADER_START + clazz->getName() + HEADER_END + "\n");
    write(BODY_START + "\n");
    write(NAME_START + clazz->getName() + NAME_END + "\n");
    write(TITLE_START + clazz->getName() + TITLE_END + "\n");
    write(MAIN_DESC_START);
    write(clazz->getDescription() + "<br>[Properties:" +
          Modifier::toString(clazz->getModifiers()) + "]<br>[Derives from:" + 
          clazz->getSuperClass()->getName() + "]" + "<br>");
   write("[Interfaces:");

   CClassVec interfaces = allInterfaces(clazz); //clazz->getInterfaces();
   if (interfaces.empty()) {
       write("None");
   }
   else {
       for (int i = 0; i < (int)interfaces.size(); i++) {
          write((i > 0 ? ", " : "") + interfaces[i]->getName());
       }
   }
   write("]<br>");
   write(MAIN_DESC_END + "\n");
}

void HTMLWriter::classEnd() {
    write(BODY_END + "\n");
    write(HEADER_END + "\n");
    write(HTML_END);
    stream->flush();
}


void HTMLWriter::document(const CClass* clazz) {
    const static string routine = "HTMLWriter::document";
    try {
        classStart(clazz);

        // construct our arrays
        IntArraySP    paramsIdx(new IntArray(0));
        StringArraySP paramName(new StringArray(0));
        StringArraySP paramType(new StringArray(0));
        StringArraySP paramOptional(new StringArray(0));
        StringArraySP paramInfo(new StringArray(0));       

        CFieldArray fields(Addin::getDataClassFields(clazz));
        int i;
        for (i = 0; i < (int)fields.size(); i++){
            if (!Modifier::isTransient(fields[i]->getModifiers())){
                (*paramsIdx).push_back(i + 1);
                (*paramName).push_back(fields[i]->getName());
                (*paramType).push_back(fields[i]->getType()->getName());
                (*paramOptional).push_back(fields[i]->isOptional()? 
                                           "Optional": "Mandatory");
                (*paramInfo).push_back(fields[i]->getDescription());
            }
        }

        if (paramsIdx->size() > 0) {
            write(TABLE_START + "\n");
            write(TABLE_ROW_START + "Input" + TABLE_ROW_END);
            for (i = 0; i < paramsIdx->size(); i++){
                write("<tr>" + COL1_START + Format::toString(i+1) + COL1_END + COL2_START);
                write((*paramName)[i] + "</b> [" + (*paramType)[i] + "]");
                write("<br>" + (*paramOptional)[i]);
                write(COL2OR3_END + COL3_START);
                write((*paramInfo)[i]);
                write(COL2OR3_END + "</tr>");
            }
            write(TABLE_END +"\n");
        }

        classEnd();
    }
    catch (exception& e){
        throw ModelException(e, routine);
    }   
}

// so we show willing with the Global DR Interface
// a "service" that documents a class
class ClassDocumentation: public CObject, virtual public ClientRunnable {
    static CClassConstSP const TYPE;

    string filename;     
    string classname;

    // EdrAction
    IObjectSP run() {
        static const string method("ClassDocumentation::run");
        try {
            CClassConstSP clazz(CClass::forName(classname));
            HTMLWriter* html = HTMLWriter::htmlFile(filename);
            try {
                html->document(clazz);
            }
            catch (exception) {
                delete html;
                throw;
            }
            delete html;
            return IObjectSP(CBool::create(true));
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** for reflection */
    ClassDocumentation():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setExport();
        clazz->setDescription("write an HTML document of a class description");
        REGISTER(ClassDocumentation, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultClassDocumentation);
        FIELD(filename, "where to write the document");
        FIELD(classname, "name of the class");
    }

    static IObject* defaultClassDocumentation(){
        return new ClassDocumentation();
    }
    
};

CClassConstSP const ClassDocumentation::TYPE = CClass::registerClassLoadMethod(
    "ClassDocumentation", typeid(ClassDocumentation), load);


DRLIB_END_NAMESPACE
