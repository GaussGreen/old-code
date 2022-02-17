//
//  Description:  ClientRunnables to read and write objects as XML.
//                (useful for enabling DRWrapper execution through the DRI).

#include "edginc/config.hpp"
#include "edginc/Addin.hpp"
#include "edginc/OutputFile.hpp"
#include "edginc/XMLReadWrite.hpp"
#include "edginc/XMLReader.hpp"
#include "edginc/XMLWriter.hpp"

DRLIB_BEGIN_NAMESPACE

XMLResultsFileReader::XMLResultsFileReader(const string& file) :
    CObject(TYPE), file(file) {
    // empty
} 

// for reflection
XMLResultsFileReader::XMLResultsFileReader() :
    CObject(TYPE), file("") {
    // empty
}

// EdrAction - Tries to read a Results object (in XML) from a file. 
IObjectSP XMLResultsFileReader::run() {
    return IObjectSP(OutputFile::xmlReadResults(file));
}

static IObjectSP addinXMLResultsFileReader(XMLResultsFileReader* xmlResultsFileReader)
{
    return xmlResultsFileReader->run();
}

class XMLResultsFileReaderHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet/DRI
        clazz->setDescription("Reads an XML results file and answers a Results object");
        REGISTER(XMLResultsFileReader, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultXMLResultsFileReader);
        FIELD(file, "Path of XML results file to read");
        // registration for addin function
        Addin::registerClassObjectMethod("XMLResultsFileReader",
                                         Addin::UTILITIES,
                                         "Reads an XML results file and answers a Results object",
                                         XMLResultsFileReader::TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)addinXMLResultsFileReader);
    }

    static IObject* defaultXMLResultsFileReader(){
        return new XMLResultsFileReader();
    }
};

CClassConstSP const XMLResultsFileReader::TYPE = CClass::registerClassLoadMethod(
    "XMLResultsFileReader", typeid(XMLResultsFileReader), XMLResultsFileReaderHelper::load); 


XMLErrorFileReader::XMLErrorFileReader(const string& file) :
    CObject(TYPE), file(file) {
    // empty
} 

// for reflection
XMLErrorFileReader::XMLErrorFileReader() :
    CObject(TYPE), file("") {
    // empty
}

// EdrAction - Tries to read error messages (in XML) from a file.
IObjectSP XMLErrorFileReader::run() {
    string err(OutputFile::xmlReadErrors(file));
    return IObjectSP(CString::create(err)); 
}

static IObjectSP addinXMLErrorFileReader(XMLErrorFileReader* xmlErrorFileReader)
{
    return xmlErrorFileReader->run();
}

class XMLErrorFileReaderHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet/DRI
        clazz->setDescription("Reads an XML error file and answers the errors as a string");
        REGISTER(XMLErrorFileReader, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultXMLErrorFileReader);
        FIELD(file, "Path of XML error file to read");
        // registration for addin function
        Addin::registerClassObjectMethod("XMLErrorFileReader",
                                         Addin::UTILITIES,
                                         "Reads an XML error file and answers the errors as a string",
                                         XMLErrorFileReader::TYPE,
					 true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)addinXMLErrorFileReader);
    }

    static IObject* defaultXMLErrorFileReader(){
        return new XMLErrorFileReader();
    }
};

CClassConstSP const XMLErrorFileReader::TYPE = CClass::registerClassLoadMethod(
    "XMLErrorFileReader", typeid(XMLErrorFileReader), XMLErrorFileReaderHelper::load); 


XMLFileReader::XMLFileReader(const string& file) :
    CObject(TYPE), file(file) {
    // empty
} 

// for reflection
XMLFileReader::XMLFileReader() :
    CObject(TYPE), file("") {
    // empty
}

// EdrAction - Tries to read an arbitrary object from a file (e.g. regression test input file)
IObjectSP XMLFileReader::run() {
        XMLReader reader(file, true);
        return IObjectSP(reader.read());
}

static IObjectSP addinXMLFileReader(XMLFileReader* xmlFileReader)
{
    return xmlFileReader->run();
}

class XMLFileReaderHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet/DRI
        clazz->setDescription("Reads an object from an XML file");
        REGISTER(XMLFileReader, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultXMLFileReader);
        FIELD(file, "Path of XML file to read");
        // registration for addin function
        Addin::registerClassObjectMethod("XMLFileReader",
                                         Addin::UTILITIES,
                                         "Reads an object from an XML file",
                                         XMLFileReader::TYPE,
					 true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)addinXMLFileReader);
    }

    static IObject* defaultXMLFileReader(){
        return new XMLFileReader();
    }
};

CClassConstSP const XMLFileReader::TYPE = CClass::registerClassLoadMethod(
    "XMLFileReader", typeid(XMLFileReader), XMLFileReaderHelper::load); 


XMLStringReader::XMLStringReader(const string& xmlString) :
    CObject(TYPE), xmlString(xmlString) {
    // empty
} 

// for reflection
XMLStringReader::XMLStringReader() :
    CObject(TYPE), xmlString("") {
    // empty
}

// EdrAction - Tries to read an arbitrary object (in XML) from a string.
IObjectSP XMLStringReader::run() {
        XMLReader reader(xmlString, false);
        return IObjectSP(reader.read());
}

static IObjectSP addinXMLStringReader(XMLStringReader* xmlStringReader)
{
    return xmlStringReader->run();
}

class XMLStringReaderHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet/DRI
        clazz->setDescription("Reads an object from a string containing XML");
        REGISTER(XMLStringReader, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultXMLStringReader);
        FIELD(xmlString, "XML string from which to read");
        // registration for addin function
        Addin::registerClassObjectMethod("XMLStringReader",
                                         Addin::UTILITIES,
                                         "Reads an object from a string containing XML",
                                         XMLStringReader::TYPE,
					 true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)addinXMLStringReader);
    }

    static IObject* defaultXMLStringReader(){
        return new XMLStringReader();
    }
};

CClassConstSP const XMLStringReader::TYPE = CClass::registerClassLoadMethod(
    "XMLStringReader", typeid(XMLStringReader), XMLStringReaderHelper::load); 


XMLFileWriter::XMLFileWriter(IObjectSP object,
                    const string& file) :
    CObject(TYPE), object(object), file(file) {
    // empty
} 

// for reflection
XMLFileWriter::XMLFileWriter() :
    CObject(TYPE), object(0), file("") {
    // empty
}

// EdrAction - Tries to write an arbitrary object to a file in XML.
IObjectSP XMLFileWriter::run() {
    XMLWriter xml(file);
    object->write("OBJECT", &xml);
    return IObjectSP(CBool::create(true));
}

static IObjectSP addinXMLFileWriter(XMLFileWriter* xmlFileWriter)
{
    return xmlFileWriter->run();
}

class XMLFileWriterHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet/DRI
        clazz->setDescription("Writes an object to a file in XML");
        REGISTER(XMLFileWriter, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultXMLFileWriter);
        FIELD(object, "Object to write as XML to file");
        FIELD(file, "Path of file to write");
        // registration for addin function
        Addin::registerClassObjectMethod("XMLFileWriter",
                                         Addin::UTILITIES,
                                         "Writes an object to a file in XML",
                                         XMLFileWriter::TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)addinXMLFileWriter);
    }

    static IObject* defaultXMLFileWriter(){
        return new XMLFileWriter();
    }
};

CClassConstSP const XMLFileWriter::TYPE = CClass::registerClassLoadMethod(
    "XMLFileWriter", typeid(XMLFileWriter), XMLFileWriterHelper::load); 


XMLStringWriter::XMLStringWriter(IObjectSP object) :
    CObject(TYPE), object(object) {
    // empty
} 

// for reflection
XMLStringWriter::XMLStringWriter() :
    CObject(TYPE), object(0) {
    // empty
}

// EdrAction - Tries to write an arbitrary object to a string in XML.
IObjectSP XMLStringWriter::run() {
    string buffer;
    XMLWriter xml(buffer, false);
    object->write("OBJECT", &xml);
    return IObjectSP(CString::create(buffer));
}

static IObjectSP addinXMLStringWriter(XMLStringWriter* xmlStringWriter)
{
    return xmlStringWriter->run();
}

class XMLStringWriterHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet/DRI
        clazz->setDescription("Returns an XML string representation of an object");
        REGISTER(XMLStringWriter, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultXMLStringWriter);
        FIELD(object, "Object to write as XML to string");
        // registration for addin function
        Addin::registerClassObjectMethod("XMLStringWriter",
                                         Addin::UTILITIES,
                                         "Returns an XML string representation of an object",
                                         XMLStringWriter::TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)addinXMLStringWriter);
    }

    static IObject* defaultXMLStringWriter(){
        return new XMLStringWriter();
    }
};

CClassConstSP const XMLStringWriter::TYPE = CClass::registerClassLoadMethod(
    "XMLStringWriter", typeid(XMLStringWriter), XMLStringWriterHelper::load); 

bool XMLReadWriteLoad(void) {
    return XMLResultsFileReader::TYPE != NULL &&
           XMLErrorFileReader::TYPE != NULL &&
           XMLFileReader::TYPE != NULL &&
           XMLStringReader::TYPE != NULL &&
           XMLFileWriter::TYPE != NULL &&
           XMLStringWriter::TYPE != NULL;
}

DRLIB_END_NAMESPACE
