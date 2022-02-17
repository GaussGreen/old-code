//
//  Description:  ClientRunnables for reading and writing objects as XML:
//
//		  * XMLResultsFileReader: Reads a Results object from an XML (DRWrapper) results file.
//		  * XMLErrorFileReader: Reads error messages from an XML (DRWrapper) error file.
//		  * XMLFileReader: Reads an arbitrary object from an XML file (e.g. regtest input file).
//		  * XMLStringReader: Reads an arbitrary object from a string containing XML.
//		  * XMLFileWriter: Writes an arbitrary object as XML to a file.
//		  * XMLStringWriter: Creates a new string containing an arbitrary object as XML.
//
//		These are useful for enabling DRWrapper execution through the DRI.

#ifndef EDG_XML_FILE_READER_H
#define EDG_XML_FILE_READER_H

#include "edginc/ClientRunnable.hpp"

DRLIB_BEGIN_NAMESPACE

/** Reads a given XML results file **/
class ADDINS_DLL XMLResultsFileReader : public CObject,
                      virtual public ClientRunnable {
public:  
    static CClassConstSP const TYPE;
    friend class XMLResultsFileReaderHelper;
    
    XMLResultsFileReader(const string& file);

    // EdrAction 
    virtual IObjectSP run();

    IObjectSP run() const;

private:
    string file; 

    XMLResultsFileReader();
    XMLResultsFileReader(const XMLResultsFileReader& rhs);
    XMLResultsFileReader& operator=(const XMLResultsFileReader& rhs);

    /** for addin - runs inputs through XMLResultsFileReader */
    static IObjectSP addinXMLResultsFileReader(XMLResultsFileReader* xmResultslFileReader);
};

/** Reads a given XML error file **/
class ADDINS_DLL XMLErrorFileReader : public CObject,
                      virtual public ClientRunnable {
public:  
    static CClassConstSP const TYPE;
    friend class XMLErrorFileReaderHelper;
    
    XMLErrorFileReader(const string& file);

    // EdrAction 
    virtual IObjectSP run();

    IObjectSP run() const;

private:
    string file; 

    XMLErrorFileReader();
    XMLErrorFileReader(const XMLErrorFileReader& rhs);
    XMLErrorFileReader& operator=(const XMLErrorFileReader& rhs);

    /** for addin - runs inputs through XMLErrorFileReader */
    static IObjectSP addinXMLErrorFileReader(XMLErrorFileReader* xmErrorlFileReader);
};

/** Reads an arbitrary object from a given XML file **/
class ADDINS_DLL XMLFileReader : public CObject,
                      virtual public ClientRunnable {
public:  
    static CClassConstSP const TYPE;
    friend class XMLFileReaderHelper;
    
    XMLFileReader(const string& file);

    // EdrAction 
    virtual IObjectSP run();

    IObjectSP run() const;

private:
    string file; 

    XMLFileReader();
    XMLFileReader(const XMLFileReader& rhs);
    XMLFileReader& operator=(const XMLFileReader& rhs);

    /** for addin - runs inputs through XMLFileReader */
    static IObjectSP addinXMLFileReader(XMLFileReader* xmlFileReader);
};

/** Reads an arbitrary object from a given string which is assumed to contain XML **/
class ADDINS_DLL XMLStringReader : public CObject,
                      virtual public ClientRunnable {
public:  
    static CClassConstSP const TYPE;
    friend class XMLStringReaderHelper;
    
    XMLStringReader(const string& xmlString);

    // EdrAction 
    virtual IObjectSP run();

    IObjectSP run() const;

private:
    string xmlString; 

    XMLStringReader();
    XMLStringReader(const XMLStringReader& rhs);
    XMLStringReader& operator=(const XMLStringReader& rhs);

    /** for addin - runs inputs through XMLStringReader */
    static IObjectSP addinXMLStringReader(XMLStringReader* xmlStringReader);
};

/** Writes a given object as XML to a file */ 
class ADDINS_DLL XMLFileWriter : public CObject,
                      virtual public ClientRunnable {
public:  
    static CClassConstSP const TYPE;
    friend class XMLFileWriterHelper;
    
    /** Does not clone object parameter */
    XMLFileWriter(IObjectSP object,
		  const string& file);

    // EdrAction 
    virtual IObjectSP run();

    IObjectSP run() const;

private:
    IObjectSP object;
    string file; 

    XMLFileWriter();
    XMLFileWriter(const XMLFileWriter& rhs);
    XMLFileWriter& operator=(const XMLFileWriter& rhs);

    /** for addin - runs inputs through XMLFileWriter */
    static IObjectSP addinXMLFileWriter(XMLFileWriter* xmlFileWriter);
};

/** Creates a new string containing a given object as XML */ 
class ADDINS_DLL XMLStringWriter : public CObject,
                      virtual public ClientRunnable {
public:  
    static CClassConstSP const TYPE;
    friend class XMLStringWriterHelper;
    
    /** Does not clone object parameter */
    XMLStringWriter(IObjectSP object);

    // EdrAction 
    virtual IObjectSP run();

    IObjectSP run() const;

private:
    IObjectSP object;

    XMLStringWriter();
    XMLStringWriter(const XMLStringWriter& rhs);
    XMLStringWriter& operator=(const XMLStringWriter& rhs);

    /** for addin - runs inputs through XMLStringWriter */
    static IObjectSP addinXMLStringWriter(XMLStringWriter* xmlStringWriter);
};

DRLIB_END_NAMESPACE

#endif
