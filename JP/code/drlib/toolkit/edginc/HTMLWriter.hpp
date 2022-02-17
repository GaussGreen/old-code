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

#ifndef _HTMLWRITER_HPP
#define _HTMLWRITER_HPP
#include <string>
#include <iostream>
#include "edginc/smartPtr.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

class CClass;

/** The HTMLWriter writes class documentation in HTML format */
class TOOLKIT_DLL HTMLWriter {
public:
    virtual ~HTMLWriter();

    /** start a comment */
    virtual void commentStart();
 
    /** end a comment */
    virtual void commentEnd();
    
    /** write free-form data */
    virtual void write(const string& data);

    /** make one */
    static HTMLWriter* htmlFile(const string& filename);
   
    /** document the class */
    void document(const CClass* clazz);

private:
    HTMLWriter();
    HTMLWriter(HTMLWriter&);

    void classStart(const CClass* clazz);
    void classEnd();

    ostream* stream;

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
};

typedef refCountPtr<HTMLWriter> HTMLWriterSP;

DRLIB_END_NAMESPACE
#endif
