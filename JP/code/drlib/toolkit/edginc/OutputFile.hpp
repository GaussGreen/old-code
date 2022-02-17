//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : OutputFile.hpp
//
//   Description : Responsible for creation of .out files
//
//   Author      : Mark A Robson
//
//   Date        : 8 Feb 2001
//
//----------------------------------------------------------------------------

#ifndef OUTPUTFILE_HPP
#define OUTPUTFILE_HPP

#include "edginc/AtomicArray.hpp"
#include "edginc/XMLWriter.hpp"

DRLIB_BEGIN_NAMESPACE

/** Controls the structure of the .out files produced on running a regression
    test */
class TOOLKIT_DLL OutputFile{
public:
    static const string OUTPUT_FILE_TAG;
    static const string OUTPUT_TAG;
    static const string SUMMARY_TAG;
    static const string RESULTS_TAG;
    static const string PERMUTATION;
    static const string LINE_PREFIX;
    static const string EXCEPTION_TAG;

    /** Open output file, ready for writing results to it */
    OutputFile(const string& filename);

    void write(ModelException& e);

    /** Writes object to the stream in the standard output format
        together with a string array which represents any exceptions
        thrown during the regression test run. The string array is
        written out within a <Exception> .. </Exception> block. If the
        string array is empty this block is skipped. Multiple objects
        can be written by repeated calls to this function. The
        supplied object can be null provided exceptions is not empty */
    void write(const IObject* object,
               const CStringArray& exceptions = CStringArray());

    /** Read an object back in from an XML serialised OutputFile */
    static IObjectSP xmlReadResults(const string& filename, int perm = 1);

    /** Read errors back in from an XML serialised OutputFile */
    static string xmlReadErrors(const string& filename, int perm = 1);

    ~OutputFile();

    /** Creates an output file name from an input file name */
    static string createOutputFileName(const string& inputFileName);

    string getFilename(void) const;

private:
    OutputFile(const OutputFile &rhs);
    OutputFile& operator=(const OutputFile& rhs);

    XMLWriter         writer;
    int               permNum;
    const string      filename; 
};

DRLIB_END_NAMESPACE

#endif
