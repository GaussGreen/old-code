//----------------------------------------------------------------------------
//
//   Group       : Global Derivatives Research
//
//   Filename    : XMatrix.cpp
//
//   Description : Wrapper for external [DR interface] objects which are
//                 matrices
//
//   Author      : Mark A Robson
//
//   Date        : 2 Dec 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XMatrix.hpp"
#include "edginc/DRAnalyticsInterface.h"
#include "edginc/Format.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/Writer.hpp"
#include "edginc/DRUtil.hpp"

DRLIB_BEGIN_NAMESPACE

// DRI_CHECK is defined in DRUtil.hpp. This simplifies the invocation of the
// DRI_CHECK macro by fixing some parameters.
#define CHECK(exec, method) \
    DRI_CHECK((exec), (method), 0, driErrorHandler, 0, ;, ;, ;)

/** tag for recording number of columns in a matrix */
static const string NUM_COLUMNS = "cols";
/** tag for recording number of rows in a matrix */
static const string NUM_ROWS = "rows";

XMatrix::~XMatrix(){}


XMatrix::XMatrix(): XObject(TYPE){}

/** Are these objects equal (ie contain the same data) */
bool XMatrix::equals(IDRObjectConstSP drObj) const{
    if (this == drObj.get()){
        return true;
    }
    if (drObj->getClass() != TYPE){
        return false;
    }
    const XMatrix* theObj = STATIC_CAST(XMatrix, drObj.get());
    int numCols, numRows;
    size(numCols, numRows);
    int numCols2, numRows2;
    theObj->size(numCols2, numRows2);
    if (numCols != numCols || numRows != numRows2){
        return false;
    }
    for (int i = 0; i < numCols; i++){
        for (int j = 0; j < numRows; j++){
            double db = get(i, j);
            if (db != theObj->get(i, j)){
                return false;
            }
        }
    }
    return false;
}

/** If we own the DRObject then nothing is done (NB May need to review this
    since an XMatrix can be modified). Otherwise
    copy the object - by building empty array and getting/setting
    (which is a bit weak but there is no clone method) */
IObject* XMatrix::clone() const{
    if (/*object.*/owns()){
        return const_cast<IObject*>((const IObject*)this);
    }
    XObjectSP newMatrix(recursiveClone());
    return newMatrix.release();
}

/** Create a clone of this object by recursively cloning every component
    of this object. This is used for regression testing and should not
    be used in general (unnecessarily slow) */
XObjectSP XMatrix::recursiveClone() const{
    // get size
    int numCols, numRows;
    size(numCols, numRows);
    XMatrixSP xMatrix(create(numCols, numRows, object->svc));
    for (int i = 0; i < numCols; i++){
        for (int j = 0; j < numRows; j++){
            xMatrix->set(i, j, get(i, j));
        }
    }
    return xMatrix;
}

/** Override default CObject implementation.  */
void XMatrix::xWrite(const string& tag, Writer* writer) const{
    try {
        // get size
        int numCols, numRows;
        size(numCols, numRows);
        // write external type as an attribute
        string attribute = XTYPE_ATTRIBUTE+"='"+ getXClassName() + "' " +
            NUM_COLUMNS+"='"+Format::toString(numCols)+"' "+
            NUM_ROWS+"='"+Format::toString(numRows)+"'";
        IObjectConstSP obj(writer->objectStart(tag, attribute, this, true));
        if (obj.get()){ // if not already written out
            char buffer[40];
            char value[512];
            for (int i = 0; i < numCols; i++){
                sprintf(buffer, "%s%d", "Column", i);
                writer->objectStart(buffer, "", 0, true);
                for (int j = 0; j < numRows; j++){
                    sprintf(value, "%.16f\n", get(i, j));
                    writer->write(value);
                }
                writer->objectEnd(buffer, 0);
            }
        }
        // all done. Just close object
        writer->objectEnd(tag, this);
    } catch (exception& e){
        throw ModelException(e, "XMatrix::xWrite");
    }
}
/** Override default CObject implementation */
void XMatrix::xImport(Reader::Node* elem, Reader* reader, DRService* svc){
    static const string method("XMatrix::xImport");
    try{
        // first get size of matrix
        string colsAsString = elem->attribute(NUM_COLUMNS);
        string rowsAsString = elem->attribute(NUM_ROWS);
        int numCols = atoi(colsAsString.c_str());
        int numRows = atoi(rowsAsString.c_str());
        XMatrixSP matrix(create(numCols, numRows, svc));
        int currentCol = 0;
        Reader::NodeListSP  nl(elem->children());
        for (unsigned int i = 0; i < nl->size(); i++) {
            Reader::Node* node = (*nl)[i].get();
            if (currentCol >= numCols){
                throw ModelException(method, "Too many columns in matrix");
            }
            // move to the next field's XML description
            string dataLines = node->value();
            const char* lines = dataLines.c_str();
            const char* unconverted = lines;
            int   idx = 0;
            bool  morenumbers;
            do {
                morenumbers = false;
                for (unsigned int j = 0;
                     j < (lines-unconverted+strlen(lines)) && !morenumbers;
                     j++) {
                    morenumbers = !isspace(*(unconverted+j));
                }

                if (morenumbers) {
                    if (idx >= numRows) {
                        throw ModelException(method,  "matrix has "+
                                             Format::toString(numRows)+
                                             " rows, but read in " +
                                             Format::toString(idx));
                    }
                    matrix->set(currentCol, idx,
                                strtod(unconverted, (char **)&unconverted));
                    idx++;
                }
            } while (morenumbers);

            if (idx < numRows) {
                throw ModelException(method,
                                     "matrix has "+
                                     Format::toString(numRows)+
                                     " rows, but read in " +
                                     Format::toString(idx));
            }
            currentCol++;
        }
        // and copy pointer over
        this->object = matrix->object;
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Override default CObject implementation */
void XMatrix::outputWrite(const string& linePrefix,
                          const string& prefix,
                          ostream&      stream) const{
    // get size
    int numCols, numRows;
    size(numCols, numRows);
    if (numCols == 0 || numRows == 0){
        stream << linePrefix << prefix << ": empty" << endl;
    } else {
        for (int i = 0; i < numCols; i++){
            for (int j = 0; j < numRows; j++){
                CDoubleSP db(CDouble::create(get(i, j)));
                db->outputWrite(linePrefix, prefix +
                                Format::toString("[%d][%d]", i, j), stream);
            }
        }
    }
}

/** Simple constructor. */
XMatrix::XMatrix(const MyDRObjectInterfaceSP& matrix):
    XObject(TYPE, matrix){}

/** Simple constructor. Takes ownership of memory */
XMatrix::XMatrix(DRMatrix matrix): XObject(TYPE, matrix){}

/** CombinableResult interface: scale by factor x */
void XMatrix::scale(double x){
    int numCols, numRows;
    size(numCols, numRows);
    for (int i = 0; i < numCols; i++){
        for (int j = 0; j < numRows; j++){
            double d = get(i, j);
            set(i, j, d * x);
        }
    }
}

/** CombinableResult interface:
    add an object (scaled by scaleFactor) to this
    result. Implementations should modify this result. If the x is
    not the same type as this then a [class cast] exception will
    be thrown */
void XMatrix::add(const CombinableResult& x, double scaleFactor){
    const XMatrix* xToAdd = DYNAMIC_CAST(XMatrix, &x);
    int numCols, numRows, numCols2, numRows2;
    size(numCols, numRows);
    xToAdd->size(numCols2, numRows2);
    if (numCols != numCols2 || numRows != numRows2){
        throw ModelException("XMatrix::add",
                             "Matrices are of different dimensions");
    }
    for (int i = 0; i < numCols; i++){
        for (int j = 0; j < numRows; j++){
            double d1 = get(i, j);
            double d2 = xToAdd->get(i, j);
            set(i, j, d1 + d2 * scaleFactor);
        }
    }
}


#if 0
CDoubleMatrixSP XMatrix::toDoubleMatrix() const{
    int numCols, numRows;
    size(numCols, numRows);
    CDoubleMatrixSP matrix(new DoubleMatrix(numCols, numRows));
    for (int i = 0; i < numCols; i++){
        for (int j = 0; j < numRows; j++){
            (*matrix)[i][j] = get(i, j);
        }
    }
    return matrix;
}

XMatrixSP XMatrix::fromDoubleMatrix(const CDoubleMatrixSP& matrix,
                                    DRService*             svc){
    int numCols = matrix->numCols();
    int numRows = matrix->numRows();
    XMatrixSP xMatrix(create(numCols, numRows, svc));
    for (int i = 0; i < numCols; i++){
        for (int j = 0; j < numRows; j++){
            xMatrix->set(i, j, (*matrix)[i][j]);
        }
    }
    return xMatrix;
}
#endif

/** Determines size of matrix */
void XMatrix::size(int& numCols, int& numRows) const{
    DRService* svc = object->svc;
    CHECK((*svc->fptr->matrixSize)(svc, (DRObject)object.get(),
                                   &numCols, &numRows),
          "XMatrix::matrixSize");
}


/** Retrieves element of matrix */
double XMatrix::get(int col, int row) const{
    DRService* svc = object->svc;
    double value;
    CHECK((*svc->fptr->matrixGet)(svc, (DRObject)object.get(),
                                  col, row, &value),
          "XMatrix::matrixGet");
    return value;
}

/** Sets element of matrix */
void XMatrix::set(int col, int row, double value){
    DRService* svc = object->svc;
    CHECK((*svc->fptr->matrixSet)(svc, (DRObject)object.get(),
                                  col, row, value),
          "XMatrix::matrixSet");
}

/** Creates empty matrix of specified size */
XMatrixSP XMatrix::create(int numCols, int numRows, DRService* svc){
    // create Array
    DRMatrix theMatrix;
    CHECK((*svc->fptr->matrixNew)(svc, numCols, numRows, &theMatrix),
          "XMatrix::createMatrix");
    return XMatrixSP(new XMatrix(theMatrix));
}

void XMatrix::load(CClassSP& clazz){
    REGISTER(XMatrix, clazz);
    SUPERCLASS(XObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

IObject* XMatrix::defaultConstructor(){
    return new XMatrix();
}

CClassConstSP const XMatrix::TYPE =
CClass::registerClassLoadMethod("XMatrix", typeid(XMatrix), load);



DRLIB_END_NAMESPACE
