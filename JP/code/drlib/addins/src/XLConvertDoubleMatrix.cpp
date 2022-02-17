//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLConvertDoubleMatrix.cpp
//
//   Description : Class for converting double matrix to/from Excel
//
//   Author      : Mark A Robson
//
//   Date        : 26 Feb 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XLConvertFactory.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/Format.hpp"
DRLIB_BEGIN_NAMESPACE


/** Specialisation of XLConvert for handling doubleMatrices */
class XLConvertDoubleMatrix: public XLConvert{
public:

    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP  desiredType,        // what type we want
        const XL_OPER&       oper,         // oper to convert
        const CFieldConstSP& field,        // where to store the converted data
        jobject&             output) const{ // store value in here using field

        static const string routine = "XLConvertDoubleMatrix::convertInput";
        int             idx, jdx;
        int             rows, cols;
        const XL_OPER*     operArray;
        if (oper.type & xltypeMulti){
            rows  = oper.val.xlarray.rows;
            cols  = oper.val.xlarray.columns;
            operArray = oper.val.xlarray.lparray;
        } else {
            rows = 1;
            cols = 1;
            operArray = &oper;
        }
        int             actRows, actCols;

        /* determine size of matrix - use top row, left most column */
        const XL_OPER* cell;
        idx = 0;
        do {
            cell = operArray + offset(idx, 0, cols, rows);
        } while (cell->type != xltypeNil && ++idx < cols);
        actCols = idx;
        jdx = 0;
        do {
            cell = operArray + offset(0, jdx, cols, rows);
        } while (cell->type != xltypeNil && ++jdx < rows);
        actRows = jdx;
        if (actCols < cols) {
            /* check that all cells in col are nil */
            for (idx = 1; idx < actRows; idx++) {
                if (!(operArray[offset(actCols, idx, 
                                       cols, rows)].type & xltypeNil)){
                    throw ModelException(
                        routine, "Matrix is not well formed\n"
                        "Expected cell ("+Format::toString(idx)+
                        ", "+Format::toString(actCols)+") to be empty");
                }
            }
        }
        if (actRows < rows)
        {
            /* check that all cells in row are nil */
            for (idx = 1; idx < actCols; idx++){
                if (!(operArray[offset(idx, actRows,
                                       cols, rows)].type & xltypeNil)){
                    throw ModelException(
                        routine, "Matrix is not well formed\n"
                        "Expected cell ("+Format::toString(idx)+
                        ", "+Format::toString(actCols)+") to be"
                        " empty");
                }
            }
        }
        /* Build matrix so that vols for the same strike form a column */
        CDoubleMatrixSP matrix(new DoubleMatrix(actCols, actRows));
        for (idx = 0; idx < actCols; idx++) {
            double* col = (*matrix)[idx];
            for (jdx = 0; jdx < actRows; jdx++) {
                cell = operArray + offset(idx, jdx, cols, rows);
                if (cell->type == xltypeNil) {
                    throw ModelException(routine, "Empty cell found at ("+
                                         Format::toString(idx)+", "+
                                         Format::toString(jdx)+") in matrix");
                }
                col[jdx] = operToDouble(*cell);
            }
        }
        
        field->set(output, matrix);
    }

    /** Take relevant value from data and fill column idx in output 
        accordingly. Default is to create a handle */
    virtual void convertObjectOutput(
        const string&  handleName,    // for handles
        const string&  name,          /* name for this bit of data
                                         (for handles/error messages) */
        const jobject& object,        // populate OPER with the object
        int            idx,           // column index number */
        XL_OPER&       oper) const{   // converted output

        if (oper.type & xltypeMulti){
            // if there are already multiple columns then something has gone
            // wrong
            throw ModelException("XLConvertDoubleMatrix::convertObjectOutput",
                                 "Internal error - trying to output matrix in"
                                 " expandMulti mode");
        }
        CDoubleMatrixSP matrix = CDoubleMatrixSP::dynamicCast(object);
        int numRows = 0;
        int numCols = 0;
        if (matrix.get()){
            // find out how big the matrix is
            numRows = matrix->numRows();
            numCols = matrix->numCols();
        }
        if (numCols > 0 && numRows > 0){
            // make sure the output oper is big enough
            setNumColsInOper(numCols, oper);
            setNumberOfRowsInOper(numRows, idx, oper);
            for (int kdx = 0; kdx < numCols; kdx++) {
                for (int jdx = 0; jdx < numRows; jdx++) {
                    XL_OPER *cell = oper.val.xlarray.lparray + 
                        offset(kdx, jdx, oper.val.xlarray.columns, 
                               oper.val.xlarray.rows);
                    cell->val.num = (*matrix)[kdx][jdx];
                    cell->type = xltypeNum;
                }
            }
        } else {
            oper.type = xltypeErr;
            oper.val.err = xlerrNA;
        }
    }

    /** Returns the width needed to dsiplay this type of object in
        XL. A return value of -1 indicates a variable width. This is
        used for addins that use Addin:expandMulti */
    virtual int width() const{
        // variable width
        return -1;
    }

    bool isSimpleType() const {
        return true;
    }
};

/** non static method which is used by XLConvertFactory to ensure that
    this class gets linked in and the method is called after
    CDoubleMatrix::TYPE is initialised. Outside of class to avoid necessity
    of header file */
void XLConvertDoubleMatrixRegister(){
    XLConvertFactory::registerXLObjectConvert(DoubleMatrix::TYPE,
                                              new XLConvertDoubleMatrix());
}

DRLIB_END_NAMESPACE
