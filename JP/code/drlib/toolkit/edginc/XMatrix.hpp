//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research
//
//   Filename    : XMatrix.hpp
//
//   Description : Wrapper for external [DR interface] objects which are
//                 matrices
//
//   Author      : Mark A Robson
//
//   Date        : 12 Nov 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDR_XMATRIX_HPP
#define EDR_XMATRIX_HPP
#include "edginc/XObject.hpp"

DRLIB_BEGIN_NAMESPACE
#if 0
class CDoubleMatrix;
typedef smartPtr<CDoubleMatrix> CDoubleMatrixSP;
#endif

//// Wrapper for external [DR interface] objects which are matrices
class TOOLKIT_DLL XMatrix: public XObject{
public:
    static CClassConstSP const TYPE;
    virtual ~XMatrix();

    /** Creates empty matrix of specified size */
    static XMatrixSP create(int numCols, int numRows, DRService* svc);

    /** Simple constructor. */
    explicit XMatrix(const MyDRObjectInterfaceSP& matrix);

    /** Simple constructor. Takes ownership of memory */
    explicit XMatrix(DRMatrix matrix);

    /** Are these objects equal (ie contain the same data) */
    virtual bool equals(IDRObjectConstSP drObj) const;

    /** If we own the DRObject then nothing is done (NB May need to review this
        since an XMatrix can be modified). Otherwise
        copy the object - by building empty array and getting/setting
        (which is a bit weak but there is no clone method) */
    IObject* clone() const;

    /** Override default CObject implementation */
    virtual void xWrite(const string& tag, Writer* writer) const;

    /** Override default CObject implementation */
    virtual void xImport(Reader::Node* elem, Reader* reader, DRService* svc);

    /** Override default CObject implementation */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix,
                             ostream&      stream) const;
    /** CombinableResult interface: scale by factor x */
    virtual void scale(double x);

    /** CombinableResult interface: 
        add an object (scaled by scaleFactor) to this
        result. Implementations should modify this result. If the x is
        not the same type as this then a [class cast] exception will
        be thrown */
    virtual void add(const CombinableResult& x, double scaleFactor);

    /** Create a clone of this object by recursively cloning every component
        of this object. This is used for regression testing and should not
        be used in general (unnecessarily slow) */
    virtual XObjectSP recursiveClone() const;

    /** Determines size of matrix */
    void size(int& numCols, int& numRows) const;

    /** Retrieves element of matrix  */
    double get(int col, int row) const;

    /** Sets element of matrix  */
    void set(int col, int row, double value);


private:
    XMatrix(const XMatrix& rhs); // don't use
    XMatrix& operator=(const XMatrix& rhs); // don't use
    static void load(CClassSP& clazz);
    XMatrix();
    static IObject* defaultConstructor();
#if 0
    static XMatrixSP fromDoubleMatrix(const CDoubleMatrixSP& matrix,
                                      DRService*             svc);
    
    CDoubleMatrixSP toDoubleMatrix() const;
#endif
};


DRLIB_END_NAMESPACE
#endif
