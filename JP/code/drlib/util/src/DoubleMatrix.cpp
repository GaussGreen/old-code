
#include "edginc/config.hpp"
#define QLIB_DOUBLEMATRIX_CPP
#include "edginc/DoubleMatrix.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Writer.hpp"

DRLIB_BEGIN_NAMESPACE
/** EigenVector analysis */
CClassConstSP const EigenVectorAnalysis::TYPE = CClass::registerClassLoadMethod(
    "EigenVectorAnalysis", typeid(EigenVectorAnalysis), load);

EigenVectorAnalysis::EigenVectorAnalysis(): CObject(TYPE), dimension(0) {};

EigenVectorAnalysis::EigenVectorAnalysis(int n): CObject(TYPE), dimension(n) {
    eigenValues  = DoubleArraySP(new DoubleArray(n));
    eigenVectors = CDoubleMatrixSP(new DoubleMatrix(n, n));
};

IObject* EigenVectorAnalysis::defaultEigenVectorAnalysis(){
    return new EigenVectorAnalysis();
}

void EigenVectorAnalysis::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(EigenVectorAnalysis, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultEigenVectorAnalysis);
    FIELD(eigenVectors, "Eigen vectors");
    FIELD(eigenValues, "Eigen values");
    FIELD(dimension, "Dimension");
}


/** SVD analysis */
CClassConstSP const SVDAnalysis::TYPE = CClass::registerClassLoadMethod(
    "SVDAnalysis", typeid(SVDAnalysis), load);

SVDAnalysis::SVDAnalysis(): CObject(TYPE) {};

SVDAnalysis::SVDAnalysis(int numCols, int numRows):
CObject(TYPE) {
    matrixU          = CDoubleMatrixSP(new DoubleMatrix(numCols, numRows));
    matrixVTranspose = CDoubleMatrixSP(new DoubleMatrix(numCols, numCols));
    sValues          = CDoubleArraySP(new CDoubleArray(numCols));
};

IObject* SVDAnalysis::defaultSVDAnalysis(){
    return new SVDAnalysis();
}

void SVDAnalysis::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(SVDAnalysis, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultSVDAnalysis);
    FIELD(matrixU, "matrixU");
    FIELD(matrixVTranspose, "matrixV transpose");
    FIELD(sValues, "Singular Values");
}

/** Beta factorization */
CClassConstSP const BetaFactorization::TYPE = CClass::registerClassLoadMethod(
    "BetaFactorization", typeid(BetaFactorization), load);

BetaFactorization::BetaFactorization(): CObject(TYPE) {};

BetaFactorization::BetaFactorization(double             avgDistance, 
                                     CDoubleMatrixSP    betaFactors, 
                                     CDoubleMatrixSP    betaMatrix,
                                     double             squeezeLimitLow,
                                     double             squeezeLimitHigh) :
CObject(TYPE), 
avgDistance(avgDistance), 
betaFactors(betaFactors), 
betaMatrix(betaMatrix),
squeezeLimitLow(squeezeLimitLow),
squeezeLimitHigh(squeezeLimitHigh) {};

IObject* BetaFactorization::defaultBetaFactorization(){
    return new BetaFactorization();
}

void BetaFactorization::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(BetaFactorization, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultBetaFactorization);
    FIELD(avgDistance, "average distance");
    FIELD(betaFactors, "beta factors");
    FIELD(betaMatrix, "beta matrix");
    FIELD(squeezeLimitLow, "lower squeeze limit");
    FIELD(squeezeLimitHigh, "upper squeeze limit");
}

/** tag for recording number of columns in a matrix */
const string CDoubleMatrix::NUM_COLUMNS = "cols";
/** tag for recording number of rows in a matrix */
const string CDoubleMatrix::NUM_ROWS = "rows";

CoeffsPQ4Pade CDoubleMatrix::PADEcoeffs = CoeffsPQ4Pade(14);

/* should possibly make this matrix inherit from CArray - so it would appear
   an an array of arrays */

CDoubleMatrix::CDoubleMatrix():
    CObject(TYPE), NumCols(0),  NumRows(0), data(0){
    // empty
}

CDoubleMatrix::CDoubleMatrix(int numCols, int numRows):
    CObject(TYPE), NumCols(numCols),  NumRows(numRows), data(0){
    allocateMem();
    fill(0.0);
}

CDoubleMatrix::CDoubleMatrix(
    int numCols, int numRows, const double** data):
    CObject(TYPE), NumCols(numCols), NumRows(numRows), data(0){
    allocateMem();
    populateMatrix(data);
}

/** Creates a 1xn matrix (i.e a curve) from a double array */
CDoubleMatrix::CDoubleMatrix(const DoubleArray& theArray):
    CObject(TYPE), NumCols(1), NumRows(theArray.size()), data(0)
{
    allocateMem();
    populateMatrix(theArray);
}

/** Creates a 1xn matrix (i.e a curve) from a double array */
CDoubleMatrix::CDoubleMatrix(int numCols, int numRows, const DoubleArray& theArray)
    : CObject(TYPE), NumCols(numCols), NumRows(numRows), data(0)
{
    allocateMem();
    populateMatrix(theArray);
}
/** Creates a n x m matrix from a 2-dimensional [n][m] double array */
CDoubleMatrix::CDoubleMatrix(const DoubleArrayArray& theArray):
CObject(TYPE), NumCols(theArray.size()), NumRows(0), data(0){
    static const string method = "CDoubleMatrix::CDoubleMatrix(const DoubleArrayArray&)";
    try{
        // check the array have columns of equal sizes
        if (NumCols > 0){
            NumRows = theArray[0].size();
            int iCol = 1;
            for (; iCol < NumCols; ++iCol){
                if (theArray[iCol].size() != NumRows){
                    throw ModelException(method,
                                         "the columns of the 2-dimensional array "
                                         "don't have equal sizes:\nthe size ("
                                         + Format::toString(theArray[iCol].size())
                                         + ") of the "
                                         + Format::toString(iCol + 1)
                                         + "-th column is different from the size ("
                                         + Format::toString(NumRows)
                                         + ") of the previous colums");
                }
            }
        }
        allocateMem();
        populateMatrix(theArray);
    } 
    catch (exception& e){
        this->~CDoubleMatrix();
        throw ModelException(e, method);
    }
}

CDoubleMatrix::CDoubleMatrix(const CDoubleMatrix& rhs):
    CObject(TYPE), NumCols(rhs.NumCols), NumRows(rhs.NumRows)
{
    allocateMem();
    populateMatrix(const_cast<const double**>(rhs.data));
}

void CDoubleMatrix::allocateMem(){
    data = new double*[NumCols];
    try{
        for (int i = 0; i < NumCols; i++){
            data[i] = new double[NumRows];
        }
    } catch (exception& e){
        this->~CDoubleMatrix();
        throw ModelException(&e, "CDoubleMatrix::allocateMem");
    }
}

void CDoubleMatrix::deallocateMem(){
    if (data){
        for (int i = 0; i < NumCols; i++){
            delete[] data[i];
        }
        delete[] data;
        data = 0;
    }
}


void CDoubleMatrix::populateMatrix(const double** rawData){
    if (!rawData && !empty()){
        throw ModelException("CDoubleMatrix::populateMatrix","NULL data");
    }
    for (int i = 0; i < NumCols; i++){
        if (!rawData[i]){
            throw ModelException("CDoubleMatrix::populateMatrix",
                                 "NULL column");
        }
        for (int j = 0; j < NumRows; j++)
            data[i][j] = rawData[i][j];
    }
}

void CDoubleMatrix::populateMatrix(const DoubleArray& theArray)
{
    if (theArray.size() == 0)
    {
        throw ModelException("CDoubleMatrix::populateMatrix",
                             "array is of zero length");
    }

    for (int iCol = 0; iCol < NumCols; iCol++){
        for (int iRow = 0; iRow < NumRows; iRow++){
            data[iCol][iRow] = theArray[iCol * NumRows + iRow];
        }
    }
}

void CDoubleMatrix::populateMatrix(const DoubleArrayArray& theArray){
    for (int iCol = 0; iCol < NumCols; iCol++){
        for (int iRow = 0; iRow < NumRows; iRow++){
            data[iCol][iRow] = theArray[iCol][iRow];
        }
    }
}

CDoubleMatrix::~CDoubleMatrix(){
    deallocateMem();
}

CDoubleMatrix& CDoubleMatrix::operator=(const CDoubleMatrix& rhs){
    if (NumRows != rhs.NumRows || NumCols != rhs.NumCols){
        deallocateMem();
        NumRows = rhs.NumRows;
        NumCols = rhs.NumCols;
        allocateMem();
    }
    populateMatrix(const_cast<const double**>(rhs.data));
    return *this;
}

IObject* CDoubleMatrix::clone() const{
    return new CDoubleMatrix(NumCols,NumRows, 
                             const_cast<const double**>(data));
}

/** adds given value to all elements of matrix */
void CDoubleMatrix::scalarAdd(double val){
    for (int i = 0; i < NumCols; i++) {
        double* col = data[i];
        for (int j = 0; j < NumRows; j++) {
            col[j] += val;
        }
    }
}

/** adds given value to all elements of given row in matrix */
void CDoubleMatrix::rowAdd(int rowIndex, double val){
    if (rowIndex < 0 || rowIndex >= NumRows){
        throw ModelException("CDoubleMatrix::rowAdd", 
                             "Index "+Format::toString(rowIndex)+" out "
                             " of bounds");
    }
    for (int j = 0; j < NumCols; j++) {
        data[j][rowIndex] += val;
    }
}

/** multiplies all elements of given row in matrix by val */
void CDoubleMatrix::rowMultiply(int rowIndex, double val){
    if (rowIndex < 0 || rowIndex >= NumRows){
        throw ModelException("CDoubleMatrix::rowMultiply", 
                             "Index "+Format::toString(rowIndex)+" out "
                             " of bounds");
    }
    for (int j = 0; j < NumCols; j++) {
        data[j][rowIndex] *= val;
    }
}

/** adds given value to all elements of given column in matrix */
void CDoubleMatrix::colAdd(int colIndex, double val){
    if (colIndex < 0 || colIndex >= NumCols){
        throw ModelException("CDoubleMatrix::colAdd", 
                             "Index "+Format::toString(colIndex)+" out "
                             " of bounds");
    }
    double* col = data[colIndex];
    for (int j = 0; j < NumRows; j++) {
        col[j] += val;
    }
}

/** adds given value to all elements of given column in matrix */
void CDoubleMatrix::colMultiply(int colIndex, double val){
    if (colIndex < 0 || colIndex >= NumCols){
        throw ModelException("CDoubleMatrix::colMultiply", 
                             "Index "+Format::toString(colIndex)+" out "
                             " of bounds");
    }
    double* col = data[colIndex];
    for (int j = 0; j < NumRows; j++) {
        col[j] *= val;
    }
}

/** Hash code function */
int CDoubleMatrix::hashCode() const {
    int hcode = 0;
    for (int i = 0; i < NumCols; i++){
        double* col = data[i];
        for (int j = 0; j < NumRows; j++){
            hcode ^= CDouble::hashCode(col[j]);
        }
    }
    return hcode;
}

/** Comparison function */
bool CDoubleMatrix::equalTo(const IObject* obj) const {
    if (this == obj) { // Obvious first test
        return true;
    }
    if (!obj || obj->getClass() != getClass()) {
        return false;
    }

    const CDoubleMatrix* dMatrix = STATIC_CAST(CDoubleMatrix, obj);

    if ((NumCols != dMatrix->NumCols) || (NumRows != dMatrix->NumRows)) {
        return false;
    }

    for (int i = 0; i < NumCols; i++){
        double* col1 = data[i];
        double* col2 = dMatrix->data[i];
        for (int j = 0; j < NumRows; j++){
            if (col1[j] != col2[j]) {
                return false;
            }
        }
    }
    return true;
}

// build from an array of double arrays
IObjectSP CDoubleMatrix::fromArray(const IObjectSP& object, 
                                   CClassConstSP    requiredType) {
    static const string method("CDoubleMatrix::fromArray");
    try {
        DoubleArrayArray& dbls = dynamic_cast<DoubleArrayArray&>(*object);

        // first check that all the arrays are the same length
        int i;
        if (dbls.size() == 0) {
            throw ModelException(method,
                                 "input array has zero size");
        }
    
        for (i = 1; i < dbls.size(); i++) {
            if (dbls[i].size() != dbls[i-1].size()) {
                throw ModelException(method,
                                     "input arrays have differing lengths");
            }
        }

        CDoubleMatrixSP matrix(new DoubleMatrix(dbls.size(), dbls[0].size()));

        for (i = 0; i < dbls.size(); i++) {
            for (int j = 0; j < dbls[0].size(); j++) {
                (*matrix)[i][j] = dbls[i][j];
            }
        }
            
        return matrix;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

CClassConstSP const CDoubleMatrix::TYPE = CClass::registerClassLoadMethod(
    "DoubleMatrix", typeid(CDoubleMatrix), load);

static IObject* defaultCDoubleMatrix(){
    return new CDoubleMatrix();
}

void CDoubleMatrix::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CDoubleMatrix, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultCDoubleMatrix);
    registerObjectFromArrayMethod(DoubleArrayArray::TYPE,
                                  TYPE,
                                  &fromArray);
}

   

//double* CDoubleMatrix::operator[](int colIndex){// Inline for performance 
//    return data[colIndex];
//}

// Inline for performance 
//  const double* CDoubleMatrix::operator[](int colIndex) const{
//      return data[colIndex];
//  }

//int CDoubleMatrix::numCols() const{ // Inline for performance 
//    return NumCols;
//}

// int CDoubleMatrix::numRows() const{ // Inline for performance 
//     return NumRows;
// }

/** Returns true if the double matrix is empty ie 0x0 dimensions */
bool CDoubleMatrix::empty() const{
    return (NumCols == 0 && NumRows == 0);
}

/** Returns true if the double matrix is the identity matrix */
bool CDoubleMatrix::isIdentity() const {
    if(NumCols != NumRows) {
        return false;
    }    
    for (int i = 0; i < NumCols; i++) {
        for (int j = 0; j < NumCols; j++) {
            if (!Maths::isZero(data[i][j]-int(i == j))) {
                return false;
            }
        }
    }
    return true;
}

/** write object out to writer */
void CDoubleMatrix::write(const string& tag, Writer* writer) const {
    try {    
        char buffer[256];
        sprintf(buffer, "%s='%d' %s='%d'",
                NUM_COLUMNS.c_str(), NumCols,
                NUM_ROWS.c_str(), NumRows);
        IObjectConstSP obj(writer->objectStart(tag, buffer, this, true));
        if (obj.get()){
            char value[256];
            for (int i = 0; i < NumCols; i++){
                sprintf(buffer, "%s%d", "Column", i);
                writer->objectStart(buffer, "", 0, true);
                for (int j = 0; j < NumRows; j++){
                    sprintf(value, "%.16f\n", data[i][j]);
                    writer->write(value);           
                }
                writer->objectEnd(buffer, 0);
            }
        }
        writer->objectEnd(tag, this);
    }
    catch (exception& e){
        throw ModelException(&e, "CDoubleMatrix::write");
    }


}

/** populate an empty object from reader */
void CDoubleMatrix::import(Reader::Node* elem, Reader* reader) {
    const static string routine = "CDoubleMatrix::import";
    try{
        // first get size of matrix
        string colsAsString = elem->attribute(NUM_COLUMNS);
        string rowsAsString = elem->attribute(NUM_ROWS);
        NumCols = atoi(colsAsString.c_str());
        NumRows = atoi(rowsAsString.c_str());
        allocateMem();
        int currentCol = 0;
        Reader::NodeListSP  nl(elem->children()); 
        for (unsigned int i = 0; i < nl->size(); i++) {
            Reader::Node* node = (*nl)[i].get();
            if (currentCol >= NumCols){
                throw ModelException(routine, 
                                     "Too many columns in matrix");
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
                    if (idx >= NumRows) {
                        throw ModelException(routine,  "matrix has "+
                                             Format::toString(NumRows)+
                                             " rows, but read in " + 
                                             Format::toString(idx));
                    }
                    data[currentCol][idx] = 
                        strtod(unconverted, (char **)&unconverted);
                    idx++;
                }
            } while (morenumbers);
            
            if (idx < NumRows) {
                throw ModelException(routine,
                                     "matrix has "+
                                     Format::toString(NumRows)+
                                     " rows, but read in " + 
                                     Format::toString(idx));                    
            }
            currentCol++;
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void DoubleMatrix::outputWrite(const string& linePrefix,
                               const string& prefix, ostream& stream) const{
    if (empty()){
        stream << linePrefix << prefix << ": empty" << endl;
    } else {
        ios::fmtflags oldFlags = stream.flags(); // save settings
        long oldPrecision = stream.precision(); // save settings
        for (int i = 0; i < NumCols; i++){
            for (int j = 0; j < NumRows; j++){
                double db = data[i][j];
                if (db < 1 && db > -1){
                    // 8 dp
                    stream.flags(ios::fixed);
                } else {
                    // 8 sf
                    stream.unsetf(ios::fixed);
                }        
                stream.precision(8);
        
                stream << linePrefix << prefix << 
                    "[" << i << "][" << j << "]: " << db << endl;
            }
        }
        stream.flags(oldFlags); // restore settings
        stream.precision(oldPrecision);
    }
}

/** resize the matrix to a new number of columns and rows.
    do nothing with the values. */   
void CDoubleMatrix::resizeOnly(int numCols, int numRows)
{
    if (NumRows != numRows || NumCols != numCols)
    {
        deallocateMem();
        NumRows = numRows;
        NumCols = numCols;
        allocateMem();
    }

}

/** Resize the matrix to a new number of columns and rows. 
    Fill it with 0s if the size has been changed */
void CDoubleMatrix::resize(int numCols, int numRows)
{
    if (NumRows != numRows || NumCols != numCols)
	{
        deallocateMem();
        NumRows = numRows;
        NumCols = numCols;
        allocateMem();
		fill(0.0);
    }
}

/** set all the element of the matrix to a given number */
void CDoubleMatrix::fill(double val)
{
   for (int i = 0; i < NumCols; i++) {
        double* col = data[i];
        for (int j = 0; j < NumRows; j++) {
            col[j] = val;
        }
    }
}

/** resize the matrix to a new number of columns and rows. 
    If the size was initially the same, just set all elements to val*/
void CDoubleMatrix::resizeAndFill(int numCols, int numRows, double val)
{
	resizeOnly(numCols, numRows);
    fill(val);
}

/** checks all values are >= 0 */
void CDoubleMatrix::checkNonNegative() const{
    for (int i = 0; i < NumCols; i++){
        for (int j = 0; j < NumRows; j++){
            if (data[i][j] < -DBL_EPSILON)
            {
                string m("element ["+Format::toString(j)+
                         "]["+Format::toString(i)+"] is < 0.0 ");
                throw ModelException("CDoubleMatrix::checkNonNegative", m);
                                     
            }
        }
    }
}

/** checks all values are > 0 */
void CDoubleMatrix::checkPositive() const{
    for (int i = 0; i < NumCols; i++){
        for (int j = 0; j < NumRows; j++){
            if (data[i][j] < DBL_EPSILON)
            {
                string m("element ["+Format::toString(j)+
                         "]["+Format::toString(i)+"] is < 0.0 ");
                throw ModelException("CDoubleMatrix::checkNonNegative", m);
                                     
            }
        }
    }
}

/** checks matrix is symmetric */
void CDoubleMatrix::checkSymmetric() const{
    static const string routine("CDoubleMatrix::checkSymmetric");
    if (NumCols != NumRows){
        throw ModelException(routine, "Matrix is "+Format::toString(NumCols)+
                             " x "+Format::toString(NumRows));
    }
    for (int i = 0; i < NumCols; i++){
        for (int j = i+1; j < NumRows; j++){
            // somewhat arbitrary tolerance
            if (!Maths::areEqualWithinTol(data[i][j], data[j][i], 
                                          0.00000000001)){
                throw ModelException( routine,
                                     "Matrix (element ["+Format::toString(i)+
                                     "]["+Format::toString(j)+"]) is not "
                                     "symmetric");
            }
        }
    }
}

/** transpose (even for m*n matrix) */
void CDoubleMatrix::transpose() {
    CDoubleMatrix trans(NumRows, NumCols);

    for (int i = 0; i < NumCols; i++) {
        for (int j = 0; j < NumRows; j++) {
            trans[j][i] = data[i][j];
        }
    }
    *this = trans;
}

/** squeeze, ie compute rho + (1-rho)*s, s>0, resp rho(1+s), s<0 */
void CDoubleMatrix::squeeze(double s) {
    /** some validation */
    static const char* method = "CDoubleMatrix::squeeze";
    if (Maths::isPositive(abs(s)-1.0)) {
        throw ModelException(method, "Squeeze must not exceed one in absolute value"
                                     "but equals " + Format::toString(s));
    }
    
    if (Maths::isZero(s)) {
        return;
    }

    CDoubleMatrix squeezed(NumRows, NumCols);      
    if (Maths::isPositive(s)) {
        for (int i = 0; i < NumCols; i++) {
            for (int j = 0; j < NumRows; j++) {
                squeezed[i][j] = data[i][j] + (1.0-data[i][j])*s;
            }
        }
    } else if (Maths::isNegative(s)) {
        for (int i = 0; i < NumCols; i++) {
            for (int j = 0; j < NumRows; j++) {
                squeezed[i][j] = data[i][j]*(1.0+s);
            }
        }
    }
    *this = squeezed;
}

/** scale by factor x */
void CDoubleMatrix::scale(double x){
    for (int i = 0; i < NumCols; i++) {
        double* col = data[i];
        for (int j = 0; j < NumRows; j++) {
            col[j] *= x;
        }
    }
}

/** Negate all entries (equivalent to scale(-1.0) */
void DoubleMatrix::negate(){
    for (int i = 0; i < NumCols; i++) {
        double* col = data[i];
        for (int j = 0; j < NumRows; j++) {
            col[j] = -col[j];
        }
    }
}

/** add a CDoubleMatrix (scaled by scaleFactor) to this
    result.  */
void CDoubleMatrix::add(const CombinableResult& x, double scaleFactor){
    static const char* method = "CDoubleMatrix::add";
    if (scaleFactor==0.0) return;
    //throw ModelException("CDoubleMatrix::add", "Not yet done");
    const CDoubleMatrix* m = dynamic_cast<const CDoubleMatrix*>(&x);
    if (!m) {
        throw ModelException(method, "The object to be added must be castable to a CDoubleMatrix.");
    }
    if (m->NumRows!=NumRows || m->NumCols!=NumCols) {
        throw ModelException(method,  
            "To add two matrices, they must have the same row and column sizes.");
    }
    for (int r=0; r<NumRows; r++) {
        for (int c=0; c<NumCols; c++) {
            data[c][r] += scaleFactor * m->data[c][r];
        }
    }
    return;
}

/** Removes last column of matrix */
void CDoubleMatrix::removeLastCol(){
    if (NumCols > 0){
        NumCols--;
        delete[] data[NumCols];
    }
}

void CDoubleMatrix::insertCols(int start, int num) {
    if (start < 0 || NumCols < start) {
        throw ModelException(__FUNCTION__,
            "Cannot insert new cols before col " + Format::toString(start) + ": "
            "matrix only has " + Format::toString(NumCols) + " cols");
    }

    ASSERT(num >= 0);

    double** data_ = new double*[NumCols + num];
    for (int c = 0; c < start; ++c) data_[c] = data[c];
    for (int c = start; c < NumCols; ++c) data_[c + num] = data[c];
    for (int c = start; c < start + num; ++c) {
        data_[c] = new double[NumRows];
        for (int r = 0; r < NumRows; ++r) data_[c][r] = 0.;
    }

    delete[] data;
    data = data_;
    NumCols += num;
}

void CDoubleMatrix::insertRows(int start, int num) {
    if (start < 0 || NumRows < start) {
        throw ModelException(__FUNCTION__,
            "Cannot insert new rows before row " + Format::toString(start) + ": "
            "matrix only has " + Format::toString(NumRows) + " rows");
    }

    ASSERT(num >= 0);

    for (int c = 0; c < NumCols; ++c) {
        double* row = new double[NumRows + num];
        for (int r = 0; r < start; ++r) row[r] = data[c][r];
        for (int r = start; r < NumRows; ++r) row[r + num] = data[c][r];
        for (int r = start; r < start + num; ++r) row[r] = 0.;
        delete[] data[c];
        data[c] = row;
    }

    NumRows += num;
}

/** Removes numRowsToRemove rows from matrix starting at startRow */
void CDoubleMatrix::removeRows(int startRow, int numRowsToRemove) {
    if (startRow + numRowsToRemove > NumRows) {
        throw ModelException("CDoubleMatrix::removeRows",
            "Internal Error: Cannot remove " + Format::toString(numRowsToRemove) +
            " from matrix starting at row " +  Format::toString(startRow) +
            ". Insuficient rows in matrix");
    }
    if (numRowsToRemove > 0) {
        NumRows -= numRowsToRemove;
        for (int i = 0; i < NumCols; i++) {
            for (int j = startRow; j < NumRows; j++) {
                data[i][j] = data[i][j + numRowsToRemove];
            }
        }
    }
}

void CDoubleMatrix::mult(const CDoubleMatrix& x, 
                         const CDoubleMatrix& y, CDoubleMatrix * gResult){
    QLIB_VERIFY((x.NumCols == y.NumRows), 
                string("CDoubleMatrix::mult :") + 
                             "NumCols (" + Format::toString(x.NumCols) +
                             ") should equal multiplier NumRows (" +
                             Format::toString(y.NumRows) + ")");
    QLIB_VERIFY(!!gResult, "CDoubleMatrix::mult : Result points to NULL");
    gResult->resizeOnly(y.NumCols, x.NumRows); // all elements will be set to 0.0 below
    for (int iRow = 0; iRow < x.NumRows; iRow++) {
        for (int iCol = 0; iCol < y.NumCols; iCol++) {
            double & _result = gResult->data[iCol][iRow];
            _result = 0.0;
            for (int i = 0; i < x.NumCols; i++) { // x.NumCols == y.NumRows
                _result += x.data[i][iRow] * y.data[iCol][i];
            }
        }
    }
}

void CDoubleMatrix::scale(const CDoubleMatrix& x, 
                         double gFactor, CDoubleMatrix * gResult)
{
    QLIB_VERIFY(!!gResult, "CDoubleMatrix::scale : Result points to NULL");
    gResult->resizeOnly(x.NumCols, x.NumRows); // all elements will be set to 0.0 below
    for (int iRow = 0; iRow < x.NumRows; iRow++) {
        for (int iCol = 0; iCol < x.NumCols; iCol++) {
            gResult->data[iCol][iRow] = gFactor * x.data[iCol][iRow];
        }
    }
}
void CDoubleMatrix::scale(double gFactor, 
                          const CDoubleMatrix & x, CDoubleMatrix * gResult)
{
    CDoubleMatrix::scale(x, gFactor, gResult);
    // because x * gFactor = gFactor * x when gFactor is a double
}

CDoubleMatrix CDoubleMatrix::operator*(const CDoubleMatrix& x) const {
    if (NumCols != x.NumRows) {
        throw ModelException("CDoubleMatrix::mult", 
                             "NumCols (" + Format::toString(NumCols) +
                             ") should equal multiplier NumRows (" +
                             Format::toString(x.NumRows) + ")");
    }
    CDoubleMatrix prod(x.NumCols, NumRows);
    for (int iRow = 0; iRow < prod.NumRows; iRow++) {
        for (int iCol = 0; iCol < prod.NumCols; iCol++) {
            prod[iCol][iRow] = 0.0;
            for (int i = 0; i < NumCols; i++) {
                prod[iCol][iRow] += data[i][iRow] * x[iCol][i];
            }
        }
    }
    return prod;
}
CDoubleMatrix CDoubleMatrix::mult(const CDoubleMatrix& x) const {
    return (*this) * x;
}

CDoubleArray CDoubleMatrix::mult(const CDoubleArray &a) const {
    if (NumCols != a.size()) {
        throw ModelException("CDoubleMatrix::mult", 
                             "NumCols (" + Format::toString(NumCols) +
                             ") should equal multiplier array size (" +
                             Format::toString(a.size()) + ")");
    }
    CDoubleArray prod(NumRows, 0.0);
    for (int iRow = 0; iRow < NumRows; iRow++) {
        prod[iRow] = 0.0;
        for (int iCol = 0; iCol < NumCols; iCol++) {
            prod[iRow] += data[iCol][iRow] * a[iCol];
        }
    }
    return prod;
}

CDoubleArray CDoubleMatrix::mult(const CIntArray &a) const {
    if (NumCols != a.size()) {
        throw ModelException("CDoubleMatrix::mult", 
                             "NumCols (" + Format::toString(NumCols) +
                             ") should equal multiplier array size (" +
                             Format::toString(a.size()) + ")");
    }
    CDoubleArray prod(NumRows, 0.0);
    for (int iRow = 0; iRow < NumRows; iRow++) {
        prod[iRow] = 0.0;
        for (int iCol = 0; iCol < NumCols; iCol++) {
            prod[iRow] += data[iCol][iRow] * a[iCol];
        }
    }
    return prod;
}

CDoubleArray CDoubleMatrix::mult(const vector<int> &a) const {
    if (NumCols != int(a.size())) {
        throw ModelException("CDoubleMatrix::mult", 
                             "NumCols (" + Format::toString(NumCols) +
                             ") should equal multiplier array size (" +
                             Format::toString(a.size()) + ")");
    }
    CDoubleArray prod(NumRows, 0.0);
    for (int iRow = 0; iRow < NumRows; iRow++) {
        prod[iRow] = 0.0;
        for (int iCol = 0; iCol < NumCols; iCol++) {
            prod[iRow] += data[iCol][iRow] * a[iCol];
        }
    }
    return prod;
}

// Returns this - x. Dimensions must match.
CDoubleMatrix CDoubleMatrix::operator-(const CDoubleMatrix& x) const {
    if (NumCols != x.NumCols ||
        NumRows != x.NumRows) {
        throw ModelException("CDoubleMatrix::operator-", 
                             "Cannot add matrix dim (#cols=" +
                             Format::toString(NumCols) + ", #rows=" +
                             Format::toString(NumRows) + ") to matrix dim (#cols=" +
                             Format::toString(x.NumCols) + ", #rows=" +
                             Format::toString(x.NumRows) +")");
    }
    CDoubleMatrix diff(NumCols, NumRows);
    for (int iRow = 0; iRow < diff.NumRows; iRow++) {
        for (int iCol = 0; iCol < diff.NumCols; iCol++) {
            diff[iCol][iRow] = data[iCol][iRow] - x[iCol][iRow];
        }
    }
    return diff;
}

CDoubleMatrix CDoubleMatrix::minus(const CDoubleMatrix& x) const {
    return (*this) - x;
}

// Returns this + x. Dimensions must match.
CDoubleMatrix CDoubleMatrix::operator+(const CDoubleMatrix& x) const {
    if (NumCols != x.NumCols ||
        NumRows != x.NumRows) {
        throw ModelException("CDoubleMatrix::operator+", 
                             "Cannot add matrix dim (#cols=" +
                             Format::toString(NumCols) + ", #rows=" +
                             Format::toString(NumRows) + ") to matrix dim (#cols=" +
                             Format::toString(x.NumCols) + ", #rows=" +
                             Format::toString(x.NumRows) +")");
    }
    CDoubleMatrix sum(NumCols, NumRows);
    for (int iRow = 0; iRow < sum.NumRows; iRow++) {
        for (int iCol = 0; iCol < sum.NumCols; iCol++) {
                sum[iCol][iRow] = data[iCol][iRow] + x[iCol][iRow];
        }
    }
    return sum;
}
CDoubleMatrix CDoubleMatrix::plus(const CDoubleMatrix& x) const {
    return (*this) + x;
}


/** computes the square root of a matrix aka Choleski / Cholesky etc. */
// returns lower triangular matrix, i.e. root[iCol][iRow]=0 for iRow>iCol
CDoubleMatrix CDoubleMatrix::computeSquareRoot() const {
    static const string routine = "CDoubleMatrix::computeSquareRoot";

    try {
        checkSymmetric();
    
        int n = NumRows;
        CDoubleMatrix root(*this);

        vector<double*> ptrs(n);
        int i;
        for(i = 0; i < n; i++) {
            ptrs[i] = root[i] - 1;
        }
    
        DoubleArray diagonal(n);
        // NR computes the root of the matrix in its lower triangle
        // Because of the implementation of [] it will compute it in the upper triangle
        try {
            choldc(&ptrs[0] - 1, 
                   n, 
                   &diagonal[0] - 1);
        } catch(exception& e) {
            string message = "Failed to compute the Cholesky root of the matrix.";
            throw ModelException::addTextToException(e, message);
        }
    
        // kill the lower triangular part and fill in the diagonal
        for(int iCol = 0; iCol < n; iCol++) {
            root[iCol][iCol] = diagonal[iCol];
            for(int iRow = iCol + 1; iRow < n; iRow++) {
                root[iCol][iRow] = 0.0;
            }
        }
        return root;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

/** computes the inverse of a square matrix using LU */
CDoubleMatrix CDoubleMatrix::computeInverse() const
{
    CDoubleMatrix       yWrapper, aWrapper;
    DoublePointerVector yVector, a;
    DoubleVector        colWrapper;
    IntVector           indxWrapper;

    this->computeInverse(&yWrapper, &aWrapper, 
                         &yVector,  &a, &colWrapper, &indxWrapper);

    return yWrapper;
}

/** computes the inverse of a square matrix using LU. No allocation */
void CDoubleMatrix::computeInverse(CDoubleMatrix * yWrapper,  // this is the result
                                   CDoubleMatrix * aWrapper,
                                   DoublePointerVector * a,
                                   DoublePointerVector * yVector,
                                   DoubleVector        * colWrapper,
                                   IntVector           * indxWrapper) const
{
    static const string routine = "CDoubleMatrix::computeInverse";

    try {
        int N = NumRows;
        if( N<1 ) 
            throw ModelException("Matrix is empty");
        if( N!=NumCols ) 
            throw ModelException("Matrix is not square");

        a->resize(N);
        double **y, d, *col;
        int i,j,*indx;

        (*aWrapper) = (*this);
        yWrapper->resize(N,N);
        
        yVector->resize(N);
        for(i=0; i<N; i++) {
            (*a)[i] = (*aWrapper)[i] - 1;
            (*yVector)[i] = (*yWrapper)[i] - 1;
        }
        y = &(*yVector)[0] - 1;

        colWrapper->resize(N);
        col = &(*colWrapper)[0]-1;
        indxWrapper->resize(N);
        indx = &(*indxWrapper)[0]-1;

        ludcmp(&(*a)[0]-1,N,indx,&d); // Decompose the matrix just once.

        for(j=0;j<N;j++) {  // Determinant
            d *= (*aWrapper)[j][j];
        }
        
        /* This was bad: some matrices with very small determinants are
           totally invertible. Consider 1E-50 times the identity matrix.
           The determinant is only a problem for sure
           if it is actually zero, not just if it is small. 
           When it is very small it should be
           up to the users to determine if inversion is safe or not. 
           I have changed this as it was causing problems with Pade approximation
           coefficient computation.
           [Charles Morcom 5/5/2006] */
        //if( Maths::isZero(d) )
        //    throw ModelException("Matrix has zero determinant");
        if (d==0.0) throw ModelException("Matrix has zero determinant");

        for(j=1;j<=N;j++) {  // Find inverse by columns.
            for(i=1;i<=N;i++) col[i]=0.0;
            col[j]=1.0;
            lubksb(&(*a)[0]-1,N,indx,col);
            for(i=1;i<=N;i++) y[i][j]=col[i];
        }

    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

/** solves for x in Ax=b, with A square matrix and b=rhs */
CDoubleArray CDoubleMatrix::solve(CDoubleArray& rhs) const {
    static const string routine = "CDoubleMatrix::solve";

    try {
        int n = NumRows;
        if( (n!=NumCols) || (n<1) ) throw;

        vector<double*>  a(n);
        double *b,d;
        int *indx;

        CDoubleMatrix   aWrapper(*this);
        for(int i=0; i<n; i++) {
            a[i] = aWrapper[i] - 1;
        }
        CDoubleArray    bWrapper(rhs);
        b = &bWrapper[0]-1;
        vector<int>     indxWrapper(n);
        indx = &indxWrapper[0]-1;

        ludcmp(&a[0]-1,n,indx,&d);
        lubksb(&a[0]-1,n,indx,b);
 
        return bWrapper;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

/** input:  symmetric matrix 
    output: correlation matrix, average squared difference between symmetric and covariance matrices
    algorithm:
        1) floor the (real) eigenvalues at 0 -> matrix is covariance
        2) normalise covariance matrix -> matrix is correlation */
CDoubleMatrix CDoubleMatrix::symmToCorrel(double* squareErr, double eigenValueFloor) const {
    static const string routine = "CDoubleMatrix::symmToCorrel";

    try {
        if(!NumRows) {
            throw ModelException("Input matrix has no rows.");
        }

        if(!NumCols) {
            throw ModelException("Input matrix has no columns.");
        }
        
        checkSymmetric();
    
        int iAsset,jAsset,kAsset;
        const CDoubleMatrix& cov = *this; 
        CDoubleMatrix covNew(NumRows,NumRows); // in order to avoid confusion

        EigenVectorAnalysisSP symmEigen = cov.computeEigenVectors();
        bool negativeEigenValue = (Maths::isNegative((*symmEigen->eigenValues)[NumRows-1]-eigenValueFloor));
        double squareErrTmp = 0.0;

        // NEGATIVE EIGENVALUE => floor by eigenValueFloor and normalize new matrix
        if( negativeEigenValue ) {
            CDoubleMatrix temp(NumRows,NumRows);
            // 1) floor the (real) eigenvalues at eigenValueFloor
            for( iAsset=0; iAsset<NumRows; iAsset++ ) {
                for( jAsset=iAsset; jAsset<NumRows; jAsset++ ) {
                    temp[iAsset][jAsset] = 0.0;
                    for( kAsset=0; kAsset<NumRows; kAsset++ ) {
                        temp[iAsset][jAsset] 
                            +=  (*symmEigen->eigenVectors)[kAsset][iAsset]
                                * (*symmEigen->eigenVectors)[kAsset][jAsset]
                                * Maths::max(eigenValueFloor,(*symmEigen->eigenValues)[kAsset]);
                    }
                    if( jAsset>iAsset ) {
                        temp[jAsset][iAsset] = temp[iAsset][jAsset];
                    }
                }
                // check for non-positive variance
                if( !Maths::isPositive(temp[iAsset][iAsset]) ) {
                    string message = "Asset " + Format::toString(iAsset) 
                                     + " has variance " + Format::toString(temp[iAsset][iAsset]);
                    throw ModelException(message);
                }
            }
            // 2) normalise covariance matrix
            for( iAsset=0; iAsset<NumRows; iAsset++ ) {
                for( jAsset=iAsset+1; jAsset<NumRows; jAsset++ ) {
                    covNew[iAsset][jAsset] 
                        = temp[iAsset][jAsset] / sqrt(temp[iAsset][iAsset] * temp[jAsset][jAsset]);
                    covNew[jAsset][iAsset] = covNew[iAsset][jAsset];
                    squareErrTmp += (cov[iAsset][jAsset] - covNew[iAsset][jAsset])
                        * (cov[iAsset][jAsset] - covNew[iAsset][jAsset]);
                }
                covNew[iAsset][iAsset] = 1.0;
            }
            squareErrTmp /= ((NumRows*NumRows - NumRows) / 2.0);

        // NO NEGATIVE EIGENVALUE => just normalize input matrix (since potentially covar matrix)
        } else {
            for( iAsset=0; iAsset<NumRows; iAsset++ ) {
                // check for non-positive variance
                if( !Maths::isPositive(cov[iAsset][iAsset]) ) {
                    string message = "Asset " + Format::toString(iAsset) 
                                     + " has variance " + Format::toString(cov[iAsset][iAsset]);
                    throw ModelException(message);
                }
                for( jAsset=iAsset+1; jAsset<NumRows; jAsset++ ) {
                    covNew[iAsset][jAsset] 
                        = cov[iAsset][jAsset] / sqrt(cov[iAsset][iAsset] * cov[jAsset][jAsset]);
                    covNew[jAsset][iAsset] = cov[iAsset][jAsset];
                    squareErrTmp += (cov[iAsset][jAsset] - covNew[iAsset][jAsset])
                        * (cov[iAsset][jAsset] - covNew[iAsset][jAsset]);
                }
                covNew[iAsset][iAsset] = 1.0;
            }
            squareErrTmp /= ((NumRows*NumRows - NumRows) / 2.0);
        }

        (*squareErr) = squareErrTmp;
        return covNew;

    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}   

/**Frobenius norm of a matrix is sqrt(Tr(A'A)).*/
double CDoubleMatrix::frobeniusNorm() const {
    double norm = 0.0;
    for (int r=0; r<this->NumRows; r++) {
        for (int c=0; c<this->NumCols; c++) {
            double e = (*this)[c][r];
            norm += e*e;
        }
    }
    return sqrt(norm);
}

/**Spectral norm of a matrix is the square root of the modulus of the biggest eigenvalue of A'A*/
double CDoubleMatrix::spectralNorm() const{
    CDoubleMatrix thisSquared(*this);
    thisSquared.transpose();
    thisSquared = thisSquared*(*this);

    EigenVectorAnalysisSP ev = thisSquared.computeEigenVectors();
    double norm = 0.0;
    for (int i=0; i<ev->eigenValues->size(); i++) {
        double absEV = fabs((*(ev->eigenValues))[i]);
        if (absEV>norm) norm = absEV;
    }
    return sqrt(norm);
}

/** Same as expPade but result is allocated inside the method **/
CDoubleMatrix CDoubleMatrix::exp(double          gT, 
                        double          gSpectralNorm, 
                        int             gPadeOrder) const
{
    CDoubleMatrix _result, _aux4exp1, aux4exp11, aux4exp2, aux4exp21, aux4exp3;
    this->expPade(gT, gSpectralNorm, gPadeOrder, &_result, 
                        &_aux4exp1, &aux4exp11, &aux4exp2, &aux4exp21, &aux4exp3);
    return _result;
}

/**Compute the exponential exp(t M) of this (M) times a constant, t, using a Pade approximation 
       of specified order (numerator and denominator are the same order). If |Mt|>=0.5, then t is halved
       repeatedly until the norm is less than 0.5 to ensure that the Pade approximation is
       likely to be in range, and the result is squared the same number of times for the
       final answer. The argument for the spectral norm is there for efficiency: computing it
       is expensive enough that, for repeated exponentiations with different t values, it
       is better to compute it once only. */
void CDoubleMatrix::expPade(     double          gT, 
                                 double          gSpectralNorm, 
                                 int             gPadeOrder,
                                 CDoubleMatrix * gResult,
                                 CDoubleMatrix * gAux1,
                                 CDoubleMatrix * gAux11,
                                 CDoubleMatrix * gAux2,
                                 CDoubleMatrix * gAux21,
                                 CDoubleMatrix * gAux3) const 
{
    ///////////////////////////////////////////////////////////////////////////
    // CHECK INPUTS
    static const string method = "CDoubleMatrix::exp :";
    QLIB_VERIFY(gPadeOrder >= 0,
                method + "The Pade approximation order may not be negative.");
    QLIB_VERIFY(NumCols == NumRows,
                method + "NumRows!=NumCols: you may only compute "
                       + "the exponential of a square matrix.");
    QLIB_VERIFY(!!gResult,
                method + "gResult points to NULL");

    ///////////////////////////////////////////////////////////////////////////
    // DO ALL ALLOCATION- and RESIZING-RELATED OPERATIONS

    gResult->resizeOnly(NumCols, NumCols);

    // If user provides gAux*, we will use these and  
    // hopefully resize will do nothing as the sizes will be the same
    // If not, alas, we will need to allocate local stuff

    // Matrices instantiation is the only price
    // we pay for the flexibility described above.
    // Should be very cheap since there is no memory allocation.
    gAux1->resizeOnly(NumCols, NumCols); DoubleMatrix * _aux1 = gAux1;
    gAux11->resizeOnly(NumCols,NumCols); DoubleMatrix * _aux11= gAux11;
    gAux2->resizeOnly(NumCols, NumCols); DoubleMatrix * _aux2 = gAux2;
    gAux21->resizeOnly(NumCols,NumCols); DoubleMatrix * _aux21= gAux21;
    gAux3->resizeOnly(NumCols, NumCols);
    
    ///////////////////////////////////////////////////////////////////////////
    static DoublePointerVector _YVector;
    static DoublePointerVector _AVector;
    static DoubleVector        _ColWrapper;
    static IntVector           _IndxWrapper;

    ///////////////////////////////////////////////////////////////////////////
    // If the spectral norm is negative, assume that you should recompute it 
    if (gSpectralNorm < 0) {
        gSpectralNorm = this->spectralNorm();
    }

    ///////////////////////////////////////////////////////////////////////////
    // If |tM|>=0.5, then halve t until it's not. 
    // Later, square result nSquaring times.
    int nSquaring = 0;
    while (fabs(gT*gSpectralNorm)>=0.5) {
        nSquaring++;
        gT /= 2.0;
    }

    // we temporarily use _aux3 to store (*this)*gT
    DoubleMatrix::scale(gT, *this, gAux3); 

    // Compute Pade coefficients 
    pairDoubleVectorPointers pqCoeffs=PADEcoeffs.calsPQ(gPadeOrder,gPadeOrder);

    // Compute P(mt) and Q(mt)
    _aux1->fill(0.0);
    _aux2->fill(0.0);

    for (int i = gPadeOrder; i>=0; i--) {
        if (i != gPadeOrder)
        {
            DoubleMatrix::mult(*gAux3, *_aux1, _aux11); //num
            DoubleMatrix * _buffer1 = _aux1; 
                            _aux1 = _aux11; 
                            _aux11 = _buffer1;
                           // swap _aux1 and _aux11 
            DoubleMatrix::mult(*gAux3, *_aux2, _aux21); //den
            DoubleMatrix * _buffer2 = _aux2; 
                           _aux2 = _aux21;
                           _aux21 = _buffer2;
                           // swap _aux2 and _aux21 
        }
        const double _pCoeff = (*pqCoeffs.first)[i]; //num
        const double _qCoeff = (*pqCoeffs.second)[i];//den
        for(int ii = 0; ii < NumCols; ii++)
        {    
            (*_aux1)[ii][ii] += _pCoeff; //num
            (*_aux2)[ii][ii] += _qCoeff; //den
       }
    }

    // invert denominator
    _aux2->computeInverse(gAux3, _aux11, &_YVector, &_AVector, 
                                          &_ColWrapper, &_IndxWrapper); 
    // now aux3 = inverse of den


    if (nSquaring == 0)
    {
        DoubleMatrix::mult(*_aux1, *gAux3, gResult); // Result = _aux1 * _aux3 = num * inv(den)  
    }
    else
    {
        DoubleMatrix::mult(*_aux1, *gAux3, _aux2); // _aux2 = _aux1 * _aux3 = num * inv(den)
        while (nSquaring>1) { // the aim of this loop is aux2 := aux2 * aux2
            DoubleMatrix::mult(*_aux2, *_aux2, _aux1);
            DoubleMatrix * _buffer = _aux1; 
                           _aux1 = _aux2; 
                           _aux2 = _buffer;
            nSquaring--;
        }
        DoubleMatrix::mult(*_aux2, *_aux2, gResult); // Result = _aux2 * _aux2  
    }
}

/** Computes eigenvalues and eigenvectors of a given matrix */
EigenVectorAnalysisSP CDoubleMatrix::computeEigenVectors() const {
    static const string routine = "CDoubleMatrix::computeEigenVectors";
    
    try {
        checkSymmetric();
        int n = NumRows;

        CDoubleMatrix matrixCopy(*this);
        CDoubleMatrixSP eigenVectors(new CDoubleMatrix(n, n));
        
        vector<double*> inputPtrs(n);
        vector<double*> eigenVectorsPtrs(n);
        int i;
        for(i = 0; i < n; i++) {
            inputPtrs[i] = matrixCopy[i] - 1;
            eigenVectorsPtrs[i] = (*eigenVectors)[i] - 1;
        }

        DoubleArraySP eigenValues(new DoubleArray(n));
        
        int nrot = 0;
        int* nrotPtr = &nrot;
        
        try {
            jacobi(&inputPtrs[0] - 1,
                   n,
                   &(*eigenValues)[0] - 1,
                   &eigenVectorsPtrs[0] - 1,
                   nrotPtr);
        
            eigsrt(&(*eigenValues)[0] - 1,
                   &eigenVectorsPtrs[0] - 1,
                   n);
        } catch(exception& e) {
            string message = "Failed to compute the eigen values of the matrix.";
            throw ModelException::addTextToException(e, message);
        }
        
        eigenVectors->transpose(); // when we call [] we should get an eigenvector back
        
        EigenVectorAnalysisSP results(new EigenVectorAnalysis(n));
        
        results->eigenVectors = eigenVectors;
        results->eigenValues  = eigenValues;
        results->dimension    = n;
        
        return results;
    
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

/** Computes the Singular Value Decomposition of a given matrix */
SVDAnalysisSP CDoubleMatrix::computeSVD() const {
    static const string routine = "CDoubleMatrix::computeSVD";

    try {
    
        CDoubleMatrixSP matrixU(new CDoubleMatrix(*this));
        matrixU->transpose();
        vector<double*> uPtrs(NumRows);
        int i;
        for(i = 0; i < NumRows; i++) {
            uPtrs[i] = (*matrixU)[i] - 1;
        }

        CDoubleMatrixSP matrixVTranspose(new CDoubleMatrix(NumCols, NumCols));
        vector<double*> vPtrs(NumCols);
        for(i = 0; i < NumCols; i++) {
            vPtrs[i] = (*matrixVTranspose)[i] - 1;
        }

        CDoubleArraySP sValues(new CDoubleArray(NumCols));
    
        try {
            svdcmp(&uPtrs[0] - 1,
                   NumRows,
                   NumCols,
                   &(*sValues)[0] - 1,
                   &vPtrs[0] - 1);
        } catch(exception& e) {
            string message = "Failed to compute the singular values of the matrix.";
            throw ModelException::addTextToException(e, message);
        }

        SVDAnalysisSP results(new SVDAnalysis(NumCols, NumRows));
        matrixU->transpose();
        results->matrixU          = matrixU;
        results->matrixVTranspose = matrixVTranspose;
        results->sValues          = sValues;

        return results;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

/** Computes the Beta Factorization of a given matrix */
BetaFactorizationSP CDoubleMatrix::computeBetaFactorization(int    nbFactors,
                                                            double precision,
                                                            int    maxIter) const {
    static const string method = "CDoubleMatrix::computeBetaFactorization";
    try {    
        /** symmetry required */
        checkSymmetric();
        
        /** positive definiteness required */
        EigenVectorAnalysisSP eigen = this->computeEigenVectors(); 
        if (Maths::isNegative(eigen->eigenValues->back())) {
            throw ModelException(method, 
                "Beta factorization only applicable to positive definite matrices");
        }
        /** dimension checking */
        int dimension = this->numCols();
        if (nbFactors > dimension) {
            throw ModelException(method, 
                "Nb of beta factors must not exceed dimension of correlation matrix. \n"
                "But nb of beta factors equals " + Format::toString(nbFactors) + 
                " and dimension of correlation matrix equals " + Format::toString(dimension) + ".");
        }
        /** initial guess is zero*/
        double average = 0.0;
        int i,j,m;
        for (i=0; i<dimension; i++) {
            for (j=i+1; j<dimension; j++) {
                average += (*this)[i][j];
            }
        }
        average /= ((double)(dimension*dimension-dimension) / 2.0);
        average = sqrt(average);
        CDoubleMatrixSP initialGuess(new DoubleMatrix(nbFactors, dimension));
        for (i=0; i<dimension; i++) {
            for (m=0; m<nbFactors; m++) {
                (*initialGuess)[m][i] = 0.0; // average;
            }
        }

        /* define F */        
        DoubleMatrix c(*initialGuess.get());        
        DoubleMatrix F(dimension, dimension), saveF(dimension, dimension);
        
        for (int iter=0; iter<maxIter; iter++) {
            DoubleMatrix ccT(dimension, dimension);
            for (i=0; i<dimension; i++) {
                for (j=0; j<dimension; j++) {
                    ccT[i][j] = 0.0;
                    for (m=0; m<nbFactors; m++) {
                        ccT[i][j] += c[m][i] * c[m][j]; // attention: row = col & col = row
                    }
                }
            }

            if (iter>0) {
                saveF = F; // save F and update F
                for (i=0; i<dimension; i++) {
                    F[i][i] = 1.0 - ccT[i][i];            
                }
                double diff = 0.0;
                for (i=0; i<dimension; i++) {
                    diff += (F[i][i] - saveF[i][i])*(F[i][i] - saveF[i][i]);
                }
                if (diff<precision) {
                    break;
                }            
            } else {
                for (i=0;i<dimension;i++) {
                    F[i][i] = 1.0 - ccT[i][i];
                }
            }
            
            DoubleMatrix E(dimension, dimension);
            for (i=0;i<dimension;i++) {
                for (j=0;j<dimension;j++) {
                    E[i][j] = (*this)[i][j] - F[i][j];
                }
            }
            EigenVectorAnalysisSP eigen = E.computeEigenVectors(); 
            /** update c */
            for (m=0; m<nbFactors; m++) {
                for (i=0; i<dimension; i++) {
                    c[m][i] = (*eigen->eigenVectors)[m][i] * sqrt(max((*eigen->eigenValues)[m],0.0));
                }
            }
        }

        /** compute beta correlation matrix */
        CDoubleMatrixSP factorMatrix(new CDoubleMatrix(c));
        CDoubleMatrixSP betaCorrelationMatrix(new CDoubleMatrix(dimension,dimension));
        for (i=0; i<dimension; i++) {
            for (j=i+1; j<dimension; j++) {
                (*betaCorrelationMatrix)[i][j] = 0.0;
                for (m=0; m<nbFactors; m++) {
                    (*betaCorrelationMatrix)[i][j] += c[m][i] * c[m][j];
                }
                (*betaCorrelationMatrix)[j][i] = (*betaCorrelationMatrix)[i][j] ;
            }
            (*betaCorrelationMatrix)[i][i]  = 1.0;
        }

        /** compute average distance */
        double averageDistance = 0.0;
        for (i=0; i<dimension; i++) {
            for (j=i+1; j<dimension; j++) {
                averageDistance += ( (*this)[i][j] - (*betaCorrelationMatrix)[i][j] )
                    * ( (*this)[i][j] - (*betaCorrelationMatrix)[i][j] );
            }
        }
        averageDistance /= ((double)(dimension*dimension-dimension) / 2.0);

        /** compute lower and upper squeeze limits */
        double squeezeLimitLow = -1.0;
        double squeezeLimitHigh = 1.0;
        if (nbFactors > 1) {
            squeezeLimitLow = -1.0;
            squeezeLimitHigh = 1.0;
            for (int i=0; i<factorMatrix->numRows(); i++) {
                double tmp = 1.0;
                for (int j=1; j<nbFactors; j++) {
                    tmp -= (*factorMatrix)[j][i]*(*factorMatrix)[j][i];                
                }
                if (Maths::isNegative(tmp)) {
                    squeezeLimitLow = 0.0;
                    squeezeLimitHigh = 0.0;
                    break;
                }
                double thisLimitLow = (-sqrt(tmp)-(*factorMatrix)[0][i])/(1.0-(*factorMatrix)[0][i]);
                double thisLimitHigh = (sqrt(tmp)-(*factorMatrix)[0][i])/(1.0-(*factorMatrix)[0][i]);
                if (Maths::isPositive(thisLimitLow - squeezeLimitLow)) {
                    squeezeLimitLow = thisLimitLow;
                }
                if (Maths::isNegative(thisLimitHigh - squeezeLimitHigh)) {
                    squeezeLimitHigh = thisLimitHigh;
                }
            }
        }
        
        return BetaFactorizationSP(new BetaFactorization(averageDistance,
                                                         factorMatrix,
                                                         betaCorrelationMatrix,
                                                         squeezeLimitLow,
                                                         squeezeLimitHigh));    
        } catch(exception& e) {
        throw ModelException(e, method);
    }

}

/** Addin for creating a valid correlation matrix */
class DoubleMatrixCorrAddin: public CObject{
public:
    /** Computes a real symmetric matrix correlation matrix from a given 
        input correlation matrix by flooring its eigenvalues. When the
        floor is positive the resulting matrix is positive definite i.e. a
        valid correlation matrix. */
    static IObjectSP symmToCorrel(DoubleMatrixCorrAddin* params){
        double sumSqr = 0.0;
        ObjectArraySP result(new ObjectArray(2));
        CDoubleMatrixSP newMatrix(new DoubleMatrix(params->matrix->symmToCorrel(&sumSqr, params->eigenValueFloor)));
        (*result)[0] = newMatrix;
        IObjectSP sqError(CDouble::create(sumSqr));
        (*result)[1] = sqError;
        return result;
    }

    /** for reflection */
    DoubleMatrixCorrAddin():  CObject(TYPE), eigenValueFloor(0.0000000001) {}

private:
/** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DoubleMatrixCorrAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDoubleMatrixCorrAddin);
        FIELD(matrix, "A matrix of doubles");
        FIELD(eigenValueFloor, "Floor for eigenvalues. Positive floor ensures positive definite.");
        FIELD_MAKE_OPTIONAL(eigenValueFloor);

        Addin::registerClassObjectMethod("MATRIX_CORRS_MAKE_VALID",
                                         Addin::UTILITIES,
                                         "Creates a valid correlation matrix",
                                         TYPE,
                                         false,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)symmToCorrel);
    }
        
    static IObject* defaultDoubleMatrixCorrAddin(){
        return new DoubleMatrixCorrAddin();
    }
    
    static CClassConstSP const TYPE;

    /** the single parameter that our addin takes */
    CDoubleMatrixSP matrix;
    double          eigenValueFloor;
};


CClassConstSP const DoubleMatrixCorrAddin::TYPE = CClass::registerClassLoadMethod(
    "DoubleMatrixCorrAddin", typeid(DoubleMatrixCorrAddin), load);


/** Addin for building handles to a double matrix */
class DoubleMatrixAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    /** the single parameter that our addin takes */
    CDoubleMatrixSP  matrix;

    /** the 'addin function' - doesn't do anything clever but just allows us
        to use the infrastructure */
    static IObjectSP createDoubleMatrix(DoubleMatrixAddin* params){
        return params->matrix;
    }

    /** the square root aka Choleski / Cholesky etc. */
    static IObjectSP squareRoot(DoubleMatrixAddin* params){
        CDoubleMatrix root = params->matrix->computeSquareRoot();
        CDoubleMatrixSP rootSP(new CDoubleMatrix(root));
        return IObjectSP(rootSP.clone());
    }

    /** the inverse of a square matrix using LU */
    static IObjectSP inverse(DoubleMatrixAddin* params){
        CDoubleMatrix inverse = params->matrix->computeInverse();
        CDoubleMatrixSP inverseSP(new CDoubleMatrix(inverse));
        return IObjectSP(inverseSP.clone());
    }

    /** computes eigen vectors and eigenvalues */
    static IObjectSP eigenVectors(DoubleMatrixAddin* params){
        EigenVectorAnalysisSP eigenResults = params->matrix->computeEigenVectors();
        return IObjectSP(eigenResults.clone());
    }

    /** computes SVD */
    static IObjectSP svd(DoubleMatrixAddin* params){
        SVDAnalysisSP svdResults = params->matrix->computeSVD();
        return IObjectSP(svdResults.clone());
    }

    /** computes the transpose matrix and returns an handle*/
    static IObjectSP transpose(DoubleMatrixAddin* params){
        CDoubleMatrixSP transSP(params->matrix);
        transSP->transpose();
        return IObjectSP(transSP.clone());
    }

    /** computes the transpose matrix and returns an handle*/
    double spectralNorm(){       
        return matrix->spectralNorm();
    }

    /** computes the transpose matrix and returns an handle*/
    double frobeniusNorm(){
        return matrix->frobeniusNorm();
    }

    /** for reflection */
    DoubleMatrixAddin():  CObject(TYPE){
        // parameter is optional so default to empty matrix
        matrix = CDoubleMatrixSP(new CDoubleMatrix());
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DoubleMatrixAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDoubleMatrixAddin);
        FIELD(matrix, "A matrix of doubles");
        FIELD_MAKE_OPTIONAL(matrix);
        Addin::registerClassObjectMethod(
            "DOUBLE_MATRIX",
            Addin::UTILITIES,
            "Constructs a handle to a matrix "
            "of doubles",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)createDoubleMatrix);

        Addin::registerClassObjectMethod("MATRIX_SQUARE_ROOT",
                                         Addin::UTILITIES,
                                         "Computes the square root of a matrix",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)squareRoot);

        Addin::registerClassObjectMethod("MATRIX_EIGENVECTORS",
                                         Addin::UTILITIES,
                                         "Computes the eigenvectors and eigenvalues of a matrix",
                                         TYPE,
                                         false,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)eigenVectors);

        Addin::registerClassObjectMethod("MATRIX_SVD",
                                         Addin::UTILITIES,
                                         "Computes the SVD of a matrix",
                                         TYPE,
                                         false,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)svd);
    
        Addin::registerClassObjectMethod("MATRIX_INVERSE",
                                         Addin::UTILITIES,
                                         "Computes the inverse of a matrix",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)inverse);

        Addin::registerClassObjectMethod("MATRIX_TRANSPOSE",
                                         Addin::UTILITIES,
                                         "Returns the tranposed matrix",
                                         TYPE,
                                         false,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)transpose);

        Addin::registerDoubleMethod("MATRIX_SPECTRAL_NORM",
            Addin::UTILITIES,
            "Returns the spectral norm of a matrix: the square root of the biggest magnitude eigen-value of A'A.",
            &DoubleMatrixAddin::spectralNorm);

        Addin::registerDoubleMethod("MATRIX_FROBENIUS_NORM",
            Addin::UTILITIES,
            "Returns the Frobenius norm of a matrix: sqrt(Tr(A'A)).",
            &DoubleMatrixAddin::frobeniusNorm);
    }

    static IObject* defaultDoubleMatrixAddin(){
        return new DoubleMatrixAddin();
    }
    
};

CClassConstSP const DoubleMatrixAddin::TYPE = CClass::registerClassLoadMethod(
    "DoubleMatrixAddin", typeid(DoubleMatrixAddin), load);

/** Addin for building handles to a double matrix */
class SqueezeMatrixAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    /** the single parameter that our addin takes */
    CDoubleMatrixSP  matrix;
    double           squeezeInput;

    /** computes the transpose matrix and returns an handle*/
    static IObjectSP squeezeMatrix(SqueezeMatrixAddin* params){
        CDoubleMatrixSP squeezedSP(params->matrix);
        squeezedSP->squeeze(params->squeezeInput);
        return IObjectSP(squeezedSP.clone());
    }

    /** for reflection */
    SqueezeMatrixAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(SqueezeMatrixAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultSqueezeMatrixAddin);
        FIELD(matrix, "A matrix of doubles");        
        FIELD(squeezeInput, "a squeeze");        

        Addin::registerClassObjectMethod("MATRIX_SQUEEZE",
                                         Addin::UTILITIES,
                                         "Returns the tranposed matrix",
                                         TYPE,
                                         false,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)squeezeMatrix);                
    }

    static IObject* defaultSqueezeMatrixAddin(){
        return new SqueezeMatrixAddin();
    }
    
};

CClassConstSP const SqueezeMatrixAddin::TYPE = CClass::registerClassLoadMethod(
    "SqueezeMatrixAddinn", typeid(SqueezeMatrixAddin), load);

/** Addin for building handles to a double matrix */
class EmptyDoubleMatrixAddin: public CObject{
    static CClassConstSP const TYPE;

    int numCols;
    int numRows;

    static IObjectSP createDoubleMatrix(EmptyDoubleMatrixAddin* params){
        return IObjectSP(new CDoubleMatrix(params->numCols, params->numRows));
    }

    /** for reflection */
    EmptyDoubleMatrixAddin():  CObject(TYPE){ }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(EmptyDoubleMatrixAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultEmptyDoubleMatrixAddin);
        FIELD(numCols, "Number of Columns");
        FIELD(numRows, "Number of Rows");
        Addin::registerClassObjectMethod(
            "EMPTY_DOUBLE_MATRIX",
            Addin::UTILITIES,
            "Constructs a handle to a matrix populated with zeros",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)createDoubleMatrix);
    }

    static IObject* defaultEmptyDoubleMatrixAddin(){
        return new EmptyDoubleMatrixAddin();
    }

};

CClassConstSP const EmptyDoubleMatrixAddin::TYPE =
CClass::registerClassLoadMethod(
    "EmptyDoubleMatrixAddin", typeid(EmptyDoubleMatrixAddin), load);

/** Addin for building handles to a double matrix from a double array */
class CreateDoubleMatrixAddin: public CObject{
    static CClassConstSP const TYPE;

    int              numCols;
    int              numRows;
    CDoubleArraySP   doubleArray;

    static IObjectSP createDoubleMatrix(CreateDoubleMatrixAddin* params){
        return IObjectSP(new CDoubleMatrix(params->numCols, params->numRows));
    }

    /** initializes the matrix with an array (in 2D) */
    static IObjectSP initFromArray(CreateDoubleMatrixAddin* params)
    {
        CDoubleMatrixSP matrix(new CDoubleMatrix(params->numCols, params->numRows,
                                                 * params->doubleArray.get()));
        return IObjectSP(matrix.get());
    }

    /** for reflection */
    CreateDoubleMatrixAddin():  CObject(TYPE){ }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(CreateDoubleMatrixAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCreateDoubleMatrixAddin);
        FIELD(numCols, "Number of Columns");
        FIELD(numRows, "Number of Rows");
        FIELD(doubleArray, "Array used to initialize the matrix");

        Addin::registerClassObjectMethod(
            "CREATE_MATRIX_FROM_1D_ARRAY",
            Addin::UTILITIES,
            "Constructs a handle to a matrix of doubles using a DoubleArray",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*) initFromArray
            );
    }

    static IObject* defaultCreateDoubleMatrixAddin(){
        return new CreateDoubleMatrixAddin();
    }
};

CClassConstSP const CreateDoubleMatrixAddin::TYPE =
CClass::registerClassLoadMethod(
    "CreateDoubleMatrixAddin", typeid(CreateDoubleMatrixAddin), load);
/** Addin for computing matrix exponents */
class DoubleMatrixAddinExp: public CObject{
public:
    static CClassConstSP const TYPE;

    /** the parameters that our addin takes */
    CDoubleMatrixSP  matrix;
    int              padeOrder;

    /** computes the transpose matrix and returns an handle*/
    IObjectSP matrixExp() {   
        CDoubleMatrixSP expm(new CDoubleMatrix());
        (*expm) = matrix->exp(1.0, -1.0, padeOrder);
        QLIB_VERIFY(!!expm, "expm is NULL");
        return IObjectSP(expm.clone());
    }

    /** for reflection */
    DoubleMatrixAddinExp():  CObject(TYPE){
        // parameter is optional so default to empty matrix
        matrix = CDoubleMatrixSP(new CDoubleMatrix());
        padeOrder = 0;
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DoubleMatrixAddinExp, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDoubleMatrixAddinExp);
        FIELD(matrix, "A matrix of doubles");
        FIELD(padeOrder, "Order of Pade approximation to use for exp(x): numerator and denominator orders are the same.");

        Addin::registerObjectMethod(
            "MATRIX_EXP",
            Addin::UTILITIES,
            "Computes the exponential of a matrix using a specified order Pade approximation.",
            false,
            Addin::returnHandle,
            &DoubleMatrixAddinExp::matrixExp);
    }

    static IObject* defaultDoubleMatrixAddinExp(){
        return new DoubleMatrixAddinExp();
    }
    
};

CClassConstSP const DoubleMatrixAddinExp::TYPE = CClass::registerClassLoadMethod(
    "DoubleMatrixAddinExp", typeid(DoubleMatrixAddinExp), load);

/** Addin for beta factorization of correlation matrix */
class BetaFactorizationAddin: public CObject{
    static CClassConstSP const TYPE;

    CDoubleMatrixSP correlationMatrix;
    int             nbFactors;
    double          precision;
    int             maxIter;

    static IObjectSP compBetas(BetaFactorizationAddin* params){
        return params->correlationMatrix->computeBetaFactorization(params->nbFactors,
                                                                   params->precision,
                                                                   params->maxIter);        
    }

    /** for reflection */
    BetaFactorizationAddin():  CObject(TYPE){ }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(BetaFactorizationAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultBetaFactorizationAddin);
        FIELD(correlationMatrix, "number of factors");
        FIELD(nbFactors, "number of factors");
        FIELD(precision, "precision");
        FIELD(maxIter, "max nb of iterations");

        Addin::registerClassObjectMethod(
            "MATRIX_BETA_FACTORIZATION",
            Addin::UTILITIES,
            "Computes beta factorization of correlation matrix",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)compBetas);
    }

    static IObject* defaultBetaFactorizationAddin(){
        return new BetaFactorizationAddin();
    }
};

CClassConstSP const BetaFactorizationAddin::TYPE =
CClass::registerClassLoadMethod(
    "BetaFactorizationAddin", typeid(BetaFactorizationAddin), load);

DEFINE_TEMPLATE_TYPE(DoubleMatrixArray);
    
DRLIB_END_NAMESPACE

