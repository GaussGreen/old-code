#ifndef __DKMAILLE2D_H
#define __DKMAILLE2D_H

#include <iostream>
#include <strstream>

#pragma warning(disable : 4251)
#pragma warning(default : 4251)

#include "DKMaille.h"


// --------------------------------------------------------------------------------
template <class Type>
class DKMaille2D 
{
public:
  // Constructors
  DKMaille2D(int r=0,int c=0);           // initial size - for performance
  DKMaille2D(int r, int c, Type **t);
  DKMaille2D(int r, int c, Type *t);     // This takes ownership of the memory
  DKMaille2D(int r,int c,const Type& t); // Initialisation value
  // DKMaille2D(const RWCString& filename); // Initialisation from a file
  DKMaille2D& operator=(const DKMaille2D<Type>&);
  DKMaille2D(const DKMaille2D<Type>&);
  virtual ~DKMaille2D();

  Type *rawAlloc(unsigned N) { return new Type[N]; }
  void rawDealloc(Type *p) { delete [] p; } 

  Type *newNagArray() const;

  bool operator==(const DKMaille2D<Type>& rhs) const { return this == &rhs; }

  void resize(int r, int c) { growTo(r,c); 	myNoRows=r;	myNoColumns=c;} 
  // Access functions
  int rows() const { return myNoRows; };		// Number of rows
  int columns() const { return myNoColumns; };	// Number of columns
	const Type *raw() const;  // Nasty Hacked function 


  // Access rows and columns - const and Non const like Rogue

#ifdef DG_FULL_CHECKING
  const Type&	operator()(int r, int c) const 
  { 
    if ((r >= rows()) || (r < 0) || (c >= columns() || (c < 0)))
    {
      throw("ERROR: Accessing outside declared DKMaille2D");
    }
    else
      return myArray[r*myMaxColumn+c]; 
  }

  Type&	operator()(int r, int c)		  
  { 
    if ((r >= rows()) || (r < 0) || (c >= columns() || (c < 0)))
    {
      throw("ERROR: Accessing outside declared DKMaille2D");
    }
    else
      return myArray[r*myMaxColumn+c]; 
  }

  const Type&	at(int r, int c)	const	   
  { 
    if ((r >= rows()) || (r < 0) || (c >= columns() || (c < 0)))
    {
      throw("ERROR: Accessing outside declared DKMaille2D");
    }
    else
      return myArray[r*myMaxColumn+c]; 
  }

  Type&	    at(int r, int c)		       
  { 
    if ((r >= rows()) || (r < 0) || (c >= columns() || (c < 0)))
    {
      throw("ERROR: Accessing outside declared DKMaille2D");
    }
    else
      return myArray[r*myMaxColumn+c]; 
  }


#else // RELEASE
  const Type&	operator()(int r, int c) const 
  { 
    return myArray[r*myMaxColumn+c]; 
  }
  Type&	operator()(int r, int c)		  
  { 
    return myArray[r*myMaxColumn+c]; 
  }
  const Type&	at(int r, int c)	const	   
  { 
    return myArray[r*myMaxColumn+c]; 
  }
  Type&	    at(int r, int c)		       
  { 
    return myArray[r*myMaxColumn+c]; 
  }
#endif // DG_FULL_CHECKING or release

	// Backwards compatibility functions
	Type **newCArray() const; // Warning this function allocates HEAP memory
  inline static void deleteCArray(int r, int c, Type **array);

  // Row & Column manipulations:
  void checkRange(int r, int c) const
  {
    if(r<0 || r>=rows() || c<0 || c>=columns())
    {
      MESSAGE("Position (");
      msg.stream() << r << ", " << c << ") is out of range.";
      throw msg;
    }
  }

  void print(FILE* fic);

  void insertRow(int r,int count=1);		// Insert 'count' new rows @ r, shift the rows >=r
  void insertColumn(int c,int count=1);
  
  void copyRow(int source,int dest);	// Copy a row
  void copyColumn(int source,int dest);
  
  void appendRow();			// Add an extra row to the end
  void appendColumn();

  void swapRows(int r1,int r2);	// Permute two rows
  void swapColumns(int c1,int c2);
  
  void removeRow(int r, int count=1);	// Remove a given row/column.
  void removeColumn(int c, int count=1);
  
  void removeAll();		// Removes all rows & columns
  void reset();		// Deletes everything and sets the size to 0,0.

  void fillColumn(int c,const Type& t);
  void fillRow(int r,const Type& t);
  void fill(const Type& t);	// Fill with a single value

  void appendRowsColumns(int r,int c);	// Add several rows & columns in one swoop
  
  // Convenience functions
  int lastRow() { return rows()-1; };
  int lastColumn() { return columns()-1; };
  int firstRow() { return 0; };
  int firstColumn() { return 0; };

  
  int boundCheck(int r,int c) const { return (r>=0 && r<myNoRows && c>=0 && c<myNoColumns); };
  void growTo(int r,int c);

  // the array:
  void coords(const Type* t,int& r, int& c) 
  {
	int d=t-myArray;
	r=d/myMaxColumn;
	c=d-r*myMaxColumn;
  };

  int checkedCoords(const Type* t,int& r, int& c) 
  {
	coords(t,r,c);
	if(r<0 || c<0 || r>=myNoRows || c>=myNoColumns) 
    {
	  r= -1;
	  c= -1;
	  return 0;
    }
	else
    {
	  return 1;
    }
  };

  DKMaille2D<Type>& invert(); // to get the inverse of the "this"

  DKMaille2D<Type>& choleskyDecompose();
	DKMaille2D<Type>& transpose();
	DKMaille2D<Type>& operator*=(double scalar);
	DKMaille2D<Type>& operator+=(double scalar);
	DKMaille2D<Type>& operator-=(double scalar);
  DKMaille<Type> diagElement() const;
  DKMaille<Type> columnVector(int col) const;
  DKMaille<Type> rowVector(int row) const;


  DKMaille2D<Type>& byElementMultiply(const Type& scalar);
	DKMaille2D<Type>& byElementAdd(const Type& scalar);
	DKMaille2D<Type>& byElementSubtract(const Type& scalar);
	DKMaille2D<Type>& byElementDivide(const Type& scalar);

	DKMaille2D<Type>& operator*=(const DKMaille2D<Type>& rhs);
	DKMaille2D<Type>& operator+=(const DKMaille2D<Type>& rhs);
	DKMaille2D<Type>& operator-=(const DKMaille2D<Type>& rhs);



  friend DKMaille2D<Type>  operator*(const DKMaille2D<Type>& lhs, const DKMaille2D<Type>& rhs);
	friend DKMaille2D<Type>  operator+(const DKMaille2D<Type>& lhs, const DKMaille2D<Type>& rhs);
	friend DKMaille2D<Type>  operator-(const DKMaille2D<Type>& lhs, const DKMaille2D<Type>& rhs);



  friend DKMaille<Type>    operator*(const DKMaille<Type>& lhs, const DKMaille2D<Type>& rhs);
  friend DKMaille<Type>    operator*(const DKMaille2D<Type>& lhs, const DKMaille<Type>& rhs);
    
private:
	//DKMaille2D<Type>& inverert();
  int myMaxRow;
  int myMaxColumn;
  int myNoRows;
  int myNoColumns;
  Type* myArray;

  void reallocate(int r2,int c2,int f=0);

};

// RWDECLARE_PERSIST_TEMPLATE_IO(DKMaille2D,RWvistream,RWvostream) // we don't want the RWFile stuff

// --------------------------------------------------------------------------------
template <class Type>
DKMaille2D<Type>::~DKMaille2D()
{
  if(myArray) rawDealloc(myArray);
  myArray = NULL;   // Stomp
}

// --------------------------------------------------------------------------------
template <class Type>
DKMaille2D<Type>::DKMaille2D(int r,int c):
    myNoRows(r),
    myNoColumns(c),
    myMaxRow(r),
    myMaxColumn(c),
    myArray(NULL)
{
  if(r && c)
	reallocate(r,c,1);
}

// --------------------------------------------------------------------------------
template <class Type>
DKMaille2D<Type>::DKMaille2D(int r,int c,const Type& t):
    myNoRows(r),
    myNoColumns(c),
    myMaxRow(r),
    myMaxColumn(c),
    myArray(NULL)
{
  reallocate(r,c,1);
  fill(t);
}

// --------------------------------------------------------------------------------
template <class Type>
DKMaille2D<Type>::DKMaille2D(const DKMaille2D<Type>& x) :
    myNoRows(0),
    myNoColumns(0),
    myMaxRow(0),
    myMaxColumn(0),
    myArray(NULL)
{
    *this = x;
}

// --------------------------------------------------------------------------------
template <class Type>
DKMaille2D<Type>::DKMaille2D(int r, int c, Type **t) :
	myNoRows(0),
    myNoColumns(0),
    myMaxRow(0),
    myMaxColumn(0),
    myArray(NULL)
{
  // size the array
  appendRowsColumns(r, c);
  for (int i = 0; i < r; ++i)
  {
	  for (int j = 0; j < c; ++j)
	  {
		  at(i, j) = t[i][j];
	  }
  }
}

// --------------------------------------------------------------------------------
template <class Type>
DKMaille2D<Type>::DKMaille2D(int r, int c, Type *t) :   // This takes ownership of the memory
	myNoRows(r),
  myNoColumns(c),
  myMaxRow(r),
  myMaxColumn(c),
  myArray(t)
{}

// --------------------------------------------------------------------------------
// Allocate a nag array
template <class Type>
Type *DKMaille2D<Type>::newNagArray() const
{
  const_cast<DKMaille2D<Type>*>(this)->reallocate(myNoRows, myNoColumns, true);
  int size = myMaxRow * myMaxColumn;

  Type* ret = new Type[size];
  for(int i = 0; i < size; ++i)
  {
    ret[i] = myArray[i];   // cannot use memcpy cause we wan't to call copy constructors
  }
  return ret;
}


// --------------------------------------------------------------------------------
// Reallocate storage space for the array.
// Resizing to a smaller size will not happen unless force is set to True.
// If the array is shrunk, then values may be lost = not copied into new array.

template <class Type>
void DKMaille2D<Type>::reallocate(int r2,int c2,int force)
{
  // Don't lose elements unless forced
  if(!force) 
  {
	if(r2<myNoRows)    r2=myNoRows;
	if(c2<myNoColumns) c2=myNoColumns;
  }
  
  // Don't resize to less than 1x1
  if(r2<1) r2=1;
  if(c2<1) c2=1;

  // Re-allocate if necessary or forced to:
  if(force || r2>myMaxRow || c2>myMaxColumn) 
  {
    Type *oldArray= myArray;
    int r1=myMaxRow;
    int c1=myMaxColumn;

    myArray=rawAlloc((r2)*(c2));

    // Store the new dimensions BEFORE doing any copying.  This is
    // important because some copy constructors may be aware of the
    // object's position in the array (grid cells for example).
    myMaxRow=r2;
    myMaxColumn=c2;

    // Determine which values to copy & shrink the real number of used rows/cols if needed.
    int rcopy=myNoRows; if(r2<myNoRows) { rcopy=myNoRows=r2; };
    int ccopy=myNoColumns; if(c2<myNoColumns) { ccopy=myNoColumns=c2; };
  
    // Copy old values if any
    if(oldArray) 
    {
      int r,c;
      for(r=0;r<rcopy;r++)
      for(c=0;c<ccopy;c++)
	    myArray[r*c2+c]=oldArray[r*c1+c];

      rawDealloc(oldArray);
    }
  }
}

// --------------------------------------------------------------------------------
template <class Type>
void DKMaille2D<Type>::growTo(int r,int c)
{
  reallocate(r,c); 
  if ( r > myNoRows)
    myNoRows = r;
  if ( c > myNoColumns)
    myNoColumns = c;
}

// --------------------------------------------------------------------------------
template <class Type>
void DKMaille2D<Type>::appendRow()
{
    if(myNoRows==myMaxRow)
	reallocate((myMaxRow+1)*2,myMaxColumn);
    myNoRows++;
}

// --------------------------------------------------------------------------------
template <class Type>
void DKMaille2D<Type>::appendColumn()
{
    if(myNoColumns==myMaxColumn)
	reallocate(myMaxRow,(myMaxColumn+1)*2);
    myNoColumns++;
}

// --------------------------------------------------------------------------------
template <class Type>
void DKMaille2D<Type>::appendRowsColumns(int r,int c)
{
    if(r>=0 && c>=0) {
	reallocate(rows()+r,columns()+c);
	myNoRows+=r;
	myNoColumns+=c;
    }
}

// --------------------------------------------------------------------------------
template <class Type>
void DKMaille2D<Type>::copyRow(int source,int dest)
{
    if(source!=dest) {
	int i;
	for(i=0;i<myNoColumns;i++)
	    myArray[dest*myMaxColumn+i]=myArray[source*myMaxColumn+i];

    }
}

// --------------------------------------------------------------------------------
template <class Type>
void DKMaille2D<Type>::copyColumn(int source,int dest)
{
    if(source!=dest) {
	int i;
	for(i=0;i<myNoRows;i++)
	    myArray[i*myMaxColumn+dest]=myArray[i*myMaxColumn+source];
    }
}

// --------------------------------------------------------------------------------
template <class Type>
void DKMaille2D<Type>::insertRow(int r, int count)
{
    int i;
    appendRowsColumns(count,0);
    for(i=rows()-1;i>=r+count;i--)
	copyRow(i-count,i);
}

// --------------------------------------------------------------------------------
template <class Type>
void DKMaille2D<Type>::insertColumn(int c, int count)
{
    int i;
    appendRowsColumns(0,count);
    for(i=columns()-1;i>=c+count;i--)
	copyColumn(i-count,i);
}

// --------------------------------------------------------------------------------
template <class Type>
void DKMaille2D<Type>::removeRow(int r, int count)
{
    int i;
    for(i=r+count;i<rows();i++)
	copyRow(i,i-count);
    myNoRows -= count;
}

// --------------------------------------------------------------------------------
template <class Type>
void DKMaille2D<Type>::removeColumn(int c, int count)
{
    int i;
    for(i=c+count;i<columns();i++)
	copyColumn(i,i-count);
    myNoColumns -= count;
}

// --------------------------------------------------------------------------------
template <class Type>
void DKMaille2D<Type>::removeAll()
{
    myNoColumns=0;
    myNoRows=0;
}

// --------------------------------------------------------------------------------
template <class Type>
void DKMaille2D<Type>::reset()
{
  if(myArray)
	  rawDealloc(myArray);

  myArray=NULL;
  myNoRows=0;
  myNoColumns=0;
  myMaxRow=0;
  myMaxColumn=0;
}

// --------------------------------------------------------------------------------
template <class Type>
void DKMaille2D<Type>::swapRows(int a,int b)
{
  Type temp;
  int i;
  DKMaille2D<Type>& t= *this;
  
  for(i=0;i<columns();i++) 
  {
	  temp= t(a,i);
	  t(a,i)= t(b,i);
	  t(b,i)=temp;
  }
}

// --------------------------------------------------------------------------------
template <class Type>
void DKMaille2D<Type>::swapColumns(int a,int b)
{
  Type temp;
  int i;
  DKMaille2D<Type>& t= *this;
  
  for(i=0;i<rows();i++) 
  {
	  temp= t(i,a);
	  t(i,a)= t(i,b);
	  t(i,b)=temp;
  }
}

// --------------------------------------------------------------------------------
template <class Type>
void DKMaille2D<Type>::fillColumn(int c,const Type& t)
{
  int i;
  for(i=0;i<rows();i++)
  {
  	myArray[i*myMaxColumn+c]=t;
  }
}

// --------------------------------------------------------------------------------
template <class Type>
void DKMaille2D<Type>::fillRow(int c,const Type& t)
{
  int i;
  for(i=0;i<columns();i++)
  {
    myArray[c*myMaxColumn+i]=t;
  }
}

// --------------------------------------------------------------------------------
template <class Type>
void DKMaille2D<Type>::fill(const Type& t)
{
  int i;
  for(i=0;i<myMaxColumn*myMaxRow;i++)
  {
	  myArray[i]=t;
  }
}

// --------------------------------------------------------------------------------
template <class Type>
DKMaille2D<Type>& DKMaille2D<Type>::operator=(const DKMaille2D<Type>& x)
{
  reset();
  appendRowsColumns(x.rows(), x.columns());
  for(int r=0; r<rows(); r++)
  {
    for(int c=0; c<columns(); c++)
    {
	    at(r,c)=x(r,c);
    }
  }
  return *this;
}


// --------------------------------------------------------------------------------
template <class Type>
std::ostream& operator<<(std::ostream&o ,const DKMaille2D<Type>& a)
{
  o << "Array 2D: rows=" << a.rows() << ", columns=" << a.columns() << '\n';
  int i,j;
  for(i=0;i<a.rows();i++) 
  {
	  for(j=0;j<a.columns();j++)
    {
	    o << '\t' << a(i,j);
    }
	  o << '\n';
  }
  return o;
}

// --------------------------------------------------------------------------------
// Backwards compatibility functions
template <class Type>
inline Type ** DKMaille2D<Type>::newCArray() const 
{
	// allocate memory and copy the data
	Type **ret = new Type* [rows()];
	for (int i = 0; i < rows(); ++i)
	{
		ret[i] = new Type [columns()];
		for (int j = 0; j < columns(); ++j)
		{
			ret[i][j] = at(i, j);
		}
	}
	return ret;
}

// --------------------------------------------------------------------------------
#ifdef WIN32
template <class Type>
inline void DKMaille2D<Type>::deleteCArray(int r, int, Type **array)
{
	for (int i=0; i < r; ++i) 
	{
		delete [] array[i];		
	}
	delete [] array;
	array = NULL;
}

// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// Global Helper Functions

template <class Ret, class In>
inline DKMaille2D<Ret> convertArray2D(const DKMaille2D<In>& inArray, const Ret *r = NULL)
{
	DKMaille2D<Ret> retArray(inArray.rows(), inArray.columns());
	for (int i=0; i < inArray.rows(); ++i) 
	{
		for (int j=0; j < inArray.columns(); ++j) 
		{
			retArray.at(i, j) = inArray.at(i, j);
		}
	}
	return retArray;
}

// --------------------------------------------------------------------------------
// Global Helper Functions

template <class In, class Out>
inline void convertArrayToRows(const DKMaille<In>& in, int cols, DKMaille2D<Out>& out)
{
  // Find the number of rows to fit the size of the input vector 
  int r = (double)(in.entries()) / (double)cols + 0.999;
  out.appendRowsColumns(r, cols);
  for (int i = 0; i < r; ++i)
  {
	  for (int j = 0; j < cols; ++j)
	  {
		  out.at(i, j) = in[j+cols*i];
	  }
  }
}



// Retrevice a row slice out of a two d array
template <class Ret, class In>
inline DKMaille<Ret> rowSlice(const DKMaille2D<In>& inArray, int row = 0)
{
	DKMaille<Ret> retArray(inArray.columns());

	for (int j=0; j < inArray.columns(); ++j) 
	{
		retArray.at(j) = inArray.at(row, j);
	}
	return retArray;
}

// Retrevice a column slice out of a two d array
template <class Ret, class In>
inline DKMaille<Ret> columnSlice(const DKMaille2D<In>& inArray, int col = 0)
{
	DKMaille<Ret> retArray(inArray.rows());

	for (int i=0; i < inArray.rows(); ++i) 
	{
		retArray.at(i) = inArray.at(i, col);
	}
	return retArray;
}


// Retreve a column slice out of a two d array
template <class Ret, class In>
inline void fillRow(DKMaille2D<Ret>& array, DKMaille<In> row,int row_no = 0)
{
  array.growTo(row_no + 1, row.entries());  // grow if necesary 

	for (unsigned j=0; j < row.entries(); ++j) 
	{
		array(row_no, j) = row(j);
	}
}


// Retreve a column slice out of a two d array
template <class Ret, class In>
inline void fillColumn(DKMaille2D<Ret>& array, DKMaille<In> column, int col_no = 0)
{
  array.growTo(column.entries(), col_no + 1);  // grow if necesary 

	for (unsigned i=0; i < column.entries(); ++i) 
	{
		array(i, col_no) = column(i);
	}
}


#endif

//==================================================================================
template <class Type>
DKMaille2D<Type> operator*(const DKMaille2D<Type>& lhs,
											  const DKMaille2D<Type>& rhs)
{
	DKMaille2D<Type> tmp = lhs;
	return tmp*=rhs;
}

//==================================================================================
template <class Type>
DKMaille2D<Type> operator+(const DKMaille2D<Type>& lhs,
											  const DKMaille2D<Type>& rhs)
{
	DKMaille2D<Type> tmp = lhs;
	return tmp+=rhs;
}

//==================================================================================
template <class Type>
DKMaille2D<Type> operator-(const DKMaille2D<Type>& lhs,
											  const DKMaille2D<Type>& rhs)
{
	DKMaille2D<Type> tmp = lhs;
	return tmp-=rhs;
}

//==================================================================================
template <class Type>
DKMaille<Type> operator*(const DKMaille<Type>& lhs, 
						 const DKMaille2D<Type>& rhs)
{
	if(lhs.entries() != rhs.rows())
	{
		throw("Multiplication of matrix and diagonal matrix of incompatibile sizes.");
	}

	DKMaille<Type> result(rhs.columns());

	for (int c = 0; c < rhs.columns(); ++c)
	{
		result(c) = lhs(0) * rhs(0,c);
		for (int r = 1; r < lhs.entries(); ++r)
		{
			result(c) += lhs(r) * rhs(r,c);
		}
	}

	return result;
}

//==================================================================================
template <class Type>
DKMaille<Type> operator*(const DKMaille2D<Type>& lhs,
										   const DKMaille<Type>& rhs)
{
	if(rhs.entries() != (unsigned)lhs.columns())
	{
		throw("Multiplication of matrix by vector of incompatible sizes.");
	}

	DKMaille<Type> result(lhs.rows());

	for (int r = 0; r < lhs.rows(); ++r)
	{
		result(r) = lhs(r,0) * rhs(0);
		for (unsigned c = 1; c < rhs.entries(); ++c)
		{
			result(r) += lhs(r,c) * rhs(c);
		}
	}

	return result;
}

//==================================================================================
/*template <class Type>
DKMaille2D<Type>& DKMaille2D<Type>::inverert()
{
   // Not implemented yet
   return *this;
}
*/
//==================================================================================
// Performs a cholesky decomposition
template <class Type>
DKMaille2D<Type>& DKMaille2D<Type>::choleskyDecompose()
{
	// check for a square matrix
	if(rows() != columns())
	{
		throw("Unable to perform a Cholesky Decompostion on a non square matrix.");
	}

	DKMaille<Type> diagonal(rows());
	double sum = 0.0;

  int i =0;

	for(i=0; i < rows(); ++i)
	{
		for (int j=0; j <= i; ++j)
		{
			int k = 0;
			for (sum = at(j,i); k < j; ++k)
			{
				sum -=  (double)at(i,k) * (double)at(j,k);
			}
			if (i == j)
			{
				if ( sum <= 0.0 )
				{
					throw("Cholesky Decompostion Matrix is not positive definite");
				}

				diagonal(i) = sqrt(sum);
			}
			else
			{
				at(i,j) = sum/((double)diagonal(j));
			}
		}
	}

	// copy in the diagonal and upper triangle
  i = 0;
	for (i = 0; i<rows(); ++i)
	{
		at(i,i) = diagonal(i);
		for (int j = i+1; j < columns(); ++j)
		{
			at(i,j) = 0.0;
		}
	}

	return *this;
}

//==================================================================================
template <class Type>
DKMaille2D<Type>& DKMaille2D<Type>::transpose()
{
	if (rows() == columns()) // if this is a square matirx
	{
	  for (int r = 0; r < rows(); ++r)
	  {
		for (int c = r+1; c < columns(); ++c)
		{
			// swap the 2 entries
			Type tmp = at(r,c);
			at(r,c) = at(c,r);
			at(c,r) = tmp;
		}
	  }
	}
	else   // Non square matrix
	{
		// save and resize
		DKMaille2D<Type> tmp = *this;
		resize(columns(), rows());

		// copy
		for (int r = 0; r < rows(); ++r)
		{
			for (int c = 0; c < columns(); ++c)
			{
				// swap the 2 entries
				at(r,c) = tmp(c,r);
			}
		}
	}

	return *this;
}

//==================================================================================
template <class Type>
DKMaille<Type> DKMaille2D<Type>::diagElement() const
// Copy Diagonal Elements of the 2D matrix to a new DKMaille.
{  
    int newsize = rows();
    if (columns() < newsize) newsize = columns();

    DKMaille<Type> results(newsize, 0.0);
    for (unsigned r = 0; r < results.entries(); ++r)
	{
        results[r] = at(r,r);
	}
    return results;
}

//==================================================================================
template <class Type>
DKMaille<Type> DKMaille2D<Type>::columnVector(int col) const
// Copy column vector of the 2D matrix to a new DKMaille.
{  
    int newsize = rows();

    DKMaille<Type> results(newsize, 0.0);
    for (unsigned r = 0; r < results.entries(); ++r)
	{
        results[r] = at(r,col);
	}
    return results;
}

//==================================================================================
template <class Type>
DKMaille<Type> DKMaille2D<Type>::rowVector(int row) const
// Copy row vector of the 2D matrix to a new DKMaille.
{  
    int newsize = columns();

    DKMaille<Type> results(newsize, 0.0);
    for (unsigned c = 0; c < results.entries(); ++c)
	{
        results[c] = at(row,c);
	}
    return results;
}

//==================================================================================
template <class Type>
const Type *DKMaille2D<Type>::raw() const
{
    if (myMaxColumn != myNoColumns || myMaxRow != myNoRows)
    {

        const_cast<DKMaille2D<Type>*>(this)->reallocate(myNoRows, myNoColumns, true);
    }
    
    return myArray;
}  

template <class Type>
void DKMaille2D<Type>::print(FILE* fic)
{
	
	int i,j;

	fprintf(fic, "\n");

	for(i=0;i<rows();i++)
	{
  		for(j=0;j<columns();j++)
		{
            fprintf(fic, "%lf ", at(i,j));
		
		}
		fprintf(fic, "\n");
	}
}

    
/*
//===================================================================================
template <class Type>
DKMaille2D<Type>::DKMaille2D(const RWCString& filename) :
    myNoRows(0),
    myNoColumns(0),
    myMaxRow(0),
    myMaxColumn(0),
    myArray(NULL)
{
  RWCString line, tok;
  int r = 0;
  int c, cBig, cSmall;

  // First readn to see row, and columns
  {
    std::ifstream ifs(filename);

	while(line.readLine(ifs))
	{
	  if(line.isNull())break;

	  int l = line.length();
	  if (l<5) continue;

	  c = 0;
	  RWCTokenizer next(line);
	  while(!(tok = next("{},\t ")).isNull())
	  {
		double val = atof(tok);		
		c++;
	  }

	  switch(r)
	  {
	  case 0:
		cSmall = c;
		cBig   = c;
		break;
	  default:
	    if (c < cSmall) cSmall = c;
		if (c > cBig)   cBig   = c;
	  };
	  
	  r++;
	}
  }

  DKMaille2D<Type> result(r,c);
  // Second read to really read the matrix
  {
	std::ifstream ifs(filename);
	r = 0;
	
	while(line.readLine(ifs))
	{ 
	  if(line.isNull())break;

	  int l = line.length();
	  if (l<5) continue;

	  c = 0;
	  RWCTokenizer next(line);
	  while(!(tok = next("{},\t ")).isNull())
	  {
		double val = atof(tok);
   		result.at(r,c) = val;
		
		c++;
	  }
	  r++;
	}
  }
  *this = result;
}


*/



#endif

