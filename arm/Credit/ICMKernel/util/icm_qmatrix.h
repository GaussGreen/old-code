
#ifndef _ICM_QMATRIX_H_
#define _ICM_QMATRIX_H_

/*********************************************************************************/
/*! \class  ICM_QMatrix icm_qmatrix.h "icm_qmatrix.h"
 *  \author D Pouponneau
 *	\version 1.0
 *	\date   June 2003
 *	\file   icm_qmatrix.h
 *		\brief Creates a column matrix elements are referenced by a[j][i], 
			j= is the column's index
			i= is the row's index 
/***********************************************************************************/

#include "ICMKernel/glob/icm_enums.h" 
#include "ICMKernel/util/icm_macro.h" 

#include <memory.h>
// #include <sstream>

#include "ARMKernel/util/refvalue.h"
#include <map>


template <class T> class ICM_QCubix ; 
template <class T> class ICM_QMatrix : public ARM_Object
{

private:
	//	Internal representation 
	//
	//	row major order.
	//	elements of the same line occupy successive storage locations.

	unsigned long itsnbcols;
	unsigned long itsnbrows;
	T* itsvalues;
private:
	T& _elt(unsigned long row,unsigned long col)  { return itsvalues[col+row*itsnbcols]; } 
	const T& _elt(unsigned long row,unsigned long col)  const { return itsvalues[col+row*itsnbcols]; } 
	friend class ICM_QCubix<T> ;
public:
	//
	ICM_QMatrix() { itsnbcols=itsnbrows=0; itsvalues=0; }		

	//JLA: copy const uses const ref 
	ICM_QMatrix(const ICM_QMatrix& MatIn) : itsnbcols(MatIn.itsnbcols), itsnbrows(MatIn.itsnbrows), itsvalues(0)
	{	
		if (MatIn.itsvalues==0) return; 
		itsvalues=new T[itsnbcols*itsnbrows ]; 
		if (!itsvalues) 
			// throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ICM_QMatrix Copy Ctor: Cann't allocate !! ");
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_QMatrix Copy Ctor: Cann't allocate !! "<<itsnbcols<<"x"<<itsnbrows); 
		memcpy(itsvalues,MatIn.itsvalues,itsnbcols*itsnbrows*sizeof(T)); 
	}

	// JLA : 
	ICM_QMatrix(const unsigned long& NbRows, const unsigned long& NbCols) :  itsnbrows(0), itsnbcols(0), itsvalues(0)
	{
		if (NbRows*NbCols==0) return; 
		itsvalues=new T[NbRows*NbCols] ; 
		if (!itsvalues) 
			// throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ICM_QMatrix (rows,cols): Cann't allocate !! ");
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_QMatrix (rows,cols): Cann't allocate !! "<<NbRows<<"x"<<NbCols); 
		itsnbrows=NbRows; 
		itsnbcols=NbCols; 
		memset(itsvalues,0,itsnbrows*itsnbcols*sizeof(T)); 
	}
	
	ICM_QMatrix(const unsigned long& NbRows, 
				const unsigned long& NbCols,
				const T& value)  :  itsnbrows(0), itsnbcols(0), itsvalues(0)
	{
		if (NbRows*NbCols==0) return; 
		itsvalues=new T[NbRows*NbCols] ; 
		if (!itsvalues) 
			// throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ICM_QMatrix (rows,cols): Cann't allocate !! ");
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_QMatrix (rows,cols): Cann't allocate !! "<<NbRows<<"x"<<NbCols); 
		itsnbrows=NbRows; 
		itsnbcols=NbCols; 
		// if the value is 0, we use memset (...) instead of looping.
		if (value==0) memset(itsvalues,0,itsnbrows*itsnbcols*sizeof(T)); 
		else for(unsigned int i=0;i<itsnbrows*itsnbcols;i++) itsvalues[i]=value; 
	}
	
	//CC vector Matrix : X rows, 1 col
	ICM_QMatrix(const ARM_Vector& aARMvector): itsnbrows(0), itsnbcols(0), itsvalues(0)
	{
		if (aARMvector.GetSize()==0) return; 
		itsvalues=new T[aARMvector.GetSize()] ; 
		if (!itsvalues) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_QMatrix (rows,1): Cann't allocate !! "<<aARMvector.GetSize()); 
		itsnbrows=aARMvector.GetSize(); 
		itsnbcols= 1 ; 
		// if the value is 0, we use memset (...) instead of looping.
		for(unsigned int i=0;i<aARMvector.GetSize();i++) itsvalues[i]=aARMvector[i]; 
	}

	~ICM_QMatrix() 
	{
		if (itsvalues) delete[] itsvalues; 
		itsvalues=0; 
	}


	ICM_QMatrix& operator= (const ICM_QMatrix& ref) 
	{
		if (this!=&ref) 
		{
			this->~ICM_QMatrix(); 
			new(this)ICM_QMatrix(ref); 
		}
		return *this; 
	}

	// JLA: for compatibility 
	ICM_QMatrix(T** mat,unsigned long NbRows, unsigned long NbCols) : itsnbrows(0), itsnbcols(0), itsvalues(0)
	{
		if (NbRows*NbCols==0) return; 
		if (mat==0) return ; 
		itsvalues=new T[NbRows*NbCols] ; 
		if (!itsvalues) 
			// throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ICM_QMatrix (rows,cols): Cann't allocate !! ");
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_QMatrix (rows,cols): Cann't allocate !! "<<NbRows<<"x"<<NbCols);
		itsnbrows=NbRows; 
		itsnbcols=NbCols; 
		for(unsigned int i=0;i<itsnbrows;i++) 
			for(unsigned int j=0;j<itsnbcols;j++) 
				_elt(i,j)=mat[i][j]; 
	}

	bool IsEmpty() const { return (itsvalues==0) ; } 
	bool IsSquare() const { return (itsnbcols==itsnbrows) ; } 

	inline bool IsSymmetric() const 
	{	bool res = true;
		if(!IsSquare()) return false;
		for ( int i=0; i< itsnbcols; i++){
			for (int j=0; j < i ; j++){
				if (_elt(i,j) != _elt(j,i))
					return false;
			}
		}
		return res;
	}
	unsigned long GetSize(void) {return itsnbrows;}//for compatibilty with ARM_Vector
	T Elt(const unsigned long& i) {return _elt(i,0);}

	void SetCol(const unsigned long& Col, const ARM_Vector* vector) 
	{
		if ((vector==0) || (vector->GetSize() != itsnbrows) ||  (Col>=itsnbcols))
	        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "SetCol(ARM_Vector) : Incompatible arguments");
		for(long i=0; i< itsnbrows; i++)	_elt(i,Col) = vector->Elt(i);
	}

	inline void SetCol(const unsigned long& Col, const T* vector) 
	{
		if ( (Col>=itsnbcols) || (vector==0 ) )
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_QMatrix::SetCol("<<Col<<","<<vector<<"): Bad Argument"); 
		for(long i=0; i< itsnbrows; i++) _elt(i,Col)=vector[i]; 
	}
	inline void SetCol(const unsigned long& Col,const T& value)
	{
		if (Col>=itsnbcols)
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_QMatrix::SetCol("<<Col<<"): Bad Argument"); 
		for(long i=0; i< itsnbrows; i++) _elt(i,Col)=value; 
	}
	inline void SetRow(const unsigned long& Row,const ARM_Vector* vector) 
	{
		if ((vector==0) ||(vector->GetSize() != itsnbcols) || (Row >= itsnbrows) )
	        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "SetRow(ARM_Vector) : Incompatible argument");

		const T* prow = vector->GetElt();
		for(unsigned int j=0;j<itsnbcols;j++) 
			_elt(Row,j)=prow[j]; 
		// memcpy(itsvalues[Row],prow,sizeof(T)*itsnbcols);
	}
	inline void SetRow(const unsigned long& Row,const T* vector) 
	{
		if (Row >= itsnbrows || (vector==0))
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "SetRow(T*) : Incompatible Tab");

		for(unsigned int j=0;j<itsnbcols;j++) 
			_elt(Row,j)=vector[j]; 
	}
	inline void SetRow(const unsigned long& Row,const T& value) 
	{
		if (Row >= itsnbrows )
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_QMatrix::SetRow("<<Row<<"): Bad Argument"); 

		for(unsigned int j=0;j<itsnbcols;j++) 
			_elt(Row,j)=value; 
	}

	// JLA : resize and set to 0. 
	inline void Resize(const unsigned long& nbrows,const unsigned long& nbcols) 
	{
		//	if any of the size is null, we get the empty object.
		if (nbrows*nbcols==0) 
		{ 
			itsnbcols=itsnbrows=0; 
			if (itsvalues) delete[]itsvalues; 
			itsvalues=0; 
			// Init(); return; 
		}

		//  if same size: realloc is not needed
		if (nbrows*nbcols != itsnbcols*itsnbrows) 
		{
			if (itsvalues) delete[] itsvalues ; 
			itsvalues=0; itsnbcols=0; itsnbrows=0; 
			itsvalues=new T[nbrows*nbcols]; 
			if (!itsvalues) 
				// throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ICM_QMatrix: Cann't allocate !! ");
				ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_QMatrix: Cann't allocate !! "<<nbrows<<"x"<<nbcols); 
		}
		itsnbcols = nbcols;
		itsnbrows = nbrows;
		//	this step can be movec upper
		memset(itsvalues,0,itsnbcols*itsnbrows*sizeof(T));		
	}
	
	

	ICM_QMatrix& operator +=(const ICM_QMatrix& MatIn) 
	{
		if ((itsnbrows!=MatIn.itsnbrows) || (itsnbcols!=MatIn.itsnbcols))
	        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		        "Incompatible Matrix for +=");
		for(unsigned int i=0;i<itsnbrows*itsnbcols;i++) 
			itsvalues[i]+=MatIt.itsvalues[i]; 
		return (*this); 
	}

	ICM_QMatrix& operator /=(const double& dDenom) 
	{
		if (dDenom == 0.0)
	        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		        "Division by 0 !!");
		for(unsigned int i=0;i<itsnbrows*itsnbcols;i++) 
			itsvalues[i]/= dDenom; 
		return (*this); 
	}
		
	T& operator()(unsigned long row,unsigned long col) 
	{
		if ( (row>= itsnbrows) || (col >= itsnbcols) ) 
		//	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        //    "operator(): Incompatible parameters");
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_QMatrix::op("<<row<<","<<col<<"): out of size "<<itsnbrows<<"x"<<itsnbcols);
		return itsvalues[col+row*itsnbcols] ; 
		
	}
	const T& operator()(unsigned long row,unsigned long col) const 
	{ 
		if ( (row>= itsnbrows) || (col >= itsnbcols) ) 
			// throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            // "operator(): Incompatible parameters");
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_QMatrix::op("<<row<<","<<col<<"): out of size "<<itsnbrows<<"x"<<itsnbcols);
		return itsvalues[col+row*itsnbcols] ; 
	}

	unsigned long  Getnbcols()	const{ return itsnbcols;	 }	
	unsigned long  Getnbrows()	const{ return itsnbrows;	 }	

	// JLA : this is for compatibility. Do not use 
	void SetValue(unsigned long row,unsigned long col, const T& value)	
	{ 
		(*this)(row,col)=value; 
	}
	const T& Getvalue(unsigned long row,unsigned long col) const		
	{ 
		return (*this)(row,col); 
	} 

	//JLA : Fill all content with the same value.
	void Fill(double v) 
	{ 
		for(unsigned int i=0;i<itsnbrows*itsnbcols;i++) itsvalues[i]=v; 
	}
	// JLA 
	inline ARM_Vector* RowAsVector(const unsigned long& row) const 
	{
		return TruncRowAsVector(row,itsnbcols); 
	}
	// JLA 
	inline ARM_Vector* TruncRowAsVector(const unsigned long& row,const unsigned long& size) const 
	{
		long i=0;
		if ((row>=itsnbrows) || size==0) return NULL;
		unsigned long actual_size = MIN(itsnbcols,size); 
		ARM_Vector* vector = new ARM_Vector(actual_size,0.);	
		for (i=0;i<actual_size;i++)	vector->Elt(i)=_elt(row,i); 
		return (vector);
	}
	// JLA: added const , changed to unsigned arg
	ARM_Vector* ColAsVector(const unsigned long& col) const 
	{
		long i=0;
		if ( col>=itsnbcols)	return NULL;
		ARM_Vector* vector = new ARM_Vector(itsnbrows,0.);
		for (i=0;i<itsnbrows;i++)	vector->Elt(i)=_elt(i,col); 
		return (vector);
	}

	// Onin
	vector<T> ColAsStdVector(const unsigned long& col) const 
	{
		long i=0;
// FIXMEFRED: mig.vc8 (28/05/2007 10:40:48):cast
		if ( col>=itsnbcols)	return static_cast<vector<T> >(NULL);
		vector<T> vector; vector.resize(itsnbrows);
		for (i=0;i<itsnbrows;i++)	vector[i]=_elt(i,col); 
		return (vector);
	}

	vector<T> RowAsStdVector(const unsigned long& row) const 
	{
		long i=0;
// FIXMEFRED: mig.vc8 (28/05/2007 10:40:56):cast
		if ( row>=itsnbrows)	return static_cast<vector<T> >(NULL);
		vector<T> vector; vector.resize(itsnbcols);
		for (i=0;i<itsnbcols;i++)	vector[i]=_elt(row,i); 
		return (vector);
	}

	// Returning copy of the specified row 
	T* DupRowAsTab(const unsigned long& row) const 
	{
		long i=0;
		if (row>=itsnbrows)	return NULL;
		T* vector = new T[itsnbcols];
		for (i=0;i<itsnbcols;i++)	vector[i] = _elt(row,i); 
		return (vector);
	}

	// Returning copy of the specified column
	//	JLA: added const, changed arg to unsigned
	T* ColAsTab(const unsigned long& col) const
	{
		long i=0;
		if (col>=itsnbcols)	return 0;
		T* vector = new T[itsnbrows];
		for (i=0;i<itsnbrows;i++)	vector[i] = _elt(i,col); 
		return (vector);
	}

	//JLA 
	void TransposeMe(void)
	{
		ICM_QMatrix copy(*this); 
		unsigned long tmp; 
		tmp=itsnbcols; 
		itsnbcols=itsnbrows;
		itsnbrows=tmp; 
		for(unsigned int i=0;i<itsnbrows;i++) 
			for(unsigned int j=0;j<itsnbcols;j++) 
				_elt(i,j)=copy._elt(j,i); 
	}

	//CC 
	ICM_QMatrix Transpose(void) const
	{
		ICM_QMatrix Qres(itsnbcols,itsnbrows,0.0);  	
		for(unsigned int i=0;i<itsnbcols ;i++) 
			for(unsigned int j=0;j<itsnbrows;j++) 
				Qres._elt(i,j)=_elt(j,i); 
		return Qres;
	}

	//CC (*this) = (m)*(m.Transpose())
	inline void MultTranspose(const ICM_QMatrix& m)
	{
		unsigned int   i=0, j=0, k=0;
		T s = 0.0;
		Resize(m.Getnbrows(), m.Getnbrows());
 		if (m.Getnbrows() == 0 || m.Getnbcols()==0 )
		{
			ICMTHROW( ERR_SQUARE_OR_SIZE_PB, "empty matrix : size =0");
		}
	   //    compute product
 		for (i = 0; i < m.Getnbrows(); i++)
		{
			for (j = 0; j < m.Getnbrows(); j++)
			{
				s = 0.0;
				for (k=0; k<m.Getnbcols(); k++)
	                s += m(i,k) * m(j,k);
			    _elt(i,j) = s;
			}
		}
 	}

	//CC  ICM_QMatrix = (*this)*(this->Transpose())
	inline ICM_QMatrix  MultTranspose() const
	{
		ICM_QMatrix Qres(itsnbrows,itsnbrows,0.0);
		unsigned int    i=0, j=0, k=0;
		T s=0.0;
 		if (itsnbcols == 0 || itsnbrows==0 )
		{
        throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                         "empty matrix : size =0");
		}
	   //    compute product
 		for (i = 0; i <itsnbrows; i++)
		{
			for (j = 0; j < itsnbrows; j++)
			{
				s = 0.0;
				for (k=0; k<itsnbcols; k++)
	                s += _elt(i,k) * _elt(j,k);
			    Qres._elt(i, j) = s;
			}
		}
		return Qres;
 	}

	void RAZ(const T& value)
	{
		for(unsigned int i=0;i<itsnbrows*itsnbcols;i++) itsvalues[i]=value;
	}

	// ----------------------------
	//	Copy of members data
	// ----------------------------
	void BitwiseCopy(const ARM_Object* src)
	{
	    const ICM_QMatrix* matrix = dynamic_cast<const ICM_QMatrix *>(src);
		(*this)=(*matrix); 
	}

	// -------------
	//	Copy Method 
	// -------------
	void Copy(const ARM_Object* src)
	{
		ARM_Object::Copy(src);
		BitwiseCopy(src);
	}

	// --------------
	//	Clone Method
	// --------------
	ARM_Object* Clone(void)				{return new ICM_QMatrix(*this); }
// 	ICM_QMatrix<T> Clone(void) const	{ return new ICM_QMatrix(*this); }


	/* ----------------------------------------------------------------------------------------
	Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
	permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
	indx[1..n] is an output vector that records the row permutation eected by the partial
	pivoting; d is output as .1 depending on whether the number of row interchanges was even
	or odd, respectively. This routine is used in combination with lubksb to solve linear equations
	or invert a matrix.
	 ---------------------------------------------------------------------------------------- */

	void LD_NR(ARM_Vector *indx, double &d)
	{
	    int i,imax,j,k;
		double big = 0.0, dum = 0.0, sum = 0.0, temp = 0.0;
		ARM_Vector    vv;
        
		// Check input

		if (!IsSquare() || indx->GetSize() != itsnbrows)
			ICMTHROW(ERR_INVALID_ARGUMENT,"Matrix is not square or Size != Number of lines");
            
		// Resize vv
		vv.Resize(itsnbrows);
		d = 1.0;
    
		for (i = 0; i < itsnbrows; i++)
		{
			big = 0.0;

			for (j = 0; j < itsnbrows; j++)
				if ((temp = fabs(_elt(i,j))) > big) big = temp;
				
			if ( big == 0.0 ) ICMTHROW(ERR_INVALID_ARGUMENT,"Problem in Ludcmp()");
			vv[i] = 1.0/big;
		}

		for	(j = 0; j < itsnbrows; j++)
		{
			for (i = 0; i < j; i++)
			{
				sum = _elt(i, j);

				for (k = 0; k < i; k++) 
					sum -= _elt(i, k) * _elt(k, j);
			

				_elt(i, j) = sum;
			}

	        big = 0.0;

		    for (i = j; i < itsnbrows; i++)
			{
				sum = _elt(i, j);

				for ( k = 0; k < j; k++) 
					sum -= _elt(i, k) * _elt(k, j);
				

				_elt(i, j) = sum;

				if ( (dum = vv[i] * fabs(sum)) >= big)
				{big=dum;imax=i;}
			}

			if ( j != imax )
			{
				for (k = 0; k < itsnbrows; k++)
				{
					dum = _elt(imax, k);
					_elt(imax, k) = _elt(j, k);
					_elt(j, k) = dum;
				}

				d = -d;
				vv[imax] = vv[j];
			}

			(*indx)[j] = (double) imax;

			if ( _elt(j, j) == 0.0 ) _elt(j, j) = K_TINY;

			if ( j != itsnbrows-1 )
			{	
				dum=1.0/_elt(j, j);
				for (i = j+1; i < itsnbrows; i++) 
					_elt(i, j) *=  dum;
			}
		}
	}	


	/*-----------------------------------------------------------------------------------------------
	Solves the set of n linear equations A.X = B. Here a[1..n][1..n] is input, not as the matrix
	A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
	as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
	B, and returns with the solution vector X. a, n, and indx are not modied by this routine
	and can be left in place for successive calls with dierent right-hand sides b. This routine takes
	into account the possibility that b will begin with many zero elements, so it is e.cient for use
	in matrix inversion.
	------------------------------------------------------------------------------------------------*/

	void LU_NR(ARM_Vector *indx, ARM_Vector *b)
	{
		int i,ii=-1,ip,j;
		double sum;

	    if (!IsSquare() || ( indx->GetSize() != itsnbrows ) || ( b->GetSize() != itsnbrows )) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"Matrix is not square or Size != Number of lines");
    
		for (i = 0; i < itsnbrows; i++) 
		{
			ip = (int) (*indx)[i];
			sum = (*b)[ip];
	        (*b)[ip] = (*b)[i];

		    if ( ii != -1 ) for ( j = ii; j <= i-1; j++) sum -= _elt(i, j) * (*b)[j];
			else if (sum) ii = i;

		    (*b)[i] = sum;
		}
    
		for (i=itsnbrows-1;i>=0;i--) 
		{
			sum = (*b)[i];
			for (j = i+1; j < itsnbrows;j++) sum -= _elt(i, j) * (*b)[j];
			(*b)[i] = sum/_elt(i, i);
		}
	}

	/*----------------------------------------------------------------------------*
	returns invert of matrix computed from LU decomposition
	*----------------------------------------------------------------------------*/
 
	inline ICM_QMatrix& Invert(double &det)
	{
		int    i, j;
		ICM_QMatrix tmp;
		ARM_Vector indx(itsnbrows);
		ARM_Vector col(itsnbrows);
  
	   //    check input
 
	    if (!IsSquare()) ICMTHROW(ERR_INVALID_ARGUMENT,"Matrix is not square in operator ^");
		ICM_QMatrix invert(itsnbrows,itsnbcols);
  
		//    copy this into tmp
	    tmp = *this;
 
		//    perform LU decomposition
	    try
		{
			tmp.LD_NR(&indx, det);
		}
	    catch(MathException& m)
		{
        m.DebugPrint();
         ICMTHROW(ERR_INVALID_ARGUMENT,"Problem in Matrix inversion, Invert()");
		}
 
	   for ( j = 0; j < itsnbrows; j++)
        det *= tmp(j, j);
	
 	   if ( fabs(det) < K_DOUBLE_TOL )
        ICMTHROW(ERR_INVALID_ARGUMENT,"Problem in Matrix inversion, Invert()");
 
		//    perform backsubstitution
 
		for ( j = 0; j < itsnbrows; j++)
		{
			for (i=0; i<itsnbrows; i++)
				col[i] = 0.0;
 
        col[j] = 1.0;
 
        try
        {
            tmp.LU_NR(&indx, &col);
        }
 
        catch(MathException& m)
        {
            m.DebugPrint();
 
            ICMTHROW(ERR_INVALID_ARGUMENT,"Problem in Matrix inversion, Invert()");
        }
 
        for (i=0; i < itsnbrows; i++)
            invert(i, j) = col[i];

		}
    
		*this = invert;

		return (*this);
    }


	/*----------------------------------------------------------------------------*
	Solves linear system a x = y
    where a is the current matrix. Uses LU decomposition.
    The result x is returned into the vector X. 
    Also returns the matrix determinant. Returns 0 if failed.
	*----------------------------------------------------------------------------*/
            
	void LinearSolve(ARM_Vector* X,double& det)
	{
    int        j;
    ICM_QMatrix tmp;
    ARM_Vector indx(itsnbrows);
 
    if ((!IsSquare()) && (X->GetSize() != itsnbrows))
       ICMTHROW(ERR_INVALID_ARGUMENT,"Matrix is not square or Size != Number of lines");
        
    // Perform LU decomposition and compute determinant

	tmp = *this;

    try
    {
        tmp.LD_NR(&indx, det);
    }

    catch(MathException& m)
    {
        m.DebugPrint(); 
 
        throw MathException(__LINE__, __FILE__, ERR_MATRIX_LIN_SOLVE_PB,
                 "Problem in LinSolve()");
    }

    for (j = 0; j < tmp.itsnbrows; j++) 
        det *= tmp(j, j);
    
    // Perform backsubstitution

    try 
    {
        tmp.LU_NR(&indx, X);
    }
    catch(MathException& m)
    {
        m.DebugPrint();
 
        throw MathException(__LINE__, __FILE__, ERR_MATRIX_LIN_SOLVE_PB,
                 "Problem in LinSolve()");
    }
	}

	/* ---------------------------------------------------------------------------
	Matricial Product
	------------------------------------------------------------------------------*/
	inline void Mult(const ICM_QMatrix &m1, const ICM_QMatrix &m2)
	{
		unsigned int  i=0, j=0, k=0;
		T s=0.0;
		Resize(m1.Getnbrows(), m2.Getnbcols());
 		if ( m1.Getnbcols() != m2.Getnbrows() )
		{
        throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                         "Inconsistent matrix sizes in operator *");
		}
 
	   //    compute product
 		for (i = 0; i < m1.Getnbrows(); i++)
		{
			for (j = 0; j < m2.Getnbcols(); j++)
			{
				s = 0.0;
 
				for (k=0; k<m1.Getnbcols(); k++)
	                s += m1(i,k) * m2(k,j);
 
			    _elt(i, j) = s;
			}
		}
 	}

	/*----------------------------------------------------------------------------*/
	/*    CholeskyDcmp(ICM_QMatrix& M (output))                                             */
	/*    Performs Cholesky decomposition. This is from Numerical Recipes in C.   */
	/*----------------------------------------------------------------------------*/

	/** Removed JLA 
	inline CholeskyDcmp(ICM_QMatrix& M)
	{

		if (!IsSquare())
			ICMTHROW(ERR_SQUARE_OR_SIZE_PB, "Matrix must be Square");
            
		int    i, j, k;
	    int numLines = Getnbrows();
		M.Resize(numLines,numLines);
		memset(M.itsvalues,'\0',numLines*numLines*sizeof(double));
	    double sum;

		for (i = 0; i<numLines; i++)
			for (j = i; j<numLines; j++) 
			{
				for (sum = _elt(i, j), k=i-1; k>=0; k--) 
                 sum -= M(i, k)*M(j, k);  
            
				if (i == j)
				{
					if (sum<=0.0)
					{
					throw MathException(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
						"Cholesky Decomp failed: Matrix is not positive definite");
					}
					M(i, i) = sqrt(sum);
				}
				else
	                M(j, i) = sum/_elt(i, i);
		        
	        }
    
	}
	**/
	/*----------------------------------------------------------------------------*/
	/*    Sort Matrix by ascending or descending order on a given column					                                          */
	/*----------------------------------------------------------------------------*/
	inline void Sort(unsigned long ColToSort, bool ascending = true)
	{
		if (ColToSort >  itsnbcols)
			ICMTHROW(ERR_INVALID_ARGUMENT,"Sort: Column "<<ColToSort<<" is not in the range"); 

		std::map<double,int> maptemp;

		ICM_QMatrix<double> MatTemp (*this);
		std::vector<int> Indice;
		Indice.resize(itsnbrows);

		for (int i =0;i<itsnbrows;i++)
			maptemp[MatTemp._elt(i,ColToSort)]=i;
		
		map<double,int>::iterator it;
		i = 0;
		for (it=maptemp.begin();it != maptemp.end(); ++it)
		{
			Indice[i] = it->second;
			i++;
		}

		if (ascending)
		{
			for (i =0;i<itsnbrows;i++)
				for (int j=0;j<itsnbcols;j++)
					(*this)(i,j) = MatTemp._elt(Indice[i],j);
		}
		else 
		{
			for (i =0;i<itsnbrows;i++)
				for (int j=0;j<itsnbcols;j++)
					(*this)(i,j) = MatTemp._elt(Indice[itsnbrows - 1 - i],j);
		}
	}

	/*----------------------------------------------------------------------------------------------*/
	/*    Permut row and column of a square and symmetric Matrix, according to a vector of new order				                                          */
	/*----------------------------------------------------------------------------------------------*/
	inline void Permut(const ICM_QMatrix<T>& QMatrix, const ARM_Vector& vNewOrder){
		// Check or the square and symmetric Matrix :
		if (!QMatrix.IsSymmetric())
			ICMTHROW(ERR_SQUARE_OR_SIZE_PB,"can't use this Permut for a no symmetric Matrix ");
		if (vNewOrder.size() != QMatrix.Getnbrows() )
			ICMTHROW(ERR_SQUARE_OR_SIZE_PB,"can't use this Permut for a vectorPremut size different from Matrix size ");
		for ( int i=0; i<itsnbcols ; i++){
			for ( int j=0; j<i ; j++){
				_elt(i,j) = QMatrix._elt(vNewOrder[i],vNewOrder[j]);
				_elt(j,i) = _elt(i,j);
			}
			_elt(i,i) = QMatrix._elt(i,i);
		}
	}

	/*----------------------------------------------------------------------------*/
	/*    View							                                          */
	/*----------------------------------------------------------------------------*/
	void View(char* id, FILE* ficOut)
	{
		FILE* fOut;
		char  fOutName[200];
		int i = 0,k = 0;

		if ( ficOut == NULL )
		{
			ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
			(void) unlink(fOutName);
			fOut = fopen(fOutName, "w"); 
		}
		else
		{
			fOut = ficOut;
		} 

		fprintf(fOut, "\n ======> QMatrix :\n\n");
		fprintf(fOut, "\n");

		for (i = 0; i<itsnbrows; i++)
		{
			for (k = 0; k<itsnbcols; k++) 
			{	
				std::stringstream sstr; sstr<<itsvalues[k+itsnbcols*i]; 
				// fprintf(fOut, "\t %f ", itsvalues[k+itsnbcols*i]  ); 
				fprintf(fOut, "\t %s ", sstr.str().c_str()  ); 
			}

		fprintf(fOut, "\n");
		}

		if ( ficOut == NULL ) fclose(fOut);
	}

	inline void GenCombinations()
	{
		if ( (itsnbrows!=itsnbcols) || (itsnbrows==0) )
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_QMatrix::GenCombinations("<<itsnbrows<<","<<itsnbcols<<"): wrong size")

		int k,n;
		_elt(0,0)=1;

		for(n=1;n<itsnbrows;n++)
		{
			_elt(n,0)=1;
			_elt(n,n)=1;
			for(k=1;k<n;k++)
			{_elt(n,k)=_elt(n-1,k-1)+_elt(n-1,k); 
			if (_elt(n,k)<0) ICMTHROW(ERR_INVALID_ARGUMENT,"combination>Max_Size_Long");}
		}
	}

};


/*********************************************************************************/
/*! \class  ICM_QCubix icm_qmatrix.h "icm_qmatrix.h"
 *  \author D Pouponneau
 *	\version 1.0
 *	\date   June 2003
 *	\file   icm_qmatrix.h
 *	\brief Creates a Cube 
/***********************************************************************************/

template <class T> class ICM_QCubix : public ARM_Object
{

private:
	unsigned long itsDim;
	unsigned long itsRows,itsCols; 
	ICM_QMatrix<T>** itsVMatrix;

	inline void Init(void)
	{
		itsDim = itsRows=itsCols=0; 
		itsVMatrix = NULL;
	}

public :

	ICM_QCubix(const int& i,const int& j,const int& k,const T& value)
	{
		Init();

		Set(i,j,k,value);

	}


	//	ATTENTION il ne s'agit pas de setter l'élément (i,jk)
	//	mais de dimensionner le cube à K matrices de dim (I,J)
	inline void Set(const unsigned long & i,const unsigned long & j,const unsigned long & k, const T& value)
	{
		int il =0;
		itsDim = k;
		itsRows=i; 
		itsCols=j; 
		if (itsVMatrix)
		{
			for (int il2=0;il2<itsDim;il2++)
				delete itsVMatrix[il2];

			delete[] itsVMatrix;
			itsVMatrix=0; 
		}

		itsVMatrix = new ICM_QMatrix<T>*[itsDim];	
		if (!itsVMatrix) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_QCubix:Set: can't allocate "); 

		for (il=0; il<itsDim; il++)
			itsVMatrix[il] = new ICM_QMatrix<T>(itsRows,itsCols,value);

	}

	/**
	const ICM_QMatrix<T>* operator[] (unsigned long i) const
	{
		if (i >= itsDim)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Qubix::()"); 
		return(itsVMatrix[i]);
	}
	**/ 

	void SetElt(const unsigned long& i,const unsigned long& j,const unsigned long& k,const T& value)
	{
		if (k >= itsDim || i>=itsRows || j>=itsCols)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Qubix::()"); 
		itsVMatrix[k]->_elt(i,j)=value; 
	}

	const T& Elt(const unsigned long& i,const unsigned long& j,const unsigned long& k) const
	{
		if (k >= itsDim || i>=itsRows || j>=itsCols)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Qubix::()"); 
		// return (*itsVMatrix[k] )(i,j);
		return itsVMatrix[k]->_elt(i,j); 
	}
	T& Elt(const unsigned long& i,const unsigned long& j,const unsigned long& k) 
	{
		if (k >= itsDim || i>=itsRows || j>=itsCols)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Qubix::()"); 
		return itsVMatrix[k]->_elt(i,j);
	}

	T& operator()(unsigned long i, unsigned long j, unsigned long k) 
	{
		if (k >= itsDim || i>=itsRows || j>=itsCols) 
			// ICMTHROW(ERR_INVALID_ARGUMENT,"Elt :Out of size for ICM_QCubix");
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Qubix::()"); 
		// if (k >= itsDim) 
		// 	ICMTHROW(ERR_INVALID_ARGUMENT,"operator("<<k<<"):Out of size for ICM_QCubix");
		// ICM_QMatrix<double>* mat = itsVMatrix[k] ;
		// if (i>=mat->Getnbrows() || j>=mat->Getnbcols()) 
		// 	ICMTHROW(ERR_INVALID_ARGUMENT,"operator("<<k<<"):Out of size for ICM_QCubix");
		return itsVMatrix[k]->_elt(i,j); 
	}
	const T& operator()(unsigned long i, unsigned long j, unsigned long k) const
	{
		if (k >= itsDim || i>=itsRows || j>=itsCols)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Qubix::()"); 
		return mat->_elt(i,j); 
	}

	~ICM_QCubix()
	{
		if (itsVMatrix)
		{
			for (int i=0;i<itsDim;i++) delete itsVMatrix[i];
			delete[] itsVMatrix; 
			itsVMatrix=0; 
		}
	}

	inline void ResizeWithCopy(const unsigned long& size, const T& value)
	{
		// We do not resize if the size is the same 
		if (size==itsDim) return ; 

		//	If size==0 , we get the empty matrix
		if (size==0) 
		{
			for (int i=0;i<itsDim;i++) delete itsVMatrix[i];
			delete[] itsVMatrix;
			itsVMatrix=0; 
			itsRows=itsCols=itsDim=0; 
			return ;
		}
		
		//	Otherwise we properly resize 
		ICM_QMatrix<T>** Matrix = new ICM_QMatrix<T>*[size];
		int i =0;

		for (i=0; i<MIN(size,itsDim); i++)
			Matrix[i] = itsVMatrix[i]; 

		if (size>itsDim)
		{
			// int size1 = (itsVMatrix[0])->Getnbrows();
			// int size2 = (itsVMatrix[0])->Getnbcols();
			for (i=itsDim; i<size; i++)
				Matrix[i] = new ICM_QMatrix<T>(itsRows,itsCols,value); 
		}
		else
		{
			for (i=size; i<itsDim; i++)
				{
				delete itsVMatrix[i]; 
				itsVMatrix[i] = NULL; 
				}
		}

		delete[] itsVMatrix;

		itsVMatrix = Matrix;

		itsDim = size;
	}

	inline void ResizeWithInitialize(const unsigned long& size, const T& value)
	{
		//	If size==0 , we get the empty matrix
		if (size==0) 
		{
			for (int i=0;i<itsDim;i++) 
				{if (itsVMatrix[i]) delete itsVMatrix[i];
				itsVMatrix[i]=NULL;
				}
			delete[] itsVMatrix; 
			itsVMatrix=NULL; 
			itsRows=itsCols=itsDim=0; 
			return ;
		}
		
		//	Otherwise we properly resize 
		ICM_QMatrix<T>** Matrix = new ICM_QMatrix<T>*[size];
		int i =0;

		for (i=0; i<size; i++)
			Matrix[i] = new ICM_QMatrix<T>(itsRows,itsCols,value); 

		if (itsVMatrix)
		{
			for (int i=0;i<itsDim;i++) 
				{if (itsVMatrix[i]) delete itsVMatrix[i];
				itsVMatrix[i]=NULL;
				}
			delete[] itsVMatrix; 
			itsVMatrix=NULL; 
		}
		itsVMatrix = Matrix;

		itsDim = size;
	}


	// CN : return an axe of the cubix (for numerical integration)
	inline ARM_Vector* GetRowVector(const unsigned long& col,const unsigned long& dim) const 
	{
		//unsigned long itsDim;
		//unsigned long itsRows,itsCols; 

		unsigned long i=0;
		if (col>=itsCols || dim>=itsDim || itsRows==0) 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Qubix::()"); 
		
		ARM_Vector* vector = new ARM_Vector(itsRows,0.);	
		for (i=0;i<itsRows;i++)	vector->Elt(i)= this->Elt(i,col,dim); 
		return (vector);
	}
	
	ARM_Vector* GetColVector(const unsigned long& row, const unsigned long& dim) const 
	{
		unsigned long i=0;
		if (row>=itsRows || dim>=itsDim || itsCols==0) 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Qubix::()"); 
		
		ARM_Vector* vector = new ARM_Vector(itsCols,0.);	
		for (i=0;i<itsCols;i++)	vector->Elt(i)= this->Elt(row,i,dim); 
		return (vector);
	}

	ARM_Vector* GetDimVector(const unsigned long& row, const unsigned long& col) const 
	{
		unsigned long i=0;
		if (row>=itsRows || col>=itsCols || itsDim==0) 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Qubix::()"); 
		
		ARM_Vector* vector = new ARM_Vector(itsDim,0.);	
		for (i=0;i<itsDim;i++)	vector->Elt(i)= this->Elt(row,col,i); 
		return (vector);
		}

	// CN : return an axe of the cubix (for numerical integration) in a VectorType
	inline vector<double> GetRowVectorV(const unsigned long& col,const unsigned long& dim) const 
	{
		unsigned long i=0;
		if (col>=itsCols || dim>=itsDim || itsRows==0) 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Qubix::()"); 
		
		std::vector<double> res(itsRows);	
		for (i=0;i<itsRows;i++)	res[i]= this->Elt(i,col,dim); 
		return (res);
	}
	
	vector<double> GetColVectorV(const unsigned long& row, const unsigned long& dim) const 
	{
		unsigned long i=0;
		if (row>=itsRows || dim>=itsDim || itsCols==0) 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Qubix::()"); 
		
		std::vector<double> res(itsCols);	
		for (i=0;i<itsCols;i++)	res[i]= this->Elt(row,i,dim); 
		return (res);
	}

	std::vector<double> GetDimVectorV(const unsigned long& row, const unsigned long& col) const 
	{
		unsigned long i=0;
		if (row>=itsRows || col>=itsCols || itsDim==0) 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Qubix::()"); 
		
		std::vector<double> res(itsDim);	
		for (i=0;i<itsDim;i++)	res[i]= this->Elt(row,col,i); 
		return (res);
	}

	// LJ: Set all elements to a single value
	void SetValue(const T& value)
	{
		for(unsigned int i=0;i<itsnbrows;i++) 
			for(unsigned int j=0;j<itsnbcols;j++) 
				_elt(i,j)=value; 
	}

	/*----------------------------------------------------------------------------*/
	/*    View							                                          */
	/*----------------------------------------------------------------------------*/
	void View(char* id, FILE* ficOut)
	{
		FILE* fOut;
		char  fOutName[200];

		if ( ficOut == NULL )
		{
			ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
			(void) unlink(fOutName);
			fOut = fopen(fOutName, "w"); 
		}
		else
		{
			fOut = ficOut;
		} 

		fprintf(fOut, "\n ======> QCubix:\n\n");
		fprintf(fOut, "\n");

		int	i, j, k;
		int	itsnbrows;
		int	itsnbcols;

		ICM_QMatrix<double>*	Current_Matrix;

		for (j=0; j<itsDim; j++)
		{
			fprintf(fOut, "\n\t==>\tQMatrix:\t%u\n", j);

			Current_Matrix	=	itsVMatrix[j];
			
			itsnbrows	=	Current_Matrix->Getnbrows();
			itsnbcols	=	Current_Matrix->Getnbcols();

			for (i=0; i<itsnbrows; i++)
			{
				for (k=0; k<itsnbcols; k++) 
					fprintf(fOut, "\t %f ", (*Current_Matrix)(i,k));

				fprintf(fOut, "\n");
			}
			
		}

		if ( ficOut == NULL ) fclose(fOut);
	}

};


#endif 
