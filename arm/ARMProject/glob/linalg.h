#ifndef _LINALG_H
#define _LINALG_H


#include <math.h>
#include <vector>

#include "armglob.h" 
#include "expt.h"
#include "matrixpp.h"

#ifdef WIN32
	#define unlink _unlink
#endif

#ifdef unix
#    include <unistd.h>
#    include <fcntl.h>
#    include <sys/stat.h>
#endif


using std::vector;



/*----------------------------------------------------------------------------*/


class ARM_Matrix; 
class ARM_TDiag; 


typedef enum {
	ARM_LIGNE,
	ARM_COL
} ARM_LignOrCol;



class ARM_GenMatrix : public ARM_Object
{
    public:

        ARM_GenMatrix(void) 
        {
            SetName(ARM_GEN_MATRIX);
        }

        ARM_GenMatrix(const ARM_GenMatrix& m) : ARM_Object(m)
        {
            SetName(ARM_GEN_MATRIX);
        }

        ARM_GenMatrix& operator = (const ARM_GenMatrix& m) 
        { 
            (*this).ARM_Object::operator=(m);

            return(*this); 
        }

       ~ARM_GenMatrix(void)
        {
        }

        void Copy(const ARM_Object* src)
        {
            ARM_Object::Copy(src);
        }
 
        ARM_Object* Clone(void)
        {
            ARM_GenMatrix* theClone = new ARM_GenMatrix();
 
 
            theClone->Copy(this);
 
            return(theClone);
        }
    
        virtual int IsSquare(void) 
        { 
            return(GetNumLines() == GetNumCols()); 
        }
    
        int IsSymmetric(void) 
        { 
            return(GetNumLines() == GetNumCols()); 
        }
    
    
    virtual int GetNumLines(void) const;

    virtual int GetNumCols(void) const;

    virtual double& Elt(int, int) const 
    {
         static double v = 0;

         printf("\n\n ==> #####????*** : BUG ARM_GenMatrix::Elt(i,j) should never be called\n\n");

         return(v);
    }

    virtual int SetToId(void);

    virtual double Trace(void);

    // Operators 

    ARM_GenMatrix& operator *= (double);
    ARM_GenMatrix& operator /= (double);
    ARM_GenMatrix& operator += (double);
    ARM_GenMatrix& operator -= (double);
    ARM_GenMatrix& operator += (const ARM_GenMatrix &);
    ARM_GenMatrix& operator -= (const ARM_GenMatrix &);


    friend int operator == (const ARM_GenMatrix &, const ARM_GenMatrix &);
    friend int operator != (const ARM_GenMatrix &, const ARM_GenMatrix &);
    friend int operator > (const ARM_GenMatrix &, const ARM_GenMatrix &);
    friend int operator < (const ARM_GenMatrix &, const ARM_GenMatrix &);
    friend int operator >= (const ARM_GenMatrix &, const ARM_GenMatrix &);
    friend int operator <= (const ARM_GenMatrix &, const ARM_GenMatrix &);


    friend ARM_Matrix& operator * (const ARM_GenMatrix &, ARM_GenMatrix &);
    friend ARM_Matrix& operator ^ (const ARM_GenMatrix &, int);

    friend ARM_Matrix operator | (const ARM_GenMatrix &, const ARM_GenMatrix &);
    friend ARM_Matrix operator & (const ARM_GenMatrix &, const ARM_GenMatrix &);

    friend ARM_Matrix log(const ARM_GenMatrix &);
    friend ARM_Matrix exp(const ARM_GenMatrix &);
    friend ARM_Matrix fabs(const ARM_GenMatrix &);
    friend ARM_Matrix sqrt(const ARM_GenMatrix &);
    friend ARM_Matrix floor(const ARM_GenMatrix &);

    friend int isequal(const ARM_GenMatrix &, double);
    friend int isdiff(const ARM_GenMatrix &, double);
    friend int issup(const ARM_GenMatrix &, double);
    friend int issupequal(const ARM_GenMatrix &, double);
    friend int isinf(const ARM_GenMatrix &, double);
    friend int isinfequal(const ARM_GenMatrix &, double);
};



/*----------------------------------------------------------------------------*
    += operator
*----------------------------------------------------------------------------*/
 
inline ARM_GenMatrix& ARM_GenMatrix::operator += (const ARM_GenMatrix & m)
{
    int    i, j;
 
 
 
    if (( m.GetNumLines() != GetNumLines() )
        ||
        ( m.GetNumCols() != GetNumCols() )
       )
    {
        throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
             "Inconsistent matrix sizes in operator +");
 
        return(*this);
    }
 
    for (i = 0; i < GetNumLines(); i++)
    {
        for (j = 0; j < GetNumCols(); j++)
        {
            Elt(i,j) += m.Elt(i,j);
        }
    }
 
    return(*this);
}


/*----------------------------------------------------------------------------*
    -= operator
*----------------------------------------------------------------------------*/
 
inline ARM_GenMatrix& ARM_GenMatrix::operator -= (const ARM_GenMatrix & m)
{
    int    i, j;
 
 
    if (m.GetNumLines() != GetNumLines() || m.GetNumCols() != GetNumCols())
    {
        throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                         "Inconsistent matrix sizes in operator +");
 
        return(*this);
    }
 
    for (i = 0; i < GetNumLines(); i++)
    {
        for (j = 0; j < GetNumCols(); j++)
        {
            Elt(i,j) -= m.Elt(i,j);
        }
    }
 
    return(*this);
}



/*----------------------------------------------------------------------------*
    *= (double) operator
*----------------------------------------------------------------------------*/
 
inline ARM_GenMatrix& ARM_GenMatrix::operator *= (double x)
{
    int    i, j;
 

    for (i = 0; i < GetNumLines(); i++)
    {
        for (j = 0; j < GetNumCols(); j++)
        {
            Elt(i,j) *= x;
        }
    }
 
    return(*this);
}



/*----------------------------------------------------------------------------*
    /= (double) operator
*----------------------------------------------------------------------------*/
 
inline ARM_GenMatrix& ARM_GenMatrix::operator /= (double x)
{
    int    i, j;
 
  
    for (i = 0; i < GetNumLines(); i++)
    {
        for (j = 0; j < GetNumCols(); j++)
        {
            Elt(i,j) /= x;
        }
    }
 
    return(*this);
}


//   += double operator  (and -= double)

inline ARM_GenMatrix& ARM_GenMatrix::operator += (double x)
{
    int    i, j;


    for (i = 0; i < GetNumLines(); i++)
    {
        for (j = 0; j < GetNumCols(); j++)
        {
            Elt(i,j) += x;
        }
    }

   return(*this);
}


inline ARM_GenMatrix& ARM_GenMatrix::operator -= (double x)
{
    (*this) += (-x);

    return(*this);
}



class ARM_Vector : public ARM_GenMatrix 
{
    protected:

        double* itsElt;
        int     itsSize;
    
    public:

        ARM_Vector(void) 
        { 
            SetName(ARM_VECTOR);

            itsSize = 0; 
            itsElt  = (double *) NULL;
        }

        ARM_Vector(int size);
        ARM_Vector(int size, double);
        ARM_Vector(int size, double*);		
		ARM_Vector(double, double, double);
        ARM_Vector(const ARM_Vector &m);
        ARM_Vector(const ARM_Matrix &m);
        ARM_Vector(const ARM_Matrix &m, int nCol);
        ARM_Vector(const ARM_Matrix &m, ARM_LignOrCol type, int n) ;

        explicit ARM_Vector(const std::vector<double>& inVect);
        ARM_Vector(std::vector<std::vector< double > >& inVect);


        explicit ARM_Vector(const ARM_Vector* vect);

        ARM_Vector(const ARM_Vector* vect, int beg, int end);

        ARM_Vector(const ARM_Vector* vect, int size, 
                   int beg, int end, int newStart = 0);


        void Init(void)
        {
            SetName(ARM_VECTOR);
 
            itsSize = 0;
            itsElt  = (double *) NULL;
        }

        double* GetElt(void);
		const double* GetElt() const; 

        void InitElt(int, double);

       ~ARM_Vector() 
        { 
            if (itsElt)
            {
               delete [] itsElt; 
               
               itsElt = NULL;
            }

            itsSize = 0;
        }
            
        void BitwiseCopy(const ARM_Object* srcVect)
        {
            ARM_Vector* src = (ARM_Vector *) srcVect; 

            if (itsElt)
            {
               delete [] itsElt;
 
               itsElt = NULL;

               itsSize = 0;
            }

            if (( src == NULL ) || ( src->itsSize == 0 ))
            {
                return;
            }
            else
            {
                itsSize = src->itsSize;

                itsElt = new double[itsSize];

               /*
                int    i;
 
                for (i = 0; i < itsSize; i++)
                {
                    itsElt[i] = (*src)[i];
                }
               */

               memcpy(itsElt, src->GetElt(), itsSize * sizeof(double));
               
            }
        }
 
        void Copy(const ARM_Object* srcVect)
        {
            ARM_GenMatrix::Copy(srcVect);

            BitwiseCopy(srcVect);
        }

        ARM_Object* Clone(void)
        {
            ARM_Vector* theClone = new ARM_Vector();
 
 
            theClone->Copy(this);
 
            return(theClone);
        }

        void View(char* id, FILE* ficOut)
        {
            FILE* fOut;
            char fOutName[200]; 
            if ( ficOut == NULL )
            {
               ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName); 
               (void) unlink(fOutName);
               fOut = fopen(fOutName, "w");
            }
            else            
               fOut = ficOut;            
            for (int i = 0; i < itsSize; i++)
            {          
               fprintf(fOut,"%14.4lf\t",itsElt[i]);               
               fprintf(fOut,"\n");
            }
        }


        ARM_Vector& operator = (const ARM_Vector &);
        ARM_Vector& operator = (const ARM_Matrix &);
        ARM_Vector& operator += (const ARM_Vector& Data);
		ARM_Vector& operator -=(const ARM_Vector& Data);
        double operator * (const ARM_Vector& Data) const;

        ARM_Vector & operator -();

        double& Elt(int) const;
        double& Elt(int, int) const;

        double& operator[] (int) const;

        void Resize(int);
		void Resize(int,double v); 
		void clear() { Resize(0); }
		/// very inefficient I repeat VERY INEFFICIENT as 
		/// ARM_Vector is not like vector<T> a dynamically growing
		/// structure
        void insert(double anElt);
		
		/// very inefficient I repeat VERY INEFFICIENT as 
		/// ARM_Vector is not like vector<T> a dynamically growing
		/// structure (STL interface)
		void push_back(double anElt);

		void pop_back(void);

        int find(double anElt,double tol=K_NEW_DOUBLE_TOL)  const ;
        void replace(size_t i,double anElt);
   
        ARM_Vector Sort(int incOrdec = K_INCREASING) const ;
		ARM_Vector Sort(ARM_Vector& vPermut, const int incOrdec) const ;
		ARM_Vector sort() const { return Sort(); }

        ARM_Vector* Sort_Compact(void);        

        ARM_Vector* Compact(void);

        inline int GetSize(void) const 
        { 
            return(itsSize); 
        }
		double last() const 
		{
			if (itsSize==0) 
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"last(): empty ARM_Vector ");
			return itsElt[itsSize-1] ;
		}

		/// for STL like interface!
		inline int size() const		{            return itsSize; 		}
		bool empty() const { return itsSize==0; }

        int GetNumLines(void) const { return(itsSize); }
        int GetNumCols(void)  const { return(1); }
        double Mean(void) const ;
       
        inline double sum() const{ double summ=0.0;for( int i =0; i<size(); ++i )	summ+=itsElt[i]; return summ;}

        inline double StdNorm(void) const 
        {
             int    i;
            double norm = 0.0;

            for (i = 0; i < itsSize; i++) 
            {
                norm += itsElt[i]*itsElt[i];
            }     

             return(sqrt(norm));        
        }


        // returns the index idx such that x in [Elt(idx), Elt(idx+1) [
        int LookupPrevIndex(double x) const ;          

        // returns idx such that |Elt(idx) - x| < tol, -1 otherwise
        // (default tol < 1/365.0)
        int LookupExactIndex(double x, double tol = 2e-3) const ;


        ARM_Vector & operator += (double); 
        ARM_Vector & operator -= (double);
        ARM_Vector & operator *= (double);
        ARM_Vector & operator /= (double);

        // pour reduire le nb de destructions
        friend void ResizeOrCreate(ARM_Vector*&, int, int = 1);   
		
		//
		friend ARM_Vector operator + (const ARM_Vector &, const ARM_Vector &);
		friend ARM_Vector operator * (const double &, const ARM_Vector &);


		/// support for mathematical function
        friend ARM_Vector log(const ARM_Vector& arg1);
        friend ARM_Vector exp(const ARM_Vector& arg1);
        friend ARM_Vector sqrt(const ARM_Vector& arg1);
        friend ARM_Vector pow(const ARM_Vector& arg1, const ARM_Vector& arg2);
        friend ARM_Vector abs(const ARM_Vector& ag1);
        friend ARM_Vector floor(const ARM_Vector& arg1);
		friend ARM_Vector VectorIf( const ARM_Vector& cond, const ARM_Vector & argTrue, const ARM_Vector & argFalse );
		friend ARM_Vector VectorMin( const ARM_Vector& arg1, const ARM_Vector& arg2 );
		friend ARM_Vector VectorMax( const ARM_Vector& arg1, const ARM_Vector& arg2);
        friend ARM_Vector Suite(double, double, double = 1.0);
        friend double Scalar(const ARM_Vector &, const ARM_Vector &);
        int GetFirstLast(double& firstValue,double& lastValue) const;

		/// typedef for iterator support
		typedef double			value_type;
		typedef double*			iterator;
		typedef const double*	const_iterator;
		typedef double&			reference;
		typedef const double&	const_reference;
		typedef size_t			size_type;
		typedef ptrdiff_t		difference_type;

		/// iterator support
		iterator begin() { return &itsElt[0]; }
		const_iterator begin() const { return &itsElt[0]; }
		iterator end() { return &itsElt[0]+itsSize; }
		const_iterator end() const { return &itsElt[0]+itsSize; }

		// transition to STL 
		ARM_Vector& operator=(const std::vector<double>&arg); 
		void populate(std::vector<double>&dest) const ; 
};

///////////////////////////////////////////////
/// function mapping for ARM_Vector
///////////////////////////////////////////////

ARM_Vector VectorExp( const ARM_Vector& v);
ARM_Vector VectorLog( const ARM_Vector& v);
ARM_Vector VectorSqrt( const ARM_Vector& v);
ARM_Vector VectorPow( const ARM_Vector& v1, const ARM_Vector& v2);



/*----------------------------------------------------------------------------*
    Initializes a vector with size elements
*----------------------------------------------------------------------------*/
 
inline ARM_Vector::ARM_Vector(int size)
{
    SetName(ARM_VECTOR);
 
    itsSize = size;
 
    if ( itsSize == 0 )
    {
        itsElt = NULL;

        return;
    }
    else
    {
        itsElt = new double[itsSize];

        MEMSET(itsElt, 0, sizeof(double)*itsSize);
    }
}



/*----------------------------------------------------------------------------*
    Initializes a vector with size elements with x0 value
*----------------------------------------------------------------------------*/
 
inline ARM_Vector::ARM_Vector(int size, double x0)
{
    int i;
 
 
 
    SetName(ARM_VECTOR);
 
    itsSize = size;
 
    if ( itsSize == 0 )
    {
        itsElt = NULL;
 
        return;
    }
    else
    {
        itsElt = new double[itsSize];
 
        for (i = 0; i < size; i++)
        {
            itsElt[i] = x0;
        }
    }
}



/*----------------------------------------------------------------------------*
    Initializes a vector with size elements with x0 value
*----------------------------------------------------------------------------*/
 
inline ARM_Vector::ARM_Vector(int size, double* x0)
{
    int i;
 
 
 
    SetName(ARM_VECTOR);
 
    itsSize = size;
 
 
    if (( x0 == NULL ) || ( itsSize == 0 ))
    {
        itsElt = NULL;
 
        return;
    }
    else
    { 
        itsElt = new double[itsSize];
 
        for (i = 0; i < itsSize; i++)
        {
            itsElt[i] = x0[i];
        }
    }
}

inline ARM_Vector::ARM_Vector(const ARM_Vector* vect)
{
    Init();
 
    SetName(ARM_VECTOR);
 
    if (vect)
        this->Copy(vect);
}



/*----------------------------------------------------------------------------*
    Initializes a vector (copy) of v.
*----------------------------------------------------------------------------*/
 
inline ARM_Vector::ARM_Vector(const ARM_Vector& v) : ARM_GenMatrix(v)
{
    Init();
 
    SetName(ARM_VECTOR);
 
    this->BitwiseCopy(&v);
}



/*----------------------------------------------------------------------------*
    Export itsElt into a tableau of double
*----------------------------------------------------------------------------*/
 
inline double* ARM_Vector::GetElt(void)
{
    return(itsElt);
}

inline const double* ARM_Vector::GetElt(void) const
{
    return(itsElt);
}


/*----------------------------------------------------------------------------*
    Initializes a vector by copy of start-th(included)
    to end-th(included) elements of v.
*----------------------------------------------------------------------------*/
 
inline ARM_Vector::ARM_Vector(const ARM_Vector *v, int start, int end)
{
    int    i;
 
    SetName(ARM_VECTOR);
 
    if (v->GetSize()<end || start>end)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "Invalid argument in ARM_Vector(ARM_Vector*) constructor");
    }
 
    itsSize = end-start+1;
 
    itsElt = new double[itsSize];
 
    for (i = 0; i < itsSize; i++)
    {
        itsElt[i] = v->Elt(start+i);
    }
}

/*----------------------------------------------------------------------------*
    Initializes a vector of size size and where elements  
    from newStart are obtained by copy of start-th(included)
    to end-th(included) elements of v.
*----------------------------------------------------------------------------*/
 
inline ARM_Vector::ARM_Vector(const ARM_Vector* v, int size, int start, int end, 
                              int newStart)
{
    int    i;
 
    SetName(ARM_VECTOR);
 
    if (v->GetSize()<end || start>end || size < end-start+newStart+1)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "Invalid argument in ARM_Vector(ARM_Vector*) constructor");
    }
 
    itsSize = size;
 
    itsElt = new double[size];
 
    for (i = 0; i < end-start+1; i++)
    {
        itsElt[i+newStart] = v->Elt(start+i);
    }
}





/*----------------------------------------------------------------------------*
    Resize the ARM_Vector to size. Elements are reset to
    zero. Returns 0 if failed.
*----------------------------------------------------------------------------*/
 
inline void ARM_Vector::Resize(int size)
{
    if (itsElt)
    {
       delete [] itsElt;

       itsElt = NULL;
    }
 
    itsSize = size;

    if ( itsSize == 0 )
       return;
 
    itsElt = new double[itsSize];

    memset(itsElt, '\0', size * sizeof(double));
}

inline void ARM_Vector::Resize(int size,double v)
{
    if (itsElt)delete [] itsElt;
    itsElt = NULL;
    itsSize = size;
	if ( itsSize == 0 )return;
    itsElt = new double[itsSize];
	for(int i=0;i<itsSize;i++) itsElt[i]=v; 
}



inline void ARM_Vector::insert(double Elt)
{
    double* newElt = new double[itsSize+1];

    for (int i = 0; i < itsSize; i++)
    {
        newElt[i] = itsElt[i];
    }

    newElt[itsSize] = Elt;

    itsSize++;

    delete [] itsElt;

    itsElt = newElt;
}


inline void ARM_Vector::pop_back(void)
{
    if (itsSize > 0)
	{
		double* newElt = new double[itsSize-1];

		for (int i = 0; i < itsSize-1; i++)
		{
			newElt[i] = itsElt[i];
		}

		itsSize--;

		delete [] itsElt;

		itsElt = newElt;
	}
	else
	{
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "no element in vector");
	}
}



/*----------------------------------------------------------------------------*
    SYNOPSIS    double & ARM_Vector::Elt(int i)
 
    Access i-th element of vector
*----------------------------------------------------------------------------*/
 
inline double& ARM_Vector::Elt(int i) const
{
    if (i<0 || i>itsSize-1)
    {
		char msg[255];
		sprintf( msg, "Index %d  not in the range 0 %d in Elt(int) method", i, itsSize-1 );
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, msg );
    }
    return itsElt[i];
}


/*----------------------------------------------------------------------------*
    Access i-th element of vector (for viewing a vector as a matrix)
*----------------------------------------------------------------------------*/
 
inline double& ARM_Vector::Elt(int i, int j) const
{
    if (i<0 || i>itsSize-1 )
    {
		char msg[255];
		sprintf( msg, "Index %d  not in the range 0 %d in Elt(int,int) method", i, itsSize-1 );
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, msg );
    }

	if( j != 0 )
	{
		char msg[255];
		sprintf( msg, "Trying to access column %d  of a vector in Elt(int,int) method, Please advise", j );
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, msg );
	}
    return itsElt[i];
}


/*----------------------------------------------------------------------------*
    SYNOPSIS    double & ARM_Vector::operator[] (int i)
 
    Access i-th element of vector. Same than Elt method
*----------------------------------------------------------------------------*/
 
inline double& ARM_Vector::operator[] (int i) const
{
    if (i<0 || i>itsSize-1)
    {
		char msg[255];
		sprintf( msg, "Index %d  not in the range 0 %d in operator[](int) method", i, itsSize-1 );
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, msg );
    }
    return itsElt[i];
}


inline int ARM_Vector::GetFirstLast(double& firstValue,double& lastValue) const
{
    if(itsSize>0)
    {
        firstValue=itsElt[0];
        lastValue=itsElt[itsSize-1];
    }
    return itsSize;
}


/*----------------------------------------------------------------------------*
    SYNOPSIS

    ARM_Vector log(ARM_Vector &v)
    ARM_Vector sqrt(ARM_Vector &v)
    ARM_Vector exp(ARM_Vector &v)
    ARM_Vector abs(ARM_Vector &v)
    ARM_Vector floor(ARM_Vector &v)

    operations element by element on a vector
        
*----------------------------------------------------------------------------*/

inline ARM_Vector log(const ARM_Vector &v)
{
    int i;
    ARM_Vector tmp(v.GetNumLines());


    if (!(issup(v, 0.0))) 
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
          "Some vector elements are non-positive in log(ARM_Vector) method");
        return(tmp);
    }

    for (i=0; i<tmp.GetNumLines(); i++) 
    {
        tmp.Elt(i) = log(v.Elt(i));
    }

    return(tmp);
}


inline ARM_Vector sqrt(const ARM_Vector &v)
{
    int i;
    ARM_Vector tmp(v.GetNumLines());

    if (!(issupequal(v, 0.0))) 
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
         "Some vector elements are non-positive in sqrt(ARM_Vector) method");

        return(tmp);
    }

    for (i=0; i<tmp.GetNumLines(); i++) 
    {
        tmp.Elt(i) = sqrt(v.Elt(i));
    }

    return(tmp);
}



inline ARM_Vector exp(const ARM_Vector &v)
{
    int i;
    ARM_Vector tmp(v.GetNumLines());


    for (i=0; i<tmp.GetNumLines(); i++) 
    {
        tmp.Elt(i) = exp(v.Elt(i));
    }

    return(tmp);
}



inline ARM_Vector abs(const ARM_Vector& v)
{
    int i;
    ARM_Vector tmp(v.GetNumLines());


    for (i=0; i<tmp.GetNumLines(); i++) 
    {
        tmp.Elt(i) = fabs(v.Elt(i));
    }

    return(tmp);
}



inline ARM_Vector floor(const ARM_Vector& v)
{
    int i;
    ARM_Vector tmp(v.GetNumLines());

    for (i=0; i<tmp.GetNumLines(); i++) 
    {
        tmp.Elt(i) = floor(v.Elt(i));
    }

    return(tmp);
}




//// a handmade rounding function
double roundsimple(double num, int places);

inline ARM_Vector pow(const ARM_Vector& bases, const ARM_Vector& exponents)
{
	if( bases.size() != exponents.size() )
        throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
             "Base and exponent matrices are not the same size... please advise!" );

	ARM_Vector results( bases.size() );

	ARM_Vector::const_iterator
			base 	= bases.begin(),
			baseEnd	= bases.end(),
			exp  	= exponents.begin();

	ARM_Vector::iterator
			result = results.begin();

	for( ; base != baseEnd; ++base, ++result, ++exp )
		*result = pow( *base, *exp );

	return results;
}




inline ARM_Vector vif( const ARM_Vector& cond, const ARM_Vector& argTrue, const ARM_Vector& argFalse )
{
	int size = argTrue.size() < argFalse.size()? argTrue.size(): argFalse.size();

	if( cond.size() < size )
		size = cond.size();

	ARM_Vector result( size );

	ARM_Vector::const_iterator
			elemCond = cond.begin(),
			elemA 	 = argTrue.begin(),
			elemB    = argFalse.begin();

	ARM_Vector ::iterator
			elemResult = result.begin(),
			endResult  = result.end();

	for( ; elemResult != endResult; ++elemResult, ++elemCond, ++elemA, ++elemB)
		*elemResult == *elemCond? *elemA : *elemB; 

	return result;
}


inline ARM_Vector vmin( const ARM_Vector& A, const ARM_Vector& B )
{
	int size =A.size() < B.size()? A.size(): B.size();
	ARM_Vector result( size );

	ARM_Vector::const_iterator
			elemA = A.begin(),
			elemB = B.begin();

	ARM_Vector::iterator
			elemResult = result.begin(),
			endResult  = result.end();

	for( ; elemResult != endResult; ++elemResult, ++elemA, ++elemB )
		*elemResult = *elemA < *elemB? *elemA: *elemB;

	return result;
}



inline ARM_Vector vmax( const ARM_Vector& A, const ARM_Vector& B)
{
	int size =A.size() < B.size()? A.size(): B.size();
	ARM_Vector result( size );

	ARM_Vector::const_iterator
			elemA = A.begin(),
			elemB = B.begin();

	ARM_Vector ::iterator
			elemResult = result.begin(),
			endResult  = result.end();

	for( ; elemResult != endResult; ++elemResult, ++elemA, ++elemB )
		*elemResult = *elemA < *elemB? *elemB: *elemA;

	return result;
}





/////////////////////////////////////////////////////////

class ARM_Matrix: public ARM_GenMatrix 
{
    protected:

        double* itsElt;
        int     itsNumLin;
        int     itsNumCol;
    
    public:
       
        
        ARM_Matrix(void) 
        { 
            itsNumLin = itsNumCol = 0; 

            itsElt = (double *) NULL;

            SetName(ARM_MATRIX);
        }

        ARM_Matrix(int, int);
        ARM_Matrix(int , int , double*);
        ARM_Matrix(int , int , double**);
        ARM_Matrix(double**, int , int);
        ARM_Matrix(int , int , double);
        ARM_Matrix(ARM_Matrix* mat);
        ARM_Matrix(ARM_Matrix* mat, int fstLine, int lastLine);
        ARM_Matrix(int , int , std::vector<std::vector< double > >& inVect);
        ARM_Matrix(const ARM_Matrix &);
        ARM_Matrix(const ARM_TDiag &);
        ARM_Matrix(const ARM_Vector &);

       ~ARM_Matrix() 
        { 
            if (itsElt)
               delete [] itsElt;

            itsNumLin = itsNumCol = 0;
        }

        void Init(void)
        {
            itsNumLin = itsNumCol = 0;
 
            itsElt = (double *) NULL;
 
            SetName(ARM_MATRIX);
        }

        void BitwiseCopy(const ARM_Object* srcMat)
        {
            ARM_Matrix* src = (ARM_Matrix *) srcMat;
            itsNumLin = src->itsNumLin;
            itsNumCol = src->itsNumCol;

            if (itsElt)
            {
                delete [] itsElt;
                itsElt = NULL;
            }

            itsElt = new double [itsNumLin*itsNumCol];

            memcpy(itsElt, src->GetElt(),itsNumLin*itsNumCol*sizeof(double));
        }
 
 
        void Copy(const ARM_Object* srcMat)
        {
            ARM_GenMatrix::Copy(srcMat);
            BitwiseCopy(srcMat);
        }

        ARM_Object* Clone(void)
        {
            ARM_Matrix* theClone = new ARM_Matrix(); 
            theClone->Copy(this); 
            return(theClone);
        }

        ARM_Matrix& operator = (const ARM_Matrix &);
        ARM_Matrix& operator = (const ARM_TDiag &);
        ARM_Matrix& operator = (const ARM_Vector &);

        ARM_Matrix& operator -();
        void    Resize(int, int);
        int     GetNumLines(void) const { return(itsNumLin); }
        int     GetNumCols(void) const {return(itsNumCol);}
        double& Elt(int, int) const;

        double& RawElt(int, int) const; //pas de test, pour accelerer!!!

        ARM_Matrix Sort(int col);
        virtual void   Transpose(void);
        virtual void   Invert(ARM_Matrix *, double &);
        virtual void   LinSolve(ARM_Vector *, double &);
        virtual double Det(void);

        void   QRSolve(ARM_Vector* , double & );

    
        inline double* GetDbleLine(int nLine)
        {
            return (itsElt + nLine*itsNumCol);
        }


        double* GetElt(void);
    	ARM_pLine	GetLine(int n);
	    ARM_pCol	GetCol (int n); 
    	ARM_pUpDiag	GetUpDiag(int n); 

        ARM_Vector* GetRow(int n);
        ARM_Vector* GetColumn (int n);
        ARM_Matrix* GetMatrix(int f_row, int l_row,
                              int f_col, int l_col);


        // pour reduire le nb de destructions        
        friend void ResizeOrCreate(ARM_Matrix*&, int, int, int = 1);

        // Swap two lines

        void SwapLines(int k, int h)
        {
            double tmpVal;

            for (int j = 0; j < itsNumCol; j++)
            {
                tmpVal = itsElt[k*itsNumCol+j];
                itsElt[k*itsNumCol+j] = itsElt[h*itsNumCol+j];
                itsElt[h*itsNumCol+j] = tmpVal; 
            }
        }

        void SortLines(int criteriaCol)
        {
            for (int i = 0; i < itsNumLin; i++)
            {
                double minValLin = itsElt[i*itsNumCol+criteriaCol];

                for (int j = i+1; j < itsNumLin; j++)
                {
                    double curVal = itsElt[j*itsNumCol+criteriaCol];

                    if ( curVal < minValLin )
                    {
                       SwapLines(i, j);

                       minValLin = curVal;
                    }
                }
            }
        }

        void Ludcmp(ARM_Vector *, double &);
        void Lubksb(ARM_Vector *, ARM_Vector *);  

        static ARM_Matrix* Jacobi(ARM_Matrix*a,ARM_Vector& d,int& nrot);

        ARM_Matrix* eigsrt(ARM_Vector& d);
        static ARM_Matrix* ACP(ARM_Matrix*a,ARM_Vector& d);
};


ARM_Matrix* MergeMatrix(ARM_Matrix* mat1, ARM_Matrix* mat2);


/*----------------------------------------------------------------------------*
    Product operator for matrixes
 
    WARNING: result is temporary stored in static ARM_Matrix tmp.
    Do not use more than once in a single expression as new product
    result would override old product result before it is assigned
    to another ARM_Matrix.
*----------------------------------------------------------------------------*/
 
inline ARM_Matrix& operator * (ARM_GenMatrix &m1, ARM_GenMatrix &m2)
{
    static ARM_Matrix tmp;
    int    i, j, k;
    double s;
 
 
 
    if ( m1.GetNumCols() != m2.GetNumLines() )
    {
        throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                         "Inconsistent matrix sizes in operator *");
        return(tmp);
    }
 
    tmp.Resize(m1.GetNumLines(), m2.GetNumCols());
 
    //    compute product
 
    for (i = 0; i < m1.GetNumLines(); i++)
    {
        for (j = 0; j < m2.GetNumCols(); j++)
        {
            s = 0.0;
 
            for (k=0; k<m1.GetNumCols(); k++)
            {
                s += m1.Elt(i,k) * m2.Elt(k,j);
            }
 
            tmp.Elt(i, j) = s;
        }
    }
 
    return(tmp);
}


/*----------------------------------------------------------------------------*
    Returns m^n.
*----------------------------------------------------------------------------*/
 
inline ARM_Matrix& operator ^ (ARM_GenMatrix& mat, int n)
{
    static    ARM_Matrix    np(mat.GetNumLines(), mat.GetNumLines());
    static    ARM_Matrix    p(mat.GetNumLines(), mat.GetNumLines());
    int       i, j;
 
 
 
    if (!mat.IsSquare())
    {
        throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                         "Matrix is not square in operator ^");
          return(np);
    }
 
    for (i=0; i<mat.GetNumLines(); i++)
    {
         for (j=0; j<mat.GetNumCols(); j++)
         {
             p.Elt(i,j) = mat.Elt(i,j);
         }
    }
 
    np.SetToId();
 
    while (n)
    {
        if (n % 2)
        {
           n--;
 
           np = np * p;
        }
        else
        {
            n /= 2;
 
            p = p * p;
        }
    }
 
    return(np);
}


/*----------------------------------------------------------------------------*/
/*    CholeskyDcmp(ARM_Matrix M)                                             */
/*    Performs Cholesky decomposition. This is from Numerical Recipes in C.   */
/*----------------------------------------------------------------------------*/

inline ARM_Matrix CholeskyDcmp(ARM_Matrix& M)
{
    //    check input

    if (!M.IsSquare() || !M.IsSymmetric())
    {
       throw MathException(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB, 
                        "Matrix must be symmetric");
    }
            
    int    i, j, k;

    int numLines = M.GetNumLines();

    ARM_Matrix tmp(numLines, numLines, 0.0);

    double sum;

    for (i = 0; i<numLines; i++)
    {
        for (j = i; j<numLines; j++) 
        {
            for (sum = M.Elt(i, j), k=i-1; k>=0; k--) 
            {
                // MA Le 01/07/2002 : (voir M.Picot sum -= M.Elt(i, k)*M.Elt(j, k);

                 sum -= tmp.Elt(i, k)*tmp.Elt(j, k);
            }
            
            if (i == j)
            {
                if (sum<=0.0)
                {
                  throw MathException(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                     "Cholesky Decomp failed: Matrix is not positive definite");
                }
                tmp.Elt(i, i) = sqrt(sum);
            }
            else
            {
                tmp.Elt(j, i) = sum/tmp.Elt(i, i);
            }
        }
    }
    
    return(tmp);
}
/*-------------------------------------------------------------------------*
      Jakobi decomposation
*--------------------------------------------------------------------------*/
   
inline void Rot(ARM_Matrix* a,const double s,
                const double tau,const int i, 
                const int j,const int k,const int l)
{
    double g,h;
    
    g = a->Elt(i,j);
    h = a->Elt(k,l);

    a->Elt(i,j) = g - s*(h+g*tau);
    a->Elt(k,l) = h + s*(g-h*tau);

}

//_______________________________________________________________________________________________________________________

/*----------------------------------------------------------------------------*
    returns the log of the matrix
*----------------------------------------------------------------------------*/
 
inline ARM_Matrix log(ARM_GenMatrix &m)
{
    ARM_Matrix tmp(m.GetNumLines(), m.GetNumCols());
    int        i, j;
 
 
 
 
    if (!(issup(m,0.0)))
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
             "Some matrix elements are non-positive in log method");
 
        return(tmp);
    }
 
    for (i=0; i<tmp.GetNumLines(); i++)
    {
        for (j=0; j<tmp.GetNumCols(); j++)
        {
            tmp.Elt(i,j) = log(m.Elt(i,j));
        }
    }
 
    return(tmp);
}


/*----------------------------------------------------------------------------*
    returns the sqrt of the matrix
*----------------------------------------------------------------------------*/
 
inline ARM_Matrix sqrt(ARM_GenMatrix &m)
{
    ARM_Matrix tmp(m.GetNumLines(), m.GetNumCols());
    int        i, j;
 
 
 
    if (!(issupequal(m, 0.0)))
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
             "Some matrix elements are negative in sqrt method");
 
        return(tmp);
    }
 
    for (i=0; i<tmp.GetNumLines(); i++)
    {
        for(j=0; j<tmp.GetNumCols(); j++)
        {
            tmp.Elt(i,j) =  sqrt(m.Elt(i,j));
        }
    }
 
    return(tmp);
}



/*----------------------------------------------------------------------------*
    SYNOPSIS    ARM_Matrix exp(ARM_GenMatrix &m)
 
    returns the exp of the matrix
*----------------------------------------------------------------------------*/
 
inline ARM_Matrix exp(ARM_GenMatrix &m)
{
    ARM_Matrix tmp(m.GetNumLines(), m.GetNumCols());
    int        i, j;
 
 
    for (i=0; i<tmp.GetNumLines(); i++)
    {
        for(j=0; j<tmp.GetNumCols(); j++)
        {
            tmp.Elt(i,j) =  exp(m.Elt(i,j));
        }
    }
 
    return(tmp);
}



/*----------------------------------------------------------------------------*
    returns the fabs of the matrix
*----------------------------------------------------------------------------*/
 
inline ARM_Matrix fabs(ARM_GenMatrix &m)
{
    ARM_Matrix tmp(m.GetNumLines(), m.GetNumCols());
    int        i, j;
 
 
    for (i=0; i<tmp.GetNumLines(); i++)
    {
        for(j=0; j<tmp.GetNumCols(); j++)
        {
            tmp.Elt(i,j) =  fabs(m.Elt(i,j));
        }
    }
 
    return(tmp);
}



/*----------------------------------------------------------------------------*
    returns the abs of the matrix
*----------------------------------------------------------------------------*/
 
inline ARM_Matrix floor(ARM_GenMatrix &m)
{
    ARM_Matrix tmp(m.GetNumLines(), m.GetNumCols());
    int        i, j;
 
 
 
    for (i=0; i<tmp.GetNumLines(); i++)
    {
        for(j=0; j<tmp.GetNumCols(); j++)
        {
            tmp.Elt(i,j) =  floor(m.Elt(i,j));
        }
    }
 
    return(tmp);
}


inline ARM_Matrix expMat(ARM_Matrix A)
//developpement en serie limitee
{
    if (A.GetNumLines()!=A.GetNumLines()) 
    {
        throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
         "the matrix must be square");
        return(A);
    }
	ARM_Matrix loc1, loc2;
	int i,n;
    n = A.GetNumLines();
	loc1 = ARM_Matrix(n,n,0.0);
    for(i = 0 ; i < n ; i++)
    {
        loc1.Elt(i,i) = 1.0;
    }

    loc2 = A;

    for(i = 0 ; i < 15 ; i++)
    {
        loc1 += loc2;
        loc2 = loc2 * A;
        loc2 *= 1.0/(2.0+i);
    }

    return(loc1);
}


inline ARM_Matrix expMatTrig(ARM_Matrix A)
//exponentielle d'une matrice triangulaire.
{
    if (A.GetNumLines()!=2  || A.GetNumLines()!=2) 
    {
        throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
         "the matrix must be 2x2");
        return(A);
    }

    if (A.Elt(1,0)*A.Elt(1,0) > 0.001) 
    {
        throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
         "the matrix must be trigonal");
        return(A);
    }

    ARM_Matrix loc1;

    loc1 = ARM_Matrix(2,2,0.0);
    loc1.Elt(0,0) = exp(A.Elt(0,0));
    loc1.Elt(1,1) = exp(A.Elt(1,1));
    if ((A.Elt(1,1)-A.Elt(0,0))*(A.Elt(1,1)-A.Elt(0,0)) < 0.000001) 
    {
        loc1.Elt(0,1) = A.Elt(0,1)*loc1.Elt(0,0);
    }
    else
    {
        loc1.Elt(0,1) = A.Elt(0,1)*(loc1.Elt(1,1)-loc1.Elt(0,0))
                        /(A.Elt(1,1)-A.Elt(0,0));
    }

    return(loc1);
}



inline void trig2(ARM_Matrix A, ARM_Matrix& P, ARM_Matrix& T, ARM_Matrix& Pm)
{
    double rien;

    ARM_Matrix locm1;

    P = ARM_Matrix(2, 2, 0.0);
   
    Pm = ARM_Matrix(2, 2, 0.0);

    T = ARM_Matrix(2, 2, 0.0);
    //force de rappel A=Pm.T.P ; Pm = inverse de P.
    //T triangulaire superieure de valeurs propres non-nulles.
    if(A.Elt(1,0)*A.Elt(1,0)>0.00000001)
    {
        rien = (A.Elt(0,0)+A.Elt(1,1))*(A.Elt(0,0)+A.Elt(1,1))
               -4.0*(A.Elt(0,0)*A.Elt(1,1)-A.Elt(0,1)*A.Elt(1,0));
        //tester rien >0
        rien = sqrt(rien);
        T.Elt(0,0)=0.5*((A.Elt(0,0)+A.Elt(1,1))+rien);
        T.Elt(1,1)=0.5*((A.Elt(0,0)+A.Elt(1,1))-rien);

        Pm.Elt(1,0)=1.0/sqrt(1+pow((T.Elt(0,0)-A.Elt(1,1))/A.Elt(1,0),2.0));
        Pm.Elt(0,0)=(T.Elt(0,0)-A.Elt(1,1))/A.Elt(1,0)*Pm.Elt(1,0);
        Pm.Elt(1,1)=1.0/sqrt(1+pow((T.Elt(1,1)-A.Elt(1,1))/A.Elt(1,0),2.0));
        Pm.Elt(0,1)=(T.Elt(1,1)-A.Elt(1,1))/A.Elt(1,0)*Pm.Elt(1,1);

        Pm.Invert(&P,rien);

        locm1 = A * Pm;
        T = P * locm1;
    }
    else
    {
        T = A;
        T.Elt(1,0) = 0.0;
        P.Elt(0,0) = 1.0;
        P.Elt(1,1) = 1.0;
        Pm = P;
    }
}


/*----------------------------------------------------------------------------*
    Initialize an nl * nc matrix
*----------------------------------------------------------------------------*/
 
inline ARM_Matrix::ARM_Matrix(int nl, int nc)
{
    SetName(ARM_MATRIX);
 
    itsElt = (double *) NULL;
 
 
    try
    {
        Resize(nl, nc);
    }
 
    catch(Exception & x)
    {
        x.DebugPrint();
 
        throw Exception(__LINE__, __FILE__, ERR_MEMORY_ALLOCATION,
            "Memory Allocation Pb in ARM_Matrix(int nl, int nc) constructor");
    }
}



/*----------------------------------------------------------------------------*
    Initializes a matrix with size elements with x00 value
*----------------------------------------------------------------------------*/
 
inline ARM_Matrix::ARM_Matrix(int szlg, int szcol, double* x00)
{
//    int i, j;

    SetName(ARM_MATRIX);

    itsNumLin = szlg;

    itsNumCol = szcol;

    if ( x00 == NULL )
       return;

    itsElt = new double[itsNumLin*itsNumCol];

    memcpy(itsElt, x00, itsNumLin * itsNumCol * sizeof(double));
}



inline ARM_Matrix::ARM_Matrix(ARM_Matrix* mat, int firstLine, int lastLine)
{
    int i, j;
 
    if (( firstLine > lastLine ) || ( lastLine > mat->GetNumLines() ))
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
               "1st line > last line or input matrix dont have enough lines");
    }

    SetName(ARM_MATRIX);
 
    itsNumLin = lastLine - firstLine + 1;
 
    itsNumCol = mat->GetNumCols();
 
    itsElt = new double[itsNumLin*itsNumCol];
 
    for (i=0; i < itsNumLin; i++)
    {
        for (j=0; j < itsNumCol; j++)
        {
             itsElt[i * itsNumCol + j] = mat->Elt(i+firstLine, j);
        }
    }
}




/*----------------------------------------------------------------------------*
    Initializes a matrix with size elements with x00 value
*----------------------------------------------------------------------------*/
 
inline ARM_Matrix::ARM_Matrix(int szlg, int szcol, double x0)
{
    int i, j;

    SetName(ARM_MATRIX);

    itsNumLin = szlg;

    itsNumCol = szcol;

    itsElt = new double[itsNumLin*itsNumCol];

    for (i=0; i < itsNumLin; i++)
    {
        for (j=0; j < itsNumCol; j++)
        {

             itsElt[i * itsNumCol + j] = x0;
        }
    }
}



/*----------------------------------------------------------------------------*
   Built a matrix from an other matrix
*----------------------------------------------------------------------------*/
inline ARM_Matrix::ARM_Matrix(ARM_Matrix* mat)
{
    Init();
 
    SetName(ARM_MATRIX);
 
    if (mat)
       this->Copy(mat);
}
 
 
 
/*----------------------------------------------------------------------------*
    Initialize matrix by copying m.
*----------------------------------------------------------------------------*/
 
inline ARM_Matrix::ARM_Matrix(const ARM_Matrix& m) : ARM_GenMatrix(m)
{
    Init();
 
    SetName(ARM_MATRIX);
 
    this->BitwiseCopy(&m);
}



/*----------------------------------------------------------------------------*
    operator = (ARM_Matrix &m)
*----------------------------------------------------------------------------*/
 
inline ARM_Matrix& ARM_Matrix::operator = (const ARM_Matrix& m)
{
 
    (*this).ARM_GenMatrix::operator = (m);
 
    this->BitwiseCopy(&m);
 
    return(*this);
}



/*----------------------------------------------------------------------------*
    Resize matrix to nl lines and nc columns. Elements are
    reset to zero. Returns 0 if failed.
*----------------------------------------------------------------------------*/
 
inline void ARM_Matrix::Resize(int nl, int nc)
{
    if (( nl < 0 ) || ( nc < 0 ))
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                 "Column or line size is negative in Resize method");
    }
 
    if (itsElt)
    {
       delete [] itsElt;
    }
 
    itsNumLin = nl;
    itsNumCol = nc;
 
 
    itsElt = new double [nl*nc];
 
    memset(itsElt, '\0', nl * nc * sizeof(double));
}



/*----------------------------------------------------------------------------*
    Returns reference to element (i,j)
*----------------------------------------------------------------------------*/
 
inline double& ARM_Matrix::Elt(int i, int j) const
{
    if ( i < 0 || j < 0 || i>itsNumLin-1 || j>itsNumCol-1)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                         "Negative column or line index in Elt method");
    }
    return itsElt[i*itsNumCol+j];
}


/*----------------------------------------------------------------------------*
    Returns reference to element (i,j) without Check, be careful!!
*----------------------------------------------------------------------------*/
 
inline double& ARM_Matrix::RawElt(int i, int j) const
{
    return(itsElt[i*itsNumCol+j]);
}

/*----------------------------------------------------------------------------*
    Export itsElt into a tableau of double
*----------------------------------------------------------------------------*/
 
inline double* ARM_Matrix::GetElt(void)
{
    return(itsElt);
}


/*----------------------------------------------------------------------------*
    Returns determinant of matrix, computed from LU
    decomposition
*----------------------------------------------------------------------------*/
 
inline double ARM_Matrix::Det(void)
{
    double    d;
    int       j;
    static    ARM_Matrix tmp;
    static    ARM_Vector indx(itsNumLin);
 
 
 
    //    check square matrix
 
    if (!IsSquare())
    {
        throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                        "Matrix is not square in operator ^");
        return(K_HUGE_DOUBLE);
    }
 
    //    copy this into tmp
    tmp = *this;
 
    //    Perform LU decomposition and compute determinant
 
    try
    {
        tmp.Ludcmp(&indx, d);
    }
 
    catch(MathException& m)
    {
        m.DebugPrint();
 
        throw MathException(__LINE__, __FILE__, ERR_CALC_DET_PB,
                            "Problem in Ludcmp() to compute determinant");
    }
 
    for (j=0; j<itsNumLin; j++)
    {
        d *= tmp.Elt(j,j);
    }
 
    return(d);
}


/*----------------------------------------------------------------------------*
    returns invert of matrix computed from LU decomposition
*----------------------------------------------------------------------------*/
 
inline void ARM_Matrix::Invert(ARM_Matrix *invert, double &det)
{
    int    i, j;
    ARM_Matrix tmp;
    ARM_Vector indx(itsNumLin);
    ARM_Vector col(itsNumLin);
 
 
 
    //    check input
 
    if (!IsSquare()
        || invert->itsNumLin != itsNumLin || invert->itsNumCol != itsNumCol)
    {
        throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                         "Matrix is not square in operator ^");
    }
 
 
 
    //    copy this into tmp
 
    tmp = *this;
 
    //    perform LU decomposition
 
    try
    {
       tmp.Ludcmp(&indx, det);
    }
 
    catch(MathException& m)
    {
        m.DebugPrint();
 
        throw MathException(__LINE__, __FILE__, ERR_MATRIX_INVER_PB,
                               "Problem in Matrix inversion, Invert()");
    }
 
    for ( j = 0; j < itsNumLin; j++)
    {
        det *= tmp.Elt(j, j);
    }

   for ( j = 0; j < itsNumLin; j++)
    {
        det *= tmp.Elt(j, j);
    }
 
    if ( fabs(det) < K_DOUBLE_TOL )
    {
        throw MathException(__LINE__, __FILE__, ERR_TOL_PB,
                              "Problem in Matrix inversion, Invert()");
    }
 
    //    perform backsubstitution
 
    for ( j = 0; j < itsNumLin; j++)
    {
        for (i=0; i<itsNumLin; i++)
        {
            col[i] = 0.0;
        }
 
        col[j] = 1.0;
 
        try
        {
            tmp.Lubksb(&indx, &col);
        }
 
        catch(MathException& m)
        {
            m.DebugPrint();
 
            throw MathException(__LINE__, __FILE__, ERR_MATRIX_INVER_PB,
                              "Problem in Matrix inversion, Invert()");
        }
 
        for (i=0; i < itsNumLin; i++)
        {
            invert->Elt(i, j) = col[i];
        }
    }
}



class ARM_TDiag : public ARM_GenMatrix 
{
    protected:

        int        itsSize;
        ARM_Vector itsLower;
        ARM_Vector itsDiag;
        ARM_Vector itsUpper;
    
    public:

        ARM_TDiag(void) 
        { 
            SetName(ARM_TDIAG);

            itsSize = 0; 
        }

        ARM_TDiag(int size)
        {
            itsSize = size;

            itsLower.Resize(size);
            itsDiag.Resize(size);
            itsUpper.Resize(size);
        }

       ~ARM_TDiag(void)
        {
        }

        void BitwiseCopy(const ARM_Object* srcTDiag)
        {
            ARM_TDiag* src = (ARM_TDiag *) srcTDiag;

            itsSize = src->itsSize;

            itsLower = src->itsLower;
            itsDiag  = src->itsDiag;
            itsUpper = src->itsUpper;
        }

        void Copy(const ARM_Object* srcMat)
        {
            ARM_GenMatrix::Copy(srcMat);
 
            BitwiseCopy(srcMat);
        }
 
        ARM_Object* Clone(void)
        {
            ARM_TDiag* theClone = new ARM_TDiag();
 
 
            theClone->Copy(this);
 
            return(theClone);
        }

        int     GetNumLines(void) const { return(itsSize); }
        int     GetNumCols(void) const {return(itsSize);}
        int     GetSize(void) const { return(itsSize); }

        ARM_Vector GetUpper(void) {return(itsUpper);}
        ARM_Vector GetDiag(void) {return(itsDiag);}
        ARM_Vector GetLower(void) {return(itsLower);}

    
        double& LowerElt(int i) {return(itsLower[i]); }
        double& DiagElt(int i) {return(itsDiag[i]); }
        double& UpperElt(int i) {return(itsUpper[i]); }
    

        ARM_Vector *LinSolve(ARM_Vector *SndMb);
};



/*******************************************************************
	Function to deal with allocation and desallocation int** 
	and double**
/*******************************************************************/

inline	int** InitializeIntMatrix(int nb1, int nb2)
{
    int** matrix;
    matrix = (int **) calloc(nb1, sizeof( int *));
    for ( long i = 0; i < nb1; i++ )
        *(matrix+i) = (int *) calloc(nb2, sizeof( int ));
    return matrix;
}


inline	void FreeDoubleMatrix( double** matrix, int nb1)
{
	for ( long i = 0; i < nb1; i++ )
    {
		free (*(matrix+i));
        *(matrix+i) = NULL;
    }
		
	free (matrix);
    matrix = NULL;
}


inline	double** InitializeDoubleMatrix(int nb1, int nb2)
{
	double** matrix;
	matrix = (double **) calloc(nb1, sizeof( double *));
	for ( long i = 0; i < nb1; i++ )
	{
		*(matrix+i) = (double *) calloc(nb2, sizeof( double ));

		memset(*(matrix+i), 0, sizeof(double)*nb2);
	}
		
	
	return matrix;
}

inline	void FreeIntMatrix(int**  matrix, int nb1)
{
	for ( long i = 0; i < nb1; i++ )
    {
		free (*(matrix+i));
        *(matrix+i) = NULL;
    }
		
	free (matrix);
    matrix = NULL;
}





inline double& ARM_pLine::operator[] (int n)
{
    return itsMatrix->Elt(itsLine,n);
}

inline  double& ARM_pLine::Elt(int n)
{
    return itsMatrix->Elt(itsLine,n);
}
	
inline int ARM_pLine::size()
{
    return itsMatrix->GetNumCols();
}

inline  double& ARM_pCol::operator[] (int n)
{
    return itsMatrix->Elt(n,itsCol);
}

inline  double& ARM_pCol::Elt(int n)
{
    return itsMatrix->Elt(n,itsCol);
}


inline int ARM_pCol::size()
{
    return itsMatrix->GetNumCols();
}


// indicage: 0 premiere colonne ligne itsDiag et on remonte le long de la diag vers la droite
inline double& ARM_pUpDiag::operator[] (int n)
{
    if (n>itsDiag)
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
        "Index is out of range");
    }
    return itsMatrix->Elt(itsDiag-n, n);
}



inline  double& ARM_pUpDiag::Elt(int n)
{
    return itsMatrix->Elt(itsDiag-n, n);
}


inline int ARM_pUpDiag::size()
{
    return itsDiag+1;
}


void VectorCopy(ARM_Vector* in, ARM_Vector*& out);
void VectorFill(ARM_Vector*& vec, int size, double x);


/*----------------------------------------------------------------------------*
    Resizes ARM_Vector to size if Vector.GetSize != size or Creates ARM_Vector.
    Elements are reset to zero if init != 0. Returns 0 if failed.
*----------------------------------------------------------------------------*/

inline void ResizeOrCreate(ARM_Vector*& Vector, int size, int init)
{

    if (Vector)
    {
        if (Vector->GetSize() != size)
        {
            Vector->Resize(size);
        }
        else if (init)
        {
            memset(Vector->GetElt(), '\0', size * sizeof(double));
        }
    }
    else
    {
        Vector = new ARM_Vector(size);
    }
}


inline void ResizeOrCreate(ARM_Matrix*& Matrix, int sizeLin, int sizeCol,
                           int init)
{
    if (Matrix)
    {
        if ((Matrix->GetNumCols()!=sizeCol)||(Matrix->GetNumLines()!=sizeLin))
        {
            Matrix->Resize(sizeLin, sizeCol);
        }
        else if (init)
        {
            memset(Matrix->GetElt(), '\0', sizeLin * sizeCol * sizeof(double));
        }
    }
    else
    {
        Matrix = new ARM_Matrix(sizeLin, sizeCol);
    }
}
inline ARM_Vector& ARM_Vector::operator=(const std::vector<double>&arg)
{
	this->~ARM_Vector(); 
	new(this)ARM_Vector(arg); 
	return *this; 
}
inline void ARM_Vector::populate(std::vector<double>&dest) const
{
	dest.resize(itsSize); 
	for(int i=0;i<itsSize;i++) dest[i]=itsElt[i]; 
}

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
