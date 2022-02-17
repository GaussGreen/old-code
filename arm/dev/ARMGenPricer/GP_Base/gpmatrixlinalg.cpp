/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gpmatrix.cpp
 *
 *  \brief file foe the generic matrices
 *	\author  main = E. Benhamou, secondary = JM Prie, E Ezzine
 *	\version 1.0
 *	\date October 2003
 */


#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/gpmatrix.h"
#include <cmath>        /// for fabs
#include <algorithm>    /// for sort
#include <functional>   /// for greater


CC_BEGIN_NAMESPACE( ARM )

#define SIGN(a,b) (b>=0 ? fabs(a):-fabs(a))
#define MAXD(a,b) (a>b ? a : b)
#define MINI(a,b) (a<b ? a : b)

////////////////////////////////////////////////
///	Utility function
///	Routine: ACPTransformation
///	Returns: ARM_GP_Matrix
///	Action : Use Analysis of composantes principlas
////////////////////////////////////////////////////
ARM_GP_Matrix* ACPTransformation(ARM_GP_Matrix* matrix,ARM_GP_Vector& eigenvalues, int nbFactors)
{
    /// local variable to store the number 
    /// of jacobi rotations that were required
	//if (!matrix->IsDiagonal())
	{
		ARM_GP_Matrix* copyMatrix = static_cast<ARM_GP_Matrix*> (matrix->Clone());
		/// check that matrix is available
		copyMatrix->CheckSymmetricMatrix();

		int nrot;
		ARM_GP_Matrix* jacobimatrix = JacobiTransformation(*copyMatrix,eigenvalues,nrot);
		ARM_GP_Matrix* matrixsorted = SortedWEigenValues(jacobimatrix, eigenvalues);

		if (nbFactors != -1)
		{
			ARM_GP_Matrix* tmpMatrix = new ARM_GP_Matrix(*matrixsorted,0,0,matrixsorted->rows(),(size_t)nbFactors);
			delete matrixsorted;
			jacobimatrix = tmpMatrix;
		}

		delete copyMatrix;
		return jacobimatrix;
	}
	// We don't do anything when the matrix is diagonal
    // Just output diagonal elts in decreasing order else non consistent with code above !
	//else
	// MARCHE PAS ...
	{
		int nbDims = matrix->GetRowsNb();

		ARM_GP_Matrix* identity = new ARM_GP_Matrix(nbDims,(nbFactors==-1?nbDims:nbFactors),0.0);

		size_t i;

        eigenvalues.resize(nbDims);
        ARM_GP_Vector oldEigenValues(nbDims);
		for(i = 0; i < nbDims; ++i)
        {
			eigenvalues[i]    = (*matrix)(i,i);
			oldEigenValues[i]   = eigenvalues[i];
        }

	    CC_NS(std,stable_sort)( eigenvalues.begin(), eigenvalues.end(), CC_NS(std,greater)<double>() );

        ARM_GP_Vector::iterator pos,prevPos=oldEigenValues.end(),beginPos=oldEigenValues.begin();
        size_t offset=0;
		for (i = 0; i < (nbFactors==-1?nbDims:nbFactors); ++i)
		{
            pos = oldEigenValues.find(eigenvalues[i]);
            if(pos==prevPos)
                /// Eigen values are identical !
                ++offset;
            else
            {
                /// Reset offset
                offset=0;
                prevPos=pos;
            }
			(*identity)(pos-beginPos+offset,i) = 1.0;
		}

		return identity;
	}
}

double ACPTransformationWithRescalling(ARM_GP_Matrix* matrix, ARM_GP_Vector& eigenvalues, ARM_GP_Matrix& ACPMatrix)
{
#if defined(__GP_STRICT_VALIDATION)
	if (eigenvalues.size() != ACPMatrix.cols())
		ARM_THROW( ERR_INVALID_ARGUMENT, "ACPTransformationWithRescalling : eigenvalues and ACPMatrix should have the same size!" );
#endif

	size_t initFactorsNb = matrix->cols();
	size_t factorsNb = ACPMatrix.cols();

#if defined(__GP_STRICT_VALIDATION)
	if (initFactorsNb < factorsNb)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ACPTransformationWithRescalling : the initial number of factors should be greater than the final one!" );
#endif

	ARM_GP_Vector tmpEigenValues(initFactorsNb);
	ARM_GP_Matrix *tmpACPMatrix = ACPTransformation(matrix,tmpEigenValues);

	size_t i,j,k;

	int maxRank = -1;
	for(i=0;i<factorsNb;++i)
		if(fabs(tmpEigenValues[i]) < K_NEW_DOUBLE_TOL)
		{
			tmpEigenValues[i] = 0.;
			if (maxRank == -1)
				maxRank = i;
		}
		else
		{
			if (tmpEigenValues[i] < 0.)
			{
#if defined(__GP_STRICT_VALIDATION)
//			ARM_THROW( ERR_INVALID_ARGUMENT, "ACPTransformationWithRescalling: negative eigenvalue" );
#endif
				tmpEigenValues[i] = 0.;
			}
		}
		
		if (maxRank == -1) maxRank = factorsNb;
		size_t effectiveRank = ( maxRank < factorsNb ? maxRank : factorsNb);

		// rescaling
		for (i=0;i<initFactorsNb;i++)
		{
			double sum=0.;
			for(j=0;j<effectiveRank-1;++j){
				sum+=(*tmpACPMatrix)(i,j)*(*tmpACPMatrix)(i,j)*tmpEigenValues[j];
			}
			double sgn = ((*tmpACPMatrix)(i,j)>0.?1.:-1.);
			double res = (*matrix)(i,i)-sum;
			if (res>K_NEW_DOUBLE_TOL)
				(*tmpACPMatrix)(i,j)=sgn*sqrt(res/tmpEigenValues[j]);
		}

		for(i=0;i<factorsNb;++i)
		{
			eigenvalues[i] =	tmpEigenValues[i];
			for(j=0; j<initFactorsNb; ++j)
			{	
				ACPMatrix(j,i) =	(*tmpACPMatrix)(j,i);
			}
		}

		double err = 0.0;
		double cov = 0.0;

		for(i=0;i<initFactorsNb;++i)
			for(j=0; j<initFactorsNb; ++j)
			{
				cov = 0.0;
				for (k = 0; k < factorsNb; ++k)
					cov +=	(*tmpACPMatrix)(i,k)*(*tmpACPMatrix)(j,k)*eigenvalues[k];

				err += (cov-(*matrix)(i,j))*(cov-(*matrix)(i,j));
			}
	err = sqrt(err);

	delete tmpACPMatrix;

	return err;
}

////////////////////////////////////////////////////
///	Utility function
///	Routine: JacobiTransformation
///	Returns: ARM_GP_Matrix
///	Action : Diagonalisation of matrix using jacobi method
////////////////////////////////////////////////////
ARM_GP_Matrix* JacobiTransformation(ARM_GP_Matrix& matrix,ARM_GP_Vector& eigenvalues,int& nrot)
{
    matrix.CheckSymmetricMatrix();
    size_t i,j,ip,iq;
	double tresh,theta,tau,t,sm,s,h,g,c;
    size_t n = matrix.GetRowsNb();    

	ARM_GP_Vector b(n);
	ARM_GP_Vector z(n);

    ARM_GP_Matrix* v = new ARM_GP_Matrix(n,n);
    /// Initialize v to Identity
    /// Initialize b and eigenvalues to the diagonal of matrix
	for (ip = 0; ip < n; ip++)
    {
        for (iq=0; iq<n; iq++)
            (*v)(ip,iq)=0.0;
		(*v)(ip,ip)=1.0;

        b[ip] = matrix(ip,ip);
        eigenvalues[ip] = matrix(ip,ip);
        /// The z vector will accumulate terms of
        /// the form ta(p,q) as in equation 11.1.14
        /// see Numerical Recipes second edition page 470
        z[ip]=0.0;
	}

	nrot=0;
	for (i=1; i<=50; i++)
    {
        /// sum magnitude of off-diagonal elements
		sm=0.0;
		for (ip=0; ip<n-1; ip++)
        {
			for (iq=ip+1; iq<n; iq++)
				sm += fabs(matrix(ip,iq));
		}
        /// the normal return , which relies on quadratic 
        /// convergence to machine underflow.
		if (sm == 0.0)        			
			return v;
        /// ... on the first three sweeps
		if (i < 4)
            tresh=0.2*sm/(n*n);
		else        
            tresh=0.0;
        /// thereafter
        for (ip=0; ip<n-1;ip++)
        {
			for (iq=ip+1; iq<n; iq++)
            {
				g = 100.0*fabs(matrix(ip,iq));
                // After four sweeps, skip the retation if the off-diagonal element is small
				if ((i > 4) && (fabs(eigenvalues[ip]) + g) == fabs(eigenvalues[ip])
					        && (fabs(eigenvalues[iq])+g) == fabs(eigenvalues[iq]))
                
					matrix(ip,iq) = 0.0;

				else if (fabs(matrix(ip,iq)) > tresh)
                {
					h = eigenvalues[iq] - eigenvalues[ip];
					if ((fabs(h)+g) == fabs(h))                    
						t= (matrix(ip,iq))/h;
					else 
                    {
						theta = 0.5*h/(matrix(ip,iq));
						t = 1.0/(fabs(theta) + sqrt(1.0+theta*theta));
						if (theta < 0.0)                       
                            t = -t;                        
					}

					c   = 1.0/sqrt(1+t*t);
					s   = t*c;
					tau = s/(1.0+c);
					h   = t*matrix(ip,iq);
					z[ip] -= h;
					z[iq] += h;
					eigenvalues[ip] -= h;
					eigenvalues[iq] += h;
					matrix(ip,iq)=0.0;

					for (j = 0;j < ip; j++)
						matrix.Rot(s,tau,j,ip,j,iq);
					for (j = ip+1; j < iq; j++)
						matrix.Rot(s,tau,ip,j,j,iq);
					for (j = iq+1; j < n; j++)
						matrix.Rot(s,tau,ip,j,iq,j);
					for (j = 0; j < n; j++)
						v->Rot(s,tau,j,ip,j,iq);
					++nrot;
				}
			}
		}
        /// Update eigenvaluess with the sum of ta(q,p)
        /// and reinitialise z
		for (ip=0 ; ip<n; ip++)
        {
			b[ip] += z[ip];
			eigenvalues[ip]  = b[ip];
			z[ip]  = 0.0;
		}
	}
    throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                                        "Too many iterations in routine jacobi");
}

///////////////////////////////////////////////////
///	Utility function
///	Routine: SortedWEigenValues
///	Returns: ARM_GP_Matrix
///	Action : this routine sorts the eigenvaluess into descending order
///          and rearranges the columns of matrix correspondingly
////////////////////////////////////////////////////
ARM_GP_Matrix* SortedWEigenValues(ARM_GP_Matrix* matrix,ARM_GP_Vector& eigenvalues)
{
    size_t k,j,i;
	double p;
    size_t n = eigenvalues.size();
    for (i=0; i<n-1; i++)
    {
        k = i;
		p = eigenvalues[k];
		for (j=i; j<n; j++)
        {
            if (eigenvalues[j] >= p)
            {
                k = j;
                p = eigenvalues[k];
            }
        }

        if (k != i)
        {
			eigenvalues[k] = eigenvalues[i];
			eigenvalues[i] = p;

			for (j=0; j<n; j++)
            {
				p = (*matrix)(j,i);
				(*matrix)(j,i) = (*matrix)(j,k);
				(*matrix)(j,k) = p;
			}
        }
	}
    return matrix;
}



///////////////////////////////////////////////////
///	Utility function
///	Routine: Ludcmp
///	Returns: void
///	Action : performs LU decomposition on a; 
/// this is from numerical recipies in C. 
///////////////////////////////////////////////////
void LUDecomposition( ARM_GP_Matrix & a, vector<int> * indx, double &d )
{
    size_t i,imax,j,k;
    double big = 0.0, dum = 0.0, sum = 0.0, temp = 0.0;
	size_t itsCols = a.cols();
	size_t itsRows = a.rows();



    // Check input

    if ( itsCols != itsRows || indx->size() != itsRows )
    {
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                        "Matrix is not square or Size != Number of lines" );
    }
            
    ARM_GP_Vector vv( itsRows );
    d = 1.0;

    for ( i = 0; i < itsRows ; i++ )
    {
        big = 0.0;
        for ( j = 0; j < itsRows ; j++ )
        {
            if (( temp = fabs( a(i,j) ) ) > big) { big = temp; }
        }

        if ( big == 0.0 )
        {
           throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                        "Problem in LUDecomposition" );
        }
        vv[i] = 1.0/big;
    }

    for ( j = 0; j < itsRows ; j++ )
    {
        for ( i = 0; i < j; i++ )
        {
            sum = a(i,j);
            for ( k = 0; k < i; k++ ) 
            {
                sum -= a(i,k) * a(k,j);
            }
            a(i,j) = sum;
        }
        big = 0.0;
        for ( i = j; i < itsRows ; i++ )
        {
            sum = a(i,j);
            for ( k = 0; k < j; k++ ) 
            {
                sum -= a(i, k) * a(k,j);
            }

            a(i,j) = sum;
            if ( (dum = vv[i] * fabs(sum)) >= big)
            {
                big=dum;

                imax=i;
            }
        }
        if ( j != imax )
        {
            for (k = 0; k < itsRows ; k++)
            {
                dum = a(imax,k);

                a(imax,k) = a(j,k);
                a(j,k) = dum;
            }
            d = -d;
            vv[imax] = vv[j];
        }

        (*indx)[j] = imax;
        if ( a(j,j) == 0.0 ) { a(j,j) = K_TINY; }

        if ( j != itsRows-1 )
        {
            dum=1.0/a(j,j);
            for ( i=j+1; i<itsRows ; i++) 
            {
                a(i,j)*= dum;
            }
        }
    }
}


///////////////////////////////////////////////////
///	Utility function
///	Routine: LUBackSubstitution
///	Returns: void
///	Action : performs LU backsubstitution on a
/// this is from numerical recipies in C. 
///////////////////////////////////////////////////

void LUBackSubstitution( ARM_GP_Matrix& a, vector<int>* indx, ARM_GP_Vector* b )
{
    size_t ii,ip,j;
	signed int i;
	size_t itsCols = a.cols();
	size_t itsRows = a.rows();
	bool iii = false;
    double sum;

    if ( a.cols() != a.rows() || indx->size() != itsRows )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                        "Matrix is not square or Size != Number of lines" );

    for ( i = 0; i < itsRows; i++ ) 
    {
        ip = (*indx)[i];
        sum = (*b)[ip];
        (*b)[ip] = (*b)[i];

        if ( iii ) //was : if ( ii != -1 ) in numerical recipies. 
        {
           for ( j = ii; j <= i-1; j++ ) 
           {
               sum -= a(i,j) * (*b)[j];
           }
        }
        else if ( sum ) 
        {
           ii = i;
		   iii = true;
        }
        (*b)[i] = sum;
    }
	for ( i=(signed int) itsRows-1 ; i>= 0 ; i-- )
    {
		sum = (*b)[i];
		for ( j = i+1 ; j < itsRows ; j++ ) 
		{ 
            sum -= a(i,j) * (*b)[j];
		}
		(*b)[i] = sum/a(i,i);
    }
}

///////////////////////////////////////////////////
///	Utility function
///	Routine: LUBackSubstitution
///	Returns: void
///	Action : performs LU backsubstitution on a
/// this is from numerical recipies in C. 
///////////////////////////////////////////////////

void LUBackSubstitution( ARM_GP_Matrix& a, vector<int>* indx, ARM_GP_Matrix* b )
{
    size_t ii,ip,j;
	signed int i;
	size_t itsCols = a.cols();
	size_t itsRows = a.rows();
	bool iii = false;
    double sum;

    if ( a.cols() != a.rows() || indx->size() != itsRows )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                        "Matrix is not square or Size != Number of lines" );

    for ( i = 0; i < itsRows; i++ ) 
    {
        ip = (*indx)[i];
        sum = (*b)(ip,0);
        (*b)(ip,0 ) = (*b)(i,0);

        if ( iii ) //was : if ( ii != -1 ) in numerical recipies. 
        {
           for ( j = ii; j <= i-1; j++ ) 
           {
               sum -= a(i,j) * (*b)(j,0);
           }
        }
        else if ( sum ) 
        {
           ii = i;
		   iii = true;
        }
        (*b)(i,0) = sum;
    }
	for ( i=(signed int) itsRows-1 ; i>= 0 ; i-- )
    {
		sum = (*b)(i,0);
		for ( j = i+1 ; j < itsRows ; j++ ) 
		{ 
            sum -= a(i,j) * (*b)(j,0);
		}
		(*b)(i,0) = sum/a(i,i);
    }
}


///////////////////////////////////////////////////
///	Utility function
///	Routine: LinSolve
///	Returns: void
///	Action : solves a x = y
/// the result x is returned in the vector y. Also gives det( x )
////////////////////////////////////////////////////
void LinSolve( ARM_GP_Matrix* a, ARM_GP_Vector* y )
{
	double det;
    size_t        j;
    ARM_GP_Matrix tmp;
    vector<int> indx( a->rows() );
       
    if ( a->rows() != a->cols() || a->rows() != y->size() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
            "problem in linsolve : matrix is not square or does not have the good size");
        
    // Copy this into tmp
    tmp = (*a);

    // Perform LU decomposition and compute determinant
    LUDecomposition( tmp, &indx, det);

    for (j = 0; j < a->rows() ; j++) det *= tmp(j,j);
        
    // Perform backsubstitution
    LUBackSubstitution( tmp, &indx, y );

	if ( det == 0 ) 
		throw MathException(__LINE__, __FILE__, ERR_MATRIX_LIN_SOLVE_PB,
			"Problem in LinSolve : determinant is 0");
}

///////////////////////////////////////////////////
///	Utility function
///	Routine: LinSolve
///	Returns: void
///	Action : solves a x = y
/// the result x is returned in the vector y. Also gives det( x )
////////////////////////////////////////////////////
void LinSolve( ARM_GP_Matrix* a, ARM_GP_Matrix* y )
{
	double det;
    size_t        j;
    ARM_GP_Matrix tmp;
    vector<int> indx( a->rows() );
       
    if ( a->rows() != a->cols() || a->rows() != y->rows() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
            "problem in linsolve : matrix is not square or does not have the good size");
        
    // Copy this into tmp
    tmp = (*a);

    // Perform LU decomposition and compute determinant
    LUDecomposition( tmp, &indx, det);

    for (j = 0; j < a->rows() ; j++) det *= tmp(j,j);
        
    // Perform backsubstitution
    LUBackSubstitution( tmp, &indx, y );

	if ( det == 0 ) 
		throw MathException(__LINE__, __FILE__, ERR_MATRIX_LIN_SOLVE_PB,
			"Problem in LinSolve : determinant is 0");
}




////////////////////////////////////////////////////
///	Utility function
///	Routine: LeastSquareRegression
///	Returns: void
///	Action : solves X'X B = X'Y. (what we want to find 
/// is B ;) it is return as an ARM_GP_Vector * 
////////////////////////////////////////////////////
ARM_GP_Vector * LeastSquareRegression( const ARM_GP_Matrix& X , const ARM_GP_Vector& Y )
{
	ARM_GP_Matrix XPrime( X );
	XPrime.transpose();

	ARM_GP_Vector * B = XPrime*Y;
	XPrime *= X;

	LinSolve( & XPrime, B );

	return B;
}



////////////////////////////////////////////////////
///	Utility function
///	Routine: LeastSquareRegressionSVD
///	Returns: void
///	Action : solves X'X B = X'Y. (what we want to find 
/// is B ;) it is return as an ARM_GP_Vector * 
////////////////////////////////////////////////////
ARM_GP_Vector * LeastSquareRegressionSVD(const ARM_GP_Matrix& X , const ARM_GP_Vector& Y )
{
	int j;
	int n = X.cols();

	ARM_GP_Matrix XPrime( X );
	ARM_GP_Matrix V(n,n,0.);
	ARM_GP_Vector sv(n);

	SingularValuesDecomposition(XPrime,sv,V);

	double wmax=0.;
	for(j=0;j<n;j++)
	{
		if (sv[j]>wmax)
			wmax=sv[j];
	}
	double wmin=wmax*1.e-12;
	for(j=0;j<n;j++)
	{
		if (sv[j]<wmin)
			sv[j]=0.;
	}

	ARM_GP_Vector* B = new ARM_GP_Vector(n);
	SingularValuesRecomposition(XPrime,sv,V,Y,B);
	return B;
}

void SingularValuesDecomposition(ARM_GP_Matrix& a,ARM_GP_Vector& w,ARM_GP_Matrix& v)
/*Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A =
U·W·V T. The matrix U replaces a on output. The diagonal matrix of singular values W is output
as a vector w[1..n]. The matrix V (not the transpose V T ) is output as v[1..n][1..n].*/
{
	int m = a.rows();
	int n = a.cols();

	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z;
	ARM_GP_Vector rv1(n);

	g=scale=anorm=0.0; //Householder reduction to bidiagonal form.
	for (i=0;i<n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs(a(k,i));
			if (scale) {
				for (k=i;k<m;k++) {
					a(k,i) /= scale;
					s += a(k,i)*a(k,i);
				}
				f=a(i,i);
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a(i,i)=f-g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += a(k,i)*a(k,j);
					f=s/h;
					for (k=i;k<m;k++) a(k,j) += f*a(k,i);
				}
				for (k=i;k<m;k++) a(k,i) *= scale;
			}
		}

		w[i]=scale *g;
		g=s=scale=0.0;
		if (i < m && i != n-1) {
			for (k=l;k<n;k++) scale += fabs(a(i,k));
			if (scale) {
				for (k=l;k<n;k++) {
					a(i,k) /= scale;
					s += a(i,k)*a(i,k);
				}
				f=a(i,l);
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a(i,l)=f-g;
				for (k=l;k<n;k++) rv1[k]=a(i,k)/h;
				for (j=l;j<m;j++) {
					for (s=0.0,k=l;k<n;k++) s += a(j,k)*a(i,k);
					for (k=l;k<n;k++) a(j,k) += s*rv1[k];
				}
				for (k=l;k<n;k++) a(i,k) *= scale;
			}
		}

		anorm=MAXD(anorm,(fabs(w[i])+fabs(rv1[i])));

	}

	for (i=n-1;i>-1;i--) { //Accumulation of right-hand transformations.
		if (i < n-1) {
			if (g) {
				for (j=l;j<n;j++) //Double division to avoid possible under.ow.
					v(j,i)=(a(i,j)/a(i,l))/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += a(i,k)*v(k,j);
					for (k=l;k<n;k++) v(k,j) += s*v(k,i);
				}
			}
			for (j=l;j<n;j++) v(i,j)=v(j,i)=0.0;
		}
		v(i,i)=1.0;
		g=rv1[i];
		l=i;
	}

	for (i=MINI(m,n)-1;i>-1;i--) { //Accumulation of left-hand transformations.
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) a(i,j)=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += a(k,i)*a(k,j);
				f=(s/a(i,i))*g;
				for (k=i;k<m;k++) a(k,j) += f*a(k,i);
			}
			for (j=i;j<m;j++) a(j,i) *= g;
		} else for (j=i;j<m;j++) a(j,i)=0.0;
		++a(i,i);
	}

	for (k=n-1;k>-1;k--) { //Diagonalization of the bidiagonal form: Loop over singular values, and over allowed iterations. 
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>-1;l--) { //Test for splitting.
				nm=l-1; //Note that rv1[0] is always zero.
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0; //Cancellation of rv1[l], if l > 1.
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=sqrt(f*f+g*g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=a(j,nm);
						z=a(j,i);
						a(j,nm)=y*c+z*s;
						a(j,i)=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) { //Convergence.
				if (z < 0.0) { //Singular value is made nonnegative.
					w[k] = -z;
					for (j=0;j<n;j++) v(j,k) = -v(j,k);
				}
				break;
			}
			if (its == 30) printf("no convergence in 30 svdcmp iterations");
			x=w[l]; //Shift from bottom 2-by-2 minor.
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=sqrt(f*f+1.);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0; //Next QR transformation:
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=sqrt(f*f+h*h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v(jj,j);
					z=v(jj,i);
					v(jj,j)=x*c+z*s;
					v(jj,i)=z*c-x*s;
				}
				z=sqrt(f*f+h*h);
				w[j]=z; //Rotation can be arbitrary if z = 0.
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=a(jj,j);
					z=a(jj,i);
					a(jj,j)=y*c+z*s;
					a(jj,i)=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
}

void SingularValuesRecomposition(ARM_GP_Matrix& u, ARM_GP_Vector& w, ARM_GP_Matrix& v, const ARM_GP_Vector& b, ARM_GP_Vector* x)
/*Solves A·X = B for a vector X, where A is specified by the arrays u[1..m][1..n], w[1..n],
v[1..n][1..n] as returned by svdcmp. m and n are the dimensions of a and will be equal for
square matrices. b[1..m] is the input right-hand side. x[1..n] is the output solution vector.
No input quantities are destroyed, so the routine may be called sequentially with different b’s.*/
{
	int m = u.rows();
	int n = u.cols();

	int jj,j,i;
	double s;
	ARM_GP_Vector tmp(n);
	for (j=0;j<n;j++) { //Calculate UTB.
		s=0.0;
		if (w[j]) { //Nonzero result only if wj is nonzero.
			for (i=0;i<m;i++) s += u(i,j)*b[i];
			s /= w[j]; //This is the divide by wj .
		}
		tmp[j]=s;
	}
	for (j=0;j<n;j++) { //Matrix multiply by V to get answer.
		s=0.0;
		for (jj=0;jj<n;jj++) s += v(j,jj)*tmp[jj];
		(*x)[j]=s;
	}
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/