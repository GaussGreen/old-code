/*
 * $Log: linalg.cpp,v $
 * Revision 1.25  2004/06/02 12:09:12  emezzine
 * Added new constructor for ARM_Matrix and ARM_Vector from sTL vector.
 *
 * Revision 1.24  2004/01/16 08:20:59  ebenhamou
 * added comment for  warning on insert and push_back
 *
 * Revision 1.22  2003/12/15 13:51:52  emezzine
 * added find() and replace()
 *
 * Revision 1.21  2003/11/25 09:11:26  mcampet
 *  MC add MergeMatrix
 *
 * Revision 1.20  2003/10/27 09:24:59  emezzine
 * Improvment
 *
 * Revision 1.19  2003/09/24 06:34:01  ebenhamou
 * iterator used built-in double* like iterator
 *
 * Revision 1.18  2003/09/23 06:56:25  ebenhamou
 * added VectorExp/Sqrt/pow ...
 *
 * Revision 1.17  2003/09/22 15:29:55  ebenhamou
 * move another function to cpp as not inline!
 *
 * Revision 1.16  2003/09/22 15:21:06  ebenhamou
 * move the function roundsimple to cpp file
 *
 * Revision 1.15  2003/09/22 13:36:32  ebenhamou
 * more operator, const and iterator
 *
 * Revision 1.14  2003/09/09 11:32:39  ebenhamou
 * added push_back for ARM_Vector
 *
 * Revision 1.13  2003/06/17 15:27:20  emezzine
 * Correct a bug
 *
 * Revision 1.12  2003/06/11 08:31:02  emezzine
 * Added jacobi methods
 *
 * Revision 1.11  2003/02/10 17:35:35  mab
 * Added :  GetRow(int n); GetColumn (int n);
 *
 * Revision 1.10  2001/07/30 09:14:25  smysona
 * ajout du constructeur vector<doubl>
 * ajout des utilitaires de compactage
 *
 * Revision 1.9  2001/03/12 19:18:55  smysona
 * Ajout VectorCopy et VectorFill
 *
 * Revision 1.8  2001/03/12 13:09:02  mab
 *  Suppression de l'operateur double ()
 *
 * Revision 1.7  2001/01/30 09:54:31  smysona
 * Rajout LookupPrev et eact index
 *
 * Revision 1.6  2001/01/24 11:11:58  sgasquet
 * Ajout GetUpDiag
 *
 * Revision 1.5  2000/11/13 10:05:15  sgasquet
 *  Ajout fonctions alloc desalloc double*** et int*** et Ajout fonction GetLine et GetCol
 *
 * Revision 1.4  2000/07/24 17:03:37  mab
 * Rajout commets. RCS
 *
 */


/*----------------------------------------------------------------------------*

    linalg.cpp

    This file implements the ARM_GenMatrix class, a generic class 
    for dealing with matrices. 
    This class is mostly used for definition of product operator for 
    matrices
    

*----------------------------------------------------------------------------*/


#include "expt.h"
#include "linalg.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <functional>





/*---------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*
  Returns the number of lines. MUST be overriden by subclasses.
*----------------------------------------------------------------------------*/

int ARM_GenMatrix::GetNumLines(void) const
{
    throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
             "Unimplemented <GetNumLines> method");

    return(-1);
}


    
/*----------------------------------------------------------------------------*
    Returns the number of columns. MUST be overriden by subclasses.
*----------------------------------------------------------------------------*/

int ARM_GenMatrix::GetNumCols(void) const
{
    throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                         "Unimplemented <GetNumCols> method");

    return(-1);
}



/*----------------------------------------------------------------------------*
    Set matrix to identity. Must be a square matrix. Returns 0 if failed.
*----------------------------------------------------------------------------*/

int ARM_GenMatrix::SetToId(void)
{
    int    i, j;



    if (!(IsSquare())) 
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                         "Matrix is not square in <SetToId> method");
       return(0);
    }
    
    for (i = 0; i < GetNumLines(); i++) 
    {
        for (j=0; j<GetNumCols(); j++) 
         {
            Elt(i, j) = 0.0;
        }

        Elt(i, i) = 1.0;
    }
    
    return(1);
}



/*----------------------------------------------------------------------------*
    returns matrix trace.
*----------------------------------------------------------------------------*/

double ARM_GenMatrix::Trace(void)
{
    int    i;
    double tr=0.0;



    if (!IsSquare()) 
    {
        throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                         "Matrix is not square in <Trace> method");
        return(K_HUGE_DOUBLE);
    }

    for (i = 0; i < GetNumLines(); i++) 
    {
        tr += Elt(i,i);
    }

    return(tr);
}



/*----------------------------------------------------------------------------*
 operator == (ARM_GenMatrix &m1, ARM_GenMatrix &m2)
 operator != (ARM_GenMatrix &m1, ARM_GenMatrix &m2)
 operator >  (ARM_GenMatrix &m1, ARM_GenMatrix &m2)
 operator <  (ARM_GenMatrix &m1, ARM_GenMatrix &m2)
 operator >= (ARM_GenMatrix &m1, ARM_GenMatrix &m2)
 operator <= (ARM_GenMatrix &m1, ARM_GenMatrix &m2)

 Implementation of comparison operators for matrices
*----------------------------------------------------------------------------*/
    
int operator == (const ARM_GenMatrix& m1, const ARM_GenMatrix& m2)
{
    int    i, j;


    if (( m1.GetNumLines() != m2.GetNumLines() ) 
        ||  
        ( m1.GetNumCols() != m2.GetNumCols() ) 
       )
    {
       return(0);    
    }     
    
    for (i = 0; i < m1.GetNumLines(); i++) 
    {
        for (j = 0; j < m1.GetNumCols(); j++) 
        {
            if ( m1.Elt(i,j) != m2.Elt(i,j) )
            {
               return(0);
            }
        }
    }

    return(1);
}
   


int operator != (const ARM_GenMatrix& m1, const ARM_GenMatrix& m2)
{
    return(!(m1 == m2));
}



int operator > (const ARM_GenMatrix& m1, const ARM_GenMatrix& m2)
{    
    int    i, j;


    if ( ( m1.GetNumLines() != m2.GetNumLines() ) 
         || 
         ( m1.GetNumCols() != m2.GetNumCols() ) 
       )
    {
       return(0);    
    }

    for (i = 0; i < m1.GetNumLines(); i++) 
    {
        for (j = 0; j < m1.GetNumCols(); j++) 
        {
            if ( m1.Elt(i,j) <= m2.Elt(i,j) )
            {
               return(0);
            }
        }
    }

    return(1);
}



int operator >= (const ARM_GenMatrix& m1, const ARM_GenMatrix& m2)
{    
    int    i, j;


    if (( m1.GetNumLines() != m2.GetNumLines() )
        || 
        ( m1.GetNumCols() != m2.GetNumCols() ) 
       )
    {
       return(0);    
    }

    for (i = 0; i < m1.GetNumLines(); i++) 
    {
        for (j = 0; j < m1.GetNumCols(); j++) 
        {
            if ( m1.Elt(i,j) < m2.Elt(i,j) )
            {
               return(0);
            }
        }
    }

    return(1);
}



int operator < (const ARM_GenMatrix& m1, const ARM_GenMatrix& m2)
{    
    int    i, j;


    if (( m1.GetNumLines() != m2.GetNumLines() ) 
       || 
        ( m1.GetNumCols() != m2.GetNumCols() )
       )
    {
       return(0);
    }

    for (i = 0; i < m1.GetNumLines(); i++) 
    {
        for (j = 0; j < m1.GetNumCols(); j++) 
        {
            if ( m1.Elt(i,j) >= m2.Elt(i,j) )
            {
               return(0);
            }
        }
    }

    return(1);
}



int operator <= (ARM_GenMatrix &m1, ARM_GenMatrix &m2)
{    
    int    i, j;



    if (( m1.GetNumLines() != m2.GetNumLines() )
        || 
        ( m1.GetNumCols() != m2.GetNumCols() )
       )
    {
       return(0);
    }

    for (i = 0; i < m1.GetNumLines(); i++) 
    {
        for (j = 0; j < m1.GetNumCols(); j++) 
        {
            if ( m1.Elt(i,j) > m2.Elt(i,j) )
            {
               return(0);
            }
        }
    }

    return(1);
}




/*----------------------------------------------------------------------------*
    ARM_Matrix operator | (ARM_GenMatrix &m1, ARM_GenMatrix &m2)

    horizontal concatenation
*----------------------------------------------------------------------------*/

ARM_Matrix operator | (const ARM_GenMatrix &m1, const ARM_GenMatrix &m2)
{
    ARM_Matrix tmp(m1.GetNumLines(), m1.GetNumCols()+m2.GetNumCols());
    int        i, j, decal;


    
    if ( m1.GetNumLines() != m2.GetNumLines() )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
           "Inconsistent matrix sizes in <operator | (horizontal concat)>");

        return(tmp);
    }
    
    decal = m1.GetNumCols();

    for (j = 0; j < decal; j++)
    {
        for (i = 0; i < tmp.GetNumLines(); i++)
        {
            tmp.Elt(i, j) = m1.Elt(i, j);
        }
    }


    for (j = decal; j < tmp.GetNumCols(); j++)
    {
        for(i = 0; i < tmp.GetNumLines(); i++)
        {
            tmp.Elt(i, j) = m2.Elt(i, j-decal);
        }
    }

    return(tmp);
}



/*----------------------------------------------------------------------------*
    ARM_Matrix operator & (ARM_GenMatrix &m1, ARM_GenMatrix &m2)

    vertical concatenation
*----------------------------------------------------------------------------*/

ARM_Matrix operator & (const ARM_GenMatrix &m1, const ARM_GenMatrix &m2)
{
    ARM_Matrix tmp(m1.GetNumLines()+m2.GetNumLines(), m1.GetNumCols());

    int i, j, decal;



    if ( m1.GetNumCols() != m2.GetNumCols() )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
         "Inconsistent matrix sizes in <operator & (vertical concat)>");

        return(tmp);
    }

    decal = m1.GetNumLines();

    for (i = 0; i < decal; i++)
    {
        for (j = 0; j < tmp.GetNumCols(); j++)
        {
            tmp.Elt(i,j) = m1.Elt(i,j );
        }
    }

    for (i = decal; i < tmp.GetNumLines(); i++)
    {
        for (j = 0; j < tmp.GetNumCols(); j++)
        {
            tmp.Elt(i,j) = m2.Elt(i-decal,j);
        }
    }

    return(tmp);
}



/*----------------------------------------------------------------------------*
    SYNOPSIS    int isequal (ARM_GenMatrix &m, double val )
    SYNOPSIS    int isdiffM_GenMatrix &m, double val )
    SYNOPSIS    int issup(ARM_GenMatrix &m, double val )
    SYNOPSIS    int isinf(ARM_GenMatrix &m, double val )
    SYNOPSIS    int issupequal(ARM_GenMatrix &m, double val )
    SYNOPSIS    int isinfequal(ARM_GenMatrix &m, double val )

    functions comparing a matrix to a real
*----------------------------------------------------------------------------*/

int isequal(const ARM_GenMatrix &m, double val)
{
    int        i, j;


    for (i = 0; i < m.GetNumLines(); i++)
    {
        for (j = 0; j < m.GetNumCols(); j++) 
        {
            if ( m.Elt(i,j) != val )
            {
               return(0);
            }
        }
    }

    return(1);
}



int isdiff(const ARM_GenMatrix &m, double val)
{
    int        i, j;



    for (i = 0; i < m.GetNumLines(); i++)
    {
        for (j = 0; j < m.GetNumCols(); j++) 
        {
            if ( m.Elt(i,j) == val )
            {
               return(0);
            }
        }
    }

    return(1);
}



int issup(const ARM_GenMatrix &m, double val)
{
    int        i, j;


    for (i = 0; i < m.GetNumLines(); i++)
    {
        for (j = 0; j < m.GetNumCols(); j++) 
        {
            if ( m.Elt(i,j) <= val )
            {
               return(0);
            }
        }
    }

    return(1);
}



int isinf(const ARM_GenMatrix &m, double val)
{
    int        i, j;



    for (i = 0; i < m.GetNumLines(); i++)
    {
        for (j = 0; j < m.GetNumCols(); j++) 
        {
            if ( m.Elt(i,j) >= val )
            {
               return(0);
            }
        }
    }

    return(1);
}



int issupequal(const ARM_GenMatrix &m, double val)
{
    int        i, j;


    for (i = 0; i < m.GetNumLines(); i++)
    {
        for (j = 0; j < m.GetNumCols(); j++) 
        {
            if ( m.Elt(i,j) < val ) 
            {
               return(0);
            }
        }
    }

    return(1);
}



int isinfequal(const ARM_GenMatrix &m, double val)
{
    int        i, j;


    for (i = 0; i < m.GetNumLines(); i++)
    {
        for (j = 0; j < m.GetNumCols(); j++) 
        {
            if ( m.Elt(i,j) > val )
            {
               return(0);
            }
        }
    }

    return(1);
}


// Be carful the Matrix is organized here by COLUMNS

ARM_Matrix::ARM_Matrix(int nl , int nc, double** mat)
{
    int i,j;

    itsNumLin = nl;
    itsNumCol = nc;
 
    itsElt = new double [itsNumLin*itsNumCol];

    for (i = 0; i < itsNumLin; i++)
    {
        for (j = 0; j < itsNumCol; j++)
        {
            Elt(i,j) = mat[j][i];
        }
    }
}



// Be carful the Matrix is organized here by LINES

ARM_Matrix::ARM_Matrix(double** mat, int nl , int nc)
{
    int i,j;
 
    itsNumLin = nl;
    itsNumCol = nc;
 
    itsElt = new double [itsNumLin*itsNumCol];
 
    for (i = 0; i < itsNumLin; i++)
    {
        for (j = 0; j < itsNumCol; j++)
        {
            Elt(i,j) = mat[i][j];
        }
    }
}


/*----------------------------------------------------------------------------*
    Initialize square matrix by copying t.
*----------------------------------------------------------------------------*/

ARM_Matrix::ARM_Matrix(const ARM_TDiag& t) : ARM_GenMatrix(t)
{    
    int    i, j;


    SetName(ARM_MATRIX);

    itsElt = (double *) NULL;

    try
    {
        Resize(t.GetSize(), t.GetSize());
    }
    
    catch(Exception & x)
    {
        x.DebugPrint();
        throw Exception(__LINE__, __FILE__, ERR_MEMORY_ALLOCATION,
   "Out of memory or negative size in <ARM_Matrix(ARM_TDiag &t)> constructor");
    }    

    for (i = 0; i < itsNumLin; i++)
    {
        for (j = 0; j < itsNumCol; j++) 
        {
            Elt(i, j) = t.Elt(i, j);
        }
    }
}



/*----------------------------------------------------------------------------*
    Initialize 1 column matrix by copying v.
*----------------------------------------------------------------------------*/

ARM_Matrix::ARM_Matrix(const ARM_Vector &v) : ARM_GenMatrix(v)
{    
    int    i;


    SetName(ARM_MATRIX);

    itsElt = (double *) NULL;

    try
    {
        Resize(v.GetSize(), 1);
    }
    
    catch(Exception & x)
    {
        x.DebugPrint();

        throw Exception(__LINE__, __FILE__, ERR_MEMORY_ALLOCATION,
 "Out of memory or negative size in <ARM_Matrix(ARM_Vector &v)> constructor");
    }

    for (i = 0; i < itsNumLin; i++) 
    {
        Elt(i, 0) = v.Elt(i);
    }
}
    





/*----------------------------------------------------------------------------*
    operator = (ARM_TDiag &t)
*----------------------------------------------------------------------------*/

ARM_Matrix& ARM_Matrix::operator = (const ARM_TDiag &t)
{    
    int    i, j;


    (*this).ARM_GenMatrix::operator = (t);

    try
    {
        Resize(t.GetSize(), t.GetSize());
    }
    
    catch(Exception & x)
    {
        x.DebugPrint();

        throw Exception(__LINE__, __FILE__, ERR_MEMORY_ALLOCATION,
     "Out of memory or negative size in <operator = ARM_Matrix> constructor");
    }    
    
    for (i = 0; i < itsNumLin; i++)
    {
        for (j = 0; j < itsNumCol; j++) 
        {
            Elt(i, j) = t.Elt(i, j);
        }
    }    

    return(*this);
}



/*----------------------------------------------------------------------------*
    operator = (ARM_Vector &v)
*----------------------------------------------------------------------------*/

ARM_Matrix& ARM_Matrix::operator = (const ARM_Vector &v)
{    
    int    i, j;



    (*this).ARM_GenMatrix::operator = (v);

    try
    {
        Resize(v.GetSize(), 1);
    }
    
    catch(Exception & x)
    {
        x.DebugPrint();

        throw Exception(__LINE__, __FILE__, ERR_MEMORY_ALLOCATION,
     "Out of memory or negative size in <operator = ARM_Vector> constructor");
    }    
    
    for (i = 0; i < itsNumLin; i++)
    {
        for (j = 0; j < itsNumCol; j++) 
        {
            Elt(i, j) = v.Elt(i, j);
        }
    }
     
    return(*this);
}



/*----------------------------------------------------------------------------*
    Returns opposite of matrix
*----------------------------------------------------------------------------*/

ARM_Matrix& ARM_Matrix::operator -()
{
    int    i;


    for (i = 0; i < itsNumLin*itsNumCol; i++) 
    {
        itsElt[i] = -itsElt[i];
    }  

    return(*this);
}



/*----------------------------------------------------------------------------*
    Transpose the matrix
*----------------------------------------------------------------------------*/

void ARM_Matrix::Transpose(void)
{
    int        i, j;
    ARM_Matrix tmp(itsNumCol, itsNumLin);


    tmp = *this;

    Resize(itsNumCol, itsNumLin);

    for (i = 0; i < itsNumLin; i++)
    {
        for (j = 0; j < itsNumCol; j++) 
        {
            Elt(i, j) = tmp.Elt(j, i);
        }
    }
}

    
void ARM_Matrix::QRSolve(ARM_Vector* b, double & sing)
{
    int n;
// FIXMEFRED: mig.vc8 (22/05/2007 15:50:19): int
	register int i,j,k;


    n = b->GetSize();
    
    if (!IsSquare() || b->GetSize() != itsNumLin)
    {
       throw MathException(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                        "Matrix is not square or Size != Number of lines");
    }
    
    double scale, sigma, sum, tau;
    ARM_Matrix TMP(*this);
    ARM_Vector c(n),d(n);

    sing = 0;

    for (k = 0; k < n-1; k++)
    {
        scale = 0.0;

        for (i = k; i < n; i++)
        {
            scale = MAX(scale,fabs(TMP.Elt(i,k)));
        }

        if ( scale == 0.0 )
        {
           sing=1;

           c.Elt(k)=d.Elt(k)=0.0;
        } 
        else 
        {
           for (i = k; i < n; i++) 
           {
               TMP.Elt(i,k) /= scale;
           }
        
           for (sum = 0.0, i = k; i <n; i++) 
           {
               sum += TMP.Elt(i,k)*TMP.Elt(i,k);
           }
    
           if ( TMP.Elt(k,k) >= 0.0)
           {
              sigma=sqrt(sum);
           }
           else
           {
              sigma=-sqrt(sum);
           }

           TMP.Elt(k,k) += sigma;
            
           c.Elt(k)=sigma*TMP.Elt(k,k);
            
           d.Elt(k) = -scale*sigma;
            
           for (j = k+1; j < n; j++) 
           {
               for (sum = 0.0, i = k; i < n; i++) 
               {
                   sum += TMP.Elt(i,k)*TMP.Elt(i,j);
               }

               tau = sum/c.Elt(k);
                
               for (i = k;i < n; i++) 
               {
                   TMP.Elt(i,j) -= tau*TMP.Elt(i,k);
               }
            }
        }
    }

    d.Elt(n-1)=TMP.Elt(n-1,n-1);

    if ( d.Elt(n-1) == 0.0 ) 
    {
       sing = 1;
    }

    if ( sing == 1 )
    {
       throw MathException(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
               "Problem in QRSolve : Matrix singular");
    }
    
    for (j = 0; j < n-1; j++)
    {
        for (sum=0.0, i = j; i < n; i++) 
        {
            sum += TMP.Elt(i,j)*b->Elt(i);
        }

        tau = sum/c.Elt(j);

        for (i = j; i < n; i++) 
        {
            b->Elt(i) -= tau*TMP.Elt(i,j);
        }
    }

    b->Elt(n-1) /= d.Elt(n-1);

    for (i = n-2; i >= 0; i--) 
    {
        for (sum = 0.0, j = i+1; j < n; j++) 
        {
            sum += TMP.Elt(i,j)*b->Elt(j);
        }

        b->Elt(i) = (b->Elt(i)-sum)/d.Elt(i);
    }
}

ARM_Matrix* ARM_Matrix ::Jacobi(ARM_Matrix*a, ARM_Vector& d,int& nrot)
{
	int i,j,ip,iq;
	double tresh,theta,tau,t,sm,s,h,g,c;
    int n = a->GetNumCols();
    //ARM_Matrix a(*this);

    if (!a->IsSquare() || !a->IsSymmetric())
    {
       throw MathException(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB, 
                        "Matrix must be symmetric");
    }

	ARM_Vector b(n);
	ARM_Vector z(n);

    ARM_Matrix* v = new ARM_Matrix(n,n);

    //ARM_Matrix* v = new ARM_Matrix(n,n);

	for (ip = 0; ip < n; ip++)
    {
        for (iq = 0; iq < n; iq++)
            v->Elt(ip,iq)=0.0;
		v->Elt(ip,ip)=1.0;

        b.Elt(ip) = a->Elt(ip,ip);
        d.Elt(ip) = a->Elt(ip,ip);
	}

    for (ip = 0; ip < n; ip++)
    {
        b.Elt(ip) = a->Elt(ip,ip);
        d.Elt(ip) = a->Elt(ip,ip);
        z.Elt(ip)=0.0;
	}
	nrot=0;

	for (i = 1; i <= 50; i++)
    {
		sm=0.0;
		for (ip = 0; ip < n-1; ip++)
        {
			for (iq = ip+1; iq < n; iq++)
				sm += fabs(a->Elt(ip,iq));
		}

		if (sm == 0.0)        			
			return v;
		if (i < 4)
            tresh=0.2*sm/(n*n);
		else        
            tresh=0.0;
        

		for (ip = 0; ip < n-1; ip++)
        {
			for (iq = ip+1; iq < n; iq++)
            {
				g = 100.0*fabs(a->Elt(ip,iq));
                // After four sweeps, skip the retation if the off-diagonal element is small

				if ((i > 4) && (fabs(d.Elt(ip)) + g) == fabs(d.Elt(ip))
					        && (fabs(d.Elt(iq))+g) == fabs(d.Elt(iq)))
                
					a->Elt(ip,iq)=0.0;
                

				else if (fabs(a->Elt(ip,iq)) > tresh)
                {
					h = d.Elt(iq) - d.Elt(ip);
					if ((fabs(h)+g) == fabs(h))                    
						t= (a->Elt(ip,iq))/h;
					else 
                    {
						theta = 0.5*h/(a->Elt(ip,iq));
						t = 1.0/(fabs(theta) + sqrt(1.0+theta*theta));
						if (theta < 0.0)                        
                            t = -t;                        
					}

					c   = 1.0/sqrt(1+t*t);
					s   = t*c;
					tau = s/(1.0+c);
					h   = t*a->Elt(ip,iq);
					z.Elt(ip) -= h;
					z.Elt(iq) += h;
					d.Elt(ip) -= h;
					d.Elt(iq) += h;
					a->Elt(ip,iq)=0.0;

					for (j = 0;j < ip; j++)
                    {
						Rot(a,s,tau,j,ip,j,iq);
					}

					for (j = ip+1; j < iq; j++)
                    {
						Rot(a,s,tau,ip,j,j,iq);
					}

					for (j = iq+1; j < n; j++)
                    {
						Rot(a,s,tau,ip,j,iq,j);
					}
					for (j = 0; j < n; j++)
                    {
						Rot(v,s,tau,j,ip,j,iq);
					}
					++nrot;
				}
			}
		}
		for (ip = 0 ; ip < n; ip++)
        {
			b.Elt(ip) += z.Elt(ip);
			d.Elt(ip)  = b.Elt(ip);
			z.Elt(ip)  = 0.0;
		}
	}
    throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                                        "Too many iterations in routine jacobi");
}


ARM_Matrix* ARM_Matrix::eigsrt(ARM_Vector& d)
{
	int k,j,i;
	double p;
    int n = d.GetSize();
    ARM_Matrix* v = (ARM_Matrix*)this->Clone();
    for (i = 0; i < n-1; i++)
    {
        k = i;
		p = d.Elt(k);
		for (j = i; j < n; j++)
        {
            if (d.Elt(j) >= p)
            {
                k = j;
                p = d.Elt(k);
            }
        }

        if (k != i)
        {
			d.Elt(k) = d.Elt(i); 
			d.Elt(i) = p;

			for (j = 0; j < n; j++)
            {
				p = v->Elt(j,i);
				v->Elt(j,i) = v->Elt(j,k);
				v->Elt(j,k) = p;
			}
        }
        
	}

    return v;
}

ARM_Matrix* ARM_Matrix::ACP(ARM_Matrix*a,ARM_Vector& d)
{
    int nrot;

    ARM_Matrix* v = Jacobi(a,d,nrot);
    ARM_Matrix* sortjacobi = v->eigsrt(d);
    return (sortjacobi);
}


/*----------------------------------------------------------------------------*
  Returns a matrix obtained by sorting v into ascending order 

  with respect to the col specified

  Th Shell algorithm in Num Rec in C is used
*----------------------------------------------------------------------------*/

ARM_Matrix ARM_Matrix::Sort(int col)
{
    int i, j, inc;

    ARM_Matrix tab(this);

    ARM_Vector array(*this, col);

    int n = tab.GetNumLines();
    int nCols = tab.GetNumCols();

    int k;

    ARM_Vector w(nCols, 0.0);

    double v;

    inc = 1;

    do
    {
        inc *= 3; 
        inc++;
    }
    while (inc <= n);

    do
    {
        inc /= 3; 

        for (i=inc+1; i<=n; i++) 
        {
            v = array.Elt(i-1);

            for (k=0; k<nCols; k++)
            {
                w.Elt(k) = tab.Elt(i-1, k);
            }
            
            j=i;
            
            while (array.Elt(j-inc-1) > v)
            {
                array.Elt(j-1) = array.Elt(j-inc-1);

                for (k=0; k<nCols; k++)
                {
                    tab.Elt(j-1, k) = tab.Elt(j-inc-1, k);
                }

                j -= inc;

                if (j <= inc) break;
            }
            
            array.Elt(j-1) = v;

            for (k=0; k<nCols; k++)
            {
                tab.Elt(j-1, k) = w.Elt(k);
            }
        }
    }
    while (inc>1);

    return(tab);
}


/*----------------------------------------------------------------------------*
    Solves linear system a x = y
    where a is the current matrix. Uses LU decomposition.
    The result x is returned into the vector y. 
    Also returns the matrix determinant. Returns 0 if failed.
*----------------------------------------------------------------------------*/
            
void ARM_Matrix::LinSolve(ARM_Vector* y, double& det)
{
    int        j;
    ARM_Matrix tmp;
    ARM_Vector indx(itsNumLin);
       

 
    if (!IsSquare() || y->GetSize() != itsNumLin)
    {
       throw MathException(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                        "Matrix is not square or Size != Number of lines");
    }
        
    // Copy this into tmp

    tmp = *this;

    // Perform LU decomposition and compute determinant

    try
    {
        tmp.Ludcmp(&indx, det);
    }

    catch(MathException& m)
    {
        m.DebugPrint(); 
 
        throw MathException(__LINE__, __FILE__, ERR_MATRIX_LIN_SOLVE_PB,
                 "Problem in LinSolve()");
    }

    for (j = 0; j < tmp.itsNumLin; j++) 
    {
        det *= tmp.Elt(j, j);
    }
    
    // Perform backsubstitution

    try 
    {
        tmp.Lubksb(&indx, y);
    }

    catch(MathException& m)
    {
        m.DebugPrint();
 
        throw MathException(__LINE__, __FILE__, ERR_MATRIX_LIN_SOLVE_PB,
                 "Problem in LinSolve()");
    }
}



/*---------------------------------------------------------------------------*/
/*    SYNOPSIS int ARM_Matrix::Ludcmp(ARM_Vector *indx, double *d)           */
/*    Performs LU decompsition. This is from Numerical Recipes in C          */
/*---------------------------------------------------------------------------*/

void ARM_Matrix::Ludcmp(ARM_Vector *indx, double &d)
{
    int i,imax,j,k;
    double big = 0.0, dum = 0.0, sum = 0.0, temp = 0.0;
    ARM_Vector    vv;


            
    // Check input

    if (!IsSquare() || indx->GetSize() != itsNumLin)
    {
       throw MathException(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB, 
                        "Matrix is not square or Size != Number of lines");
    }
            
    // Resize vv

    vv.Resize(itsNumLin);

    d = 1.0;
    
    for (i = 0; i < itsNumLin; i++)
    {
        big = 0.0;

        for (j = 0; j < itsNumLin; j++)
        {
            if ((temp = fabs(Elt(i,j))) > big) 
            {
               big = temp;
            }
        }

        if ( big == 0.0 )
        {
           throw MathException(__LINE__, __FILE__, RET_KO, 
                        "Problem in Ludcmp()");
        }

        vv[i] = 1.0/big;
    }

    for (j = 0; j < itsNumLin; j++)
    {
        for (i = 0; i < j; i++)
        {
            sum = Elt(i, j);

            for (k = 0; k < i; k++) 
            {
                sum -= Elt(i, k) * Elt(k, j);
            }

            Elt(i, j) = sum;
        }

        big = 0.0;

        for (i = j; i < itsNumLin; i++)
        {
            sum = Elt(i, j);

            for ( k = 0; k < j; k++) 
            {
                sum -= Elt(i, k) * Elt(k, j);
            }

            Elt(i, j) = sum;

            if ( (dum = vv[i] * fabs(sum)) >= big)
            {
                big=dum;

                imax=i;
            }
        }

        if ( j != imax )
        {
            for (k = 0; k < itsNumLin; k++)
            {
                dum = Elt(imax, k);

                Elt(imax, k) = Elt(j, k);

                Elt(j, k) = dum;
            }

            d = -d;

            vv[imax] = vv[j];
        }

        (*indx)[j] = (double) imax;

        if ( Elt(j, j) == 0.0 )
        {
            Elt(j, j) = K_TINY;
        }

        if ( j != itsNumLin-1 )
        {
            dum=1.0/Elt(j, j);

            for (i = j+1; i < itsNumLin; i++) 
            {
                Elt(i, j) *=  dum;
            }
        }
    }
}



/*---------------------------------------------------------------------------*/
/*    SYNOPSIS int ARM_Matrix::Lubksb(ARM_Vector *indx, ARM_Vector *b)       */
/*    Performs LU backsubstitution. This is from Numerical Recipes in C.     */
/*---------------------------------------------------------------------------*/

void ARM_Matrix::Lubksb(ARM_Vector *indx, ARM_Vector *b)
{
    int i,ii=-1,ip,j;
    double sum;




    if (!IsSquare() || ( indx->GetSize() != itsNumLin ) 
        || ( b->GetSize() != itsNumLin ) 
       ) 
    {
       throw MathException(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB, 
                        "Matrix is not square or Size != Number of lines");
    }
    
    for (i = 0; i < itsNumLin; i++) 
    {
        ip = (int) (*indx)[i];

        sum = (*b)[ip];

        (*b)[ip] = (*b)[i];

        if ( ii != -1 )
        {
           for ( j = ii; j <= i-1; j++) 
           {
               sum -= Elt(i, j) * (*b)[j];
           }
        }
        else if (sum) 
        {
           ii = i;
        }

        (*b)[i] = sum;
    }
    
    for (i=itsNumLin-1;i>=0;i--) 
    {
        sum = (*b)[i];

        for (j = i+1; j < itsNumLin;j++) 
        { 
            sum -= Elt(i, j) * (*b)[j];
        }

        (*b)[i] = sum/Elt(i, i);
    }
}



/*----------------------------------------------------------------------------*
    Initializes a vector with an arithmetic sequence
*----------------------------------------------------------------------------*/
 
ARM_Vector::ARM_Vector(double xdeb, double xfin, double pas)
{
    itsSize = 0;
 
    int nl, i;
 
 
 
    SetName(ARM_VECTOR);
 
    if (( (xfin-xdeb)*pas < 0.0 ) || ( pas == 0.0 ))
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
       "Invalid argument in <ARM_Vector(size, double, double)> constructor");
    }
 
    nl = 1+int(floor((double) (xfin-xdeb)/pas));
 
    itsSize = nl;
 
    itsElt = new double[itsSize];
 
    for (i = 0; i<nl; i++)
    {
        itsElt[i] = xdeb+i*pas;
    }
}



/*----------------------------------------------------------------------------*
    Initializes a vector by copy of nCol-th column of m.
*----------------------------------------------------------------------------*/

ARM_Vector::ARM_Vector(const ARM_Matrix& m, int nCol) : ARM_GenMatrix(m)
{
    int    i;



    SetName(ARM_VECTOR);

    if ( m.GetNumCols() < nCol )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "Invalid argument in <ARM_Vector(ARM_Matrix)> constructor");
    }

    itsElt = (double *) NULL;

    try
    {
        Resize(m.GetNumLines());
    }

    catch(Exception & x)
    {
        x.DebugPrint();

        throw Exception(__LINE__, __FILE__, ERR_MEMORY_ALLOCATION,
                    "Out of memory in <ARM_Vector(ARM_Matrix)> constructor");
    }    

    for (i = 0; i < m.GetNumLines(); i++) 
    {
        itsElt[i] = m.Elt(i, nCol);
    }
}



/*----------------------------------------------------------------------------*
    Initializes a vector by copy of first column of m.
*----------------------------------------------------------------------------*/
 
ARM_Vector::ARM_Vector(const ARM_Matrix& m) : ARM_GenMatrix(m)
{
    int    i;
 
 
 
    SetName(ARM_VECTOR);
 
    if ( m.GetNumCols() != 1 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "Invalid argument in <ARM_Vector(ARM_Matrix)> constructor");
    }
 
    itsElt = (double *) NULL;
 
    try
    {
        Resize(m.GetNumLines());
    }
 
    catch(Exception & x)
    {
        x.DebugPrint();

        throw Exception(__LINE__, __FILE__, ERR_MEMORY_ALLOCATION,
                    "Out of memory in <ARM_Vector(ARM_Matrix)> constructor");
    }
 
    for (i = 0; i < m.GetNumLines(); i++)
    {
        itsElt[i] = m.Elt(i, 0);
    }
}



ARM_Vector::ARM_Vector(const ARM_Matrix &m, ARM_LignOrCol type, int nCorL)
{
	if (type == ARM_COL) {
		ARM_Vector(nCorL);
	}
	else {			// type == ARM_LIGNE
		int    i;



		SetName(ARM_VECTOR);

		if ( m.GetNumLines() < nCorL )
		{
	       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "Invalid argument in <ARM_Vector(ARM_Matrix)> constructor");
		}

		itsElt = (double *) NULL;

	    try
		{
			Resize(m.GetNumCols());
		}	

		catch(Exception & x)
		{
			x.DebugPrint();

			throw Exception(__LINE__, __FILE__, ERR_MEMORY_ALLOCATION,
                    "Out of memory in <ARM_Vector(ARM_Matrix)> constructor");
		}    

		for (i = 0; i < m.GetNumCols(); i++) 
		{
			itsElt[i] = m.Elt(nCorL, i);
		}
	}
}


/*----------------------------------------------------------------------------*
    Initializes an ARM_Vector with a vector<double>
*----------------------------------------------------------------------------*/

ARM_Vector::ARM_Vector(const vector<double>& inVect)
{
    int i;
    
    Init();

    itsSize = inVect.size();

    itsElt = new double[itsSize];
 
    for (i = 0; i<itsSize; i++)
    {
        itsElt[i] = inVect[i];
    }
}

/*----------------------------------------------------------------------------*
    Initializes an ARM_Vector with a vector<vector<double> >
*----------------------------------------------------------------------------*/

ARM_Vector::ARM_Vector(vector<vector< double > >& inVect)
{
    Init();   
    for(size_t i = 0; i<inVect.size(); i++)
        itsSize += inVect[i].size();

    itsElt = new double[itsSize];
    
    for(size_t i = 0, k=0; i<inVect.size(); ++k,++i)
        for(size_t j = 0; j<inVect[i].size(); j++)
            itsElt[k] = inVect[i][j];
}


/*----------------------------------------------------------------------------*
    SYNOPSIS    ARM_Vector & ARM_Vector::operator = (ARM_Matrix &m)

    Assigns first matrix column to this.
*----------------------------------------------------------------------------*/

ARM_Vector& ARM_Vector::operator = (const ARM_Matrix &m)
{
    int    i;



    if ( m.GetNumCols() != 1 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                 "Invalid argument in <operator = (ARM_Matrix &m)>");
    }

    (*this).ARM_GenMatrix::operator = (m);

    try
    {
        Resize(m.GetNumLines());
    }

    catch(Exception & x)
    {
        x.DebugPrint();

        throw Exception(__LINE__, __FILE__, ERR_MEMORY_ALLOCATION,
                    "Out of memory in <operator = (ARM_Matrix &)>");
    }    

    for (i = 0; i < m.GetNumLines(); i++) 
    {
        itsElt[i] = m.Elt(i, 0);
    }

    return(*this);
}
    


ARM_Vector& ARM_Vector::operator -()
{
    int    i;


    for (i = 0; i < itsSize; i++) 
    {
        itsElt[i] = -itsElt[i];
    }

    return(*this);
}



/*----------------------------------------------------------------------------*
    Returns the mean of vector elements
*----------------------------------------------------------------------------*/

double ARM_Vector::Mean(void) const 
{
    int    i;
    double m = 0.0;



    for (i = 0; i < itsSize; i++) 
    {
        m += itsElt[i];
    }

    m /= (double) itsSize;

    return(m);
}




/*----------------------------------------------------------------------------*
    Returns scalar product of v1 and v2.
*----------------------------------------------------------------------------*/

double Scalar(const ARM_Vector& v1, const ARM_Vector& v2)
{
    int        i;
    double    s=0.0;

    if ( v1.itsSize != v2.itsSize )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,

     "The sizes of two vectors are different in <Scalar(ARM_Vector, ARM_Vector)> method");

        return(K_HUGE_DOUBLE);
    }

    for (i = 0; i < v1.itsSize; i++) 
    {
        s += v1.itsElt[i] * v2.itsElt[i];
    }

    return(s);
}


/*----------------------------------------------------------------------------*
    Returns a vector obtained by sorting v into ascending order

  using the Shell algorithm in Num Rec in C
*----------------------------------------------------------------------------*/

ARM_Vector ARM_Vector::Sort(int incOrdec) const 
{
    int i, j, inc;

    ARM_Vector array(this);

    int n = array.GetSize();

    double    v;
    inc = 1;

    do
    {
        inc *= 3; 
        inc++;
    }
    while (inc <= n);

    do
    {
        inc /= 3; 

        for (i=inc+1; i<=n; i++) 
        {
            v = array.Elt(i-1);
            
            j=i;
            
            while (array.Elt(j-inc-1) > v)
            {
                array.Elt(j-1) = array.Elt(j-inc-1);

                j -= inc;

                if (j <= inc) break;
            }
            array.Elt(j-1) = v;
        }
    }
    while (inc>1);

	if (incOrdec == K_DECREASING)
	{
		for (i = 0; i < n / 2; i++)
		{
			double tmp = array.Elt(i);
			array.Elt(i) = array.Elt(n-i-1);
			array.Elt(n-i-1) = tmp;
		}
	}

    return(array);
}

/*----------------------------------------------------------------------------*
    Returns a vector obtained by sorting v into ascending order

  using the Shell algorithm in Num Rec in C

  and the vector of the old nom of row ( vPermut)
*----------------------------------------------------------------------------*/

ARM_Vector ARM_Vector::Sort(ARM_Vector& vPermut, int incOrdec) const 
{
    
	int i, j, inc, k=0, d=0;

    ARM_Vector array(this);

    int n = array.GetSize();
	vPermut.Resize(n);
	for (i=0;i< n; i++) 
		vPermut[i]=(int) i;
    double    v;
    inc = 1;

    do
    {
        inc *= 3; 
        inc++;
    }
    while (inc <= n);

    do
    {
        inc /= 3; 

        for (i=inc+1; i<=n; i++) 
        {
            v = array.Elt(i-1);
            d = vPermut[i-1];
            j=i;
            
            while (array.Elt(j-inc-1)  > v )
            {
                array.Elt(j-1) = array.Elt(j-inc-1);
				vPermut[j-1] = vPermut[j-inc-1];
				j -= inc;
                if (j <= inc) break;
            }
            array.Elt(j-1) = v;
			vPermut[j-1] = d;
			
        }
    }
    while (inc>1);

	if (incOrdec == K_DECREASING)
	{
		for (i = 0; i < n / 2; i++)
		{
			double tmp = array.Elt(i);
			int ti = vPermut[i];
			array.Elt(i) = array.Elt(n-i-1);
			vPermut[i] = vPermut[n-i-1];
			array.Elt(n-i-1) = tmp;
			vPermut[n-i-1] = ti;
		}
	}
	
    return(array);
}



/*----------------------------------------------------------------------------*
    Initializes i-th elt of vector 
*----------------------------------------------------------------------------*/

void ARM_Vector::InitElt(int i, double value)
{

    if (i < 0 || i > itsSize-1 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
        "Index is negative or out of range in <InitElt(int , double)> method");
    }
    else 
    {
       itsElt[i] = value;
    }
}



/*----------------------------------------------------------------------------*
    ARM_Vector suite(double xdeb, double xfin, double pas)

    initialize a vector xdeb, xdeb+pas,.... Default pas = 1
*----------------------------------------------------------------------------*/

ARM_Vector Suite(double xdeb, double xfin, double pas)
{
    int nl,    i;


    if (( (xfin-xdeb)*pas < 0.0 ) || ( pas == 0.0 ))
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "Invalid argument in <Suite(double xdeb, double xfin, double pas)>");
    }

    nl = 1+int(floor((xfin-xdeb)/pas));

    ARM_Vector tmp(nl);

    for (i = 0; i < nl; i++) 
    {
        tmp.Elt(i) = xdeb+i*pas;
    }

    return(tmp);
}




/*----------------------------------------------------------------------------*
     ARM_TDiag
*----------------------------------------------------------------------------*/

ARM_Vector* ARM_TDiag::LinSolve(ARM_Vector *SndMb)
{
    int size = GetSize();
    int i;
    double piv = 0.0;
 
    ARM_Vector newDiag(size);
    ARM_Vector newSndMb(size);
    ARM_Vector *result = new ARM_Vector(size);
 
    newDiag.Elt(0) = itsDiag.Elt(0);
    newSndMb.Elt(0) = SndMb->Elt(0);
 
    for (i = 1; i < size; i++)
    {
        piv = itsLower.Elt(i)/newDiag.Elt(i-1);
 
        newDiag.Elt(i) = itsDiag.Elt(i) - piv*itsUpper.Elt(i-1);
 
        newSndMb.Elt(i) = SndMb->Elt(i) - piv*newSndMb.Elt(i-1);
 
        if (newDiag.Elt(i) == 0)
        {
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "Non Regular matrix ");
            return(NULL);
        }
    }
 
    result->Elt(size-1) = newSndMb.Elt(size-1)/newDiag.Elt(size-1);
 
    for (i = size-2; i >= 0 ; i--)
    {
        result->Elt(i) = newSndMb.Elt(i) - itsUpper.Elt(i)*result->Elt(i+1);

        result->Elt(i) /= newDiag.Elt(i);
    }
 
    return(result);
}



/* bisection routine fom NR                                       */
/* returns the index idx such that x in [Elt(idx), Elt(idx+1) [   */
/* Previous or equal*/

int ARM_Vector::LookupPrevIndex(double x) const 
{
	long ju,jm,jl;
	int ascnd, j ;
    int n = itsSize - 1;

	jl = -1;
	ju = n;
	ascnd = (itsElt[n-1] >= itsElt[0]);
	while (ju-jl > 1) {
		jm=(ju+jl) >> 1;
		if (x >= itsElt[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	if (x == itsElt[0])  j = 0;
	else if(x == itsElt[n-1])  j = n-1;
	else  j = jl;

    return j;
}




int ARM_Vector::LookupExactIndex(double x, double tol) const         // returns idx such that Elt(idx) == x, -1 otherwise
{
    int j = LookupPrevIndex(x);

    if (j>=0)
        if (fabs (x - itsElt[j]) < tol)
            return j;

    return -1;    
}


/***************************************************/
/*                                                 */
/*   Sort and compact to suppress repetitions      */
/*                                                 */
/***************************************************/

ARM_Vector* ARM_Vector::Sort_Compact(void)
{
    std::sort( begin(), end() );
	ARM_Vector::iterator iter = std::unique_copy( begin(), end(), begin() );
	int size = iter - begin();
	return new ARM_Vector( size, begin() );
}



ARM_pLine ARM_Matrix::GetLine(int n)
{
	ARM_pLine line(this, n);
	return line;
}


ARM_pCol ARM_Matrix::GetCol(int n)
{
	ARM_pCol col(this, n);
	return col;
}


ARM_pUpDiag ARM_Matrix::GetUpDiag(int n)
{
	ARM_pUpDiag diag(this, n);
	return diag;
}


ARM_Vector* ARM_Matrix::GetRow(int n)
{
    ARM_Vector* myRow = new ARM_Vector(GetNumCols());

    for (int i=0;i<GetNumCols();i++)
        myRow->Elt(i) = this->Elt(n,i);

    return myRow;
}

// f_ :first index, l_ :last index,

ARM_Matrix* ARM_Matrix::GetMatrix(int f_row, int l_row,
                                  int f_col, int l_col)
{
    int nbrow = l_row - f_row + 1;
    int nbcol = l_col - f_col + 1;

    ARM_Matrix* newMatrix = new ARM_Matrix(nbrow,nbcol);

    for (int i= 0; i < nbrow; i++)
    {
        for (int j = 0; j < nbcol; j++)
        {
            newMatrix->Elt(i,j) = Elt(f_row+i,f_col+j);
        }
    } 
    Resize(nbrow,nbcol);
    return newMatrix;
}

ARM_Matrix::ARM_Matrix(int rowsNb , int colsNb, vector<vector< double > >& inVect)
{
    Init();
    SetName(ARM_MATRIX);
    itsNumLin = rowsNb;
    itsNumCol = colsNb;

    if ( !(inVect.size()))
       return;

    itsElt = new double[itsNumLin*itsNumCol];

    int i,j;
    int k= 0;
    for(i = 0; i<rowsNb;++i)
    {
        for(j = 0; j<colsNb; j++)
        {
            double xx = inVect[i][j];
            itsElt[k++] = inVect[i][j];
        }
    }

}

ARM_Vector* ARM_Matrix::GetColumn(int n)
{
    ARM_Vector* myCol = new ARM_Vector(GetNumLines());

    for (int i=0;i<GetNumLines();i++)
        myCol->Elt(i) = this->Elt(i,n);

    return myCol;
}


void VectorCopy(ARM_Vector* in, ARM_Vector*& out)
{
    if (in)
    {
        int i;
        int size = in->GetSize();

        ResizeOrCreate(out, size);
        double* inElt  = in->GetElt();
        double* outElt = out->GetElt();

        for (i=size-1; i>= 0; i--)
        {
            outElt[i] = inElt[i];
        }
    }
    else
    {
        if (out)
            delete out;
    }
}

void VectorFill(ARM_Vector*& vec, int size, double x)
{
    if (vec)
    {
        int i;

        ResizeOrCreate(vec, size);

        if (x !=0.0)
        {
            double* vecElt  = vec->GetElt();

            for (i=size-1; i>= 0; i--)
            {
                vecElt[i] = x;
            }
        }
    }
    else
    {
        vec = new ARM_Vector(size, x);
    }
}



/*!
 * push_back: more standard interface
 */
void ARM_Vector::push_back(double anElt)
{ 
	insert( anElt ); 
}


int ARM_Vector::find(double anElt,double tol) const 
{
    int j = -1;
    for (int i = 0; i<itsSize; ++i)
    {

        if (fabs(anElt- itsElt[i]) < tol)
        {
            j = i;
            break;
        }
    }

    return j;

}
void ARM_Vector::replace(size_t i, double anElt)
{
    itsElt[i] = anElt;

}


/*----------------------------------------------------------------------------*
    SYNOPSIS  ARM_Vector& ARM_Vector::operator = (ARM_Vector &v)
 
    Assignment.
*----------------------------------------------------------------------------*/
 
ARM_Vector& ARM_Vector::operator = (const ARM_Vector& v)
{
    (*this).ARM_GenMatrix::operator = (v);
 
    this->BitwiseCopy(&v);
 
    return(*this);
}


ARM_Vector& ARM_Vector::operator +=(const ARM_Vector& Data)
{
    int Size = Data.GetSize();

    if ( GetSize() < Data.GetSize() )
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                       "Invalid Size");
    int i;

    for (i = 0; i < Size; i++)
        itsElt[i] += Data.itsElt[i];

    return(*this);
}

ARM_Vector& ARM_Vector::operator -=(const ARM_Vector& Data)
{
	int Size = Data.GetSize();
	if (GetSize() < Data.GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                 "Invalid Size");
	int i;

	for (i = 0; i < Size; i++)
		itsElt[i]-=Data.itsElt[i];
	return(*this);
}



/*!
 * operator -= with double
 */
double ARM_Vector::operator * (const ARM_Vector& Data) const
{
    int Size = Data.GetSize();

    if ( GetSize() != Data.GetSize() )
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                       "Vector Sizes are differents!");
    int i;
    double sum = 0.0;

    for (i = 0; i < Size; i++)
        sum+= itsElt[i]*Data.itsElt[i];

    return(sum);
}





///////////////////////////////////////////////////////////
/// operator with double
///////////////////////////////////////////////////////////


/*!
 * operator += with double
 */
ARM_Vector& ARM_Vector::operator += (double x)
{
    int    i;
    for (i = 0; i < GetSize(); i++)
        itsElt[i] += x;
    return(*this);
}


/*!
 * operator -= with double
 */
ARM_Vector& ARM_Vector::operator -= (double x)
{
    int    i;
    for (i = 0; i < GetSize(); i++)
        itsElt[i] -= x;
    return(*this);
}


/*!
 * operator /= with double
 */
ARM_Vector& ARM_Vector::operator /= (double x)
{
    int    i;

    if ( x == 0.0 )
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                 "Division by zero");

    for (i = 0; i < GetSize(); i++)
        itsElt[i] /= x;

    return(*this);
}


/*!
 * operator *= with double
 */
ARM_Vector& ARM_Vector::operator *= (double x)
{
    int    i;

    for (i = 0; i < GetSize(); i++)
        itsElt[i] *= x;

    return(*this);
}



double roundsimple(double num, int places)
{
    double temp, mult;
    mult = pow(10.0, places);
    temp = floor(num * mult + 0.5);
    temp = temp / mult;
    return temp;

}

ARM_Vector round(const ARM_Vector& v, int digitNb)
{
    ARM_Vector tmp(v );
	ARM_Vector::iterator begin	= tmp.begin();
	ARM_Vector::iterator end	= tmp.end();
	std::transform( begin, end, begin, std::bind2nd( std::ptr_fun( roundsimple ), digitNb ) );
    return(tmp);
}



ARM_Vector VectorIf( const ARM_Vector& cond, const ARM_Vector& argTrue, const ARM_Vector& argFalse )
{
	int size = argTrue.size() < argFalse.size()? argTrue.size(): argFalse.size();

	if( cond.size() < size )
		size = cond.size();

	ARM_Vector result( size );

	ARM_Vector::const_iterator
			elemCond = cond.begin(),
			elemA 	 = argTrue.begin(),
			elemB    = argFalse.begin();

	ARM_Vector::iterator
			elemResult = result.begin(),
			endResult  = result.end();

	for( ; elemResult != endResult; ++elemResult, ++elemCond, ++elemA, ++elemB)
		*elemResult == *elemCond? *elemA : *elemB; 

	return result;
}



ARM_Vector VectorMin( const ARM_Vector& A, const ARM_Vector& B )
{
	int size = A.size() < B.size()? A.size(): B.size();
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



ARM_Vector VectorMax( const ARM_Vector& A, const ARM_Vector& B)
{
	int size = A.size() < B.size()? A.size(): B.size();
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



ARM_Vector VectorExp( const ARM_Vector& v)
{
	return exp( v );
}

ARM_Vector VectorLog( const ARM_Vector& v)
{
	return log( v );
}

ARM_Vector VectorSqrt( const ARM_Vector& v)
{
	return sqrt( v );
}

ARM_Vector VectorPow( const ARM_Vector& v1, const ARM_Vector& v2)
{
	return pow( v1, v2 );
}



// Merge de deux matrices
ARM_Matrix* MergeMatrix(ARM_Matrix* mat1, ARM_Matrix* mat2)
{

	if (mat1==NULL)
		return mat2;
	if (mat2==NULL)
		return mat1;
	int nbl1, nbl2;
	int nbc1, nbc2;
	nbl1= mat1->GetNumLines();
	nbl2= mat2->GetNumLines();
	nbc1= mat1->GetNumCols();
	nbc2= mat2->GetNumCols();
	ARM_Matrix* mat=NULL;

	if (nbl1 == nbl2)
		// on merge les 2 matrices
	{
		mat = new ARM_Matrix(nbl1,nbc1+nbc2);
		for (int i=0; i<nbl1; i++)
		{
			for (int j=0; j<nbc1+nbc2; j++)
			{
				if (j<nbc1)
					mat->Elt(i,j)=mat1->Elt(i,j);
				else
					mat->Elt(i,j)=mat2->Elt(i,j-nbc1);

			}
		}
	}
	else
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"Cannot merge 2 matrixs");
	return mat;
}

ARM_Vector operator + (const ARM_Vector& m1, const ARM_Vector& m2)
{
    int    i;
	ARM_Vector Output;

    if ( m1.GetSize() != m2.GetSize() ) 
    {
       return Output;    
    }     

	Output.Resize(m1.GetSize());
    
    for (i = 0; i < m1.GetSize(); i++) 
    {Output.Elt(i)=m1.Elt(i)+m2.Elt(i);}

    return (Output);
}

ARM_Vector operator * (const double & a, const ARM_Vector & v)
{
    int    i;
	ARM_Vector Output;

	Output.Resize(v.GetSize());
    
    for (i = 0; i < v.GetSize(); i++) 
    {Output.Elt(i)=a* v.Elt(i);}

    return (Output);
}
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

