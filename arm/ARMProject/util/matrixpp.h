#ifndef _MATRIXPP_H
#define _MATRIXPP_H

#include "expt.h"


class ARM_Matrix;

class ARM_1D
{
public :


        ARM_1D (void) {}

	virtual double& operator[] (int n)
	{
    throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                         "Unimplemented <operator> method");
	}; 

	virtual double& Elt(int n)
	{
    throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                         "Unimplemented <Elt> method");
	};
	
	virtual int size()
	{
		return 0;
	};
};


// pointeur sur une ligne de ARM_Matrix
class ARM_pLine : public ARM_1D
{
private :
	ARM_Matrix	*itsMatrix;
	int itsLine;


public :

	ARM_pLine (void) {
		itsMatrix = NULL;
		itsLine = 0;
	}

	ARM_pLine (ARM_Matrix *pMatrix, int n) {
		itsMatrix = pMatrix;
		itsLine = n;
	}

	
	ARM_pLine (const ARM_pLine& ligne) {
		itsMatrix = ligne.itsMatrix;
		itsLine = ligne.itsLine;
	}
	
	ARM_pLine& operator =(const ARM_pLine& ligne) {
		itsMatrix = ligne.itsMatrix;
		itsLine = ligne.itsLine;
		return *this;
	}


	inline ARM_pLine& operator++ (void) { itsLine++; return *this;}
	inline ARM_pLine& operator-- (void) { itsLine--; return *this;}
	inline ARM_pLine& operator+ (int n) { itsLine+=n; return *this;}
	inline ARM_pLine& operator- (int n) { itsLine-=n; return *this;}
	
	inline double& operator[] (int n);
	inline double& Elt(int n);
	
	int size();
};


class ARM_pCol : public ARM_1D
{
private :
	ARM_Matrix	*itsMatrix;
	int itsCol;


public :

	ARM_pCol (void) {
		itsMatrix = NULL;
		itsCol = 0;
	}
	
	ARM_pCol (const ARM_pCol& col) {
		itsMatrix = col.itsMatrix;
		itsCol = col.itsCol;
	}

	ARM_pCol (ARM_Matrix *pMatrix, int n) {
		itsMatrix = pMatrix;
		itsCol = n;
	}

	
	ARM_pCol& operator =(const ARM_pCol& col) {
		itsMatrix = col.itsMatrix;
		itsCol = col.itsCol;
		return *this;
	}


	inline ARM_pCol& operator++ (void) { itsCol++; return *this;};
	inline ARM_pCol& operator-- (void) { itsCol--; return *this;};
	inline ARM_pCol& operator+ (int n) { itsCol+=n; return *this;};
	inline ARM_pCol& operator- (int n) { itsCol-=n; return *this;};
	
	inline double& operator[] (int n);
	inline double& Elt(int n);
	
	int size();
};


	
// pointeur sur une ligne de ARM_Matrix
class ARM_pUpDiag : public ARM_1D
{
private :
	ARM_Matrix	*itsMatrix;
	int itsDiag;


public :

	ARM_pUpDiag (void) {
		itsMatrix = NULL;
		itsDiag = 0;
	}

	ARM_pUpDiag (ARM_Matrix *pMatrix, int n) {
		itsMatrix = pMatrix;
		itsDiag = n;
	}

	
	ARM_pUpDiag (const ARM_pUpDiag& ligne) {
		itsMatrix = ligne.itsMatrix;
		itsDiag = ligne.itsDiag;
	}
	
	ARM_pUpDiag& operator =(const ARM_pUpDiag& ligne) {
		itsMatrix = ligne.itsMatrix;
		itsDiag = ligne.itsDiag;
		return *this;
	}


	inline ARM_pUpDiag& operator++ (void) { itsDiag++; return *this;}
	inline ARM_pUpDiag& operator-- (void) { itsDiag--; return *this;}
	inline ARM_pUpDiag& operator+ (int n) { itsDiag+=n; return *this;}
	inline ARM_pUpDiag& operator- (int n) { itsDiag-=n; return *this;}
	
	inline double& operator[] (int n);
	inline double& Elt(int n);
	
	int size();
};



#endif


