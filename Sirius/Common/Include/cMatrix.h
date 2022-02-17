

#ifndef CMatrixH
#define CMatrixH

#include <assert.h>
#include "TContainer.h"

static int Max(int n,int m)
{
	if ( n > m )
		return n;
	else
		return m;
};





class CVectorString
{

	char** p;
	int _size;

public:

	CVectorString();

	virtual ~CVectorString();
	void resize(int n);
	int getsize()const{return _size;};
	char*& operator[](int n)const;
	CVectorString& operator=(const CVectorString & vec);


	char** ptr()const; // use only for numerical recipes

};
	
class CVector : public TContainer < double >
{	
	
	int size;
	
public:

	double* p;// should be private sos

	CVector() : TContainer<double>() {resize(0);}
	CVector(const CVector & rhs) : TContainer<double>() {*this = rhs;}
	CVector(int n) : TContainer<double>() {resize(n);}
	CVector(int n, const double & initval) : TContainer<double>() {resize(n,initval);}
	CVector(int n,const double* pVec)
	{
		resize(n);
		for ( int m = 0 ; m < n; m++ )
		{
			p[m] = pVec[m];
		};
	};
	CVector(const double** pVec,int nrows,int ncols=1)
	{

		assert(nrows == 1 || ncols == 1 );
		int n = Max(nrows,ncols);

		resize(n);

		if ( nrows > ncols )
		{
			for ( int m = 0 ; m < n; m++ )
			{
				p[m] = pVec[m][0];
			};
		}
		else
		{
			for ( int m = 0 ; m < n; m++ )
			{
				p[m] = pVec[0][m];
			};
		}
	};

	int resize(int n, const double &initval = 0.0)
	{
		int status = TContainer<double>::resize(n,initval);
		if ( status >0 )
		{
			size = _size;
			p = _p;
		}
		return status;
	}
	
	CVector& operator=(const CVector& rhs)
	{
		if ( this == &rhs)return *this;
		TContainer<double>::operator =( rhs );
		size=_size;
		p=TContainer<double>::_p;
		return *this;
	}
	
	double getsize() {return size;}
	const double getsize() const {return size;}
	double scalar_product(CVector& vec);
	double abs_max(int& i);
	operator  double*() { return _p;};

	 const double* getConstPtr()const{return p;};
	 double* getPtr(){return p;};


	
};	
	
class CMatrix : public TContainer < CVector >
{	
protected:
	typedef TContainer < CVector> dmatrix_type;
// data from SVD
	
	
protected:
	
	
	CVector* p;
	int _rows;
	int _cols;
	
public:
	
	CMatrix() : dmatrix_type() {resize(0,0);}
	CMatrix(const CMatrix & rhs) : dmatrix_type() {if(this!=&rhs){*this = rhs;};}
	CMatrix(int n , int m) : dmatrix_type() {resize(n,m,0.0);}
	CMatrix(int n, int m, const double & initval) : dmatrix_type() {resize(n,m,initval);}

	CMatrix(int n, int m, const double** pMat) : dmatrix_type() 
	{

		resize(n,m);
		for ( int i = 0 ; i < n; i++ )
			for ( int j = 0 ; j < m; j++ )
			{
				p[i][j] = pMat[i][j];
			}
	};



	int resize(int n, int m, const double &initval = 0.0)
	{
		int status = dmatrix_type::resize(n);
		if ( status >0 )
		{
			for ( int i = 0 ; i < n ; i ++ )
				status = (*this)[i].resize(m,initval);
				_rows = n;
				_cols = m;
				p = _p;
			
		}
		return status;
	}
	
	
	CMatrix& operator=(const CMatrix& rhs)
	{
		if ( this == &rhs)return *this;

		dmatrix_type::operator =( rhs );
		_rows = rhs.rows();
		_cols = rhs.cols();
		p=_p;
		return *this;
	}

	void transpose()
	{
		CMatrix transposed(_cols,_rows);
		
		for(int i=0;i<_cols;i++)
			for(int j=0;j<_rows;j++)
				transposed[i][j]=p[j][i];


		(*this)=transposed;
		
	}


	int rows(){return _rows;}
	int cols(){return _cols;}
	const int rows()const{return _rows;}
	const int cols()const{return _cols;}

	operator  double**() { return   &(p->p);};


};	

	

class iVector : public TContainer < int >
{	
	
	int* p;
	int size;
	
public:


	iVector() : TContainer<int>() {resize(0);}
	iVector(const iVector & rhs) : TContainer<int>() {if(this!=&rhs){*this = rhs;};}
	iVector(int n) : TContainer<int>() {resize(n);}
	iVector(int n, const int & initval) : TContainer<int>() {resize(n,initval);}
	iVector(int n,const int* pVec)
	{
		resize(n);
		for ( int m = 0 ; m < n; m++ )
		{
			p[m] = pVec[m];
		};
	};
	iVector(const int** pVec,int nrows,int ncols=1)
	{

		assert(nrows == 1 || ncols == 1 );
		int n = Max(nrows,ncols);
		resize(n);

		if ( nrows > ncols )
		{
			for ( int m = 0 ; m < n; m++ ){
				p[m] = pVec[m][0];
			};
		}
		else
		{
			for ( int m = 0 ; m < n; m++ ){
				p[m] = pVec[0][m];
			};
		}
	};

	int resize(int n, const int &initval = 0.0)
	{
		int status = TContainer<int>::resize(n,initval);
		if ( status >0 )
		{
			size = _size;
			p = _p;
		}
		return status;
	}
	
	iVector& operator=(const iVector& rhs){
		TContainer<int>::operator =( rhs );
		size=_size;
		p=TContainer<int>::_p;
		return *this;
	}
	
	double getsize() {return size;}
	const double getsize() const {return size;}
	double scalar_product(iVector& vec);
	double abs_max(int& i);
	operator  int*() { return _p;};

	 const int* getConstPtr(){return p;};
	 int* getPtr(){return p;};


	
};	


class iMatrix : public TContainer < iVector >
{	
protected:
	typedef TContainer < iVector> imatrix_type;
// data from SVD
	
	
protected:
	
	
	iVector* p;
	int _rows;
	int _cols;
	
public:
	
	iMatrix() : imatrix_type() {resize(0,0);}
	iMatrix(const iMatrix & rhs) : imatrix_type() {*this = rhs;}
	iMatrix(int n , int m) : imatrix_type() {resize(n,m,0.0);}
	iMatrix(int n, int m, const int & initval) : imatrix_type() {resize(n,m,initval);}

	iMatrix(int n, int m, const int** pMat) : imatrix_type() 
	{

		resize(n,m);
		for ( int i = 0 ; i < n; i++ )
			for ( int j = 0 ; j < m; j++ )
			{
				p[i][j] = pMat[i][j];
			}
	};



	int resize(int n, int m, const int &initval = 0.0)
	{
		int status = imatrix_type::resize(n);
		if ( status >0 )
		{
			for ( int i = 0 ; i < n ; i ++ )
				status = (*this)[i].resize(m,initval);
				_rows = n;
				_cols = m;
				p = _p;
			
		}
		return status;
	}
	
	iMatrix& operator=(const iMatrix& rhs){
		imatrix_type::operator =( rhs );
		_rows = rhs.rows();
		_cols = rhs.cols();
		p=_p;
		return *this;
	}

	int rows(){return _rows;}
	int cols(){return _cols;}
	const int rows()const{return _rows;}
	const int cols()const{return _cols;}



};	

void Max(double& current_max,int& i,CVector& x);
void Min(double& current_min,int& i,CVector& x);
	

class CMatrixString
{



protected:


	CVectorString* p;
	int _rows;
	int _cols;

public:

	CMatrixString();

	virtual ~CMatrixString();
//	void resize(int n,int m);
	int rows()const{return _rows;}
	int cols()const{return _cols;}
	CVectorString& operator[](int n)const;
	CMatrixString& operator =(CMatrixString& mat);

};	

	

	


class cTensor : public TContainer < CMatrix >
{	
protected:
	typedef TContainer < CMatrix> dtensor_type;
// data from SVD
	
	
protected:
	
	
	CMatrix * p;

	CVector _dim;

	
public:
	
	cTensor() : dtensor_type() {_dim.resize(3);resize(0,0,0);}
	cTensor(const cTensor & rhs) : dtensor_type() {*this = rhs;}
	cTensor(int n , int m,int k) : dtensor_type() {_dim.resize(3);;resize(n,m,k,0.0);}
	cTensor(int n, int m,int k, const double & initval) : dtensor_type() {_dim.resize(3);resize(n,m,k,initval);}

	int resize(int n, int m,int k, const double &initval = 0.0)
	{
		int status = dtensor_type::resize(n);
		if ( status >0 )
		{
			for ( int i = 0 ; i < n ; i ++ )
				status = (*this)[i].resize(m,k,initval);
				_dim[0] = n;
				_dim[1] = m;
				_dim[2] = k;

				p = _p;
			
		}
		return status;
	}
	
	cTensor& operator=(const cTensor& rhs){
		dtensor_type::operator =( rhs );
		_dim  = rhs._dim;
		p=_p;
		return *this;
	}



	int ndim(int i){return _dim[i];};
	const int ndim(int i)const{return _dim[i];};


    operator double***() { return (double***)p;};	
};	
	



	
#endif