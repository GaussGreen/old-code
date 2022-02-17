
#ifndef __TCONTAINER_H__
#define __TCONTAINER_H__



template < class T >
class TContainer
{
protected:

	T* _p;
	int _size;

public:


TContainer()
{
	_p = NULL;
	_size = 0;
}

TContainer(const TContainer < T > & rhs)
{
	_p = NULL;
	_size = 0;
}

TContainer(int n)
{
	_p = NULL;
	_size = 0;
	resize(n);
}

TContainer(int n, const T & initval)
{
	_p = NULL;
	_size = 0;
	resize(n,initval);
}

virtual ~TContainer(){resize(0);}

int resize(int n, const T & initval)
{
	int status = resize(n);

	for ( int i = 0 ; i < _size ; i ++ )
		_p[i] = initval;

	return status;
}

int resize(int n)
{
	//is there a smart way to do this when going from large to smaller arrays
	int i, status = -1;
	if ( n == _size ) return 1;
	
	if ( n == 0 && _size > 0 ){
		delete [] _p;
		_size = 0;
		return 1;
	}

	T *tp = _p;
	_p = NULL;

	if ( n > 0 )
		_p = new T[n];

	if ( _p != NULL || n == 0){
		if ( n < _size )
			for ( i = 0 ; i < n ; i ++ )
				_p[i] = tp[i];

		else
			for ( i = 0 ; i < _size ; i ++ )
				_p[i] = tp[i];
		
		delete [] tp;
		_size = n;
		status = 1;
	}
	else
	{
		_p = tp;
		tp = NULL;
		status = -1;
	}
	return status;
}

int size() { return _size; }

const int size() const { return _size; }

T& operator[](int n)
{
	if ( n >= 0 && n < _size )
		return _p[n];
	else {//don't know how to handle this yet ????
		assert(false);
		return _p[0];
	}
}

const T& operator[](int n) const
{
	if ( n >= 0 && n < _size )
		return _p[n];
	else {
		assert(false);//don't know how to handle this yet ????
		return _p[0];
		}
}

TContainer<T>& operator=(const TContainer<T>& rhs)
{
	resize(rhs.size());
	for ( int i = 0 ; i < _size ; i ++ )
		_p[i] = rhs[i];

	return *this;
}
};

template < class T >
class GVector : public TContainer < T >
{	
	int size;
	T* p;
	
public:


	GVector() : TContainer<T>() {resize(0);}
	GVector(const GVector & rhs) : TContainer<T>() {*this = rhs;}
	GVector(int n) : TContainer<T>() {resize(n);}
	GVector(int n, const T & initval) : TContainer<T>() {resize(n,initval);}
	GVector(int n,const T* pVec)
	{
		resize(n);
		for ( int m = 0 ; m < n; m++ ){
			p[m] = pVec[m];
		};
	};
	GVector(const T** pVec,int nrows,int ncols=1)
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
		else{
			for ( int m = 0 ; m < n; m++ ){
				p[m] = pVec[0][m];
			};
		}
	};

	T* getPtr(){return p;};// this function should not be used with the exception of numerical recipes



	int resize(int n)
	{
		int status = TContainer< T >::resize(n);
		if ( status >0 )
		{
			size = _size;
			p = _p;
		}
		return status;
	}
	
	GVector& operator=(const GVector& rhs)
	{
		if ( this == &rhs)return *this;
		TContainer<T>::operator =( rhs );
		size=_size;
		p=TContainer<T>::_p;
		return *this;
	}
	
	double getsize() {return size;}
	const double getsize() const {return size;}
	
};	


template < class T >
class GMatrix : public TContainer < GVector < T > >
{	
protected:
	typedef TContainer < GVector< T > > gmatrix_type;	
	
	GVector < T >* p;
	int _rows;
	int _cols;
	
public:
	
	GMatrix() : gmatrix_type() {resize(0,0);}
	GMatrix(const GMatrix & rhs) : gmatrix_type() {*this = rhs;}
	GMatrix(int n , int m) : gmatrix_type() {resize(n,m);}
	GMatrix(int n, int m, const T & initval) : gmatrix_type() {resize(n,m,initval);}
	GMatrix(int n, int m, const T** pMat) : gmatrix_type() 
	{
		resize(n,m);
		for ( int i = 0 ; i < n; i++ )
			for ( int j = 0 ; j < m; j++ ){
				p[i][j] = pMat[i][j];
			}
	};



	int resize(int n, int m)
	{
		int status = gmatrix_type::resize(n);
		if ( status >0 )
		{
			for ( int i = 0 ; i < n ; i ++ )
				status = (*this)[i].resize(m);
				_rows = n;
				_cols = m;
				p = _p;		
		}
		return status;
	}

/*	int resize(int n, int m, const T &initval )
	{
		int status = gmatrix_type::resize(n);
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
*/	
	GMatrix& operator=(const GMatrix& rhs)
	{
		if ( this == &rhs)return *this;
		gmatrix_type::operator =( rhs );
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




#endif
