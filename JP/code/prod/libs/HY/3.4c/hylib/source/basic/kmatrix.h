#ifndef _MATRIX_H_
#define	_MATRIX_H_
#include "kinput.h"
#include "kfunction.h"
#include "kvalarray.h"

template <class T>
class	KMatrix	
{
public:
	static KMatrix identity (int dim);
	KMatrix(void) : _m(NULL), _rows(0), _cols(0) {}
	KMatrix(int	row, int col, const T&t=T()) {alloc(row,col); *this = t;}
	KMatrix(KInput &);
	KMatrix (const KMatrix&);
	~KMatrix(void){destroy();}
	int		row_size(void)const{return _rows;}
	int		col_size(void)const{return _cols;}
	void	resize(int row, int col){if (row != _rows || col !=_cols) {destroy(); alloc(row, col);}} 	
	KValarray<T> get_row(int i)const;
	KValarray<T> get_column(int i)const;
	KValarray<T>& get_row_reference(int i);
	KValarray<T>& get_column_reference(int i);
	void	store_row(int i, const KValarray<T>& v);
	void	store_column(int i, const KValarray<T>& v);
	KMatrix& operator= (const KMatrix& m);
	KMatrix& operator = (const T&t);
	KMatrix&	operator += (const KMatrix &m){ApplyOp(m, mPlus); return *this;}
	KMatrix&	operator -= (const KMatrix &m){ApplyOp(m, mMinus); return *this;}
	KMatrix&	operator *= (const T &t){ApplyOp(t, mMult); return *this;}
   	KMatrix&	operator /= (const T &t){ApplyOp(t, mDiv); return *this;}
	KMatrix	operator + (const KMatrix &m)const{return Combine(*this, m, mPlus);}
	KMatrix	operator - (const KMatrix &m)const{return Combine(*this, m, mMinus);}
	KMatrix	operator * (const KMatrix &m)const;
	KMatrix	operator * (const T &t)const{return Combine(*this, t, mMult);}
   	KMatrix	operator / (const T &t)const{return Combine(*this, t, mDiv);}
	//	T determinant() const
	//	KMatrix	inverse()const
	KMatrix	transpose(void) const;
   	KInput	print()const;
	friend std::ostream & operator<<(std::ostream &out, const KMatrix &m)
	{out<<m.print(); return out;}
	friend	KMatrix	operator * (const T &t, const KMatrix &m){return m * t;}
	operator T**() {return _m;}
	operator T**() const {return _m;}
public:
	void apply_equal(const KFunction<T, T>&f);
	/*	template <class Result>
	KMatrix<Result> apply(const KFunction<T, Result>&f)const
	{
		KMatrix<Result>	ans(_rows, _cols);
		int	i, j;
		T	*tmpPtr1;
		Result	*tmpPtr2;
		for (i = 0; i < _rows; i++) 				
		{
			tmpPtr1 = _m[i];
			tmpPtr2 = ans._m[i];
			for (j = 0; j < _cols; j++)
				*tmpPtr2++ = f(*tmpPtr1++);
		}
		return ans;
		}	*/			
private:
	int _rows, _cols;
	T** _m;
	KValarray<T> _row, _col;
	void destroy()
	{
		if (_m != NULL) {
			for (int i = 0 ; i < _rows; i++) 
				if (_m[i] != NULL) 
					delete [] _m[i];
				delete [] _m;
		}
	}
	void alloc (int row, int col)
	{
		_rows = row;
		_cols = col;
		if ((_m = new T* [_rows]) == NULL)
			throw KException ("Out of space");
		int	i;
		for (i = 0; i < _rows; i++) 
			_m[i] = NULL;
		for (i = 0; i < _rows; i++)
			if ((_m[i] = new T [_cols]) == NULL)
			{
				destroy(); 
				throw KException ("Out of space");
			}
	}
	bool	CheckSameSize(const KMatrix &m)
	{return (_rows == m._rows) && (_cols == m._cols);}
	
	typedef void (*Op) (T&, const T&);
	static void mMinus (T& a, const T &b) {a-=b;}
	static void mPlus (T& a, const T &b) {a+=b;}
	static void mMult (T& a, const T &b) {a*=b;}
	static void mDiv (T& a, const T &b) {a/=b;}
	
	static KMatrix Combine (const KMatrix& m1, const KMatrix& m2, Op op)
	{
		if (!m1.CheckSameSize(m2))
			throw KException ("KMatrix size mismatch ") << m1 << " and " << m2;
		KMatrix ans(m1._rows, m1._cols);
		int	i, j;
		T	*tmpPtr1, *tmpPtr2, *tmpPtr3;
		for (i = 0; i < _rows; i++)
		{
			tmpPtr1 = ans[i];
			tmpPtr2 = m1[i];
			tmpPtr3 = m2[i];
			for (j = 0; j < _cols; j++) 
			{
				*tmpPtr1 = *tmpPtr2++;
				op (*tmpPtr1++, *tmpPtr3++);
			}
		}
		return ans;
	}
	
	void ApplyOp (const KMatrix& m, Op op)
	{
		if (!CheckSameSize(m))
			throw KException ("KMatrix size mismatch ") << *this << " and " << m;
		int	i, j;
		T	*tmpPtr1, *tmpPtr2;
		for (i = 0; i < _rows; i++)
		{
			tmpPtr1 = _m[i];
			tmpPtr2 = m._m[i];
			for (j = 0; j < _cols; j++) 
				op (*tmpPtr1++, *tmpPtr2++);
		}
	}
	
	static KMatrix Combine (const KMatrix& m1, const T& t, Op op)
	{
		KMatrix ans(m1._rows, m1._cols);
		int	i, j;
		T	*tmpPtr1, *tmpPtr2;
		for (i = 0; i < _rows; i++) 
		{
			tmpPtr1 = ans._m[i];
			tmpPtr2 = m1._m[i];
			for (j = 0; j < _cols; j++) 
			{
				*tmpPtr1 = *tmpPtr2++;
				op (*tmpPtr1++, t);
			}
		}
		return ans;
	}
	
	void ApplyOp (const T& t, Op op)
	{
		int	i, j;
		T	*tmpPtr;
		for (i = 0; i < _rows; i++) 
		{
			tmpPtr = _m[i];
			for (j = 0; j < _cols; j++) 
				op (*tmpPtr++,, t);
		}
	}
	
};

template <class T1, class T2, class Result>
KMatrix<Result> apply(const KMatrix<T1>& m1, const KMatrix<T2> &m2, KBinaryFunction<T1, T2, Result>&f)
{
	if (m1.row_size()!= m2.row_size() || m1.col_size() != m2.col_size())
		throw KException ("KMatrix size mismatch ");
	
	KMatrix<Result>	ans(m1.row_size(), m1.col_size());
	for (int i = 0; i < m1.row_size(); i++) 
		for (int j = 0; j < m1.col_size(); j++)
			ans[i][j] = f(m1[i][j], m2[i][j]);

	return ans;
}

template	<class T> inline
KValarray<T> KMatrix<T>::get_row(int i)const 
{
	KValarray<T> ans (_cols); 
	for (int j = 0; j < _cols; j++) 
		ans[j] = _m[i][j]; 
	return ans;
}
template	<class T> inline
KValarray<T> KMatrix<T>::get_column(int i)const
{
	KValarray<T> ans (_rows); 
	for (int j = 0; j < _rows; j++) 
		ans[j] = _m[j][i]; 
	return ans;
}
template	<class T> inline
KValarray<T>& KMatrix<T>::get_row_reference(int i)
{
	_row.resize(_cols); 
	for (int j = 0; j < _cols; j++) 
		_row[j] = _m[i][j]; 
	return _row;
}
template	<class T> inline
KValarray<T>& KMatrix<T>::get_column_reference(int i)
{
	_col.resize(_rows); 
	for (int j = 0; j < _rows; j++) 
		_col[j] = _m[j][i]; 
	return _col;
}

template	<class T> inline
void	KMatrix<T>::store_row(int i, const KValarray<T>& v)
{
	if (v.size() != _cols) 
		throw KException ("Size mismatch");
	for (int j = 0; j < _cols; j++) 
		_m[i][j] = v[j];
}
template	<class T> inline
void	KMatrix<T>::store_column(int i, const KValarray<T>& v)
{
	if (v.size() != _rows) 
		throw KException ("Size mismatch");
	for (int j = 0; j < _rows; j++) 
		_m[j][i] = v[j];
}
template	<class T>
void	KMatrix<T>::apply_equal(const KFunction<T, T>&f)
{
	int	i, j;
	T	*tmpPtr;
	for (i = 0; i < _rows; i++) 
	{
		tmpPtr = _m[i];
		for (j = 0; j < _cols; j++)
		{
			*tmpPtr = f(*tmpPtr);
			tmpPtr++;
		}
	}
}

template	<class T> inline
KMatrix<T>& KMatrix<T>::operator= (const KMatrix<T>& m)
{
	if (!CheckSameSize(m)) 
	{
		destroy(); 
		alloc(m._rows, m._cols);
	}
	for (int i = 0; i < _rows; i++) 
		for (int j = 0; j < _cols; j++) 
			_m[i][j] = m[i][j];
		return *this;
}
template	<class T> inline
KMatrix<T>::KMatrix(const KMatrix<T>& m)
{
	alloc(m._rows, m._cols);
	for (int i = 0; i < _rows; i++) 
		for (int j = 0; j < _cols; j++) 
			_m[i][j] = m[i][j];
}

template	<class T> inline
KMatrix<T>& KMatrix<T>::operator = (const T&t)
{
	for (int i = 0; i < _rows; i++) 
		for (int j = 0; j < _cols; j++) 
			_m[i][j] = t;
		return *this;
}

template	<class T> inline
KMatrix<T>	KMatrix<T>::transpose(void) const
{
	KMatrix<T> ans (_cols, _rows);
	for (int i = 0; i < _rows; i++) 
		for (int j = 0; j < _cols; j++) ans[j][i] = _m[i][j];
		return ans;
}

template	<class T> inline
KMatrix<T>	KMatrix<T>::operator * (const KMatrix<T> &m)const
{
	if (_cols != m._rows) 
		throw KException ("KMatrix mult mismatch") << *this << " and " << m;
	KMatrix<T> ans (_rows, m_cols);
	for (int i = 0; i < ans._rows; i++) {
		for (int j = 0; j < ans._cols; j++) {
			for (int k = 0; k < _cols; k++) ans[i][j] += _m[i][k] * m[k][j];
		}
	}
	return ans;
}

template	<class T> inline
KMatrix<T>::KMatrix<T>(KInput &input)
{
	int	i, j;
	KInputVector	&v = input.get_vector();
	
	alloc(v.size(), v[0].get_vector().size());
	for(i=0; i<_rows; i++)
	{
		KInputVector	&r = v[i].get_vector();
		for(j=0; j<_cols; j++)
			r[j] >> _m[i][j];
	}
}

template	<class T> inline
KInput	KMatrix<T>::print()const
{
	int	i, j;
	KInputVector input_v;
	for(i = 0; i< _rows; i++)
	{
		KInputVector	input_r;
		for(j = 0; j<_cols; j++)
		{
			KInput	tmp;
			_m[i][j] >> tmp;
			input_r.push_back(tmp);
		}
		input_v.push_back(KInput(input_r));
	}
	return input_v;
}

#endif



