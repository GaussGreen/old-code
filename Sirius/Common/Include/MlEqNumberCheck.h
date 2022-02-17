//	MlEqNumberCheck.h : Checks for infinities in code
//
/////////////////////////////////////////////////////////////////////////////


#undef min
#undef max
template <class T> class Number
{
protected:
	T									m_t;
	void								Check(void) const
	{
		T l = std::numeric_limits<T>::max();
		if (-m_t >= std::numeric_limits<T>::max()) ATLASSERT(false);	// note that min() defines the minimum POSITIVE value
		if (m_t >= std::numeric_limits<T>::max()) ATLASSERT(false);
	}

public:
	Number()
	{
	}
	
	Number(const T& t)
	{
		m_t = t;
		Check();
	}

	operator T(void) const
	{
		return m_t;
	}

	const T& operator=(const T& t)
	{
		m_t = t;
		Check();
		return m_t;
	}

	
	const T& operator+=(const T& t)
	{
		m_t += t;
		Check();
		return m_t;
	}

	const T& operator-=(const T& t)
	{
		m_t -= t;
		Check();
		return m_t;
	}


	const T& operator*=(const T& t)
	{
		m_t *= t;
		Check();
		return m_t;
	}

	const T& operator/=(const T& t)
	{
		m_t /= t;
		Check();
		return m_t;
	}

	
	T operator++(void)		// Prefix (++i)
	{
		// increment before return
		m_t++;
		Check();
		return m_t;
	}

	T operator++(int)		// Postfix (i++)
	{
		// increment after return
		T t = m_t;
		m_t++;
		Check();
		return t;
	}

	T operator--(void)		// Prefix (--i)
	{
		// decrement before return
		m_t--;
		Check();
		return m_t;
	}

	T operator--(int)		// Postfix (i--)
	{
		// decrement after return
		T t = m_t;
		m_t--;
		Check();
		return t;
	}
};
