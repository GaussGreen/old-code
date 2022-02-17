
#pragma once


class complex
{
public:
	complex(){x = y = 0.;}
	complex(double x, double y);
	complex(double x);
	~complex(){}

	complex operator+(complex& other);
	complex operator*(complex& other);
	complex operator-(complex& other);
	complex operator/(complex& other);
	complex operator~(void);// conjugate
	complex operator-(void);

public:
	double x;	// real part
	double y;	// imaginary part
};

complex log(complex c);
complex exp(complex c);
complex sqrt(complex c);
double norm(complex c);
double arg(complex c);
complex operator*(const complex c, double f);
complex operator+(const complex c, double f);
complex operator*(double f, const complex c);
complex operator+(double f, const complex c);
complex operator/(const complex c, double f);
complex operator-(const complex c, double f);
complex operator/(double f, const complex c);
complex operator-(double f, const complex c);

