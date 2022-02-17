
#include "stdafx.h"
#include "complex.h"


complex::complex(double x_, double y_)
{
	x = x_;
	y = y_;
}

complex::complex(double x_)
{
	x = x_;
	y = 0.;
}

complex complex::operator +( complex& c )
{
	return complex( x + c.x, y + c.y);
}

complex complex::operator -( complex& c )
{
	return complex( x - c.x, y - c.y);
}


complex complex::operator *( complex& c )
{
	return complex( x*c.x - y*c.y, x*c.y + y*c.x);
}


complex complex::operator /( complex& c )	
{
	double normm = norm(c);
	normm *= normm;
	return complex( (x*c.x + y*c.y)/normm, (-x *c.y + y*c.x)/normm);
}

complex complex::operator ~(void)
{
	return complex(x,-y);
}

complex operator*(const complex c, double f)
{
	return complex(c.x*f, c.y*f);
}


complex operator+(const complex c, double f)
{
	return complex(c.x+f, c.y);
}

complex operator*(double f, const complex c)
{
	return complex(c.x*f, c.y*f);
}


complex operator+(double f, const complex c)
{
	return complex(c.x+f, c.y);
}


complex operator/(const complex c, double f)
{
	return complex(c.x/f, c.y/f);
}


complex operator-(const complex c, double f)
{
	return complex(c.x-f, c.y);
}

complex operator/(double f, const complex c)
{
	double normm = norm(c);
	double argg	= arg(c);
	return complex( f/normm*cos(argg), -f/normm*sin(argg) );
}


complex operator-(double f, const complex c)
{
	return complex(f-c.x, -c.y);
}

complex complex::operator -(void)
{
	return complex(-x,-y);
}

double norm(complex c)
{
	double x = c.x;
	double y = c.y;
	return sqrt(x*x + y*y);
}

double arg(complex c)
{
	return atan2(c.y, c.x);	
}


complex log(complex c)	
{
	double x = log(norm(c));
	double y = arg(c);

	return complex(x, y);
}

complex exp(complex c)
{
	double normm = exp(c.x);
	double argg = c.y;

	return complex(normm*cos(argg), normm*sin(argg));
}

complex sqrt(complex c)	
{
	double rho = sqrt(norm(c));
	double theta = 0.5*arg(c);

	return complex(rho*cos(theta), rho*sin(theta));
}


