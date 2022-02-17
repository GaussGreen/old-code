// C++ wrapper for GTO bond routines.
#ifndef	_K_BOND_H_
#define	_K_BOND_H_

extern "C" {
#include "bonds.h"
#include "bondcnst.h"
}

#include "kstring.h"

class KBondType
{
public:
	KBondType(void){_bondType = GTO_GENERIC;}
	KBondType(const KString &);
	//to be only used in wrapping GTO routines
	TBondType type()const{return _bondType;}
private:
	TBondType	_bondType;
};

class KBond
{
public:
	KBond(void);
	KBond(const KBond &);
	KBond(const KBondType &s,  
			int        couponFrequency,  
			double     couponRate,       
			KDate      &maturityDate     
			){_tbond = GtoNewTBond(s.type(), couponFrequency,couponRate,maturityDate);
				if(_tbond)	KException("KBond constructor failed!");}

	~KBond(void){destroy();}
	
	KBond& operator=(const KBond&);

	//add two bonds
	KBond	operator +(const KBond &b)const{KBond	ans(*this); ans += b; return ans;}
	void	operator +=(const KBond &);
	
	double	get_bond_price(const KRateCurve	&zr)const;
	double	get_io_price(const KRateCurve	&zr)const{return _couponPayment.get_pv(zr);} 
	double	get_po_price(const KRateCurve	&zr)const{return _prinPayment.get_pv(zr);}
	//compute flat yield
	double	get_bond_yield(const KDate &valueDate, double p, double guess = 0.05);
	double	get_io_yield(const	KDate &valueDate, double p, double guess = 0.05)
			{	return _couponPayment.get_yield(valueDate, p, guess);}
	double	get_po_yield(const KDate &valueDate, double p, double guess = 0.05);
			{	return _prinPayment.get_yield(valueDate, p, guess);}

	friend	ostream & operator<<(ostream & out,  KBond& debt);

private:
	void	destroy(){if(_tbond); GtoFreeTBond(_tbond);}
	TBond	*_tbond;
};


#endif	