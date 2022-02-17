#ifndef _TREE_EQ_H_
#define _TREE_EQ_H_

#include "kplatform.h"   

#include <map>
#include "kmatrix.h"
#include "tree.h"
#include "kprocess.h"
#include "kplatdep.h"



class EqTreeParam
{
public:
	KVector(int)		_idx;
	KVector(double)	_prob;
	double			_under;

	EqTreeParam():_idx(3), _prob(3){}
	~EqTreeParam(){}
	friend	std::ostream &operator<<(std::ostream &out, const EqTreeParam&param)
	{
		out<<"Coordinates:"<<"\t"<<param._idx<<"\t";
		out<<"Probabilities:"<<"\t"<<param._prob<<"\t";
		out<<"Underlying:"<<"\t"<<param._under<<"\t";
		return out;
	}
};


class EqTree:virtual public	BaseTree
{
	TimeSlice<EqTreeParam>	_eqParam;
protected:
	KEqProcess	_process;
	virtual void	memInit(void){_eqParam.resize(this);}
	virtual	void	limits (void);
	virtual	void	treeSpotVol(void);
	virtual	void	drift(void);
	virtual void	buildLattice();
public:
	EqTree(void){}
	~EqTree(void){}
	//set stock and vol information
	void	set_process(const KEqProcess &p);
	virtual KValarray<int>	get_curr_limit(void) const;
	virtual KValarray<int>	get_next_limit(void) const;
	virtual KValarray<int>	get_max_limit(void) const;
	virtual void	dev(DTimeSlice &);
	virtual void	dev(FTimeSlice &);
	virtual	DTimeSlice	get_stock_price(void)const;
	void	print(std::ostream &out)const{out<<_eqParam;}
};




#endif

