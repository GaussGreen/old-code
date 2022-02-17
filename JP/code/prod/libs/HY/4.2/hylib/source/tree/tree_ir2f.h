#ifndef	_TREE_IR_2F_H_
#define	_TREE_IR_2F_H_

#include "kplatform.h"  
#include "tree_ir.h"
#include <iostream>


class IR2FTreeParam
{
public:
	KVector(int)	_idx1;
	KVector(int)	_idx2;
	KVector(double)	_prob[3];
	double			_idxDisc;

	IR2FTreeParam():_idx1(3), _idx2(3){for(int i=0; i<3; i++) _prob[i].resize(3);}
	~IR2FTreeParam(){}
	friend	std::ostream &operator<<(std::ostream &out, const IR2FTreeParam&param)
	{
		out<<"Coordinate1:"<<"\t"<<param._idx1<<"\t";
		out<<"Coordinate2:"<<"\t"<<param._idx2<<"\t";
		out<<"Probabilities:"<<"\t"<<param._prob[0]<<"\t";
		out<<"Probabilities:"<<"\t"<<param._prob[1]<<"\t";
		out<<"Probabilities:"<<"\t"<<param._prob[2]<<"\t";
		out<<"IndexDiscount:"<<"\t"<<param._idxDisc<<"\t";
		return out;
	}


};


class IRTree2F:public	IRTree
{
   	TimeSlice<IR2FTreeParam>	
			_irParam;
protected:
	virtual	void	memInit(){_irParam.resize(this);}
	virtual	void	drift(void);
	virtual	void	treeSpotVol(void);
	virtual	void	limits (void);
	virtual	void	buildLattice(void);
public:
	IRTree2F(void){}
	~IRTree2F(void){}
	virtual	KValarray<int>	get_curr_limit(void) const;
	virtual	KValarray<int>	get_next_limit(void) const;
	virtual	KValarray<int>	get_max_limit(void) const;
	virtual	void	dev(DTimeSlice &);
	virtual	void	dev(FTimeSlice &);
	void	print(std::ostream &out)const{out<<_irParam;}
};

#endif




