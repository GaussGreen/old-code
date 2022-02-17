
#ifndef	_TREE_IR_1F_H_
#define	_TREE_IR_1F_H_
#include "tree_ir.h" 



class IR1FTreeParam
{
public:
	KVector(int)	_idx;
	KVector(double)	_prob;
	double	_idxDisc;

	IR1FTreeParam():_idx(3), _prob(3){}
	~IR1FTreeParam(){}
	friend	std::ostream &operator<<(std::ostream &out, const IR1FTreeParam&param)
	{
		out<<"Coordinates:"<<"\t"<<param._idx<<"\t";
		out<<"Probabilities:"<<"\t"<<param._prob<<"\t";
		out<<"Index Discount:"<<"\t"<<param._idxDisc<<"\t";
		return out;
	}
};

class IRTree1F:public	IRTree
{	 
   	TimeSlice<IR1FTreeParam>	
			_irParam;
protected:
	virtual	void	memInit(){_irParam.resize(this);}
	virtual	void	drift(void);
	virtual	void	limits(void);
	virtual	void	treeSpotVol(void);
	virtual	void	buildLattice(void);
public:
	IRTree1F(void){}
	~IRTree1F(void){}
	virtual	KValarray<int>	get_curr_limit(void) const;
	virtual	KValarray<int>	get_next_limit(void) const;
	virtual	KValarray<int>	get_max_limit(void) const;
	virtual	void	dev(DTimeSlice &);
	virtual	void	dev(FTimeSlice &);
	void	print(std::ostream &out)const{out<<_irParam;}
};


#endif




