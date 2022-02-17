#ifndef _TREE_IR_H_
#define _TREE_IR_H_

#include "tree.h"
#include "kprocess.h"

class	IRTree:virtual public	BaseTree
{
protected:
	//underlying process
	KIRProcess		_process;
	virtual	void	memInit()=0;
	virtual	void	treeSpotVol(void) = 0;
	virtual	void	drift(void) = 0;
	virtual	void	limits(void) = 0;
public:
	IRTree(void){}
	~IRTree(void){}
	virtual void	dev(DTimeSlice &) = 0;
	virtual void	dev(FTimeSlice &) = 0;
};


#endif
