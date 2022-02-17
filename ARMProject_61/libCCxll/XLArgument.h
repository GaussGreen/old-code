#ifndef XLARGUMENT_H
#define XLARGUMENT_H

#include <CCString.h>

class XLArgument
{
public:
	XLArgument () {};
	XLArgument (const char* n_name, const char* n_helptext, int n_type);
	~XLArgument () {};

	inline void setName (const char* n_name)
	{
		name.Set (n_name);
	}
	inline CCString getName () const
	{
		return name;
	}

	inline void setHelptext (const char* n_helptext)
	{
		helptext.Set (n_helptext);
	}
	inline CCString getHelptext () const
	{
		return helptext;
	}

	inline void setType (int n_type)
	{
		type = n_type;
	}
	inline int getType () const
	{
		return type;
	}

private:
	CCString name;
	CCString helptext;
	int type;
};

#ifdef XLArgument_cpp

#endif	// XLArgument_cpp

#endif	// XLARGUMENT_H
