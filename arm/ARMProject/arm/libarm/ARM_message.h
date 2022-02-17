#ifndef ARM_MESSAGE_H
#define ARM_MESSAGE_H

class ARM_message
{
public:
	ARM_message (const char* n_msg) : msg (n_msg) {}
	ARM_message (const CCString& n_msg) : msg (n_msg) {}
	~ARM_message () {}

	inline void setMsg (const CCString& n_msg)
	{
		msg = n_msg;
	}

	inline void setMsg (const char* n_msg)
	{
		msg.Set (n_msg);
	}

	inline CCString getMsg ()
	{
		return msg;
	}
		
private:
	CCString msg;
};

#endif	// ARM_MESSAGE_H