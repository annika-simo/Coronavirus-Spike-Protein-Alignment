#ifndef HELLO_H
#define HELLO_H
#include <string>

class Hello {
public:
	Hello(std::string from, std::string to);
	void wave() const;
	bool is_tongue_twister() const;
	~Hello();
private:
	std::string m_hello_from;
	std::string m_hello_to;
	int* m_ptr_to_int;
};


#endif
