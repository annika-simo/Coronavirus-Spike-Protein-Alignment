#include "hello.h"
#include <string>
#include <iostream>

Hello::Hello(std::string from, std::string to) : m_hello_from(from), m_hello_to(to) {

	m_ptr_to_int = new int{ 5 };

	std::cout << "I am a " << m_hello_from << ". I'm a friendly greeter!\n";

}

bool Hello::is_tongue_twister() const {
	return (m_hello_from.length() + m_hello_to.length()) > 10;
}

void Hello::wave() const {
	std::cout << "Waving from " << m_hello_from << "\n";
}

Hello::~Hello() {
	std::cout << "Destroy! " << m_hello_from << " is saying goodbye!\n";
	delete m_ptr_to_int;
}
