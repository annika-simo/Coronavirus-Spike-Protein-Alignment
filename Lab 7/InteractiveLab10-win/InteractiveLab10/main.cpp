#include "hello.h"
#include <iostream>

void wave_a_bunch(const int &num_times_to_wave, const Hello& hello) {
	for (int i = 0; i < num_times_to_wave; i++) {
		hello.wave();
	}
}

bool check_for_tongue_twister(const Hello hello) {
	return hello.is_tongue_twister();
}

int main() {
	std::cout << "About to start main().\n";
	Hello say_hello_from_annika_to_mom("Annika","Mom");
	Hello say_hello_from_dad_to_brother("Super 7", "Taavo");
	std::cout << "About to wave a bunch.\n";
	wave_a_bunch(5, say_hello_from_annika_to_mom);
	
	if (check_for_tongue_twister(say_hello_from_dad_to_brother)) {
		std::cout << "Yes, it is a tongue twister!\n";
	}

	std::cout << "About to return.\n";
	return 0;

}