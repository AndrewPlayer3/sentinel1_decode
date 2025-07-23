#include <iostream>
#include <complex>
#include <numeric>

int main()
{
	std::complex<double> a(1.3, 2.7);
	std::cout << a * 1.0 << std::endl;
	std::cout << a * 0.0 << std::endl;
	return 0;
}
