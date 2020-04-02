#pragma once
#include <vector>
#include <string>

class Polynom{
private:
	unsigned int deg;
public:
	std::vector<int> coefs;
	Polynom();
	Polynom(unsigned int deg, std::vector<int> coefs);
	Polynom(unsigned int deg);
	Polynom(const Polynom& pol);
	Polynom(const std::string& str);

	~Polynom();

	unsigned int get_max_deg(); // получение максимальной степени полинома
	unsigned int get_deg();
	std::vector<int> get_coefs();
	void set_deg(unsigned int d);
	bool is_polynom_equals_0(); 

	const Polynom operator%(const Polynom& rhs) const; // деление полиномов по модулю
	const bool operator<(const Polynom& rhs) const;
	const Polynom operator-(const Polynom& rhs) const;
	const Polynom operator*(const Polynom& rhs) const;
	const Polynom operator+(const Polynom& rhs) const;
};