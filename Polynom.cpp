#include "Polynom.h"
#include <vector>
#include <iostream>

Polynom::Polynom(){
	deg = 0;
	coefs.resize(1);
	coefs.push_back(0);
}
Polynom::Polynom(unsigned int deg, std::vector<int> coefs){
	if(deg > 0){
		this->deg = deg;
		this->coefs.resize(deg + 1);
		this->coefs = coefs;
		deg = get_max_deg();
	}
	else{
		this->deg = 0;
		this->coefs.resize(1);
		this->coefs[0] = 0;
	}
}
Polynom::Polynom(unsigned int deg){
	this->deg = deg;
	//std::cout << deg << ' ';
	this->coefs.resize(deg + 1);
}
Polynom::Polynom(const Polynom& pol){
	this->deg = pol.deg;
	this->coefs = pol.coefs;
}
Polynom::Polynom(const std::string& str){
	deg = str.size() - 1;
	coefs.resize(str.size());
	for(int i = 0; i < str.size(); i++){
		if(str[i] == '0'){
			coefs[i] = 0;
		}

		// ���� ������ �� ��������, �� ����� ������ ����� "0" ��������� "1"
		else{
			coefs[i] = 1;
		}
	}
}

Polynom::~Polynom(){}

unsigned int Polynom::get_deg(){
	return deg;
}
std::vector<int> Polynom::get_coefs(){
	return coefs;
}

void Polynom::set_deg(unsigned int d){
	deg = d;
}

bool Polynom::is_polynom_equals_0(){
	return deg == 0 && this->coefs[0] == 0;
}

const Polynom Polynom::operator%(const Polynom& rhs) const{
	Polynom res(this->deg - rhs.deg);
	Polynom remain(*this);
	Polynom tmp_for_sub(this->deg);

	for(int i = this->deg; i >= 0; i--){
		//��������� ��������� ������� ���������� ��������, ��� �������� ������������ ��������  �������� ������� � ��������
		unsigned int leftdeg_sub_rightdeg = remain.deg - rhs.deg;
		res.coefs[leftdeg_sub_rightdeg] = 1;

		for(int j = rhs.deg; j >= 0; j--){
			// ������� ���������, ������� ���� ������� �� �������. ������ ������������ ��, ��� �������� �� ���������� ���� � ��������
			tmp_for_sub.coefs[leftdeg_sub_rightdeg + j] = rhs.coefs[j] & res.coefs[leftdeg_sub_rightdeg];
		}

		//��������� ����� �������
		remain = remain - tmp_for_sub;
		remain.deg = remain.get_max_deg();
		tmp_for_sub.coefs.clear();
		tmp_for_sub.coefs.resize(this->deg + 1);
		
		// ���� ������� ������ ��������, �� ������� ���������
		if(remain < rhs){
			break;
		}
	}

	return remain;
}

const bool Polynom::operator<(const Polynom& rhs) const{
	return this->deg < rhs.deg;
}

const Polynom Polynom::operator-(const Polynom& rhs) const{

	unsigned int min_deg;
	unsigned int max_deg;
	if(this->deg > rhs.deg){
		min_deg = rhs.deg;
		max_deg = this->deg;
	}
	else{
		min_deg = this->deg;
		max_deg = rhs.deg;
	}
	Polynom res(max_deg);

	for(int i = min_deg + 1; i <= max_deg; i++){
		if(this->deg > rhs.deg){
			res.coefs[i] = this->coefs[i];
		}
		else{
			res.coefs[i] = rhs.coefs[i];
		}
		
	}
	for(int i = 0; i <= min_deg; i++){
		res.coefs[i] = this->coefs[i] ^ rhs.coefs[i];
	}
	res.coefs.shrink_to_fit();
	return res;
}

unsigned int Polynom::get_max_deg(){
	for(int i = this->coefs.size() - 1; i >= 0; i--){
		if(this->coefs[i] == 1){
			return i;
		}
	}
	return 0;
}

// � ������ ������ ������ ���������� �������� ��������� �� ��������, ���� �������� ��������� ������ ���. ��������� ��������� �� ��������� �� �����������
const Polynom Polynom::operator*(const Polynom& rhs) const{ 
	unsigned int max_deg = this->deg + rhs.deg;
	Polynom res(this->deg + rhs.deg);
	res.coefs[res.coefs.size() - 1] = 0;
	Polynom max_pol;
	Polynom min_pol;

	// ���������� � ����� ������� ������� ������
	if(this->deg > rhs.deg){
		for(int i = max_deg; i >= 0; i--){

			if(i < rhs.deg){
				res.coefs[i] = this->coefs[i] && rhs.coefs[rhs.deg];
			}
			// ����� ��� �������� � ������� ��� �� ������ ��������� ���������� �� �������
			else{
				res.coefs[i] = this->coefs[i - rhs.deg] && rhs.coefs[rhs.deg];
			}
		}
		//for(int i = 0; i < max_deg - rhs.deg - 1; i++){

		//	// ��� ����������� ��� ���������, ��� �� ���������� ������ ���� ����� 0
		//	res.coefs[i] = 0;
		//}
	}
	else{
		for(int i = max_deg; i >= max_deg - this->deg; i--){
			res.coefs[i] = this->coefs[i - rhs.deg] && rhs.coefs[rhs.deg];
		}
		for(int i = 0; i < max_deg - this->deg - 1; i++){
			res.coefs[i] = 0;
		}
	}

	return res;
}

const Polynom Polynom::operator+(const Polynom& rhs) const{
	unsigned int max_deg;
	unsigned int min_deg;
	if(this->deg > rhs.deg){
		max_deg = this->deg;
		min_deg = rhs.deg;
	}
	else{
		max_deg = rhs.deg;
		min_deg = this->deg;
	}

	Polynom res = Polynom(max_deg);
	res.coefs[res.deg] = 0;
	for(int i = 0; i <= min_deg; i++){
		res.coefs[i] = this->coefs[i] ^ rhs.coefs[i];
	}
	for(int i = min_deg + 1; i <= max_deg; i++){
		if(this->deg > rhs.deg){
			res.coefs[i] = this->coefs[i];
		}
		else{
			res.coefs[i] = rhs.coefs[i];
		}
	}
	return res;
}

