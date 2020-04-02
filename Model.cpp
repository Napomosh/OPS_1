#include "Model.h"
#include <random>
#include <cmath>
#include <time.h>
#include <iostream>
#include <fstream>

Model::Model(){
	l = 0;
	r = 0;
	polynom = Polynom();
}
Model::Model(int l, Polynom polynom, float epsilon){
	this->l = l;
	this->polynom = Polynom(polynom);
	this->epsilon = epsilon;
	r = polynom.get_deg();
	message.resize(l);
}

void Model::generate_message(){
	for(int i = 0; i < message.size(); i++){
		message[i] = rand() % 2;
	}
	if(calculate_vector_weight(message) == 0){
		message[0] = 1;
	}
}

Model::~Model(){}

Polynom Model::calculate_cx(){
	Polynom rx(r);
	rx.coefs[r] = 1;
	Polynom cx = message_polynom * rx;
	cx = cx % polynom;
	cx.set_deg(cx.get_max_deg());
	return cx;
}

void Model::message_to_polynom(){
	message_polynom = Polynom(message.size() - 1, message);
	message_polynom.set_deg(message_polynom.get_max_deg());
}

Polynom Model::calculate_ax(){
	Polynom rx(r);
	rx.coefs[r] = 1;
	Polynom ax = message_polynom * rx + calculate_cx();
	ax.set_deg(ax.get_max_deg());
	return ax;
}

void Model::form_aVector_from_polynom(Polynom& const polynom){
	a.resize(polynom.get_deg() + 1);
	for(int i = 0; i <= polynom.get_deg(); i++){
		a[i] = polynom.get_coefs()[i];
	}
}

void Model::generate_eVector(float p){
	e.resize(a.size());

	for(int i = 0; i < e.size(); i++){
		float x = (rand() % 100) / 100.0;
		//std::cout << x << ' ';
		if(x < p){
			e[i] = 1;
		}
		else{
			e[i] = 0;
		}
	}
}

void Model::generate_bVector(){
	b.resize(a.size());
	for(int i = 0; i < a.size(); i++){
		b[i] = a[i] ^ e[i];
	}
}

Polynom Model::calculate_sindrom(){
	Polynom b(this->b.size() - 1, b);
	Polynom res = b % polynom;

	return res;
}

int Model::calculate_vector_weight(const std::vector<int>& v){
	int w = 0;
	for(int i = 0; i < v.size(); i++){
		w += v[i];
	}
	return w;
}

std::vector<float> Model::modeling(){
	float N = 9 / (4 * epsilon); // число испытаний
	float Ne = 0; // число ошибок декодера

	std::vector<float> Pe; // вектор вероятностей оишибки декодера при различных вероятностях ошибки
	for(float p = 0.1; p < 1.1; p += 0.1){
		for(int i = 0; i < N; i++){
			generate_message();
			message_to_polynom();
			Polynom ax = calculate_ax();
			form_aVector_from_polynom(ax);
			generate_eVector(p);
			generate_bVector();
			Polynom s = calculate_sindrom();

			// декодер ошибся, если он говорит, что ошибки нет (s(x) == 0), но она есть (e != 0)
			if(s.is_polynom_equals_0() != 0 && calculate_vector_weight(e) != 0){
				Ne++;
			}
			b.clear();
			e.clear();
		}
		Pe.push_back(Ne / N);
		Ne = 0;
	}
	return Pe;
}

void Model::write_data_in_file(std::string file_name, const std::vector<float>& v){
	std::ofstream file;
	file.open(file_name, std::ofstream::binary);

	for(int i = 0; i < v.size(); i++){
		file << v[i] << ' ';
	}
	file.close();
}