#pragma once
#include <vector>
#include <string>
#include "Polynom.h"

class Model{
private:
	Polynom polynom; //порождающий многочлен
	int l; // длина кодируемой последовательности
	int r; // степень порождающего полинома
	float epsilon; // точность вычислений
	std::vector<int> message; // вектор сообщения
	std::vector<int> a; // вектор а
	std::vector<int> e; // вектор ошибок
	std::vector<int> b; // вектор b
	Polynom message_polynom; // сообщение в виде полинома
public:
	Model();
	Model(int l, Polynom polynom, float epsilon);
	~Model();

	void message_to_polynom(); // перевод сообщения из вектора в полином
	Polynom calculate_cx(); // вычисление полинома с(x)
	void generate_message(); // получение случайного сообщения длины l
	Polynom calculate_ax(); // вычисление полинома а(х)
	void form_aVector_from_polynom(Polynom& const polynom); // перевод а(х) в вектор а
	void generate_eVector(float p); // получение вектора ошибок. Вероятность ошибки определяется p
	void generate_bVector(); // полуение вектора b
	Polynom calculate_sindrom(); // вычисление синдрома s(x)
	int calculate_vector_weight(const std::vector<int>& v); // вычисление веса вектора (количество единиц), нужно для проверки e = 0 or e != 0

	std::vector<float> modeling(); // реализует само моделирование
	void write_data_in_file(std::string file_name, const std::vector<float>& v); // запись значений Ре в файл
};

