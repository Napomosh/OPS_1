#pragma once
#include <vector>
#include <string>
#include "Polynom.h"

class Model{
private:
	Polynom polynom; //����������� ���������
	int l; // ����� ���������� ������������������
	int r; // ������� ������������ ��������
	float epsilon; // �������� ����������
	std::vector<int> message; // ������ ���������
	std::vector<int> a; // ������ �
	std::vector<int> e; // ������ ������
	std::vector<int> b; // ������ b
	Polynom message_polynom; // ��������� � ���� ��������
public:
	Model();
	Model(int l, Polynom polynom, float epsilon);
	~Model();

	void message_to_polynom(); // ������� ��������� �� ������� � �������
	Polynom calculate_cx(); // ���������� �������� �(x)
	void generate_message(); // ��������� ���������� ��������� ����� l
	Polynom calculate_ax(); // ���������� �������� �(�)
	void form_aVector_from_polynom(Polynom& const polynom); // ������� �(�) � ������ �
	void generate_eVector(float p); // ��������� ������� ������. ����������� ������ ������������ p
	void generate_bVector(); // �������� ������� b
	Polynom calculate_sindrom(); // ���������� �������� s(x)
	int calculate_vector_weight(const std::vector<int>& v); // ���������� ���� ������� (���������� ������), ����� ��� �������� e = 0 or e != 0

	std::vector<float> modeling(); // ��������� ���� �������������
	void write_data_in_file(std::string file_name, const std::vector<float>& v); // ������ �������� �� � ����
};

