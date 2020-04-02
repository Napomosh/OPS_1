#include <iostream>
#include "Polynom.h"
#include "Model.h"
#include <time.h>

int main(){
	srand(time(0));

	std::cout << "Please, enter coeficients of g(x) polynom in little-endian order. For example, g(x) = x^3 + x + 1, you should enter 1101. Press Enter to confirm." << std::endl;
	std::string str;
	std::cin >> str;

	std::cout << "Please enter the length of sequence: " << std::endl;
	unsigned int l = 0;
	std::cin >> l;

	std::cout << "Please enter the epsilon. For example, 0.001." << std::endl;
	float epsilon = 0.0;
	std::cin >> epsilon;

	Polynom gx = Polynom(str);

	Model model = Model(l, gx, epsilon);
	std::vector<float> Pe = model.modeling();

	model.write_data_in_file("C:\\Users\\sashu\\Documents\\MATLAB\\file.txt", Pe);


	/*float epsilon = 0.0001;

	std::string str = "1101";
	Polynom gx = Polynom(str);

	Model model = Model(4, gx, epsilon);
	std::vector<float> Pe = model.modeling();
	model.write_data_in_file("C:\\Users\\sashu\\Documents\\MATLAB\\file4.txt", Pe);
	std::cout << "Modeling 4 finished" << std::endl << std::endl;

	Model model1 = Model(6, gx, epsilon);
	std::vector<float> Pe1 = model1.modeling();
	model.write_data_in_file("C:\\Users\\sashu\\Documents\\MATLAB\\file6.txt", Pe);
	std::cout << "Modeling 6 finished" << std::endl << std::endl;

	Model model2 = Model(12, gx, epsilon);
	std::vector<float> Pe2 = model2.modeling();
	model.write_data_in_file("C:\\Users\\sashu\\Documents\\MATLAB\\file12.txt", Pe);
	std::cout << "Modeling 12 finished" << std::endl << std::endl;

	Model model3 = Model(16, gx, epsilon);
	std::vector<float> Pe3 = model3.modeling();
	model.write_data_in_file("C:\\Users\\sashu\\Documents\\MATLAB\\file16.txt", Pe);
	std::cout << "Modeling 16 finished" << std::endl;*/


	return 0;
}