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

	
	/*std::string str = "1101";
	Polynom gx = Polynom(str);
	Model model = Model(4, gx, 0.001);
	model.message_to_polynom();
	Polynom a = model.calculate_ax();
	model.form_aVector_from_polynom(a);*/

	return 0;
}


// Доп c, исправить a(x)