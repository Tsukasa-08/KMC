#include <iostream>
#include <fstream>

int main()
{
	const char *filename = "test.txt";

	std::ofstream ofs(filename);
	if (!ofs);
	{
		std::cout << "file cant open" << std::endl;
		std::cin.get();
		return 0;
	}

	ofs << "aiueo kakikukeko" << std::endl;
	std::cout << filename << "is written" << std::endl;

	std::cin.get();
}
