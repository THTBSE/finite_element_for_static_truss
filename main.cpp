#include "FiniteElemSolver.h"
#include <iostream>

int main()
{
	static_solver trussFem;
	trussFem.readTrussData("inputfile.txt");
	trussFem.Solve();
	int exitSwitch = 1;
	std::cout << "Press '2' to print results" << std::endl 
		<< "Press '0' to quit the program" << std::endl;
	std::cin >> exitSwitch;
	while (exitSwitch)
	{
		trussFem.save_results("result.txt");
		std::cout << "success output data" << std::endl 
			<< "Press '0' to quit the program" << std::endl;
		std::cin >> exitSwitch;
	}
	return 0;
}