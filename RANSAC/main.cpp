#include <iostream>
#include "Ransac.h"


int main()
{
	std::vector<Point2<double>> data =
	{
		{ 0.0, 0.0 },
		{ 9.0, 2.0 },
		{ 4.0, 2.0 },
		{ 3.0, 3.0 },
		{ 4.0, 4.0 },
		{ 5.0, 5.0 },
		{ 1.0, 2.0 },
		{ 6.0, 6.0 },
		{ 7.2, 7.2 },
		{ 9.0, 8.5 },
		{ 2.5, 3.0 },
		{ 0.7, 2.0 },
	};

	Ransac2D<double> ransac;
	auto r = ransac.Calc(data, 1000, 0.5, 2);

	std::cout << "y = " << r[0] << "x + " << r[1] << std::endl;

	return 0;
}