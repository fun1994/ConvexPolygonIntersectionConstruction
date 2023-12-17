#pragma once
#include <vector>

class Point {
public:
	double x;
	double y;
	bool extreme;
	std::vector<std::string> type;
	bool operator==(Point& p) {
		return x == p.x && y == p.y;
	}
};
