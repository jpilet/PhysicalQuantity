// A small executable to generate documentation.

#include "PhysicalQuantity.h"

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <fstream>

std::string tolower(const std::string& input) {
	std::string result;
	for (auto c: input) {
		result += tolower(c);
	}
	return result;
}

int main(int argc, char** argv) {
	std::map<std::string, std::set<std::string>> quantities;

#define ADD_QUANTITY(name, siUnit, t, l, a, m, T, I, N, J)\
	quantities[#name].insert(#siUnit);

	FOREACH_QUANTITY(ADD_QUANTITY)
#undef ADD_QUANTITY

#define ADD_UNIT(type, name, factor, offset) \
	quantities[#type].insert(#name);
	FOREACH_UNIT(ADD_UNIT)
#undef ADD_UNIT

	std::ofstream out("UNITS.md");

	out << "# All quantities\n";
	for (auto quantity: quantities) {
		out << " - [" << quantity.first << "](#"
			<< tolower(quantity.first) << ")\n";
	}

	out << "\n\n# All Units\n";
	for (auto quantity: quantities) {
		out << "\n## " << quantity.first << "\n";
		for (auto unit: quantity.second) {
			out << " - " << unit << "\n";
		}
	}
	return 0;
}
