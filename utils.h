#include "defs.h"
#pragma once

namespace sim {
	class utils {
	public:
		static void split(const std::string s, const char sep, std::vector<std::string>& ret); // Splits a string into a vector of substrings based on the given separator character.
		static bool isEmpty(std::ofstream& s);
	};
}
