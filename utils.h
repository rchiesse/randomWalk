#include "defs.h"
#pragma once

namespace sim {
	class utils {
	public:
		//Splits a string into a vector of substrings based on the given separator character.
		static void split(const std::string s, const char sep, std::vector<std::string>& ret);
		
		//Checks whether an output file is empty. Returns 'true' if either the file exists with no entries in it, or does not exist at all (in which case a new one is created).
		static bool isEmpty(std::ofstream& s);

		//Checks whether the input file contains at least one line matching the string 'content'
		static bool containsLine(const std::string& filename, const std::string& content);
	};
}


