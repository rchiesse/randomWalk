#include "defs.h"
#pragma once

namespace sim {
	class utils {
	private:
		static void cutNext(const std::string s, const char sep, const size_t level, size_t& wordCount, std::vector<std::string>& ret);

	public:
		static std::vector<std::string> split(const std::string& s, const char sep); // Splits a string into a vector of substrings based on the given separator character.
	};
}
