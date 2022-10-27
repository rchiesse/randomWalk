#include "utils.h"
void sim::utils::split(const std::string s, const char sep, std::vector<std::string>& ret) {
	size_t pos = std::min(s.find(sep), s.size());
	ret.emplace_back(s.substr(0, pos));
	if (pos < s.size())
		split(s.substr(pos + 1), sep, ret);
}

bool sim::utils::isEmpty(std::ofstream& s) {
	s.seekp(0, std::ios::end);
	return s.tellp() == 0;
}