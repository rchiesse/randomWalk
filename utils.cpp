#include "utils.h"
void sim::utils::cascadeCut(const std::string s, const char sep, const size_t level, size_t& wordCount, std::vector<std::string>& ret) {
	size_t pos = s.find(sep);
	size_t tokenSize = s.end() - s.begin();
	size_t limit = (pos < tokenSize) ? pos : tokenSize;
	if (pos < tokenSize)
		cascadeCut(std::string(s.begin() + pos + 1, s.end()), sep, level + 1, ++wordCount, ret);
	ret[level] = std::string(s.begin(), s.begin() + limit);
}

std::vector<std::string> sim::utils::split(const std::string& s, const char sep) {
	std::vector<std::string> ret(s.size());
	size_t wordCount = 1;
	cascadeCut(s, sep, 0, wordCount, ret);
	ret.resize(wordCount);
	return ret;
}

bool sim::utils::isEmpty(std::ofstream& s) {
	s.seekp(0, std::ios::end);
	return s.tellp() == 0;
}