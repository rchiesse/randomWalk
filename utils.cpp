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

bool sim::utils::containsLine(const std::string& filename, const std::string& content) {
	std::ifstream s;
	s.open(filename);
	if (!s.is_open()) {
		std::cout << "\nMethod \"containsLine():\" Error while trying to open the file " << filename << '\n';
		return false;
	}

	std::string line;
	do {
		getline(s, line);
		if (line == content) {
			s.close();
			return true;
		}
	} while (s.good());
	s.close();
	return false;
}