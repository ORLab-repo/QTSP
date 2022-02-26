#pragma once
#include "lib.h"

using namespace std;

namespace Util {
	inline void removeSpaceAround(string& s) {
		const char* t = " \t\n\r\f\v";
		s.erase(0, s.find_first_not_of(t));// remove space at the end of string
		s.erase(s.find_last_not_of(t) + 1);// remover space at the end of string
	}
	// Split a string
	inline vector<string> splitString(string& s, string sep) {
		vector<string> res;
		//removeSpaceAround(s);
		size_t current, previous = 0;
		current = s.find(sep);
		while (current != std::string::npos) {
			res.push_back(s.substr(previous, current - previous));
			previous = current + 1;
			current = s.find(sep, previous);
		}
		res.push_back(s.substr(previous, current - previous));
		for (int i = 0; i < res.size(); ++i)removeSpaceAround(res[i]);
		if (res.back() == "")res.pop_back();
		return res;
	}
	//Convert a string into number
	inline int convertStringToNum(string& s) {
		int res = 0;
		for (int i = 0; i < s.size(); ++i) {
			assert(s[i] >= '0' && s[i] <= '9');
			res = res * 10 + s[i] - '0';
		}
		return res;
	}

	inline double round2num(double var)
	{
		int value = (int)(var * 100);
		return (double)value / 100;
	}

	//hash function with 2 mod
	inline pair<int, int> getHash(int* arr, int length) {
		int* arrSort = new int[length];
		for (int i = 0; i < length; ++i)arrSort[i] = arr[i + 1];
		int mod[] = { (int)1e9 + 2277, (int)1e9 + 5277 };
		int base = 311; //should bigger than num of customers
		sort(arrSort, arrSort + length);
		pair<int, int> res = make_pair(0, 0);
		for (int i = 0; i < length; ++i) {
			res.first = ((long long)res.first * base + arrSort[i]) % mod[0];
			res.second = ((long long)res.second * base + arrSort[i]) % mod[1];
		}
		delete[] arrSort;
		return res;
	}
}
