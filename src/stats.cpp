#include "stats.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

using namespace std;

namespace vcflib
{
	double median(vector<double>& v)
	{
		size_t n = v.size() / 2;
		nth_element(v.begin(), v.begin() + n, v.end());
		return v[n];
	}

	double mean(const vector<double>& v)
	{
		double sum = accumulate(v.cbegin(), v.cend(), 0.0);
		return sum / v.size();
	}


	double mean(const vector<int>& v)
	{
		double sum = 0;

		for (vector<int>::const_iterator it = v.cbegin(); it != v.cend(); ++it) {
			sum += (*it);
		}
		return sum / v.size();
	}


	double variance(const vector <double>& data, const double mean) {
		double total = 0;
		for (vector <double>::const_iterator i = data.cbegin(); i != data.cend(); ++i) {
			total += (*i - mean) * (*i - mean);
		}
		return total / (data.size());
	}

	double standard_deviation(const vector <double>& data, const double mean) {
		return sqrt(variance(data, mean));
	}
}
