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

	double mean(const vector<double>& data)
	{
		double sum = accumulate(data.cbegin(), data.cend(), 0.0);
		return sum / data.size();
	}


	double mean(const vector<int>& data)
	{
		double sum = 0;

		for (const auto d : data) {
			sum += d;
		}
		return sum / data.size();
	}


	double variance(const vector <double>& data, const double mean) {
		double total = 0;
		for (const auto x : data) {
			total += (x - mean) * (x - mean);
		}
		return total / (data.size());
	}

	double sample_variance(const std::vector<double> &data, double mean) {
		double variance = 0;

		for(const auto& d : data){
			variance += pow(d - mean,2);
		}

		return variance / (data.size() - 1);
	}

	double standard_deviation(const vector <double>& data, const double mean) {
		return sqrt(variance(data, mean));
	}

	double sample_standard_deviation(const std::vector<double> &data, double mean) {
		return sqrt(sample_variance(data, mean));
	}
}
