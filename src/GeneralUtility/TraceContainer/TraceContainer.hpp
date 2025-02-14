/* 
 * MIT License
 * 
 * Copyright (c) 2022 Francesco Da Dalt
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
 * Check the references for the algorithm publication and eventual copyrights that apply to it.
 */

#ifndef TraceContainer_hpp
#define TraceContainer_hpp

#include <stdio.h>
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <unordered_set>
#include <map>
#include <cmath>
#include <random>
#include <algorithm>
#include <ranges>
#include "MurmurHash3.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>


/*
 A simple container class for a trace.
 Computes the list of distinct ids and the ground truth of the trace.
 */
class TraceContainer {
private:
	std::map<int, double> ground_truth;
	std::tuple<std::vector<int>, std::vector<double>> ground_truth_vectorized;
	std::string tracetype;
	std::string trace;
	
	std::vector<std::tuple<std::string, std::string>> stats;
public:
	/*
	 Read in whole trace from path
	 */
	TraceContainer(std::string path);
	
	/*
	 Read in a trace and randomly subsample keys from it
	 */
	TraceContainer(std::string path, int num_samples, int seed);
	
	/*
	 Define a trace with a number of keys. Must be followed by a synthetization step
	 */
	TraceContainer(int num_ids);
	
	/*
	 Synthesizes an exponential trace.
	 Items are i.i.d.
	 */
	void synthesizeExponential(double scale, double offset, int seed);
	
	/*
	 Synthesizes a Pareto trace.
	 Items are i.i.d.
	 */
	void synthesizePareto(double shape, double scale, double offset, int seed);
	
	/*
	 Synthesizes a Normal distributed trace of teh kind trunc(N(scale, scale^2) + low, [low, high]) if scale > 0.
	 Items are i.i.d.
	 */
	void synthesizeNormal(double scale, double low, double high, int seed);
	
	const std::map<int, double>* const get_ground_truth_mapping();
	const std::tuple<std::vector<int>, std::vector<double>>* const get_ground_truth_vector();
	
	std::string get_tracetype() {return tracetype;};
	std::string get_trace() {return trace;};
	
	/*
	 Overwrites item-ids / keys with the sequence of natural numbers 0, 1, 2, ...
	 */
	void override_item_ids();
	
	/*
	 Get information about the trace
	 */
	auto getStats() {
		std::cout << "\n-> Trace" << std::endl;
		std::cout << "--> Type: " << tracetype << std::endl;
		std::cout << "--> Name: " << trace << std::endl;
		std::cout << "--> Number of Items: " << ground_truth.size() << std::endl;
		stats.clear();
		stats.push_back(std::make_tuple("TraceType", tracetype));
		stats.push_back(std::make_tuple("TraceName", trace));
		stats.push_back(std::make_tuple("NumItems", std::to_string(ground_truth.size())));
		return &stats;
	}
	
	void writeGroundTruth (std::string namee) {
		std::ofstream dump_file("Data/Results/" + namee + "_" + tracetype + ".csv", std::ios::out | std::ios::trunc);
		int num_cols = 2;
		for (const auto& pair : ground_truth) {
			dump_file << pair.first << ", " << pair.second << std::endl;
		}
		dump_file.close();
	}
	
	/*
	 Scales the frequencies to have mean 100, for numercial stability and comparability between datasets
	 */
	double normalize();
};

#endif /* TraceContainer_hpp */
