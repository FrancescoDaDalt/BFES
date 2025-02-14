/* 
 * MIT License
 * 
 * Copyright (c) 2023 Francesco Da Dalt
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

#ifndef PRSketch_hpp
#define PRSketch_hpp

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <vector>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "MurmurHash3.h"

// CHECKED, all good
template <typename key_type, typename support_type>
class PRSketch {
	const int depth = 1;
	const int width_freq;
	const int width_filter;
	std::vector<support_type> counters_freq;
	std::vector<bool> counters_filter;
	std::vector<key_type> identified_keys;
	static constexpr auto basicHash = std::hash<key_type>{};
	std::vector<std::tuple<std::string, double>> memstats;
public:
	PRSketch (const int total_counters):
	width_freq(total_counters * (1.0 - 0.125)),
	width_filter((total_counters - width_freq) * 32) { // We allow more "memory" because a bool vector is 32 times more efficient than a vector of ints or floats.
		counters_freq.resize(depth * width_freq);
		counters_filter.resize(depth * width_filter);
		std::fill(counters_freq.begin(), counters_freq.end(), 0);
		std::fill(counters_filter.begin(), counters_filter.end(), false);
	}
	
	void insert (const std::tuple<key_type, support_type> key_value) {
		const auto key = std::get<0>(key_value);
		const auto value = std::get<1>(key_value);
		const size_t firsthash = basicHash(key);
		for (int i = 0; i < depth; i++) {
			size_t hashedkey;
			size_t extradest[2];
			MurmurHash3_x64_128(&firsthash,
								sizeof(size_t),
								7218563 * i + 1543724,
								&extradest);
			hashedkey = extradest[0] % width_freq;
			counters_freq.at(i * width_freq + hashedkey) += value;
		}
		for (int i = 0; i < depth; i++) {
			size_t hashedkey;
			size_t extradest[2];
			MurmurHash3_x64_128(&firsthash,
								sizeof(size_t),
								952249 * i + 1008924,
								&extradest);
			hashedkey = extradest[0] % width_filter;
			if (!counters_filter.at(i * width_filter + hashedkey)) {
				identified_keys.push_back(key);
				counters_filter.at(i * width_filter + hashedkey) = true;
			}
		}
	}
	
	auto batch_query (const std::vector<key_type>* const keys) {
		const int query_size = keys->size();
		auto* const index = new std::vector<size_t>(query_size);
		auto* const estimates = new std::vector<std::vector<double>>(query_size, std::vector<double>(1));
		
		std::vector<size_t> key_indexer(depth * width_freq, -1);
		
		
		int identified_keys_size = identified_keys.size();
		std::vector<Eigen::Triplet<double>> sensingmatrix_tripletlist(identified_keys_size * depth);
		for (int j = 0; j < identified_keys_size; j++) {
			const auto key = identified_keys.at(j);
			const size_t firsthash = basicHash(key);
			for (int i = 0; i < depth; i++) {
				size_t hashedkey;
				size_t extradest[2];
				MurmurHash3_x64_128(&firsthash,
									sizeof(size_t),
									7218563 * i + 1543724,
									&extradest);
				hashedkey = extradest[0] % width_freq;
				key_indexer.at(hashedkey) = j;
				sensingmatrix_tripletlist.push_back(Eigen::Triplet<double>((int) hashedkey, j, 1.0));
			}
		}
		
		
		for (int i = 0; i < depth * width_freq; i++) {
			if (key_indexer.at(i) == -1) {
				key_indexer.at(i) = identified_keys_size;
				sensingmatrix_tripletlist.push_back(Eigen::Triplet<double>((int) i, identified_keys_size, 1.0));
				identified_keys_size++;
			}
		}
		
		
		Eigen::SparseMatrix<double> sensingmatrix(depth * width_freq, identified_keys_size);
		sensingmatrix.setFromTriplets(sensingmatrix_tripletlist.begin(), sensingmatrix_tripletlist.end());
		
		auto eigen_counters = Eigen::Map<Eigen::Vector<support_type, Eigen::Dynamic>>(counters_freq.data(), depth * width_freq);
		
		Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<support_type>> solver;
		solver.compute(sensingmatrix);
		Eigen::Vector<support_type, Eigen::Dynamic> sol = solver.solve(eigen_counters);
		
#pragma omp parallel for
		for (int j = 0; j < query_size; j++) {
			const auto key = keys->at(j);
			const size_t firsthash = basicHash(key);
			double est = std::numeric_limits<double>::max();
			for (int i = 0; i < depth; i++) {
				size_t hashedkey;
				size_t extradest[2];
				MurmurHash3_x64_128(&firsthash,
									sizeof(size_t),
									7218563 * i + 1543724,
									&extradest);
				hashedkey = extradest[0] % width_freq;
				est = sol(key_indexer.at(hashedkey));
			}
			index->at(j) = j;
			estimates->at(j).at(0) = est;
		}
		return std::make_tuple(index, estimates, estimates);
	}
	
	static std::string getDescription() {
		return "PRS";
	}
	
	auto getNummemcellsFast() const {
		return depth * (width_freq + width_filter / 32);
	}
	
	auto getNummemcells() const {
		return depth * (width_freq + width_filter / 32) + identified_keys.size();
	}
	
	auto getMemStats() {
		std::cout << "\n-> " + getDescription() + " Memory" << std::endl;
		std::cout << "--> On-Line Memory: " << getNummemcells() << std::endl;
		std::cout << "--> On-Line Counters: " << getNummemcellsFast() << std::endl;
		memstats.clear();
		memstats.push_back(std::make_tuple("OnLineMemory", getNummemcells()));
		memstats.push_back(std::make_tuple("OnLineCounters", getNummemcellsFast()));
		return &memstats;
	}
};

// CHECKED, all good
template <typename key_type, typename support_type>
class PRSketch_Native {
	const int depth = 1;
	const int width_freq;
	std::vector<support_type> counters_freq;
	static constexpr auto basicHash = std::hash<key_type>{};
	std::vector<std::tuple<std::string, double>> memstats;
public:
	PRSketch_Native (const int total_counters):
	width_freq(total_counters ) {
		counters_freq.resize(depth * width_freq);
		std::fill(counters_freq.begin(), counters_freq.end(), 0);
	}
	
	void insert (const std::tuple<key_type, support_type> key_value) {
		const auto key = std::get<0>(key_value);
		const auto value = std::get<1>(key_value);
		const size_t firsthash = basicHash(key);
		for (int i = 0; i < depth; i++) {
			size_t hashedkey;
			size_t extradest[2];
			MurmurHash3_x64_128(&firsthash,
								sizeof(size_t),
								7218563 * i + 1543724,
								&extradest);
			hashedkey = extradest[0] % width_freq;
			counters_freq.at(i * width_freq + hashedkey) += value;
		}
	}
	
	auto batch_query (const std::vector<key_type>* const keys) {
		const int query_size = keys->size();
		auto* const index = new std::vector<size_t>(query_size);
		auto* const estimates = new std::vector<std::vector<double>>(query_size, std::vector<double>(1));
		
		std::vector<Eigen::Triplet<double>> sensingmatrix_tripletlist(query_size * depth);
		for (int j = 0; j < query_size; j++) {
			const auto key = keys->at(j);
			const size_t firsthash = basicHash(key);
			for (int i = 0; i < depth; i++) {
				size_t hashedkey;
				size_t extradest[2];
				MurmurHash3_x64_128(&firsthash,
									sizeof(size_t),
									7218563 * i + 1543724,
									&extradest);
				hashedkey = extradest[0] % width_freq;
				sensingmatrix_tripletlist.push_back(Eigen::Triplet<double>((int) hashedkey, j, 1.0));
			}
		}
		
		Eigen::SparseMatrix<double> sensingmatrix(depth * width_freq, query_size);
		sensingmatrix.setFromTriplets(sensingmatrix_tripletlist.begin(), sensingmatrix_tripletlist.end());
		
		auto eigen_counters = Eigen::Map<Eigen::Vector<support_type, Eigen::Dynamic>>(counters_freq.data(), depth * width_freq);
		
		Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<support_type>> solver;
		solver.compute(sensingmatrix);
		Eigen::Vector<support_type, Eigen::Dynamic> sol = solver.solve(eigen_counters);
		
#pragma omp parallel for
		for (int j = 0; j < query_size; j++) {
			index->at(j) = j;
			estimates->at(j).at(0) = sol(j);
		}
		return std::make_tuple(index, estimates, estimates);
	}
	
	static std::string getDescription() {
		return "PRSN";
	}
	
	auto getNummemcellsFast() const {
		return depth * width_freq;
	}
	
	auto getNummemcells() const {
		return depth * width_freq;
	}
	
	auto getMemStats() {
		std::cout << "\n-> " + getDescription() + " Memory" << std::endl;
		std::cout << "--> On-Line Memory: " << getNummemcells() << std::endl;
		std::cout << "--> On-Line Counters: " << getNummemcellsFast() << std::endl;
		memstats.clear();
		memstats.push_back(std::make_tuple("OnLineMemory", getNummemcells()));
		memstats.push_back(std::make_tuple("OnLineCounters", getNummemcellsFast()));
		return &memstats;
	}
};

#endif /* PRSketch_hpp */
