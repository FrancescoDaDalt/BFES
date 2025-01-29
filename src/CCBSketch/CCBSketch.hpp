//
//  CCBSketch.hpp
//  MCMC_Sketch_Aux
//
//  Created by Francesco Da Dalt on 13.06.23.
//

#ifndef CCBSketch_hpp
#define CCBSketch_hpp

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <vector>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "Murmurhash3.h"
#include "distinct_counter.h"


// CHECKED, all good
template <typename key_type, typename support_type>
class CCBSketch {
	const int num_memory_cells;
	const int freq_counters;
	const int depth_freq;
	const int width_freq;
	const int width_filter;
	support_type total_counter = 0;
	std::vector<support_type> counters_freq;
	std::vector<bool> counters_filter;
	std::vector<key_type> identified_keys;
	static constexpr auto basicHash = std::hash<key_type>{};
	hyperloglog_hip::distinct_counter<size_t> distinct_counter;
	std::vector<std::tuple<std::string, double>> memstats;
public:
	CCBSketch (const int total_counters):
	num_memory_cells(total_counters - (int) 4 * std::log2(total_counters)),
	freq_counters(num_memory_cells * (1.0 - 0.125)),
	depth_freq((int) 5),
	width_freq(freq_counters / depth_freq),
	width_filter((num_memory_cells - freq_counters) * 32) { // We allow more "memory" because a bool vector is 32 times more efficient than a vector of ints or floats.
		counters_freq.resize(depth_freq * width_freq);
		counters_filter.resize(width_filter);
		new (&distinct_counter) hyperloglog_hip::distinct_counter<size_t>((int) 4 * std::log2(total_counters));
		std::fill(counters_freq.begin(), counters_freq.end(), 0);
		std::fill(counters_filter.begin(), counters_filter.end(), false);
	}
	
	void insert (const std::tuple<key_type, support_type> key_value) {
		const auto key = std::get<0>(key_value);
		const auto value = std::get<1>(key_value);
		const size_t firsthash = basicHash(key);
		total_counter += value;
		for (int i = 0; i < depth_freq; i++) {
			size_t hashedkey;
			size_t extradest[2];
			MurmurHash3_x64_128(&firsthash,
								sizeof(size_t),
								7218563 * i + 1543724,
								&extradest);
			if (i == 0) {distinct_counter.insert(extradest[0]);}
			hashedkey = extradest[0] % width_freq;
			counters_freq.at(i * width_freq + hashedkey) += value;
		}
		for (int i = 0; i < 1; i++) {
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
		
		const double mu = total_counter / ((double) distinct_counter.count());
		const double chi_inv = (double) distinct_counter.count();
		
		std::vector<support_type> counters_card(depth_freq * width_freq, 1);
		
		
		int identified_keys_size = identified_keys.size();
		for (int j = 0; j < identified_keys_size; j++) {
			const auto key = identified_keys.at(j);
			const size_t firsthash = basicHash(key);
			for (int i = 0; i < depth_freq; i++) {
				size_t hashedkey;
				size_t extradest[2];
				MurmurHash3_x64_128(&firsthash,
									sizeof(size_t),
									7218563 * i + 1543724,
									&extradest);
				hashedkey = extradest[0] % width_freq;
				counters_card.at(i * width_freq + hashedkey) += 1;
			}
		}
#pragma omp parallel for
		for (int j = 0; j < query_size; j++) {
			const auto key = keys->at(j);
			const size_t firsthash = basicHash(key);
			double est = 0;
			double est_aux = 0;
			for (int i = 0; i < depth_freq; i++) {
				size_t hashedkey;
				size_t extradest[2];
				MurmurHash3_x64_128(&firsthash,
									sizeof(size_t),
									7218563 * i + 1543724,
									&extradest);
				hashedkey = extradest[0] % width_freq;
				est += counters_freq.at(i * width_freq + hashedkey) / std::max(counters_card.at(i * width_freq + hashedkey) - 1, 1.0);
				est_aux += (identified_keys_size - counters_card.at(i * width_freq + hashedkey)) / std::max(counters_card.at(i * width_freq + hashedkey) - 1, 1.0);
			}
			est = (mu * chi_inv + identified_keys_size * est - depth_freq * total_counter) / (chi_inv + est_aux);
			index->at(j) = j;
			estimates->at(j).at(0) = est;
		}
		return std::make_tuple(index, estimates, estimates);
	}
	
	static std::string getDescription() {
		return "CCBS";
	}
	
	auto getNummemcellsFast() const {
		return depth_freq * width_freq + width_filter / 32 + 1 + distinct_counter.size() / 32;
	}
	
	auto getNummemcells() const {
		return depth_freq * width_freq + width_filter / 32 + 1 + distinct_counter.size() / 32 + identified_keys.size();
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
class CCBSketch_Native {
	const int freq_counters;
	const int depth_freq;
	const int width_freq;
	support_type total_counter = 0;
	std::vector<support_type> counters_freq;
	static constexpr auto basicHash = std::hash<key_type>{};
	std::vector<std::tuple<std::string, double>> memstats;
public:
	CCBSketch_Native (const int total_counters):
	freq_counters(total_counters),
	depth_freq((int) std::sqrt(freq_counters)),
	width_freq(freq_counters / depth_freq) {
		counters_freq.resize(depth_freq * width_freq);
		std::fill(counters_freq.begin(), counters_freq.end(), 0);
	}
	
	void insert (const std::tuple<key_type, support_type> key_value) {
		const auto key = std::get<0>(key_value);
		const auto value = std::get<1>(key_value);
		const size_t firsthash = basicHash(key);
		total_counter += value;
		for (int i = 0; i < depth_freq; i++) {
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
		
		const double mu = total_counter / ((double) query_size);
		const double chi_inv = (double) query_size;
		
		std::vector<support_type> counters_card(depth_freq * width_freq, 1);
		
		for (int j = 0; j < query_size; j++) {
			const auto key = keys->at(j);
			const size_t firsthash = basicHash(key);
			for (int i = 0; i < depth_freq; i++) {
				size_t hashedkey;
				size_t extradest[2];
				MurmurHash3_x64_128(&firsthash,
									sizeof(size_t),
									7218563 * i + 1543724,
									&extradest);
				hashedkey = extradest[0] % width_freq;
				counters_card.at(i * width_freq + hashedkey) += 1;
			}
		}
#pragma omp parallel for
		for (int j = 0; j < query_size; j++) {
			const auto key = keys->at(j);
			const size_t firsthash = basicHash(key);
			double est = 0;
			double est_aux = 0;
			for (int i = 0; i < depth_freq; i++) {
				size_t hashedkey;
				size_t extradest[2];
				MurmurHash3_x64_128(&firsthash,
									sizeof(size_t),
									7218563 * i + 1543724,
									&extradest);
				hashedkey = extradest[0] % width_freq;
				est += counters_freq.at(i * width_freq + hashedkey) / std::max(counters_card.at(i * width_freq + hashedkey) - 1, 1.0);
				est_aux += (query_size - counters_card.at(i * width_freq + hashedkey)) / std::max(counters_card.at(i * width_freq + hashedkey) - 1, 1.0);
			}
			est = (mu * chi_inv + query_size * est - depth_freq * total_counter) / (chi_inv + est_aux);
			index->at(j) = j;
			estimates->at(j).at(0) = est;
		}
		return std::make_tuple(index, estimates, estimates);
	}
	
	static std::string getDescription() {
		return "CCBSN";
	}
	
	auto getNummemcellsFast() const {
		return depth_freq * width_freq + 1;
	}
	
	auto getNummemcells() const {
		return depth_freq * width_freq + 1;
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




#endif /* CCBSketch_hpp */
