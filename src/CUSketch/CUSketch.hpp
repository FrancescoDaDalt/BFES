//
//  CUSketch.hpp
//  MCMC_Sketch_Aux
//
//  Created by Francesco Da Dalt on 28.02.24.
//

#ifndef CUSketch_hpp
#define CUSketch_hpp

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <vector>

#include "MurmurHash3.h"


// CHECKED, all good
template <typename key_type, typename support_type>
class CUSketch {
	const int depth;
	const int width;
	std::vector<support_type> counters;
	static constexpr auto basicHash = std::hash<key_type>{};
	std::vector<std::tuple<std::string, double>> memstats;
public:
	CUSketch (const int nummemcells):
	depth((int) 5),
	width(nummemcells / depth) {
		counters.resize(depth * width);
		std::fill(counters.begin(), counters.end(), 0);
	}
	
	void insert (const std::tuple<key_type, support_type> key_value) {
		const auto key = std::get<0>(key_value);
		const auto value = std::get<1>(key_value);
		const size_t firsthash = basicHash(key);
		
		double current_est = std::numeric_limits<double>::max();
		for (int i = 0; i < depth; i++) {
			size_t hashedkey;
			size_t extradest[2];
			MurmurHash3_x64_128(&firsthash,
								sizeof(size_t),
								7218563 * i + 1543724,
								&extradest);
			hashedkey = extradest[0] % width;
			current_est = std::min(counters.at(i * width + hashedkey), current_est);
		}
		
		for (int i = 0; i < depth; i++) {
			size_t hashedkey;
			size_t extradest[2];
			MurmurHash3_x64_128(&firsthash,
								sizeof(size_t),
								7218563 * i +1543724,
								&extradest);
			hashedkey = extradest[0] % width;
			counters.at(i * width + hashedkey) = std::max(current_est + value, counters.at(i * width + hashedkey));
		}
	}
	
	auto batch_query (const std::vector<key_type>* const keys) {
		const int query_size = keys->size();
		auto* const index = new std::vector<size_t>(query_size);
		auto* const estimates = new std::vector<std::vector<double>>(query_size, std::vector<double>(1));
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
				hashedkey = extradest[0] % width;
				est = std::min(counters.at(i * width + hashedkey), est);
			}
			index->at(j) = j;
			estimates->at(j).at(0) = est;
		}
		return std::make_tuple(index, estimates, estimates);
	}
	
	static std::string getDescription() {
		return "CUS";
	}
	
	auto getNummemcellsFast() const {
		return width * depth;
	}
	
	auto getNummemcells() const {
		return width * depth;
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

#endif /* CUSketch_hpp */
