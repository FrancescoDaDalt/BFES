//
//  CSketch.hpp
//  MCMC_Sketch_Aux
//
//  Created by Francesco Da Dalt on 13.06.23.
//

#ifndef CSketch_hpp
#define CSketch_hpp

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <vector>

#include "MurmurHash3.h"


// CHECKED, all good
template <typename key_type, typename support_type>
class CSketch {
	const int depth;
	const int width;
	std::vector<support_type> counters;
	static constexpr auto basicHash = std::hash<key_type>{};
	std::vector<std::tuple<std::string, double>> memstats;
public:
	CSketch (const int nummemcells):
	depth((int) 5),
	width(nummemcells / depth) {
		counters.resize(depth * width);
		std::fill(counters.begin(), counters.end(), 0);
	}
	
	void insert (const std::tuple<key_type, support_type> key_value) {
		const auto key = std::get<0>(key_value);
		const auto value = std::get<1>(key_value);
		const size_t firsthash = basicHash(key);
		
		for (int i = 0; i < depth; i++) {
			size_t hashedkey_index;
			int hashedkey_multiplier;
			size_t extradest_index[2];
			size_t extradest_multiplier[2];
			MurmurHash3_x64_128(&firsthash,
								sizeof(size_t),
								7218563 * i + 1543724,
								&extradest_index);
			MurmurHash3_x64_128(&firsthash,
								sizeof(size_t),
								6494787 * i + 8209594,
								&extradest_multiplier);
			hashedkey_index = extradest_index[0] % width;
			hashedkey_multiplier = (extradest_multiplier[0] % 2) * 2 - 1;
			counters.at(i * width + hashedkey_index) += hashedkey_multiplier * value;
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
			std::vector<double> estvec(depth);
			for (int i = 0; i < depth; i++) {
				size_t hashedkey_index;
				int hashedkey_multiplier;
				size_t extradest_index[2];
				size_t extradest_multiplier[2];
				MurmurHash3_x64_128(&firsthash,
									sizeof(size_t),
									7218563 * i + 1543724,
									&extradest_index);
				MurmurHash3_x64_128(&firsthash,
									sizeof(size_t),
									6494787 * i + 8209594,
									&extradest_multiplier);
				hashedkey_index = extradest_index[0] % width;
				hashedkey_multiplier = (extradest_multiplier[0] % 2) * 2 - 1;
				estvec.at(i) = counters.at(i * width + hashedkey_index) * hashedkey_multiplier;
			}
			std::nth_element(estvec.begin(), estvec.begin() + depth / 2, estvec.end());
			index->at(j) = j;
			estimates->at(j).at(0) = *(estvec.begin() + depth / 2);
		}
		return std::make_tuple(index, estimates, estimates);
	}
	
	static std::string getDescription() {
		return "CS";
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

#endif /* CSketch_hpp */
