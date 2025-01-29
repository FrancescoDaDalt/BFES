//
//  CBSketch.hpp
//  MCMC_Sketch_Aux
//
//  Created by Francesco Da Dalt on 13.06.23.
//

#ifndef CBSketch_hpp
#define CBSketch_hpp

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <vector>

#include "distinct_counter.h"
#include "Murmurhash3.h"

// CHECKED, all good
template <typename key_type, typename support_type>
class CBSketch {
	const int num_memory_cells;
	const int depth;
	const int width;
	support_type total_counter = 0;
	std::vector<support_type> counters;
	static constexpr auto basicHash = std::hash<key_type>{};
	hyperloglog_hip::distinct_counter<size_t> distinct_counter;
	std::vector<std::tuple<std::string, double>> memstats;
public:
	CBSketch (const int memory_cells):
	num_memory_cells(memory_cells - ((int) 4 * std::log2(memory_cells))),
	depth((int) 5),
	width(num_memory_cells / depth) {
		counters.resize(depth * width);
		new (&distinct_counter) hyperloglog_hip::distinct_counter<size_t>((int) 4 * std::log2(memory_cells));
		std::fill(counters.begin(), counters.end(), 0);
	}
	
	void insert (const std::tuple<key_type, support_type> key_value) {
		const auto key = std::get<0>(key_value);
		const auto value = std::get<1>(key_value);
		const size_t firsthash = basicHash(key);
		total_counter += value;
		
		for (int i = 0; i < depth; i++) {
			size_t hashedkey;
			size_t extradest[2];
			MurmurHash3_x64_128(&firsthash,
								sizeof(size_t),
								7218563 * i +1543724,
								&extradest);
			hashedkey = extradest[0] % width;
			if (i == 0) {distinct_counter.insert(extradest[0]);}
			counters.at(i * width + hashedkey) += value;
		}
	}
	
	auto batch_query (const std::vector<key_type>* const keys) {
		const int query_size = keys->size();
		auto* const index = new std::vector<size_t>(query_size);
		auto* const estimates = new std::vector<std::vector<double>>(query_size, std::vector<double>(1));
		const double mu = total_counter / ((double) distinct_counter.count());
		const double chi_inv = (double) distinct_counter.count();
#pragma omp parallel for
		for (int j = 0; j < query_size; j++) {
			const auto key = keys->at(j);
			const size_t firsthash = basicHash(key);
			double est = 0;
			for (int i = 0; i < depth; i++) {
				size_t hashedkey;
				size_t extradest[2];
				MurmurHash3_x64_128(&firsthash,
									sizeof(size_t),
									7218563 * i + 1543724,
									&extradest);
				hashedkey = extradest[0] % width;
				est += counters.at(i * width + hashedkey);
			}
			est = (mu * chi_inv + width * est - depth * total_counter) / (chi_inv + depth * (width - 1.0));
			index->at(j) = j;
			estimates->at(j).at(0) = est;
		}
		return std::make_tuple(index, estimates, estimates);
	}
	
	static std::string getDescription() {
		return "CBS";
	}
	
	auto getNummemcellsFast() const {
		return width * depth + 1 + distinct_counter.size() / 32;
	}
	
	auto getNummemcells() const {
		return width * depth + 1 + distinct_counter.size() / 32;
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

#endif /* CBSketch_hpp */
