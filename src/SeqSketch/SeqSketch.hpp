//
//  SeqSketch.hpp
//  MCMC_Sketch_Aux
//
//  Created by Francesco Da Dalt on 13.06.23.
//

#ifndef SeqSketch_hpp
#define SeqSketch_hpp

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <vector>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <gsl/gsl_cdf.h>

#include "MurmurHash3.h"
#include "fusion.h"

namespace msk = mosek::fusion;
namespace mty = monty;


// CHECKED, all good
template <typename key_type, typename support_type>
class SeqSketch {
	const int depth = 1;
	const int width_freq;
	const int width_filter;
	std::vector<support_type> counters_freq;
	std::vector<bool> counters_filter;
	std::vector<key_type> identified_keys;
	static constexpr auto basicHash = std::hash<key_type>{};
	std::vector<std::tuple<std::string, double>> memstats;
public:
	SeqSketch (const int total_counters):
	width_freq(total_counters * (1.0 - 0.125)), // 12.5% of memory for bloom filter and rest for actual frequency counters
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
			size_t hashedkey_index;
			double hashedkey_multiplier;
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
			hashedkey_index = extradest_index[0] % width_freq;
			hashedkey_multiplier = gsl_cdf_ugaussian_Pinv(((extradest_multiplier[0] % 10000) + 1) / 10001.0) / std::sqrt(width_freq);
			counters_freq.at(i * width_freq + hashedkey_index) += value * hashedkey_multiplier;
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
		
		std::vector<int> system_row;
		std::vector<int> system_col;
		std::vector<double> system_val;
		std::vector<double> system_vec;
		
		system_row.reserve(identified_keys_size);
		system_col.reserve(identified_keys_size);
		system_val.reserve(identified_keys_size);
		for (int j = 0; j < identified_keys_size; j++) {
			const auto key = identified_keys.at(j);
			const size_t firsthash = basicHash(key);
			for (int i = 0; i < depth; i++) {
				size_t hashedkey_index;
				double hashedkey_multiplier;
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
				hashedkey_index = extradest_index[0] % width_freq;
				hashedkey_multiplier = gsl_cdf_ugaussian_Pinv(((extradest_multiplier[0] % 10000) + 1) / 10001.0) / std::sqrt(width_freq);
				key_indexer.at(hashedkey_index) = j;
				system_val.push_back(hashedkey_multiplier);
				system_row.push_back((int) hashedkey_index);   // row index
				system_col.push_back(j);
			}
		}
		
		for (int i = 0; i < depth * width_freq; i++) {
			if (key_indexer.at(i) == -1) {
				key_indexer.at(i) = identified_keys_size;
				system_val.push_back(1.0);
				system_row.push_back((int) i);
				system_col.push_back(identified_keys_size);
				identified_keys_size++;
			}
		}
		
		
		std::vector<double> counters_freq_double(counters_freq.begin(), counters_freq.end());
		auto system_vec_mosak = mty::new_array_ptr<double>(counters_freq_double);
		
		msk::Model::t mdl = new msk::Model("SeqSketch Solver");
		auto _mdl = mty::finally([&](){mdl->dispose();});
		msk::Variable::t x = mdl->variable("x", identified_keys_size, msk::Domain::unbounded());
		msk::Variable::t u = mdl->variable("u", identified_keys_size, msk::Domain::unbounded());
		
		auto constraint_mat = msk::Matrix::sparse((int) depth * width_freq, (int) identified_keys_size, mty::new_array_ptr<int>(system_row), mty::new_array_ptr<int>(system_col), mty::new_array_ptr<double>(system_val));
		
		mdl->constraint(msk::Expr::mul(constraint_mat, x), msk::Domain::equalsTo(system_vec_mosak));
		mdl->constraint(msk::Expr::sub(u, x), msk::Domain::greaterThan(0.0));
		mdl->constraint(msk::Expr::add(u, x), msk::Domain::greaterThan(0.0));
		
		std::vector<double> aux_obj_vec(identified_keys_size, 1.0);
		auto aux_obj_mosek = mty::new_array_ptr<double>(aux_obj_vec);
		mdl->objective(msk::ObjectiveSense::Minimize, msk::Expr::dot(aux_obj_mosek, u));
		mdl->solve();
		
		Eigen::VectorXd sol = Eigen::VectorXd::Zero(identified_keys_size);
		auto sol_msk = x->level();
		std::memcpy((void*) sol.data(),(void*) sol_msk->begin(), identified_keys_size * sizeof(double));
		
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
		return "SeqS";
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
class SeqSketch_Native {
	const int depth = 1;
	const int width_freq;
	std::vector<support_type> counters_freq;
	static constexpr auto basicHash = std::hash<key_type>{};
	std::vector<std::tuple<std::string, double>> memstats;
public:
	SeqSketch_Native (const int total_counters):
	width_freq(total_counters) {
		counters_freq.resize(depth * width_freq);
		std::fill(counters_freq.begin(), counters_freq.end(), 0);
	}
	
	void insert (const std::tuple<key_type, support_type> key_value) {
		const auto key = std::get<0>(key_value);
		const auto value = std::get<1>(key_value);
		const size_t firsthash = basicHash(key);
		for (int i = 0; i < depth; i++) {
			size_t hashedkey_index;
			double hashedkey_multiplier;
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
			hashedkey_index = extradest_index[0] % width_freq;
			hashedkey_multiplier = gsl_cdf_ugaussian_Pinv(((extradest_multiplier[0] % 10000) + 1) / 10001.0) / std::sqrt(width_freq);
			counters_freq.at(i * width_freq + hashedkey_index) += value * hashedkey_multiplier;
		}
	}
	
	auto batch_query (const std::vector<key_type>* const keys) {
		
		const int query_size = keys->size();
		auto* const index = new std::vector<size_t>(query_size);
		auto* const estimates = new std::vector<std::vector<double>>(query_size, std::vector<double>(1));
		
		std::vector<int> system_row;
		std::vector<int> system_col;
		std::vector<double> system_val;
		std::vector<double> system_vec;
		
		system_row.reserve(query_size);
		system_col.reserve(query_size);
		system_val.reserve(query_size);
		for (int j = 0; j < query_size; j++) {
			const auto key = keys->at(j);
			const size_t firsthash = basicHash(key);
			for (int i = 0; i < depth; i++) {
				size_t hashedkey_index;
				double hashedkey_multiplier;
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
				hashedkey_index = extradest_index[0] % width_freq;
				hashedkey_multiplier = gsl_cdf_ugaussian_Pinv(((extradest_multiplier[0] % 10000) + 1) / 10001.0) / std::sqrt(width_freq);
				system_val.push_back(hashedkey_multiplier);
				system_row.push_back((int) hashedkey_index);   // row index
				system_col.push_back(j);
			}
		}
		
		std::vector<double> counters_freq_double(counters_freq.begin(), counters_freq.end());
		auto system_vec_mosak = mty::new_array_ptr<double>(counters_freq_double);
		
		msk::Model::t mdl = new msk::Model("SeqSketch Solver");
		auto _mdl = mty::finally([&](){mdl->dispose();});
		msk::Variable::t x = mdl->variable("x", query_size, msk::Domain::unbounded());
		msk::Variable::t u = mdl->variable("u", query_size, msk::Domain::unbounded());
		
		auto constraint_mat = msk::Matrix::sparse((int) depth * width_freq, (int) query_size, mty::new_array_ptr<int>(system_row), mty::new_array_ptr<int>(system_col), mty::new_array_ptr<double>(system_val));
		
		mdl->constraint(msk::Expr::mul(constraint_mat, x), msk::Domain::equalsTo(system_vec_mosak));
		mdl->constraint(msk::Expr::sub(u, x), msk::Domain::greaterThan(0.0));
		mdl->constraint(msk::Expr::add(u, x), msk::Domain::greaterThan(0.0));
		
		std::vector<double> aux_obj_vec(query_size, 1.0);
		auto aux_obj_mosek = mty::new_array_ptr<double>(aux_obj_vec);
		mdl->objective(msk::ObjectiveSense::Minimize, msk::Expr::dot(aux_obj_mosek, u));
		mdl->solve();
		
		Eigen::VectorXd sol = Eigen::VectorXd::Zero(query_size);
		auto sol_msk = x->level();
		std::memcpy((void*) sol.data(),(void*) sol_msk->begin(), query_size * sizeof(double));
		
#pragma omp parallel for
		for (int j = 0; j < query_size; j++) {
			index->at(j) = j;
			estimates->at(j).at(0) = sol(j);
		}
		return std::make_tuple(index, estimates, estimates);
	}
	
	static std::string getDescription() {
		return "SeqSN";
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

#endif /* SeqSketch_hpp */
