//
//  DPSketch.hpp
//  braindust
//
//  Created by Anynomous
//

#ifndef DPSketch_hpp
#define DPSketch_hpp

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_sf_gamma.h>

#include "MurmurHash3.h"

template <typename key_type, typename support_type>
class DPSketch {
public:
	const double scale_multiplier;
	const int numsamples;
	const int depth;
	const int width;
	std::vector<support_type> counters;
	double M;
	static constexpr auto basicHash = std::hash<key_type>{};
	std::vector<std::tuple<std::string, double>> memstats;
	
	std::vector<gsl_rng*> rngs;
	
	
	static double empiricalBayes_f(double log_alpha, void *params) {
		const DPSketch* obj = (DPSketch*) params;
		double alpha = std::exp(log_alpha);
		double M = obj->M;
		double N = (double) obj->depth;
		double J = (double) obj->width;
		double dataLogLikelihood = N * (gsl_sf_lngamma(M + 1) + gsl_sf_lngamma(alpha) - gsl_sf_lngamma(M + alpha) - J * gsl_sf_lngamma(alpha / J));
		for (int i = 0; i < obj->depth * obj->width; i++) {
			dataLogLikelihood += (gsl_sf_lngamma(obj->counters[i] + alpha / J) - gsl_sf_lngamma(obj->counters[i] + 1));
		}
		
		return -dataLogLikelihood;
	}
	
	DPSketch (const int nummemcells, const int numsamples, const double scale_multiplier):
	scale_multiplier(scale_multiplier),
	numsamples(numsamples),
	depth((int) 5),
	width(nummemcells / depth) {
		counters.resize(depth * width);
		std::fill(counters.begin(), counters.end(), 0);
		rngs = std::vector<gsl_rng*>(omp_get_max_threads());
		for (int i = 0; i < omp_get_max_threads(); i++) {
			rngs[i] = gsl_rng_alloc(gsl_rng_default);
			gsl_rng_set(rngs[i], 42 * i + 1234567);
		}
	}
	
	void insert (const std::tuple<key_type, support_type> key_value) {
		const auto key = std::get<0>(key_value);
		const auto value = std::get<1>(key_value) * scale_multiplier;
		const size_t firsthash = basicHash(key);
		
		for (int i = 0; i < depth; i++) {
			size_t hashedkey;
			size_t extradest[2];
			MurmurHash3_x64_128(&firsthash,
								sizeof(size_t),
								7218563 * i +1543724,
								&extradest);
			hashedkey = extradest[0] % width;
			counters.at(i * width + hashedkey) += value;
		}
	}
	
	auto batch_query (const std::vector<key_type>* const keys) {
		const int query_size = keys->size();
		auto* const index = new std::vector<size_t>(query_size);
		auto* const estimates = new std::vector<std::vector<double>>(query_size, std::vector<double>(numsamples));
		
		M = std::accumulate(counters.begin(), counters.begin() + width, 0);
		
		gsl_min_fminimizer *s;
		gsl_function F;
		
		F.function = &this->empiricalBayes_f;
		F.params = (void*) this;
		
//		lower bound and upper bounds of log(alpha) chosen based on machine precision
		double lower_bound = -20.0;
		double upper_bound = 20.0;
		double x_initial = 0.0;
		
		if (empiricalBayes_f(lower_bound, F.params) <= empiricalBayes_f(lower_bound + 1, F.params) and
			empiricalBayes_f(upper_bound - 1, F.params) <= empiricalBayes_f(upper_bound, F.params)) {
			x_initial = lower_bound;
		} else if (empiricalBayes_f(lower_bound, F.params) >= empiricalBayes_f(lower_bound + 1, F.params) and
				   empiricalBayes_f(upper_bound - 1, F.params) >= empiricalBayes_f(upper_bound, F.params)) {
			x_initial = upper_bound;
		} else {
			
			//		Uses brent-method which is a pseudo-quadratic method which should perform similar to the Newton-Raphson method described in the original paper.
			s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
			gsl_min_fminimizer_set(s, &F, x_initial, lower_bound, upper_bound);
			
			int status;
			int iter = 0;
			int max_iter = 1000;
			
			do {
				iter++;
				status = gsl_min_fminimizer_iterate(s);
				
				x_initial = gsl_min_fminimizer_x_minimum(s);
				lower_bound = gsl_min_fminimizer_x_lower(s);
				upper_bound = gsl_min_fminimizer_x_upper(s);
				
				status = gsl_min_test_interval(lower_bound, upper_bound, 1e-5, 0.0);
				
				if (status == GSL_SUCCESS) {
					std::cout << "--> DPS Alpha Optimizer: Converged to minimum at: " << x_initial << std::endl;
					std::cout << "--> DPS Alpha Optimizer: Maximum value of the function: " << -empiricalBayes_f(x_initial, this) << std::endl;
				}
			} while (status == GSL_CONTINUE && iter < max_iter);
			
			gsl_min_fminimizer_free(s);
		}
		
		double alpha = std::exp(x_initial);
		
#pragma omp parallel for
		for (int j = 0; j < query_size; j++) {
			std::vector<unsigned int> samples_indices(numsamples);
			std::vector<double> samples_value(numsamples);
			std::vector<double> samples_loglikelihood(numsamples);
			std::vector<double> samples_likelihood(numsamples);
			
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
			const double min = 0.0;
			const double max = est;
			double max_loglikelihood = -std::numeric_limits<double>::infinity();
			for (int i = 0; i < numsamples; i++) {
				samples_value[i] = min + ((max - min) * i) / (numsamples - 1);
				samples_loglikelihood[i] = 0.0;
				for (int k = 0; k < depth; k++) {
					size_t hashedkey;
					size_t extradest[2];
					MurmurHash3_x64_128(&firsthash,
										sizeof(size_t),
										7218563 * k + 1543724,
										&extradest);
					hashedkey = extradest[0] % width;
					double counter = counters.at(k * width + hashedkey);
					samples_loglikelihood[i] += -std::log(1 + counter * width / alpha) + gsl_sf_lngamma(counter + 1) + gsl_sf_lngamma(counter - samples_value[i] + alpha / width) - gsl_sf_lngamma(counter + alpha / width) - gsl_sf_lngamma(counter - samples_value[i] + 1);
				}
				max_loglikelihood = std::max(max_loglikelihood, samples_loglikelihood[i]);
			}
			double total_likelihood = 0.0;
			for (int i = 0; i < numsamples; i++) {
				samples_loglikelihood[i] -= max_loglikelihood;
				samples_likelihood[i] = std::exp(samples_loglikelihood[i]);
				total_likelihood += samples_likelihood[i];
			}
//			double total_likelihood = std::accumulate(samples_likelihood.begin(), samples_likelihood.end(), 0);
//			for (int i = 0; i < numsamples; i++) {
//				std::cout << samples_likelihood[i] << std::endl;
//			}
//			std::cout << std::endl;
//			std::cout << total_likelihood << std::endl;
//			if (total_likelihood == 0) {
//				for (int i = 0; i < numsamples; i++) {
//					samples_likelihood[i] = 1.0 / numsamples;
//				}
//			} else if (!std::isfinite(total_likelihood)) {
//				total_likelihood = 0.0;
//				for (int i = 0; i < numsamples; i++) {
//					if (!std::isfinite(samples_likelihood[i])) {
//						samples_likelihood[i] = 1.0;
//						total_likelihood += 1.0;
//					} else {
//						samples_likelihood[i] = 0.0;
//					}
//				}
//				for (int i = 0; i < numsamples; i++) {
//					samples_likelihood[i] /= total_likelihood;
//				}
//			} else {
//				for (int i = 0; i < numsamples; i++) {
//					samples_likelihood[i] /= total_likelihood;
//				}
//			}
			for (int i = 0; i < numsamples; i++) {
				samples_likelihood[i] /= total_likelihood;
			}
//			for (int i = 0; i < numsamples; i++) {
//				std::cout << samples_likelihood[i] << std::endl;
//			}
//			std::cout << std::endl;
			gsl_ran_multinomial(rngs[omp_get_thread_num()], numsamples, numsamples, samples_likelihood.data(), samples_indices.data());
			
			index->at(j) = j;
			
			int idx = 0;
			for (int i = 0; i < numsamples; i++) {
				for (int k = 0; k < samples_indices[i]; k++) {
					estimates->at(j).at(idx) = samples_value.at(i) / scale_multiplier;
					idx++;
				}
			}
		}
		return std::make_tuple(index, estimates, estimates);
	}
	
	static std::string getDescription() {
		return "DPS";
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

#endif /* DPSketch_hpp */
