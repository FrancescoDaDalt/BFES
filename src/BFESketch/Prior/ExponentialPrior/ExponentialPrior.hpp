/* 
 * MIT License
 * 
 * Copyright (c) 2024 Francesco Da Dalt
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
 */

#ifndef AffineExponentialPrior_hpp
#define AffineExponentialPrior_hpp


#include "../../Sensor/LatentBaseVector/LatentBaseVector.hpp"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

// Prior that implements an exponential prior. Can be used both for level-1 and level-2 sampling.
// Assertion statements limit the mplmentation to cases relevant for the benhcmarking
template <typename supp_type>
class ExponentialPrior {
	
//	The distribution defined by this object is trunc(Exp(scale) + low, [low, high])
	supp_type low = 0;
	supp_type high = 0;
	supp_type scale = 0;
	supp_type supp_size = 0;
	int direction = 0;
	
	int num_batch_approximators = 0;
	std::vector<supp_type> batch_approximators_low;
	std::vector<supp_type> batch_approximators_high;
	std::vector<supp_type> batch_approximators_scales;
	
public:
	using supp_t = supp_type;
	using intermediate_t = double;
	
	ExponentialPrior () {
	}
	
	const auto get_num_batch_approximators() const {
		return num_batch_approximators;
	}
	
	const auto variance_factor() const {
		return 2.0;
	}
	
	const auto mean(const supp_t counter) const {
		assert(direction > 0);
		const double corrected_high = std::min(counter, high);
		const double corrected_supp_size = corrected_high - low;
		if (corrected_supp_size > 16 * scale) {
			return scale + low;
		} else {
			return scale + corrected_supp_size / (std::exp(-corrected_supp_size / scale) - 1) + corrected_high;
		}
	}
	
	template<typename array_type>
	void normPlusBatchMean(const supp_t counter, array_type* const ans) const {
		assert(direction > 0);
		ans->resize(num_batch_approximators + 1);
		ans->operator()(0) = mean(counter);
		for (int i = 1; i < num_batch_approximators + 1; i++) {
			if (supp_size > 16 * scale or std::abs(batch_approximators_scales[i - 1]) <= 1e-6) {
				const auto scl = batch_approximators_scales[i - 1];
				
				ans->operator()(i) = (scl > 0) ? scl + batch_approximators_low[i - 1] : scl + batch_approximators_high[i - 1];
			} else {
				assert(false);
				const auto scl = batch_approximators_scales[i - 1];
				const auto supp_siz = batch_approximators_high[i - 1] - batch_approximators_low[i - 1];
				ans->operator()(i) = scl + supp_siz / (std::exp(supp_size / scl) - 1) + batch_approximators_high[i - 1];
			}
		}
	}
	
	const auto variance(const supp_t counter) const {
		assert(direction > 0);
		const double corrected_high = std::min(counter, high);
		const double corrected_supp_size = corrected_high - low;
		if (corrected_supp_size > 16 * scale) {
			return scale * scale;
		} else {
			return scale * scale + 0.5 * corrected_supp_size * corrected_supp_size / (1.0 - std::cosh(corrected_supp_size / scale));
		}
	}
	
	template<typename array_type>
	void normPlusBatchVar(const supp_t counter, array_type* const ans) const {
		ans->resize(num_batch_approximators + 1);
		ans->operator()(0) = variance(counter);
		for (int i = 1; i < num_batch_approximators + 1; i++) {
			if (supp_size > 16 * scale or std::abs(batch_approximators_scales[i - 1]) <= 1e-6) {
				const auto scl = batch_approximators_scales[i - 1];
				ans->operator()(i) = scl * scl;
			} else {
				assert(false);
				const auto scl = batch_approximators_scales[i - 1];
				const auto supp_siz = batch_approximators_high[i - 1] - batch_approximators_low[i - 1];
				ans->operator()(i) = scl * scl + 0.5 * supp_siz * supp_siz / (1.0 - std::cosh(batch_approximators_high[i - 1] / scl));
			}
		}
	}
	
	void adjust_to_match_mean(const supp_type mean) {
		assert(supp_size > 16 * scale and mean < high and mean > low);
		if (direction > 0) {
			scale = (mean - low);
		} else {
			scale = (mean - high);
		}
	}
	
	void adjust_to_match_mean_var(const supp_type mean, const supp_type var) {
		assert(supp_size > 16 * scale and mean < high and mean > low);
		if (direction > 0) {
			scale = (mean - low);
		} else {
			scale = (mean - high);
		}
	}
	
	void set_batch_approximation(const int num_approximators, const double batchtarget_size, const supp_t batch_mean_, const supp_t batch_var_) {
		num_batch_approximators = num_approximators;
		assert(num_batch_approximators > 0);
		assert(direction > 0);
		
		supp_t batch_mean_local = batch_mean_;
		supp_t batch_std_local = std::sqrt( batch_var_);
		
		if (batch_std_local >= 1e-6) {
			
			for (int i = 0; i < 50; i++) {
				double x = (0 - batch_mean_local) / batch_std_local;
				
				if (x > 3) {
					double mean_actual = batch_mean_local + gsl_ran_ugaussian_tail_pdf(x, x) * batch_std_local;
					batch_mean_local += batch_mean_ - mean_actual;
					x = (0 - batch_mean_local) / batch_std_local;
					
					double std_actual = batch_std_local * std::sqrt(1 + x * gsl_ran_ugaussian_tail_pdf(x, x) - std::pow(gsl_ran_ugaussian_tail_pdf(x, x), 2.0));
					batch_std_local *= std::sqrt(std::sqrt( batch_var_) / std_actual);
				} else {
					double Z = gsl_cdf_ugaussian_Q(x);
					assert(Z > 1e-6);
					
					double mean_actual = batch_mean_local + gsl_ran_ugaussian_pdf(x) * batch_std_local / Z;
					batch_mean_local += batch_mean_ - mean_actual;
					x = (0 - batch_mean_local) / batch_std_local;
					Z = gsl_cdf_ugaussian_Q(x);
					assert(Z > 1e-6);
					
					double std_actual = batch_std_local * std::sqrt(1 + x * gsl_ran_ugaussian_pdf(x) / Z - std::pow(gsl_ran_ugaussian_pdf(x) / Z, 2.0));
					batch_std_local *= std::sqrt(std::sqrt( batch_var_) / std_actual);
					x = (0 - batch_mean_local) / batch_std_local;
					Z = gsl_cdf_ugaussian_Q(x);
					assert(Z > 1e-6);
				}
			}
		}
		supp_t batch_mean = batch_mean_local;
		supp_t batch_var = batch_std_local * batch_std_local;
		
		
		
		const double balancers_scale = std::sqrt(batch_var / (variance(std::numeric_limits<supp_t>::max()) * num_approximators));
		const double balancers_offset = batch_mean / num_approximators - balancers_scale * mean(std::numeric_limits<supp_t>::max());
		batch_approximators_scales = std::vector<supp_type>(num_approximators, scale * balancers_scale);
		batch_approximators_low = std::vector<supp_type>(num_approximators, low * balancers_scale + balancers_offset);
		batch_approximators_high = std::vector<supp_type>(num_approximators, high * balancers_scale + balancers_offset);
		batch_approximators_scales[0] *= -1;
		batch_approximators_high[0] = batch_mean / num_approximators + scale * balancers_scale;
		batch_approximators_low[0] = batch_approximators_high[0] - (high - low) * balancers_scale;
	}
	
	ExponentialPrior(supp_t low, supp_t high, supp_t scale): low(low), high(high), scale(scale), supp_size(high - low), direction((scale > 0) * 2 - 1) {
		assert(std::isfinite(scale) and std::abs(scale) > 1e-6);
		assert(std::isfinite(low));
		assert(std::isfinite(high));
	}
	
	const std::pair<supp_t, supp_t> get_bounds() const {
		return std::make_pair(low, high);
	}
	
	const std::pair<std::vector<supp_t>, std::vector<supp_t>> get_batchApproximator_bounds() const {
		std::pair<std::vector<supp_t>, std::vector<supp_t>> ans;
		assert(direction > 0);
		ans = std::make_pair(batch_approximators_low, batch_approximators_high);
		return ans;
	}
	
	template 	<	typename delta_type,
	typename coeff_type,
	typename index_type,
	typename state_type>
	const delta_type sampleHyperplane(const LatentBaseVectorDynamic<coeff_type, delta_type, index_type>* const hyperplane, const state_type* const state, gsl_rng* const rng_local) const {
		assert(direction > 0);
		
		using ArrayXs = Eigen::Array<supp_t, Eigen::Dynamic, 1>;
		using ArrayXh = Eigen::Array<coeff_type, Eigen::Dynamic, 1>;
		
		const auto projected_state = hyperplane->projectState(state);
		
		const auto hyperplane_coeffs = *hyperplane->get_coeff_map();
		
		assert(hyperplane_coeffs.sum() == 0);
		
		const auto plane_size = hyperplane->get_size();
		
		delta_type maximum_lower = std::numeric_limits<delta_type>::lowest();
		delta_type minimum_upper = std::numeric_limits<delta_type>::max();
		for (int i = 0; i < plane_size; i++) {
			const auto coeff = hyperplane_coeffs(i);
			const bool dir = coeff > 0;
			if (dir) {
				maximum_lower = std::max(maximum_lower, (low - projected_state(i)) / coeff);
				minimum_upper = std::min(minimum_upper, (high - projected_state(i)) / coeff);
			} else {
				minimum_upper = std::min(minimum_upper, (low - projected_state(i)) / coeff);
				maximum_lower = std::max(maximum_lower, (high - projected_state(i)) / coeff);
			}
		}
		
		if (maximum_lower >= minimum_upper) {
			return (delta_type) 0;
		}
		
		if constexpr (std::is_integral<delta_type>::value) {
			const delta_type delta = (delta_type) (maximum_lower + gsl_rng_uniform_int(rng_local, minimum_upper - maximum_lower + 1));
			return delta;
		} else {
			const delta_type delta = (delta_type) (maximum_lower + (minimum_upper - maximum_lower) * gsl_rng_uniform(rng_local));
			return delta;
		}
	}
	
	template 	<	typename delta_type,
	typename coeff_type,
	typename index_type,
	typename state_type>
	const delta_type sampleHyperplaneBalanced(const LatentBaseVectorDynamic<coeff_type, delta_type, index_type>* const hyperplane, const state_type* const state, const supp_t counter, gsl_rng* const rng_local) const {
		assert(direction > 0);
		const double corrected_high = std::min(counter, high);
		
		
		using ArrayXs = Eigen::Array<supp_t, Eigen::Dynamic, 1>;
		using ArrayXh = Eigen::Array<coeff_type, Eigen::Dynamic, 1>;
		
		const auto projected_state = hyperplane->projectState(state);
		
		const auto hyperplane_coeffs = *hyperplane->get_coeff_map();
		
		const auto plane_size = hyperplane->get_size();
		
		delta_type maximum_lower = std::numeric_limits<delta_type>::lowest();
		delta_type minimum_upper = std::numeric_limits<delta_type>::max();
		supp_t total_scale_inv = (supp_t) 0;
		
		const auto him = hyperplane->get_index_map();
		for (int i = 0; i < plane_size; i++) {
			const auto coeff = hyperplane_coeffs(i);
			const auto hyperplane_index = him->operator()(i);
			
			if (hyperplane_index == 0) {
				const bool dir = coeff > 0;
				total_scale_inv += coeff / scale;
				
				if (dir) {
					maximum_lower = std::max(maximum_lower, (low - projected_state(i)) / coeff);
					minimum_upper = std::min(minimum_upper, (corrected_high - projected_state(i)) / coeff);
				} else {
					minimum_upper = std::min(minimum_upper, (low - projected_state(i)) / coeff);
					maximum_lower = std::max(maximum_lower, (corrected_high - projected_state(i)) / coeff);
				}
			} else {
				const auto batch_approximators_scale = batch_approximators_scales[hyperplane_index - 1];
				total_scale_inv += coeff / batch_approximators_scale;
				
				const auto batch_approximators_l = batch_approximators_low[hyperplane_index - 1];
				const auto batch_approximators_h = batch_approximators_high[hyperplane_index - 1];
				
				const bool dir = coeff > 0;
				if (dir) {
					maximum_lower = std::max(maximum_lower, (batch_approximators_l - projected_state(i)) / coeff);
					minimum_upper = std::min(minimum_upper, (batch_approximators_h - projected_state(i)) / coeff);
				} else {
					minimum_upper = std::min(minimum_upper, (batch_approximators_l - projected_state(i)) / coeff);
					maximum_lower = std::max(maximum_lower, (batch_approximators_h - projected_state(i)) / coeff);
				}
			}
			
		}
		
		
		supp_t total_scale = 1.0 / total_scale_inv;
		
		if (maximum_lower >= minimum_upper) {
			return (delta_type) 0;
		}
		
		if (std::abs(total_scale_inv) < 1e-6) {
			if constexpr (std::is_integral<delta_type>::value) {
				const delta_type delta = (delta_type) (maximum_lower + gsl_rng_uniform_int(rng_local, minimum_upper - maximum_lower + 1));
				return delta;
			} else {
				const delta_type delta = (delta_type) (maximum_lower + (minimum_upper - maximum_lower) * gsl_rng_uniform(rng_local));
				return delta;
			}
		} else {
			const double exponential_sample = gsl_ran_exponential(rng_local, std::abs(total_scale));
			if constexpr (std::is_integral<delta_type>::value) {
				const delta_type aux = (delta_type) std::round(exponential_sample);
				const delta_type truncated_sample = aux % (minimum_upper - maximum_lower + 1);
				const delta_type delta = (total_scale > 0) ? truncated_sample + maximum_lower : minimum_upper - truncated_sample;
				return delta;
			} else {
				const delta_type truncated_sample = std::fmod(exponential_sample, minimum_upper - maximum_lower);
				const delta_type delta = (total_scale > 0) ? truncated_sample + maximum_lower : minimum_upper - truncated_sample;
				return delta;
			}
		}
	}
	
	
	auto randomVectorSample (size_t size, gsl_rng* const rng_local) const {
		assert(direction > 0);
		auto ans = new std::vector<supp_t>(size);
		if constexpr (std::is_integral<supp_t>::value) {
			if (direction > 0) {
				for (int i = 0; i < size; i++) {
					const supp_t exponential_sample = (supp_t) std::round(gsl_ran_exponential(rng_local, std::abs(scale))) % supp_size;
					ans->operator[](i) = (supp_t) (low + exponential_sample);
				}
			} else {
				assert(false);
				for (int i = 0; i < size; i++) {
					const supp_t exponential_sample = (supp_t) std::round(gsl_ran_exponential(rng_local, std::abs(scale)));
					ans->operator[](i) = (supp_t) (low - exponential_sample);
				}
			}
		} else {
			if (direction > 0) {
				for (int i = 0; i < size; i++) {
					const supp_t exponential_sample = std::fmod(gsl_ran_exponential(rng_local, std::abs(scale)), supp_size);
					ans->operator[](i) = (supp_t) (low + exponential_sample);
				}
			} else {
				assert(false);
				for (int i = 0; i < size; i++) {
					const supp_t exponential_sample = gsl_ran_exponential(rng_local, std::abs(scale));
					ans->operator[](i) = (supp_t) (low - exponential_sample);
				}
			}
		}
		return ans;
	}
	
	template <size_t size>
	auto batchApproximatorRandomVectorSample (const supp_t counter, gsl_rng* const rng_local) const {
		assert(direction > 0);
		assert(size == num_batch_approximators + 1);
		auto ans = new std::vector<supp_t>(size);
		const supp_t corrected_supp_size = std::min(counter, high) - low;
		if constexpr (std::is_integral<supp_t>::value) {
			if (direction > 0) {
				const supp_t exponential_sample = (supp_t) std::round(gsl_ran_exponential(rng_local, std::abs(scale))) % corrected_supp_size;
				ans->operator[](0) = (supp_t) (low + exponential_sample);
				for (int i = 1; i < size; i++) {
					const auto batch_approximators_l = batch_approximators_low[i - 1];
					const auto batch_approximators_h = batch_approximators_high[i - 1];
					const auto batch_approximators_scale = batch_approximators_scales[i - 1];
					const supp_t exponential_sample = (supp_t) std::round(gsl_ran_exponential(rng_local, std::abs(batch_approximators_scale))) % (batch_approximators_h - batch_approximators_l);
					if (batch_approximators_scale > 0) {
						ans->operator[](i) = (supp_t) (batch_approximators_l + exponential_sample);
					} else {
						ans->operator[](i) = (supp_t) (batch_approximators_h - exponential_sample);
					}
				}
			} else {
				assert(false);
			}
		} else {
			if (direction > 0) {
				const supp_t exponential_sample = std::fmod(gsl_ran_exponential(rng_local, std::abs(scale)), corrected_supp_size);
				ans->operator[](0) = (supp_t) (low + exponential_sample);
				for (int i = 1; i < size; i++) {
					const auto batch_approximators_l = batch_approximators_low[i - 1];
					const auto batch_approximators_h = batch_approximators_high[i - 1];
					const auto batch_approximators_scale = batch_approximators_scales[i - 1];
					const supp_t exponential_sample = (supp_t) std::fmod(gsl_ran_exponential(rng_local, std::abs(batch_approximators_scale)), batch_approximators_h - batch_approximators_l);
					if (batch_approximators_scale > 0) {
						ans->operator[](i) = (supp_t) (batch_approximators_l + exponential_sample);
					} else {
						ans->operator[](i) = (supp_t) (batch_approximators_h - exponential_sample);
					}
				}
			} else {
				assert(false);
			}
		}
		return ans;
	}
	
};

#endif /* AffineExponentialPrior_hpp */
