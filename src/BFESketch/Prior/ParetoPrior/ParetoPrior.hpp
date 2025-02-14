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

#ifndef ParetoPrior_hpp
#define ParetoPrior_hpp

#include "../../Sensor/LatentBaseVector/LatentBaseVector.hpp"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

// Prior that implements a pareto prior. Can be used both for level-1 and level-2 sampling.
// Assertion statements limit the mplmentation to cases relevant for the benhcmarking.
// Sampling from the density of a pdoruct fo linearly transformed pareto odstributions is by far the hardest amongst the pirors implemented in this library. In fact, there is no efficient exact way of doing this we know of so what we do is approximate the distribution to be samples by a piecewise-pareto distribution.
template <typename supp_type>
class ParetoPrior {
#define num_subdivision_intervals 3
	
//	The distribution defined by this object is equal to trunc(scale * (Pareto(1, shape) - 1) + low). Shape must be greater than 2 because otherwise variance and/or mean of teh distribution are infinite which causes issues with statistical stuff. Also scale must be greater than 0 becauzse we didnt fully implement the general case.
	double shape = 0;
	supp_type low = 0;
	supp_type high = 0;
	double scale = 0;
	supp_type supp_size = 0;
	int direction = 0;
	
	
	int num_batch_approximators = 0;
	double batch_approximators_shape = 0;
	std::vector<supp_type> batch_approximators_low;
	std::vector<supp_type> batch_approximators_high;
	std::vector<double> batch_approximators_scales;
	
	static constexpr auto node_helper_generator() {
		std::array<double, num_subdivision_intervals + 1> ans;
		for (int i = 0; i < num_subdivision_intervals + 1; i++) {
			ans[i] = i / ((double) num_subdivision_intervals);
		}
		return ans;
	}
	
	static constexpr std::array<double, num_subdivision_intervals + 1> node_helper = node_helper_generator();
	
public:
	using supp_t = supp_type;
	using intermediate_t = double;
	
	ParetoPrior () {
	}
	
	ParetoPrior(supp_t low, supp_t high, double scale, double shape): low(low), high(high), scale(scale), shape(shape), supp_size(high - low), direction((scale > 0) * 2 - 1) {
		assert(std::isfinite(scale) and std::abs(scale) > 1e-6);
		assert(std::isfinite(low));
		assert(std::isfinite(high));
		assert(shape > 2.0);
	}
	
	const auto variance_factor() const {
		return 10.0;
	}
	
	const auto mean(const supp_t counter) const {
		assert(direction > 0);
		const double corrected_high = std::min(counter, high);
		const double corrected_supp_size = corrected_high - low;
		if (direction > 0) {
			if (corrected_supp_size > 16 * scale) {
				return low + scale / (shape - 1);
			} else {
				const double aux1 = std::pow(scale / (scale + corrected_supp_size), shape);
				return low + ((scale + shape * corrected_supp_size) * aux1 - scale) / ((shape - 1) * (aux1 - 1));
			}
		} else {
			assert(false);
			return 0.0;
		}
	}
	
	const auto mean(const supp_t counter, double customShape) const {
		assert(direction > 0);
		const double corrected_high = std::min(counter, high);
		const double corrected_supp_size = corrected_high - low;
		if (direction > 0) {
			if (corrected_supp_size > 16 * scale) {
				return low + scale / (customShape - 1);
			} else {
				const double aux1 = std::pow(scale / (scale + corrected_supp_size), customShape);
				return low + ((scale + customShape * corrected_supp_size) * aux1 - scale) / ((customShape - 1) * (aux1 - 1));
			}
		} else {
			assert(false);
			return 0.0;
		}
	}
	
	template<typename array_type>
	auto normPlusBatchMean(const supp_t counter, array_type* const ans) const {
		assert(direction > 0);
		ans->resize(num_batch_approximators + 1);
		ans->operator()(0) = mean(counter);
		for (int i = 1; i < num_batch_approximators + 1; i++) {
			const auto scl = batch_approximators_scales[i - 1];
			const auto abs_scl = std::abs(scl);
			const auto supp_siz = batch_approximators_high[i - 1] - batch_approximators_low[i - 1];
			if (scl > 0) {
				if (supp_siz > 16 * scl) {
					ans->operator()(i) =  batch_approximators_low[i - 1] + scl / (batch_approximators_shape - 1);
				} else {
					const double aux1 = std::pow(scl / (scl + supp_siz), batch_approximators_shape);
					ans->operator()(i) =  batch_approximators_low[i - 1] + ((scl + batch_approximators_shape * supp_siz) * aux1 - scl) / ((batch_approximators_shape - 1) * (aux1 - 1));
				}
			} else {
				if (supp_siz > -16 * scl) {
					ans->operator()(i) =  batch_approximators_high[i - 1] + scl / (batch_approximators_shape - 1);
				} else {
					const double aux1 = std::pow(-scl / (-scl + supp_siz), batch_approximators_shape);
					ans->operator()(i) =  batch_approximators_high[i - 1] + -((-scl + batch_approximators_shape * supp_siz) * aux1 + scl) / ((batch_approximators_shape - 1) * (aux1 - 1));
				}
			}
		}
		return ans;
	}
	
	const auto variance(const supp_t counter, double customShape) const {
		assert(direction > 0);
		const double corrected_high = std::min(counter, high);
		const double corrected_supp_size = corrected_high - low;
		if (corrected_supp_size > 16 * scale) {
			return scale * scale * customShape / ((customShape - 1) * (customShape - 1) * (customShape - 2));
		} else {
			const double abs_scl = std::abs(scale);
			const double aux1 = std::pow(abs_scl / (abs_scl + corrected_supp_size), customShape);
			const double squaredexpectation = std::pow(((abs_scl + customShape * corrected_supp_size) * aux1 - abs_scl) / ((customShape - 1) * (aux1 - 1)), 2);
			const double expectationsquared = (aux1 * (2 * abs_scl * abs_scl + 2 * abs_scl * customShape * corrected_supp_size + (customShape - 1) * customShape * corrected_supp_size * corrected_supp_size) - 2 * abs_scl * abs_scl) / ((customShape - 2) * (customShape - 1) * (aux1 - 1));
			
			return expectationsquared - squaredexpectation;
		}
	}
	
	const auto variance(const supp_t counter) const {
		assert(direction > 0);
		const double corrected_high = std::min(counter, high);
		const double corrected_supp_size = corrected_high - low;
		if (corrected_supp_size > 16 * scale) {
			return scale * scale * shape / ((shape - 1) * (shape - 1) * (shape - 2));
		} else {
			const double abs_scl = std::abs(scale);
			const double aux1 = std::pow(abs_scl / (abs_scl + corrected_supp_size), shape);
			const double squaredexpectation = std::pow(((abs_scl + shape * corrected_supp_size) * aux1 - abs_scl) / ((shape - 1) * (aux1 - 1)), 2);
			const double expectationsquared = (aux1 * (2 * abs_scl * abs_scl + 2 * abs_scl * shape * corrected_supp_size + (shape - 1) * shape * corrected_supp_size * corrected_supp_size) - 2 * abs_scl * abs_scl) / ((shape - 2) * (shape - 1) * (aux1 - 1));
			
			return expectationsquared - squaredexpectation;
		}
	}
	
	template<typename array_type>
	auto normPlusBatchVar(const supp_t counter, array_type* const ans) const {
		ans->resize(num_batch_approximators + 1);
		ans->operator()(0) = variance(counter);
		for (int i = 1; i < num_batch_approximators + 1; i++) {
			const auto scl = batch_approximators_scales[i - 1];
			const auto abs_scl = std::abs(scl);
			const auto supp_siz = batch_approximators_high[i - 1] - batch_approximators_low[i - 1];
			
			if (supp_siz > 16 * abs_scl) {
				ans->operator()(i) = abs_scl * abs_scl * batch_approximators_shape / ((batch_approximators_shape - 1) * (batch_approximators_shape - 1) * (batch_approximators_shape - 2));
			} else {
				const double aux1 = std::pow(abs_scl / (abs_scl + supp_siz), batch_approximators_shape);
				const double squaredexpectation = std::pow(((abs_scl + batch_approximators_shape * supp_siz) * aux1 - abs_scl) / ((batch_approximators_shape - 1) * (aux1 - 1)), 2);
				const double expectationsquared = (aux1 * (2 * abs_scl * abs_scl + 2 * abs_scl * batch_approximators_shape * supp_siz + (batch_approximators_shape - 1) * batch_approximators_shape * supp_siz * supp_siz) - 2 * abs_scl * abs_scl) / ((batch_approximators_shape - 2) * (batch_approximators_shape - 1) * (aux1 - 1));
				ans->operator()(i) = expectationsquared - squaredexpectation;
			}
		}
		return ans;
	}
	
	void adjust_to_match_mean(const supp_type mean) {
		assert(supp_size > 16 * scale and mean < high and mean > low);
		if (direction > 0) {
			scale = (mean - low) * (shape - 1);
		} else {
			scale = (mean - high) * (shape - 1);
		}
	}
	
	void adjust_to_match_mean_var(const supp_type mean, const supp_type var) {
		assert(supp_size > 16 * scale and mean < high and mean > low);
		if (direction > 0) {
			scale = (mean - low) * (shape - 1);
		} else {
			scale = (mean - high) * (shape - 1);
		}
	}
	
	void set_batch_approximation(const int num_approximators, const double batchtarget_size, const supp_t batch_mean_, const supp_t batch_var_) {
		num_batch_approximators = num_approximators;
		assert(num_batch_approximators > 0);
		assert(direction > 0);
		
		supp_t batch_mean_local = batch_mean_;
		supp_t batch_std_local = std::sqrt( batch_var_);
		supp_t batch_mean = batch_mean_local;
		supp_t batch_var = batch_std_local * batch_std_local;
		batch_approximators_shape = std::sqrt(scale * scale / (batch_var)) + 2;
		const double balancers_scale = 1;
		const double balancers_offset = batch_mean / num_approximators - balancers_scale * mean(std::numeric_limits<supp_t>::max(), batch_approximators_shape);
		batch_approximators_scales = std::vector<supp_type>(num_approximators, scale * balancers_scale);
		batch_approximators_low = std::vector<supp_type>(num_approximators, low * balancers_scale + balancers_offset);
		batch_approximators_high = std::vector<supp_type>(num_approximators, high * balancers_scale + balancers_offset);
		batch_approximators_scales[0] *= -1;
		const double diff_low = batch_mean / num_approximators - batch_approximators_low[0];
		const double diff_high = batch_mean / num_approximators - batch_approximators_high[0];
		batch_approximators_high[0] = batch_approximators_low[0] + 2 * diff_low;
		batch_approximators_low[0] = batch_approximators_high[1] + 2 * diff_high;
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
		
		const auto support_size = minimum_upper - maximum_lower;
		const auto interval_size = support_size / (num_subdivision_intervals  * 1.0);
		
		const Eigen::Array<intermediate_t, num_subdivision_intervals + 1, 1> nodes = Eigen::Array<intermediate_t, num_subdivision_intervals + 1, 1>::LinSpaced((intermediate_t) maximum_lower, (intermediate_t) minimum_upper);
		Eigen::Array<intermediate_t, num_subdivision_intervals + 1, Eigen::Dynamic> linear_terms = Eigen::Array<intermediate_t, num_subdivision_intervals + 1, Eigen::Dynamic>::Zero(num_subdivision_intervals + 1, plane_size);
		
		for (int i = 0; i < plane_size; i++) {
			if (scale > 0) {
				linear_terms.col(i) = projected_state(i) + hyperplane_coeffs(i) * nodes + scale - low;
			} else {
				assert(false);
			}
		}
		
		linear_terms = -linear_terms.log();
		const auto individual_loglikelihoods = linear_terms * (shape + 1);
		Eigen::Array<intermediate_t, num_subdivision_intervals + 1, 1> loglikelihoods = individual_loglikelihoods.rowwise().sum();
		loglikelihoods = loglikelihoods - loglikelihoods.mean();
		auto inverse_likelihoods_depowered = (-loglikelihoods / (shape + 1)).exp();
		Eigen::Array<intermediate_t, num_subdivision_intervals, 1> linear_rates = (inverse_likelihoods_depowered.template tail<num_subdivision_intervals>() - inverse_likelihoods_depowered.template head<num_subdivision_intervals>()) / interval_size;
		Eigen::Array<intermediate_t, num_subdivision_intervals, 1> linear_offsets = inverse_likelihoods_depowered.template head<num_subdivision_intervals>() - linear_rates * nodes.template head<num_subdivision_intervals>();
		Eigen::Array<intermediate_t, num_subdivision_intervals, 1> partial_integrals = ((linear_offsets + linear_rates * nodes.head<num_subdivision_intervals>()).pow(-shape) - (linear_offsets + linear_rates * nodes.tail<num_subdivision_intervals>()).pow(-shape)) / (linear_rates * shape);
		const Eigen::Array<intermediate_t, num_subdivision_intervals, 1> partial_integrals_backup = interval_size * linear_offsets.pow(-shape - 1);
		partial_integrals = (linear_rates.abs() > 1e-6).select(partial_integrals, partial_integrals_backup);
		const auto total_volume = partial_integrals.sum();
		
		double selector_randomness = gsl_rng_uniform(rng_local) * total_volume;
		
		double sample;
		
		for (int sample_interval = 0; sample_interval < num_subdivision_intervals; sample_interval++) {
			if (selector_randomness <= partial_integrals(sample_interval) || sample_interval == num_subdivision_intervals - 1) {
				if (std::abs(linear_rates(sample_interval)) > 1e-6) {
					const double local_randomness = selector_randomness - std::pow(linear_offsets(sample_interval) + linear_rates(sample_interval) * nodes(sample_interval), -shape) / (linear_rates(sample_interval) * shape);
					sample = (std::pow(- local_randomness * linear_rates(sample_interval) * shape, - 1 / shape) - linear_offsets(sample_interval)) / linear_rates(sample_interval);
					
					assert(std::isfinite(sample));
				} else {
					sample = nodes(sample_interval) + selector_randomness * std::pow(linear_offsets(sample_interval), shape + 1);
					assert(std::isfinite(sample));
				}
				break;
			}
			selector_randomness -= partial_integrals(sample_interval);
		}
		if constexpr (std::is_integral<delta_type>::value) {
			const delta_type aux = (delta_type) std::round(sample);
			return aux;
		} else {
			return (delta_type) sample;
		}
	}
	
	template 	<	typename delta_type,
	typename coeff_type,
	typename index_type,
	typename state_type>
	const delta_type sampleHyperplaneBalanced(const LatentBaseVectorDynamic<coeff_type, delta_type, index_type>* const hyperplane,
											  const state_type* const state,
											  const supp_t counter,
											  gsl_rng* const rng_local) const {
		assert(direction > 0);
		
		const double corrected_high = std::min(counter, high);
		
		using ArrayXs = Eigen::Array<supp_t, Eigen::Dynamic, 1>;
		using ArrayXh = Eigen::Array<coeff_type, Eigen::Dynamic, 1>;
		
		const auto projected_state = hyperplane->projectState(state);

		const auto hyperplane_coeffs = *hyperplane->get_coeff_map();

		const auto plane_size = hyperplane->get_size();
		
		assert(plane_size == 2);
		delta_type maximum_lower = std::numeric_limits<delta_type>::lowest();
		delta_type minimum_upper = std::numeric_limits<delta_type>::max();
		for (int i = 0; i < plane_size; i++) {
			const auto coeff = hyperplane_coeffs(i);

			if (hyperplane->get_index_map()->operator()(i) == 0) {
				assert(i == 0);
				const bool dir = coeff > 0;
				
				if (dir) {
					maximum_lower = std::max(maximum_lower, (low - projected_state(i)) / coeff);
					minimum_upper = std::min(minimum_upper, (corrected_high - projected_state(i)) / coeff);
				} else {
					minimum_upper = std::min(minimum_upper, (low - projected_state(i)) / coeff);
					maximum_lower = std::max(maximum_lower, (corrected_high - projected_state(i)) / coeff);
				}
			} else {
				const auto batch_approximators_scale = batch_approximators_scales[hyperplane->get_index_map()->operator()(i) - 1];
				const auto batch_approximators_l = batch_approximators_low[hyperplane->get_index_map()->operator()(i) - 1];
				const auto batch_approximators_h = batch_approximators_high[hyperplane->get_index_map()->operator()(i) - 1];
				
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

		const auto support_size = minimum_upper - maximum_lower;
		const auto interval_size = support_size / (num_subdivision_intervals  * 1.0);
		
		const Eigen::Array<intermediate_t, num_subdivision_intervals + 1, 1> nodes = Eigen::Array<intermediate_t, num_subdivision_intervals + 1, 1>::LinSpaced((intermediate_t) maximum_lower, (intermediate_t) minimum_upper);
		Eigen::Array<intermediate_t, num_subdivision_intervals + 1, 2> linear_terms = Eigen::Array<intermediate_t, num_subdivision_intervals + 1, Eigen::Dynamic>::Zero(num_subdivision_intervals + 1, 2);
		
		int linear_terms_total_dir = 0;
		Eigen::Array<intermediate_t, 2, 1> shapearray;
		double localshapes_left_right[2] = {0};
		double localshapes_cumulative = 0;
		double localscales_array[2];
		for (int i = 0; i < plane_size; i++) {
			if (hyperplane->get_index_map()->operator()(i) == 0) {
				shapearray(i) = shape ;
				if (scale > 0) {
					localshapes_left_right[0] = (shape + 1) * (hyperplane_coeffs(i) > 0);
					localshapes_left_right[1] = (shape + 1) * (hyperplane_coeffs(i) < 0);
					localshapes_cumulative += shape + 1;
					localscales_array[i] = scale;
					linear_terms.col(i) = projected_state(i) + hyperplane_coeffs(i) * nodes + scale - low;
					linear_terms_total_dir += hyperplane_coeffs(i);
				} else {
					assert(false);
				}
			} else {
				shapearray(i) = batch_approximators_shape ;
				const auto batch_approximators_scale = batch_approximators_scales[hyperplane->get_index_map()->operator()(i) - 1];
				const auto batch_approximators_l = batch_approximators_low[hyperplane->get_index_map()->operator()(i) - 1];
				const auto batch_approximators_h = batch_approximators_high[hyperplane->get_index_map()->operator()(i) - 1];
				localscales_array[i] = std::abs(batch_approximators_scale);
				if (batch_approximators_scale > 0) {
					localshapes_left_right[0] += (batch_approximators_shape + 1) * (hyperplane_coeffs(i) > 0);
					localshapes_left_right[1] += (batch_approximators_shape + 1) * (hyperplane_coeffs(i) < 0);
					localshapes_cumulative += batch_approximators_shape + 1;
					linear_terms.col(i) = projected_state(i) + hyperplane_coeffs(i) * nodes + batch_approximators_scale - batch_approximators_l;
					linear_terms_total_dir += hyperplane_coeffs(i);
				} else {
					localshapes_left_right[1] += (batch_approximators_shape + 1) * (hyperplane_coeffs(i) > 0);
					localshapes_left_right[0] += (batch_approximators_shape + 1) * (hyperplane_coeffs(i) < 0);
					localshapes_cumulative += batch_approximators_shape + 1;
					linear_terms.col(i) = -(projected_state(i) + hyperplane_coeffs(i) * nodes) - batch_approximators_scale + batch_approximators_h;
					linear_terms_total_dir -= hyperplane_coeffs(i);
				}
			}
		}
		double localscale_cumulative = std::pow(localscales_array[0], (shapearray(0) + 1) / localshapes_cumulative) * std::pow(localscales_array[1], (shapearray(1) + 1) / localshapes_cumulative);
		double sample;
		Eigen::Array<intermediate_t, num_subdivision_intervals, 1> localshapes;
		Eigen::Array<intermediate_t, num_subdivision_intervals + 1, 1> localshapesExtended;
		
		if (linear_terms_total_dir != 0) {
			localshapes = Eigen::Array<intermediate_t, num_subdivision_intervals, 1>::Constant(localshapes_cumulative - 1);
			localshapesExtended = Eigen::Array<intermediate_t, num_subdivision_intervals + 1, 1>::Constant(localshapes_cumulative - 1);
		} else {
			localshapes = Eigen::Array<intermediate_t, num_subdivision_intervals, 1>::LinSpaced((intermediate_t) localshapes_left_right[0] - 1, (intermediate_t) localshapes_left_right[1] - 1);
			localshapesExtended = Eigen::Array<intermediate_t, num_subdivision_intervals + 1, 1>::LinSpaced((intermediate_t) localshapes_left_right[0] - 1, (intermediate_t) localshapes_left_right[1] - 1);
		}
		
		if (support_size > 1e+6) {
			sample = std::fmod((gsl_ran_pareto(rng_local, localshapes_cumulative - 1, 1) - 1) * localscale_cumulative, support_size);
			if (linear_terms_total_dir > 0) {
				sample += maximum_lower;
			} else {
				sample = minimum_upper - sample;
			}
		} else {
			linear_terms = -linear_terms.log();
			const auto individual_loglikelihoods = linear_terms.rowwise() * (shapearray.transpose() + 1);
			Eigen::Array<intermediate_t, num_subdivision_intervals + 1, 1> loglikelihoods = individual_loglikelihoods.rowwise().sum();
			loglikelihoods = loglikelihoods - loglikelihoods.mean();
			auto inverse_likelihoods_depowered = (-loglikelihoods / (localshapesExtended + 1)).exp();
			Eigen::Array<intermediate_t, num_subdivision_intervals, 1> linear_rates = (inverse_likelihoods_depowered.tail<num_subdivision_intervals>() - inverse_likelihoods_depowered.head<num_subdivision_intervals>()) / interval_size;
			Eigen::Array<intermediate_t, num_subdivision_intervals, 1> linear_offsets = inverse_likelihoods_depowered.head<num_subdivision_intervals>() - linear_rates * nodes.head<num_subdivision_intervals>();
			Eigen::Array<intermediate_t, num_subdivision_intervals, 1> partial_integrals = ((linear_offsets + linear_rates * nodes.head<num_subdivision_intervals>()).pow(-localshapes) - (linear_offsets + linear_rates * nodes.tail<num_subdivision_intervals>()).pow(-localshapes)) / (linear_rates * localshapes);
			const Eigen::Array<intermediate_t, num_subdivision_intervals, 1> partial_integrals_backup = interval_size * linear_offsets.pow(-localshapes - 1);
			partial_integrals = (linear_rates.abs() > 1e-6).select(partial_integrals, partial_integrals_backup);
			const auto total_volume = partial_integrals.sum();
			
			double selector_randomness = gsl_rng_uniform(rng_local) * total_volume;
			
			for (int sample_interval = 0; sample_interval < num_subdivision_intervals; sample_interval++) {
				if (selector_randomness <= partial_integrals(sample_interval) || sample_interval == num_subdivision_intervals - 1) {
					if (std::abs(linear_rates(sample_interval)) > 1e-6) {
						const double local_randomness = selector_randomness - std::pow(linear_offsets(sample_interval) + linear_rates(sample_interval) * nodes(sample_interval), -localshapes(sample_interval)) / (linear_rates(sample_interval) * localshapes(sample_interval));
						sample = (std::pow(- local_randomness * linear_rates(sample_interval) * localshapes(sample_interval), - 1 / localshapes(sample_interval)) - linear_offsets(sample_interval)) / linear_rates(sample_interval);
						
						assert(std::isfinite(sample));
					} else {
						sample = nodes(sample_interval) + selector_randomness * std::pow(linear_offsets(sample_interval), localshapes(sample_interval) + 1);
						assert(std::isfinite(sample));
					}
					break;
				}
				selector_randomness -= partial_integrals(sample_interval);
			}
		}
		
		if constexpr (std::is_integral<delta_type>::value) {
			const delta_type aux = (delta_type) std::round(sample);
			return aux;
		} else {
			return (delta_type) sample;
		}
	}
	
	
	auto randomVectorSample (size_t size, gsl_rng* const rng_local) const {
		assert(supp_size > 16 * scale);
		auto ans = new std::vector<supp_t>(size);
		if constexpr (std::is_integral<supp_t>::value) {
			if (direction > 0) {
				for (int i = 0; i < size; i++) {
					const supp_t aux_sample = (supp_t) std::round(gsl_ran_pareto(rng_local, shape, scale) - scale);
					ans->operator[](i) = (supp_t) (low + aux_sample);
				}
			} else {
				for (int i = 0; i < size; i++) {
					const supp_t aux_sample = (supp_t) std::round(gsl_ran_pareto(rng_local, shape, std::abs(scale)) - std::abs(scale));
					ans->operator[](i) = (supp_t) (high - aux_sample);
				}
			}
		} else {
			if (direction > 0) {
				for (int i = 0; i < size; i++) {
					const supp_t aux_sample = (supp_t) gsl_ran_pareto(rng_local, shape, scale) - scale;
					ans->operator[](i) = (supp_t) (low + aux_sample);
				}
			} else {
				for (int i = 0; i < size; i++) {
					const supp_t aux_sample = (supp_t) gsl_ran_pareto(rng_local, shape, std::abs(scale)) - std::abs(scale);
					ans->operator[](i) = (supp_t) (high - aux_sample);
				}
			}
		}
		return ans;
	}
};

#endif /* ParetoPrior_hpp */
