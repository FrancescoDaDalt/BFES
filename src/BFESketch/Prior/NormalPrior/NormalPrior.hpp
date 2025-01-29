//
//  NormalPrior.hpp
//  braindust
//
//  Created by Francesco Da Dalt on 16.05.24.
//

#ifndef NormalPrior_hpp
#define NormalPrior_hpp

#include "../../Sensor/LatentBaseVector/LatentBaseVector.hpp"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

// This class implements a truncated normal prior. Can be used both as level-1 or level-2 prior.
// This class is not commented as there are too many technical details which cannot be fully explained in comments
// There are assertions that limit the flexibility of the prior to for example the case where the lower truncation bound of the prior must be greater or equal to 0. Otherwise the code needs to be augmented to handle the other case. For benchmarking however we onlyever assume positive item frequencies so we have limited the implementation to this case.
template <typename supp_type>
class NormalPrior {
	
//	The distribution defined by this prior is equal to trunc(low + mu + sigma * N(0,1), [low, high]) where N(0,1) is the standard normal distrib, if sigma > 0.
//	When the object is initialized, it stores the ratio of mu to sigma and keep this ratio constant throughout the computations. This means after initialization the distribution behaves like a 1-parameter distribution
//		which allows us to fit its parameters to some mean
	double mu_to_sigma_ratio = 1;
	supp_type low = 0;
	supp_type high = 0;
	supp_type sigma = 0;
	supp_type mu = 0;
	supp_type supp_size = 0;
	
	
	int num_batch_approximators = 0;
	std::vector<supp_type> batch_approximators_low;
	std::vector<supp_type> batch_approximators_high;
	std::vector<supp_type> batch_approximators_scales;
	
public:
	using supp_t = supp_type;
	using intermediate_t = double;
	
	NormalPrior() {
	}
	
	NormalPrior(supp_t low,
				supp_t high,
				supp_t mu,
				supp_t sigma):
	low(low),
	high(high),
	sigma(sigma),
	mu(mu),
	supp_size(high - low),
	mu_to_sigma_ratio((mu - low) /sigma) {
		assert(std::isfinite(sigma) and std::abs(sigma) > 1e-6);
		assert(std::isfinite(low));
		assert(std::isfinite(high));
		assert(low >=0);
		assert(sigma >= 0);
	}
	
	
	const auto variance_factor() const {
		return 1;
	}
	
	const std::pair<supp_t, supp_t> get_bounds() const {
		return std::make_pair(low, high);
	}
	
	const auto get_num_batch_approximators() const {
		return num_batch_approximators;
	}
	
	const auto mean(const supp_t counter) const {
		const double corrected_high = std::min(counter, high);
		double a = (low - mu) / sigma;
		double b = (corrected_high - mu) / sigma;
		double Z = gsl_cdf_ugaussian_P(b) - gsl_cdf_ugaussian_P(a);
		
		double correctedmean = mu + (gsl_ran_ugaussian_pdf(a) - gsl_ran_ugaussian_pdf(b)) * sigma / Z;
		return correctedmean;
	}
	
	template<typename array_type>
	auto normPlusBatchMean(const supp_t counter, array_type* const ans) const {
		ans->resize(num_batch_approximators + 1);
		ans->operator()(0) = mean(counter);
		for (int i = 1; i < num_batch_approximators + 1; i++) {
			double sigma_local = batch_approximators_scales[i - 1];
			double mu_local = sigma_local * mu_to_sigma_ratio;
			if (sigma_local > 0) {
				mu_local += batch_approximators_low[i - 1];
			} else {
				mu_local += batch_approximators_high[i - 1];
			}
			double a = (batch_approximators_low[i - 1] - mu_local) / std::abs(sigma_local);
			double b = (batch_approximators_high[i - 1] - mu_local) / std::abs(sigma_local);
			double Z = gsl_cdf_ugaussian_P(b) - gsl_cdf_ugaussian_P(a);
			
			double correctedmean = mu_local + (gsl_ran_ugaussian_pdf(a) - gsl_ran_ugaussian_pdf(b)) * sigma_local / Z;
			ans->operator()(i) = correctedmean;
		}
	}
	
	const auto variance(const supp_t counter) const {
		const double corrected_high = std::min(counter, high);
		double a = (low - mu) / sigma;
		double b = (corrected_high - mu) / sigma;
		double Z = gsl_cdf_ugaussian_P(b) - gsl_cdf_ugaussian_P(a);
		
		double correctedvar = sigma * sigma * (1 + (a * gsl_ran_ugaussian_pdf(a) - b * gsl_ran_ugaussian_pdf(b)) / Z - std::pow((gsl_ran_ugaussian_pdf(a) - gsl_ran_ugaussian_pdf(b)) / Z, 2.0));
		return correctedvar;
	}
	
	template<typename array_type>
	auto normPlusBatchVar(const supp_t counter, array_type* const ans) const {
		ans->resize(num_batch_approximators + 1);
		ans->operator()(0) = variance(counter);
		for (int i = 1; i < num_batch_approximators + 1; i++) {
			double sigma_local = batch_approximators_scales[i - 1];
			double mu_local = sigma_local * mu_to_sigma_ratio;
			if (sigma_local > 0) {
				mu_local += batch_approximators_low[i - 1];
			} else {
				mu_local += batch_approximators_high[i - 1];
			}
			double a = (batch_approximators_low[i - 1] - mu_local) / std::abs(sigma_local);
			double b = (batch_approximators_high[i - 1] - mu_local) / std::abs(sigma_local);
			double Z = gsl_cdf_ugaussian_P(b) - gsl_cdf_ugaussian_P(a);
			
			double correctedvar = sigma_local * sigma_local * (1 + (a * gsl_ran_ugaussian_pdf(a) - b * gsl_ran_ugaussian_pdf(b)) / Z - std::pow((gsl_ran_ugaussian_pdf(a) - gsl_ran_ugaussian_pdf(b)) / Z, 2.0));
			
			ans->operator()(i) = correctedvar;
		}
	}
	
	void adjust_to_match_mean_var(const supp_type mean, const supp_type var) {
		assert(mean < high and mean > low);
		supp_t batch_mean_local = mean;
		supp_t batch_std_local = std::sqrt(var);
		for (int i = 0; i < 50; i++) {
			double a = (low - batch_mean_local) / batch_std_local;
			double b = (high - batch_mean_local) / batch_std_local;
			double Z = gsl_cdf_ugaussian_P(b) - gsl_cdf_ugaussian_P(a);
			
			if (a > 3 and b > 3) {
				double correctedmean = batch_mean_local + gsl_ran_ugaussian_tail_pdf(a, a) * batch_std_local;
				batch_mean_local += mean - correctedmean;
				a = (low - batch_mean_local) / batch_std_local;
				b = (high - batch_mean_local) / batch_std_local;
				
				double correctedvar = batch_std_local * batch_std_local * (1 + a * gsl_ran_ugaussian_tail_pdf(a, a) - std::pow(gsl_ran_ugaussian_tail_pdf(a, a), 2.0));
				batch_std_local *= std::sqrt(std::sqrt(var) / std::sqrt(correctedvar));
			} else if (a < -3 and b < -3) {
				assert(false);
			} else if (Z > 1e-6) {
				double correctedmean = batch_mean_local + (gsl_ran_ugaussian_pdf(a) - gsl_ran_ugaussian_pdf(b)) * batch_std_local / Z;
				batch_mean_local += mean - correctedmean;
				a = (low - batch_mean_local) / batch_std_local;
				b = (high - batch_mean_local) / batch_std_local;
				Z = gsl_cdf_ugaussian_P(b) - gsl_cdf_ugaussian_P(a);
				
				double correctedvar = batch_std_local * batch_std_local * (1 + (a * gsl_ran_ugaussian_pdf(a) - b * gsl_ran_ugaussian_pdf(b)) / Z - std::pow((gsl_ran_ugaussian_pdf(a) - gsl_ran_ugaussian_pdf(b)) / Z, 2.0));
				batch_std_local *= std::sqrt(std::sqrt(var) / std::sqrt(correctedvar));
			} else {
				assert(false);
			}
		}
		
		mu = batch_mean_local;
		sigma = batch_std_local;
	}
	
	void adjust_to_match_mean(const supp_type mean) {
		assert(mean < high and mean > low and sigma > 0);
		supp_t batch_std_local = (mean - low) / mu_to_sigma_ratio;
		supp_t batch_mean_local = mean;
		for (int i = 0; i < 50; i++) {
			double a = (low - batch_mean_local) / batch_std_local;
			double b = (high - batch_mean_local) / batch_std_local;
			double Z = gsl_cdf_ugaussian_P(b) - gsl_cdf_ugaussian_P(a);
			
			if (a > 3 and b > 3) {
				double correctedmean = batch_mean_local + gsl_ran_ugaussian_tail_pdf(a, a) * batch_std_local;
				batch_mean_local += mean - correctedmean;
				batch_std_local = (batch_mean_local - low) / mu_to_sigma_ratio;
			} else if (a < -3 and b < -3) {
				assert(false);
			} else if (Z > 1e-6) {
				double correctedmean = batch_mean_local + (gsl_ran_ugaussian_pdf(a) - gsl_ran_ugaussian_pdf(b)) * batch_std_local / Z;
				batch_mean_local += mean - correctedmean;
				batch_std_local = (batch_mean_local - low) / mu_to_sigma_ratio;
			} else {
				assert(false);
			}
		}
		mu = batch_mean_local;
		sigma = batch_std_local;
	}
	
	void set_batch_approximation(const int num_approximators,
								 const double batchtarget_size,
								 const supp_t batch_mean_,
								 const supp_t batch_var_) {
		num_batch_approximators = num_approximators;
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
		
		batch_approximators_scales = std::vector<supp_type>(num_approximators, sigma * balancers_scale);
		batch_approximators_low = std::vector<supp_type>(num_approximators, low * balancers_scale + balancers_offset);
		batch_approximators_high = std::vector<supp_type>(num_approximators, high * balancers_scale + balancers_offset);
		batch_approximators_scales[0] *= -1;
		const double diff_low = batch_mean / num_approximators - batch_approximators_low[0];
		const double diff_high = batch_mean / num_approximators - batch_approximators_high[0];
		batch_approximators_high[0] = batch_approximators_low[0] + 2 * diff_low;
		batch_approximators_low[0] = batch_approximators_high[1] + 2 * diff_high;
	}
	
	template 	<	typename delta_type,
	typename coeff_type,
	typename index_type,
	typename state_type>
	const delta_type sampleHyperplane(const LatentBaseVectorDynamic<coeff_type, delta_type, index_type>* const hyperplane, const state_type* const state, gsl_rng* const rng_local) const {
		using ArrayXs = Eigen::Array<supp_t, Eigen::Dynamic, 1>;
		using ArrayXh = Eigen::Array<coeff_type, Eigen::Dynamic, 1>;
		
		const auto projected_state = hyperplane->projectState(state);
		
		const auto hyperplane_coeffs = *hyperplane->get_coeff_map();
		
		const auto plane_size = hyperplane->get_size();
		
		delta_type maximum_lower = std::numeric_limits<delta_type>::lowest();
		delta_type minimum_upper = std::numeric_limits<delta_type>::max();
		supp_t total_mean = (supp_t) 0;
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
			total_mean += (mu - projected_state(i)) / coeff;
		}
		total_mean /= plane_size;
		supp_t total_var = (sigma * sigma) / plane_size;
		supp_t total_stddev = std::sqrt(total_var);
		
		double a = (maximum_lower - total_mean) / total_stddev;
		double b = (minimum_upper - total_mean) / total_stddev;
		
		if (a > 3 and b > 3) {
			double sample = (std::fmod(gsl_ran_ugaussian_tail(rng_local, a) - a, b - a) + a) * total_stddev + total_mean;
//			assert(sample >= maximum_lower and sample <= minimum_upper);
			return sample;
		} else if (a < -3 and b < -3) {
			double sample = (std::fmod(gsl_ran_ugaussian_tail(rng_local, -b) + b, b - a) - b) * -total_stddev + total_mean;
//			assert(sample >= maximum_lower and sample <= minimum_upper);
			return sample;
		} else {
			
			double random_unif_sample = gsl_rng_uniform(rng_local);
			double truncated_gaussian_sample = total_mean + total_stddev * gsl_cdf_ugaussian_Pinv(gsl_cdf_ugaussian_P(a) + random_unif_sample * (gsl_cdf_ugaussian_P(b) - gsl_cdf_ugaussian_P(a)));
//			assert(truncated_gaussian_sample >= maximum_lower and truncated_gaussian_sample <= minimum_upper);
			return truncated_gaussian_sample;
		}
	}
	
	template 	<	typename delta_type,
	typename coeff_type,
	typename index_type,
	typename state_type>
	const delta_type sampleHyperplaneBalanced(const LatentBaseVectorDynamic<coeff_type, delta_type, index_type>* const hyperplane, const state_type* const state, const supp_t counter, gsl_rng* const rng_local) const {
		assert(sigma > 0);
		const double corrected_high = std::min(counter, high);
		
		
		using ArrayXs = Eigen::Array<supp_t, Eigen::Dynamic, 1>;
		using ArrayXh = Eigen::Array<coeff_type, Eigen::Dynamic, 1>;
		
		const auto projected_state = hyperplane->projectState(state);
		
		const auto hyperplane_coeffs = *hyperplane->get_coeff_map();
		
		const auto plane_size = hyperplane->get_size();
		
		delta_type maximum_lower = std::numeric_limits<delta_type>::lowest();
		delta_type minimum_upper = std::numeric_limits<delta_type>::max();
		supp_t total_mean = (supp_t) 0;
		double total_var_inv = 0;
		const auto him = hyperplane->get_index_map();
		for (int i = 0; i < plane_size; i++) {
			const auto coeff = hyperplane_coeffs(i);
			const auto hyperplane_index = him->operator()(i);
			const bool dir = coeff > 0;
			if (hyperplane_index == 0) {
				if (dir) {
					maximum_lower = std::max(maximum_lower, (low - projected_state(i)) / coeff);
					minimum_upper = std::min(minimum_upper, (corrected_high - projected_state(i)) / coeff);
				} else {
					minimum_upper = std::min(minimum_upper, (low - projected_state(i)) / coeff);
					maximum_lower = std::max(maximum_lower, (corrected_high - projected_state(i)) / coeff);
				}
				total_mean += (mu - projected_state(i)) / (coeff * sigma * sigma);
				total_var_inv += 1.0 / (sigma * sigma);
			} else {
				const auto batch_approximators_scale = batch_approximators_scales[hyperplane_index - 1];
				double sigma_local = batch_approximators_scale;
				double mu_local = sigma_local * mu_to_sigma_ratio;
				if (sigma_local > 0) {
					mu_local += batch_approximators_low[hyperplane_index - 1];
				} else {
					mu_local += batch_approximators_high[hyperplane_index - 1];
				}
				const auto batch_approximators_l = batch_approximators_low[hyperplane_index - 1];
				const auto batch_approximators_h = batch_approximators_high[hyperplane_index - 1];
				if (dir) {
					maximum_lower = std::max(maximum_lower, (batch_approximators_l - projected_state(i)) / coeff);
					minimum_upper = std::min(minimum_upper, (batch_approximators_h - projected_state(i)) / coeff);
				} else {
					minimum_upper = std::min(minimum_upper, (batch_approximators_l - projected_state(i)) / coeff);
					maximum_lower = std::max(maximum_lower, (batch_approximators_h - projected_state(i)) / coeff);
				}
				total_mean += (mu_local - projected_state(i)) / (coeff * sigma_local * sigma_local);
				total_var_inv += 1.0 / (sigma_local * sigma_local);
			}
		}
		total_mean /= total_var_inv;
		supp_t total_var = 1 / total_var_inv;
		supp_t total_stddev = std::sqrt(total_var);
		
		double a = (maximum_lower - total_mean) / total_stddev;
		double b = (minimum_upper - total_mean) / total_stddev;
		
		if (a > 3 and b > 3) {
			double sample = std::fmod(gsl_ran_ugaussian_tail(rng_local, a), b - a) * total_stddev + total_mean;
			return sample;
		} else if (a < -3 and b < -3) {
			double sample = std::fmod(gsl_ran_ugaussian_tail(rng_local, -b), b - a) * -total_stddev + total_mean;
			return sample;
		} else {
			
			double random_unif_sample = gsl_rng_uniform(rng_local);
			double truncated_gaussian_sample = total_mean + total_stddev * gsl_cdf_ugaussian_Pinv(gsl_cdf_ugaussian_P(a) + random_unif_sample * (gsl_cdf_ugaussian_P(b) - gsl_cdf_ugaussian_P(a)));
			return truncated_gaussian_sample;
		}
	}
	
	auto randomVectorSample (size_t size, gsl_rng* const rng_local) const {
		double a = (0 - mu) / sigma;
		double b = (high - low - mu) / sigma;
		auto ans = new std::vector<supp_t>(size);
		if constexpr (std::is_integral<supp_t>::value) {
			assert(false);
		} else {
			for (int i = 0; i < size; i++) {
				if (a > 3 and b > 3) {
					double sample = std::fmod(gsl_ran_ugaussian_tail(rng_local, a), b - a) * sigma + mu + low;
					ans->operator[](i) = sample;
				} else if (a < -3 and b < -3) {
					double sample = std::fmod(gsl_ran_ugaussian_tail(rng_local, -b), b - a) * -sigma + mu + low;
					ans->operator[](i) = sample;
				} else {
					double random_unif_sample = gsl_rng_uniform(rng_local);
					double truncated_gaussian_sample = mu + sigma * gsl_cdf_ugaussian_Pinv(gsl_cdf_ugaussian_P(a) + random_unif_sample * (gsl_cdf_ugaussian_P(b) - gsl_cdf_ugaussian_P(a)));
					ans->operator[](i) = truncated_gaussian_sample + low;
				}
			}
		}
		return ans;
	}
};

#endif /* NormalPrior_hpp */
