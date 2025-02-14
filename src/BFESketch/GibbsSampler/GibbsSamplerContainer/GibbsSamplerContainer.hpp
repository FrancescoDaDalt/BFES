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



#ifndef GibbsSamplerContainer_hpp
#define GibbsSamplerContainer_hpp

#include <random>

// Container for the level-1 sampling. Multithreaded
template <size_t num_threads, typename gse_type>
class GibbsSamplerContainer {
	const size_t max_num_samples_per_thread;
	const size_t num_items;
	const gse_type sampler;
	std::array<gsl_rng*, num_threads> rng;
	std::array<Eigen::Array<typename gse_type::supp_t, Eigen::Dynamic, Eigen::Dynamic>, num_threads> buffers;
	std::array<size_t, num_threads> num_samples_in_buffers;
	
public:
	GibbsSamplerContainer (const gse_type* const sampler,
						   const size_t max_num_samples):
	max_num_samples_per_thread(max_num_samples / num_threads),
	num_items(sampler->target_dim),
	sampler(*sampler) {
		for (int i = 0; i < num_threads; i++) {
			num_samples_in_buffers[i] = 0;
			buffers[i] = Eigen::Array<typename gse_type::supp_t, Eigen::Dynamic, Eigen::Dynamic>::Zero(sampler->target_dim,
																									   max_num_samples_per_thread);
			rng[i] = gsl_rng_alloc(gsl_rng_mt19937);
			std::random_device rd;
			gsl_rng_set(rng[i],
						rd());
		}
	}
	
	~GibbsSamplerContainer () {
		for (auto r: rng) {
			gsl_rng_free(r);
		}
	}
	
	auto get_rng (const size_t idx) const {
		return rng[idx];
	}
	
	template<typename mattype>
	auto refresh (const mattype* const initial_states) {
		for (int i = 0; i < num_threads; i++) {
			num_samples_in_buffers[i] = 1;
			buffers[i].col(0) = initial_states->col(i);
		}
	}
	
	template<typename prior_type>
	auto burnin (const prior_type* const prior,
				 const int burning_num) {
		#pragma omp parallel for
		for (int i = 0; i < num_threads; i++) {
			const unsigned long idx = std::max(num_samples_in_buffers[i] - 1, (unsigned long) 0);
			for (int j = 0; j < burning_num; j++) {
				const auto src = buffers[i].col(idx);
				auto trgt = buffers[i].col(idx + 1);
				sampler.sampleNextState_external(prior,
												 &src,
												 &trgt,
												 rng[i]);
				buffers[i].col(idx) = buffers[i].col(idx + 1);
			}
		}
	}
	
	template<typename prior_type>
	auto sample (const prior_type* const prior,
				 const int interleaving,
				 const int numsamples_) {
		const int numsamples = std::min((int) (numsamples_ / num_threads), (int) max_num_samples_per_thread);
		
		#pragma omp parallel for
		for (int i = 0; i < num_threads; i++) {
			const unsigned long idx = std::max(num_samples_in_buffers[i], (unsigned long) 1);
			for (int j = idx; j < numsamples; j++) {
				for (int k = 0; k < interleaving - 1; k++) {
					const auto src = buffers[i].col(j - 1);
					auto trgt = buffers[i].col(j);
					sampler.sampleNextState_external(prior,
														 &src,
														 &trgt,
														 rng[i]);
					buffers[i].col(j - 1) = buffers[i].col(j);
				}
				const auto src = buffers[i].col(j - 1);
				auto trgt = buffers[i].col(j);
				sampler.sampleNextState_external(prior,
													 &src,
													 &trgt,
													 rng[i]);
				num_samples_in_buffers[i]++;
			}
		}
	}
	
	auto get_samples (Eigen::Array<typename gse_type::supp_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>* const ans) const {
		const size_t num_samples_total = std::accumulate(num_samples_in_buffers.cbegin(),
														 num_samples_in_buffers.cend(),
														 0);
		std::array<size_t, num_threads + 1> indices;
		indices[0] = 0;
		std::partial_sum(num_samples_in_buffers.cbegin(),
						 num_samples_in_buffers.cend(),
						 indices.begin() + 1);
		
		ans->resize(num_items, num_samples_total);
		
		#pragma omp parallel for
		for (int i = 0; i < num_threads; i++) {
			const size_t begin_idx = indices[i];
			const size_t end_idx = indices[i + 1];
			ans->operator()(Eigen::all, Eigen::seq(begin_idx, end_idx - 1, 1)) = buffers[i].leftCols(end_idx - begin_idx);
		}
		return ans;
	}
};

// Container for the level-1 sampling. Not multithreaded because we do multithreadeing at a bigger granularity to reduce thread overhead
template <typename gse_type>
class GibbsSamplerContainerSingleItem {
	const size_t max_num_samples;
	const size_t num_items;
	const gse_type sampler;
	gsl_rng* rng;
	Eigen::Array<typename gse_type::supp_t, Eigen::Dynamic, Eigen::Dynamic> buffer;
	size_t num_samples_in_buffer;
	
public:
	GibbsSamplerContainerSingleItem (const gse_type* const sampler,
									 const size_t max_num_samples):
	max_num_samples(max_num_samples),
	num_items(sampler->target_dim),
	sampler(*sampler) {
		num_samples_in_buffer = 0;
		buffer = Eigen::Array<typename gse_type::supp_t, Eigen::Dynamic, Eigen::Dynamic>::Zero(sampler->target_dim, max_num_samples);
		rng = gsl_rng_alloc(gsl_rng_mt19937);
		std::random_device rd;
		gsl_rng_set(rng, rd());
	}
	
	GibbsSamplerContainerSingleItem ():
	max_num_samples(0),
	num_items(0),
	rng(nullptr),
	num_samples_in_buffer(0) {}
	
	~GibbsSamplerContainerSingleItem () {
		gsl_rng_free(rng);
	}
	
	auto get_rng () const {
		return rng;
	}
	
	template<typename mattype>
	auto refresh (const mattype* const initial_states) {
		num_samples_in_buffer = 1;
		buffer.col(0) = *initial_states;
		return;
	}
	
	template<typename prior_type>
	auto burnin (const typename gse_type::supp_t counter,
				 const prior_type* const prior,
				 const int burnin_num) {
		for (int j = 0; j < burnin_num; j++) {
			const auto src = buffer.col(num_samples_in_buffer - 1);
			auto trgt = buffer.col(num_samples_in_buffer);
			sampler.sampleNextState_external(prior,
											 &src,
											 &trgt,
											 counter,
											 rng);
			buffer.col(num_samples_in_buffer - 1) = buffer.col(num_samples_in_buffer);
		}
	}
	
	template<typename prior_type>
	auto sample (const typename gse_type::supp_t counter,
				 const prior_type* const prior,
				 int interleaving,
				 int numsamples) {
		numsamples = std::min(numsamples, (int) max_num_samples);
		for (int j = (int) num_samples_in_buffer; j < numsamples; j++) {
			for (int k = 0; k < interleaving - 1; k++) {
				const auto src = buffer.col(j - 1);
				auto trgt = buffer.col(j);
				sampler.sampleNextState_external(prior,
												 &src,
												 &trgt,
												 counter,
												 rng);
				buffer.col(j - 1) = buffer.col(j);
			}
			const auto src = buffer.col(j - 1);
			auto trgt = buffer.col(j);
			sampler.sampleNextState_external(prior,
											 &src,
											 &trgt,
											 counter,
											 rng);
			num_samples_in_buffer++;
		}
		
	}
	
	auto get_samples (std::vector<typename gse_type::supp_t>* const ans) const {
		ans->resize(num_samples_in_buffer);
		for (int i = 0; i < num_samples_in_buffer; i++) {
			ans->operator[](i) = buffer(0, i);
		}
		return ans;
	}
};

#endif /* GibbsSamplerContainer_hpp */
