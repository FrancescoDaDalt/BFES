/* 
 * MIT License
 * 
 * Copyright (c) 2023 Francesco Da Dalt
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

#ifndef LatentBase_hpp
#define LatentBase_hpp

#include "../LatentBaseVector/LatentBaseVector.hpp"
#include <numeric>
#include <gsl/gsl_rng.h>
#include <Eigen/Dense>

// This class construct the latent base N of teh matrix H defined by the LDPC hypergrid
// N is constructed in a sparse manner and optimized for the purpose of gibbs sampling
template <typename hypergrid_type, typename delta_type>
class HyperplaneLatentBase {
public:
	
	using hypergrid_t = hypergrid_type;
	using coef_t = signed char;
	
//	Defines the LDPC code and thus the matrix H and thus the latent base N
	const hypergrid_type hg;
	
//	Takes support type and an integer and float type and returns the integer type if support is integral type and float type otherwise
	template <typename supp_t, typename integer_t, typename float_t>
	struct LatentType {
		using t = typename std::conditional<
		std::is_integral<supp_t>::value,
		integer_t,
		float_t
		>::type;
	};
	
//	Latent dim is d - m
	const int latent_dim;
//	Target dim is d
	const int target_dim;
	static constexpr int num_axes = hypergrid_type::num_axes;
	
//	Auxiliary function to compute num_axes-dimensional gray-code. Needed for generating N
	static constexpr std::array<std::array<int, num_axes>, 1 << num_axes> generate_gray_code() {
		const int gray_length = 1 << num_axes;
		std::array<std::vector<int>, num_axes> auxilliary;
		for (int i = 0; i < num_axes; i++) {
			std::vector<int> ones(1 << i, 1);
			std::vector<int> zeros(1 << i, 0);
			zeros.insert(zeros.end(), ones.begin(), ones.end());
			for (int j = 0; j < i; j++) {
				std::vector<int> current = auxilliary[j];
				std::vector<int> current_flipped = auxilliary[j];
				std::reverse(current_flipped.begin(), current_flipped.end());
				current.insert(current.end(), current_flipped.begin(), current_flipped.end());
				auxilliary[j] = current;
			}
			auxilliary[i] = zeros;
		}
		std::array<std::array<int, num_axes>, 1 << num_axes> ans;
		for (int i = 0; i < gray_length; i++) {
			std::array<int, num_axes> tmparr;
			for (int j = 0; j < num_axes; j++) {
				tmparr[j] = auxilliary[j].at(i);
			}
			ans[i] = tmparr;
		}
		return ans;
	}
	
	static constexpr std::array<int, num_axes + 1> free_plus_base_axes_coordinates(const hypergrid_type &hg, const int index) {
		std::array<int, num_axes + 1> ans;
		int running_index = index;
		auto hypergrid_base = hg.axes;
		for (int i = 0; i < num_axes; i++) {
			const size_t axis = hypergrid_base[i];
			ans[i] = running_index % axis;
			running_index /= axis;
		}
		ans[num_axes] = running_index;
		return ans;
	}
	
	static constexpr std::array<int, num_axes + 1> index_to_coordinates(const int index, std::array<int, num_axes> coordinate_base) {
		std::array<int, num_axes + 1> ans;
		int running_index = index;
		for (int i = 0; i < num_axes; i++) {
			const int axis = coordinate_base[i];
			ans[i] = running_index % axis;
			running_index /= axis;
		}
		ans[num_axes] = running_index;
		return ans;
	}
	
	static constexpr int free_plus_base_axes_index(const hypergrid_type &hg, const std::array<int, num_axes + 1> coordinates) {
		int ans = 0;
		int running_base = 1;
		auto hypergrid_base = hg.axes;
		for (int i = 0; i < num_axes; i++) {
			ans += coordinates[i] * running_base;
			running_base *= hypergrid_base[i];
		}
		ans += coordinates[num_axes] * running_base;
		return ans;
	}
	
	static constexpr void mask_iterator(std::array<bool, hypergrid_type::num_axes> &mask) {
		bool carryover = true;
		for (int i = 0; i < hypergrid_type::num_axes; i++) {
			if (mask[i] && carryover) {
				mask[i] = false;
			} else if (carryover) {
				mask[i] = true;
				carryover = false;
			} else {
				break;
			}
		}
	}
	
	const auto compute_num_subbases_vector (const hypergrid_type &hg) {
		constexpr int axis_combinations = (1 << num_axes) - 1 - num_axes + 1;
		std::array<int, axis_combinations> num_subbases_vector = {0};
		auto hypergrid_base = hg.axes;
		std::array<bool, num_axes> mask = {0};
		mask_iterator(mask);
		int it = 0;
		while (!std::all_of(mask.begin(), mask.end(), [](bool value) { return !value; })) {
			const int count = (const int) std::count(mask.begin(), mask.end(), true);
			if (count > 1) {
				std::array<int, num_axes> coordinate_base;
				for (int i = 0; i < num_axes; i++) {
					if (mask[i]) {
						coordinate_base[i] = hypergrid_base[i] - 1;
					} else {
						coordinate_base[i] = 1;
					}
				}
				const int num_subbases = std::accumulate(coordinate_base.begin(), coordinate_base.end(), 1, std::multiplies<int>());
				num_subbases_vector[it] = num_subbases;
				it++;
			}
			mask_iterator(mask);
		}
		num_subbases_vector[axis_combinations - 1] = 0;
		
		return num_subbases_vector;
	}
	
	const std::array<int, (1 << num_axes) - 1 - num_axes + 1> num_subbases_vector;
	
	const auto compute_num_subbases_acc_vector (const std::array<int,(1<<num_axes)-1-num_axes+1> &num_subbases_vector) {
		constexpr int axis_combinations = (1 << num_axes) - 1 - num_axes + 1;
		std::array<int, axis_combinations> num_subbases_acc_vector = {0};
		for (int i = 1; i < axis_combinations; i++) {
			num_subbases_acc_vector[i] = num_subbases_acc_vector[i - 1] + num_subbases_vector[i - 1];
		}
		return num_subbases_acc_vector;
	}
	
	const std::array<int, (1 << num_axes) - 1 - num_axes + 1> num_subbases_acc_vector;
	
	static constexpr auto compute_nnz_per_basis_vector (const hypergrid_type &hg, const int &latent_dim) {
		std::vector<short> nnz_per_basis_vector(latent_dim);
		auto hypergrid_base = hg.axes;
		
		std::array<bool, num_axes> mask = {0};
		mask_iterator(mask);
		int it = 0;
		while (!std::all_of(mask.begin(), mask.end(), [](bool value) { return !value; })) {
			const int count = (const int) std::count(mask.begin(), mask.end(), true);
			if (count > 1) {
				const int num_hypercube_vertices = 1 << count;
				std::array<int, num_axes> coordinate_base;
				for (int i = 0; i < num_axes; i++) {
					if (mask[i]) {
						coordinate_base[i] = hypergrid_base[i] - 1;
					} else {
						coordinate_base[i] = 1;
					}
				}
				const int num_subbases = std::accumulate(coordinate_base.begin(), coordinate_base.end(), 1, std::multiplies<int>());
				for (int i = 0; i < num_subbases; i++) {
					nnz_per_basis_vector[it] = num_hypercube_vertices;
					it++;
				}
			}
			mask_iterator(mask);
		}
		for (; it < latent_dim; it++) {
			nnz_per_basis_vector[it] = 2;
		}
		return nnz_per_basis_vector;
	}
	
	const std::vector<short> nnz_per_basis_vector;
	
	const size_t total_nnz;
	
	const auto compute_basis () {
		constexpr int axis_combinations = (1 << num_axes) - 1 - num_axes + 1;
		const int total_num_nnz = total_nnz;
		auto basis_outer_idx = std::vector<unsigned int>(latent_dim + 1);
		auto basis_inner_idx = std::vector<unsigned int>(total_num_nnz);
		auto basis_coeffs = std::vector<coef_t>(total_num_nnz);
		auto hypergrid_base = hg.axes;
		std::array<bool, num_axes> mask = {0};
		mask_iterator(mask);
		std::array<std::array<int, num_axes>, 1 << num_axes> gray_code = generate_gray_code();
		int running_outer_idx = 0;
		for (int i = 0; i < axis_combinations - 1; i++) {
			const int subbase_offset = num_subbases_acc_vector[i];
			const int subbase_num = num_subbases_vector[i];
			while (!std::all_of(mask.begin(), mask.end(), [](bool value) { return !value; })) {
				const int count = (const int) std::count(mask.begin(), mask.end(), true);
				if (count > 1) {
					std::array<int, num_axes> coordinate_base;
					for (int i = 0; i < num_axes; i++) {
						if (mask[i]) {
							coordinate_base[i] = hypergrid_base[i] - 1;
						} else {
							coordinate_base[i] = 1;
						}
					}
					for (int k = 0; k < subbase_num; k++) {
						const int it_offset = (int) k;
						const int total_offset = it_offset + subbase_offset;
						
						std::array<int, num_axes + 1> base_coordinate = index_to_coordinates(it_offset, coordinate_base);
						
						
						for (int j = 0; j < nnz_per_basis_vector[total_offset]; j++) {
							std::array<int, num_axes + 1> running_coordinate = base_coordinate;
							std::array<int, num_axes> gray_coordinate = gray_code[j];
							int k_aux = 0;
							for (int k = 0; k < num_axes; k++) {
								if (mask[k]) {
									running_coordinate[k] += gray_coordinate[k_aux++];
								}
							}
							int index = free_plus_base_axes_index(hg, running_coordinate);
							int coeff = (j % 2) * 2 - 1;
							basis_inner_idx[running_outer_idx] = index;
							basis_coeffs[running_outer_idx] = coeff;
							running_outer_idx++;
						}
						basis_outer_idx[total_offset + 1] = running_outer_idx;
					}
					mask_iterator(mask);
					break;
				}
				mask_iterator(mask);
			}
		}
		const int subbase_offset = num_subbases_acc_vector[axis_combinations - 1];
		const int subbase_num = num_subbases_vector[axis_combinations - 1];
		for (int i = 0; i < subbase_num; i++) {
			int lower_index = i;
			int upper_index = lower_index + hg.num_nodes;
			basis_inner_idx[running_outer_idx] = lower_index;
			basis_coeffs[running_outer_idx] = 1;
			running_outer_idx++;
			basis_inner_idx[running_outer_idx] = upper_index;
			basis_coeffs[running_outer_idx] = -1;
			running_outer_idx++;
			basis_outer_idx[subbase_offset + i + 1] = running_outer_idx;
		}
		auto basis = std::make_tuple(basis_outer_idx, basis_inner_idx, basis_coeffs);
		return basis;
	}
//	This contains N in a sparse format using only std-containers
	const std::tuple<std::vector<unsigned int>, std::vector<unsigned int>, std::vector<coef_t>> basis_std;
	
//	We wrap the std-container version of N into a custom object version where each column of N is stored as an object of lbv_type
	using lbv_type = LatentBaseVectorDynamic<coef_t, delta_type, unsigned int>;
//	Function to wrap the std container of N
	const std::vector<lbv_type> compute_basis_wrapped() {
		std::vector<lbv_type> ans(latent_dim);
		for (int i = 0; i < latent_dim; i++) {
			const size_t size = get<0>(basis_std)[i + 1] - get<0>(basis_std)[i];
			const size_t index_base = get<0>(basis_std)[i];
			const auto& coeff_ptr = &get<2>(basis_std)[index_base];
			const auto& index_ptr = &get<1>(basis_std)[index_base];
			new (ans.data() + i) lbv_type(size, coeff_ptr, index_ptr);
		}
		return ans;
	}
	
//	Contains N in a custom format
	const std::vector<lbv_type> basis_wrapped;
//	Pointer to N
	const lbv_type* basis_wrapped_data_ptr;
	
//	Get i-th column of N
	const auto get_col (size_t i) const {
		return (const lbv_type* const) (basis_wrapped.data() + i);
	}

	HyperplaneLatentBase (const hypergrid_type &hg):
	hg(hg),
	latent_dim(hg.num_nodes - hg.num_counters),
	target_dim(hg.num_nodes),
	num_subbases_vector(compute_num_subbases_vector(hg)),
	num_subbases_acc_vector(compute_num_subbases_acc_vector(num_subbases_vector)),
	nnz_per_basis_vector(compute_nnz_per_basis_vector(hg, latent_dim)),
	total_nnz(std::accumulate(nnz_per_basis_vector.begin(), nnz_per_basis_vector.end(), 0)),
	basis_std(compute_basis()),
	basis_wrapped(compute_basis_wrapped()),
	basis_wrapped_data_ptr(basis_wrapped.data()) {
		
	}
	
};

// Conceptually the same as HyperplaneLatentBase except that it is not based on an LDPC code. Instead the matrix "H" here is a row-vector of ones of size target_dim_.
// Delta type is the type used for the state-change during gibbs sampling
template<int target_dim_, typename delta_type>
class DegenerateLatentBase {
public:
	
	using coef_t = signed char;
	
	template <typename supp_t, typename integer_t, typename float_t>
	struct LatentType {
		using t = typename std::conditional<
		std::is_integral<supp_t>::value,
		integer_t,
		float_t
		>::type;
	};
	
	static constexpr int target_dim = target_dim_;
	static constexpr int latent_dim = target_dim - 1;
	
//	Computes the std-container version of N
	const auto compute_basis () {
		const int total_num_nnz = 2 * latent_dim;
		auto basis_outer_idx = std::vector<unsigned int>(latent_dim + 1);
		auto basis_inner_idx = std::vector<unsigned int>(total_num_nnz);
		auto basis_coeffs = std::vector<coef_t>(total_num_nnz);
		
		int running_outer_idx = 0;
		basis_outer_idx[0] = running_outer_idx;
		
		for (int i = 0; i < latent_dim; i++) {
			
			int lower_index = (i + 2) % (latent_dim + 1);
			int upper_index = (i + 1) % (latent_dim + 1);
			basis_inner_idx[running_outer_idx] = lower_index;
			basis_coeffs[running_outer_idx] = 1;
			running_outer_idx++;
			
			basis_inner_idx[running_outer_idx] = upper_index;
			basis_coeffs[running_outer_idx] = -1;
			running_outer_idx++;
			
			basis_outer_idx[i + 1] = running_outer_idx;
		}
		auto basis = std::make_tuple(basis_outer_idx, basis_inner_idx, basis_coeffs);
		return basis;
	}
	
//	Computes an augmented basis of N with more columns than necessary. Useful for alternating between different bases N in order to improve gibbs sampling mixture speed
	const auto compute_multi_basis () {
		const int total_num_nnz = 2 * latent_dim * 2;
		
		auto basis_outer_idx = std::vector<unsigned int>(latent_dim * 2 + 1);
		auto basis_inner_idx = std::vector<unsigned int>(total_num_nnz);
		auto basis_coeffs = std::vector<coef_t>(total_num_nnz);
		
		int running_outer_idx = 0;
		basis_outer_idx[0] = running_outer_idx;
//		First basis
		for (int i = 0; i < latent_dim; i++) {
			
			int lower_index = i;
			int upper_index = i + 1;
			basis_inner_idx[running_outer_idx] = lower_index;
			basis_coeffs[running_outer_idx] = 1;
			running_outer_idx++;
			
			basis_inner_idx[running_outer_idx] = upper_index;
			basis_coeffs[running_outer_idx] = -1;
			running_outer_idx++;
			
			basis_outer_idx[i + 1] = running_outer_idx;
		}
//		Second basis that is different from the first basis
		for (int i = latent_dim; i < latent_dim * 2; i++) {
			
			int lower_index = (i) % (latent_dim + 1);
			int upper_index = (i + 2) % (latent_dim + 1);
			
			if (lower_index > upper_index) {
				std::swap(upper_index, lower_index);
			}
			basis_inner_idx[running_outer_idx] = lower_index;
			basis_coeffs[running_outer_idx] = 1;
			running_outer_idx++;
			
			basis_inner_idx[running_outer_idx] = upper_index;
			basis_coeffs[running_outer_idx] = -1;
			running_outer_idx++;
			
			basis_outer_idx[i + 1] = running_outer_idx;
		}
		auto basis = std::make_tuple(basis_outer_idx, basis_inner_idx, basis_coeffs);
		return basis;
	}
	
//	Store std-container of bases
	const std::tuple<std::vector<unsigned int>, std::vector<unsigned int>, std::vector<coef_t>> basis_std;
	const std::tuple<std::vector<unsigned int>, std::vector<unsigned int>, std::vector<coef_t>> multi_basis_std;
	
//	We wrap the std-container version of N into a custom object version where each column of N is stored as an object of lbv_type
	using lbv_type = LatentBaseVectorDynamic<coef_t, delta_type, unsigned int>;
	
//	Wrap std-container
	const std::vector<lbv_type> compute_basis_wrapped() {
		std::vector<lbv_type> ans(latent_dim);
		for (int i = 0; i < latent_dim; i++) {
			const size_t size = get<0>(basis_std)[i + 1] - get<0>(basis_std)[i];
			const size_t index_base = get<0>(basis_std)[i];
			const auto& coeff_ptr = &get<2>(basis_std)[index_base];
			const auto& index_ptr = &get<1>(basis_std)[index_base];
			
			new (ans.data() + i) lbv_type(size, coeff_ptr, index_ptr);
		}
		return ans;
	}
	
	const std::vector<lbv_type> compute_multi_basis_wrapped() {
		std::vector<lbv_type> ans(latent_dim * 2);
		for (int i = 0; i < latent_dim * 2; i++) {
			const size_t size = get<0>(multi_basis_std)[i + 1] - get<0>(multi_basis_std)[i];
			const size_t index_base = get<0>(multi_basis_std)[i];
			const auto& coeff_ptr = &get<2>(multi_basis_std)[index_base];
			const auto& index_ptr = &get<1>(multi_basis_std)[index_base];
			
			new (ans.data() + i) LatentBaseVectorDynamic<coef_t, delta_type, unsigned int>(size, coeff_ptr, index_ptr);
		}
		return ans;
	}
//	Store wrapped bases
	const std::vector<lbv_type> basis_wrapped;
	const std::vector<lbv_type> multi_basis_wrapped;
	
//	Get i-th column of N
	const auto get_col (size_t i) const {
		return (const lbv_type* const) (basis_wrapped.data() + i);
	}
	
//	Get a random column from one of the two bases of N computed
	template<typename rng_type>
	const auto random_col (rng_type* const rng) const {
		const int r = gsl_rng_uniform_int(rng, latent_dim * 2);
		return (const lbv_type* const) (multi_basis_wrapped.data() + r);
	}
	
	DegenerateLatentBase ():
	basis_std(compute_basis()),
	basis_wrapped(compute_basis_wrapped()),
	multi_basis_std(compute_multi_basis()),
	multi_basis_wrapped(compute_multi_basis_wrapped()){
		
	}
	
};

#endif /* LatentBase_hpp */
