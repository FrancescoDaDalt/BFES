//
//  LatentBaseVector.hpp
//  BaInDaSt
//
//  Created by Francesco Da Dalt on 27.10.23.
//

#ifndef LatentBaseVector_hpp
#define LatentBaseVector_hpp

#include <Eigen/Dense>

// Object that stores one column vector of the latent base N in sparse format plus operations on the column
template <typename coeff_type, typename delta_type, typename index_type>
class LatentBaseVectorDynamic {
public:
//	Number of nonzero eleemnts in the column
	const size_t size;
//	Pointer to array of row-indices of nonzero elements in column
	const index_type* const index_ptr;
//	Eigen-map to indices
	const Eigen::Map<Eigen::Array<index_type, Eigen::Dynamic, 1>> index_map;
//	Map to nonzero coefficients of column
	const coeff_type* const coeff_ptr;
//	Eigen map to coeffs
	const Eigen::Map<Eigen::Array<coeff_type, Eigen::Dynamic, 1>> coeff_map;
//	Eigen array of teh soefficients but pre-cased to the currect type we require for gibbs sampling
	const Eigen::Array<delta_type, Eigen::Dynamic, 1> precasted_coeffs;
	
//	Generate Eigen mapping
	static constexpr auto compute_index_map (const size_t size,
											 const index_type* const index_ptr) {
		return Eigen::Map<Eigen::Array<index_type, Eigen::Dynamic, 1>>((index_type*) index_ptr, size);
	}
	
	//	Generate Eigen mapping
	static constexpr auto compute_coeff_map (const size_t size,
											 const coeff_type* const coeff_ptr) {
		return Eigen::Map<Eigen::Array<coeff_type, Eigen::Dynamic, 1>>((coeff_type*) coeff_ptr, size);
	}
	
	LatentBaseVectorDynamic (const size_t size,
							 const coeff_type* const coeff_ptr,
							 const index_type* const index_ptr):
	size(size),
	index_ptr(index_ptr),
	coeff_ptr(coeff_ptr),
	index_map(compute_index_map(size, index_ptr)),
	coeff_map(compute_coeff_map(size, coeff_ptr)),
	precasted_coeffs(coeff_map.template cast<delta_type>()) {
		
	}
	
	LatentBaseVectorDynamic(const LatentBaseVectorDynamic& other):
	size(other.size),
	index_ptr(other.index_ptr),
	coeff_ptr(other.coeff_ptr),
	index_map(other.index_map),
	coeff_map(other.coeff_map),
	precasted_coeffs(other.precasted_coeffs) {

	}
	
	LatentBaseVectorDynamic ():
	size(0),
	index_ptr(nullptr),
	coeff_ptr(nullptr),
	index_map(compute_index_map(0, nullptr)),
	coeff_map(compute_coeff_map(0, nullptr)),
	precasted_coeffs(0) {

	}
	
	const auto get_index_map () const {
		return &index_map;
	}
	
	const auto get_coeff_map () const {
		return &coeff_map;
	}
	
	const auto get_size() const {
		return size;
	}
	
//	Takes in a dense vector representing the current state of the gibbs sampler and discards the rows for which this column vector is zero.
	template <typename state_type>
	const auto projectState(const state_type* const state) const {
		return state->operator()(index_map);
	}
	
//	Takes in the current state and some delta d and increases the state by delta times this column vector
	template <typename state_type>
	const void updateState(state_type* const state,
						   const delta_type d) const {
		state->operator()(index_map) += d * precasted_coeffs;
	}
};

#endif /* LatentBaseVector_hpp */
