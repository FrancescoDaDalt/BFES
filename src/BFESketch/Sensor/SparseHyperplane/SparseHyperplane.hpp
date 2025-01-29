//
//  SparseHyperplane.hpp
//  BaInDaSt
//
//  Created by Francesco Da Dalt on 25.10.23.
//

#ifndef SparseHyperplane_hpp
#define SparseHyperplane_hpp

#include <vector>
#include <utility>
#include <tuple>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

//template <typename hyperplaneCoeff_type, std::pair<int, hyperplaneCoeff_type>...args>
//class SparseHyperplane {
//	static constexpr int nnz = (int) sizeof...(args);
//
//	static constexpr Eigen::Array<int, nnz, 1> compute_indices() {
//		Eigen::Array<int, nnz, 1> ans;
//		std::pair<int, hyperplaneCoeff_type> arg_array[nnz] = {args...};
//		for (int i = 0; i < nnz; i++) {
//			ans(i) = arg_array[i].first;
//		}
//		return ans;
//	}
//
//	static constexpr Eigen::Array<hyperplaneCoeff_type, nnz, 1> compute_coeffs() {
//		Eigen::Array<hyperplaneCoeff_type, nnz, 1> ans;
//		std::pair<int, hyperplaneCoeff_type> arg_array[nnz] = {args...};
//		for (int i = 0; i < nnz; i++) {
//			ans(i) = arg_array[i].second;
//		}
//		return ans;
//	}
//
//	const Eigen::Array<int, nnz, 1> indices = compute_indices();
//	const Eigen::Array<hyperplaneCoeff_type, nnz, 1> coefficients = compute_coeffs();
//
//public:
//	SparseHyperplane(Eigen::ArrayXi* indices, Eigen::Array<hyperplaneCoeff_type, Eigen::Dynamic, 1>* coefficients): indices(*indices), coefficients(*coefficients) {}
//
//	template <typename support_type>
//	void projectState(const Eigen::Array<support_type, Eigen::Dynamic, 1>* const state, Eigen::Array<support_type, Eigen::Dynamic, 1>* const target) const {
//		target->resize(indices.size());
//		*target = state->operator()(indices);
//	}
//
//	const Eigen::Array<hyperplaneCoeff_type, Eigen::Dynamic, 1>* const get_coefficients() const {
//		return &coefficients;
//	}
//
//	const Eigen::ArrayXi* const get_indices() const {
//		return &indices;
//	}
//};

template <int length, typename t>
constexpr Eigen::Array<t, length, 1> std_to_eigen_array (std::array<t, length> arr) {
	Eigen::Array<t, length, 1> ans;
	for (int i = 0; i < length; i++) {
		ans(i) = arr[i];
	}
	return ans;
}

template <typename hyperplaneCoeff_type>
class SparseHyperplane {
	const Eigen::ArrayXi indices;
	const Eigen::Array<hyperplaneCoeff_type, Eigen::Dynamic, 1> coefficients;

public:
	SparseHyperplane(Eigen::ArrayXi* indices, Eigen::Array<hyperplaneCoeff_type, Eigen::Dynamic, 1>* coefficients): indices(*indices), coefficients(*coefficients) {}

	template <typename support_type>
	void projectState(const Eigen::Array<support_type, Eigen::Dynamic, 1>* const state, Eigen::Array<support_type, Eigen::Dynamic, 1>* const target) const {
		target->resize(indices.size());
		*target = state->operator()(indices);
	}

	const Eigen::Array<hyperplaneCoeff_type, Eigen::Dynamic, 1>* const get_coefficients() const {
		return &coefficients;
	}

	const Eigen::ArrayXi* const get_indices() const {
		return &indices;
	}
};

#endif /* SparseHyperplane_hpp */
