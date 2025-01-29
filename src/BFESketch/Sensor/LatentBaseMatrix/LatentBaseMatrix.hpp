//
//  LatentBaseMatrix.hpp
//  BaInDaSt
//
//  Created by Francesco Da Dalt on 29.10.23.
//

#ifndef LatentBaseMatrix_hpp
#define LatentBaseMatrix_hpp




#include "../LatentBaseVector/LatentBaseVector.hpp"
#include "../../Utility/Utility.hpp"


template <class basis_type>
class LatentBaseMatrix {
	const basis_type basis;
	const int num_cols = std::tuple_size<basis_type>::value;
	const int num_rows;
	const std::array<int, std::tuple_size<basis_type>::value> num_nnz_per_col;
public:
	constexpr std::array<int, std::tuple_size<basis_type>::value> compute_num_nnz_per_col (const basis_type basis) {
		std::array<int, std::tuple_size<basis_type>::value> ans;
		constexpr_for<std::tuple_size<basis_type>::value>::loop ([&] (auto i) {
			ans[i] = std::get<i>(basis).get_size();
		});
		return ans;
	}

	LatentBaseMatrix (const basis_type basis, const int num_rows): basis(basis), num_rows(num_rows), num_nnz_per_col(compute_num_nnz_per_col(basis))  {}

	const basis_type* get_basis() const {
		return &basis;
	}

	const std::array<int, std::tuple_size<basis_type>::value>* get_num_nnz_per_col() const {
		return &num_nnz_per_col;
	}

	const int get_num_cols() const {
		return num_cols;
	}

	const int get_num_rows() const {
		return num_rows;
	}
};



//template <typename T1, typename T2, int N>
//struct arraySize_in_tuple;
//
//template <typename T1, typename T2, int N>
//struct arraySize_in_tuple<std::tuple<std::array<T1, N>, std::array<T2, N>>> {
//	static constexpr int value = N;
//};

template <typename basistype, const basistype &basis>
class LatentBaseMatrixWrapper {
public:
	static constexpr auto compute_nnz_per_col () {
		constexpr int size = std::tuple_size<basistype>::value;
		std::array<int, size> ans;
		constexpr_for<size>::loop([&](auto i){
			ans[i] = std::get<0>(std::get<i>(basis)).size();
		});
		return ans;
	}
	static constexpr auto nnz_per_col = compute_nnz_per_col();
	
	static constexpr auto compute_num_rows () {
		constexpr int size = std::tuple_size<basistype>::value;
		int ans = 0;
		constexpr_for<size>::loop([&](auto i){
			constexpr auto a = std::get<0>(std::get<i>(basis));
			ans = std::max(*std::max_element(a.begin(), a.end()), ans);
		});
		return ans;
	}
	static constexpr auto num_rows = compute_num_rows();
	
	static constexpr auto compute_num_cols () {
		constexpr int size = std::tuple_size<basistype>::value;
		return size;
	}
	static constexpr auto num_cols = compute_num_cols();
	
	static constexpr auto compute_LatentBaseVectorWrapper () {
		
	}
	
	
	template <int i>
	static constexpr auto compute_basis_wrapped_aux_aux () {
		//		using arrElemType = typename std::remove_all_extents<decltype(std::get<0>(std::get<i>(basis)))>::type;
		static constexpr int arrLen = std::get<0>(std::get<i>(basis)).size();
		static constexpr auto b = std::get<i>(basis);
		//		static constexpr auto b = std::tuple<std::array<int, 1>, std::array<int, 1>>(std::array<int, 1>{1}, std::array<int, 1>{2});
		return LatentBaseVectorWrapper<int, arrLen, b>();
	}
	
	template <int... vals>
	static constexpr auto compute_basis_wrapped_aux (std::integer_sequence<int, vals...>) {
		return std::make_tuple(compute_basis_wrapped_aux_aux<vals>()...);
	}
	
	static constexpr auto compute_basis_wrapped () {
		return compute_basis_wrapped_aux(typename std::make_integer_sequence<int, std::tuple_size<basistype>::value>{});
	}
	
	static const inline auto basis_wrapped = compute_basis_wrapped();
	
	const auto blablatest () const {
		return &basis_wrapped;
	}
	template <int i>
	static constexpr auto get_col () {
		return &std::get<i>(basis_wrapped);
	}
	
//	LatentBaseMatrixWrapper (): basis_wrapped(compute_basis_wrapped()) {}

	
	
};



#endif /* LatentBaseMatrix_hpp */
