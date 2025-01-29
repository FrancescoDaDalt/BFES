//
//  Utility.hpp
//  BaInDaSt
//
//  Created by Francesco Da Dalt on 02.11.23.
//

#ifndef Utility_hpp
#define Utility_hpp

#include <tuple>
#include <array>
#include <Eigen/Dense>

using size_t = std::size_t;
using size32s_t = int;

template <const auto& Arr, typename Seq = std::make_index_sequence<std::size(Arr)>>
struct make_seq;

template <typename T, std::size_t N, const T (&Arr)[N], std::size_t ... Is>
struct make_seq<Arr, std::index_sequence<Is...>>
{
	using type = std::integer_sequence<T, Arr[Is]...>;
};

template <typename T, std::size_t N, const std::array<T, N>& Arr, std::size_t ... Is>
struct make_seq<Arr, std::index_sequence<Is...>>
{
	using type = std::integer_sequence<T, Arr[Is]...>;
};

template <size_t... vals, typename Tuple>
constexpr auto select_elements(const Tuple& t, std::integer_sequence<size_t, vals...>) {
	return std::make_tuple(std::get<vals>(t)...);
}

template <size_t arrsize, typename arrtype, const std::array<arrtype, arrsize>& arr, typename Tuple>
constexpr auto select_elements_arr(Tuple t) {
	return select_elements(t, typename make_seq<arr>::type());
}

template<size_t Size, typename F>
static constexpr void constexpr_for_backend(F&& function) {
	auto unfold = [&]<size_t... Ints>(std::index_sequence<Ints...>) {
		(std::forward<F>(function)(std::integral_constant<size_t, Ints>{}), ...);
	};
	
	unfold(std::make_index_sequence<Size>());
}

template<size_t Low, size_t High, typename F>
static constexpr void constexpr_ranged_for_backend(F&& function) {
	constexpr size_t Size = High - Low;
	
	auto unfold = [&]<size_t... Ints>(std::index_sequence<Ints...>) {
		(std::forward<F>(function)(std::integral_constant<size_t, Ints>{} + Low), ...);
	};
	
	unfold(std::make_index_sequence<Size>());
}

template<size_t Size>
struct constexpr_for {
	template<typename F>
	static constexpr void loop_alt(F&& function) {
		constexpr size_t numSizes = (Size - 1) / 256 + 1;
		
		auto unfold = [&]<size_t... Ints>(std::index_sequence<Ints...>) {
			(std::forward<F>([&] (auto i) {
				constexpr size_t low = std::min(i * 256, Size);
				constexpr size_t high = std::min((i + 1) * 256, Size);
				constexpr size_t Size_Local = high - low;
				
				auto unfold2 = [&]<size_t... IntsX>(std::index_sequence<IntsX...>) {
					(std::forward<F>(function)(std::integral_constant<size_t, IntsX>{} + low), ...);
				};
				
				unfold2(std::make_index_sequence<Size_Local>());
			})(std::integral_constant<size_t, Ints>{}), ...);
		};
		
		unfold(std::make_index_sequence<numSizes>());
	}
	
	template<typename F>
	static constexpr void loop(F&& function) {
		constexpr_for_backend<Size>(function);
	}
};

template<size_t numaxes>
class hypergrid {
public:
	const std::array<size_t, numaxes> axes;
	const size_t num_counters;
	const size_t num_nodes;
	hypergrid (std::array<size_t, numaxes> axes_arg):
	axes(axes_arg),
	num_counters(compute_num_counters()),
	num_nodes(compute_num_nodes()) {
		
	}
	
	static constexpr std::string description () {
		return "hg" + std::to_string(numaxes) + "dim";
	}
	static constexpr size_t num_axes = numaxes;
private:
	const size_t compute_num_counters() const {
		size_t length = 0;
		for (auto a : axes) {
			length += a;
		}
		return length - num_axes + 1;
	}
	const size_t compute_num_nodes() const {
		size_t length = 1;
		for (auto a: axes) {
			length *= a;
		}
		return length;
	}
};

template<size_t I>
static std::array<size_t, I> array_aux (const size_t X, const size_t Y) {
	std::array<size_t, I> A;
	A.fill(Y);
	A[0] = X;
	return A;
}


#endif /* Utility_hpp */
