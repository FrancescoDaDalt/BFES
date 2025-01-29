//
//  HyperplaneStructure.hpp
//  BaInDaSt
//
//  Created by Francesco Da Dalt on 26.10.23.
//

#ifndef HyperplaneStructure_hpp
#define HyperplaneStructure_hpp

template <int... axes>
class HyperplaneStructure {
public:
	static constexpr int hypergrid_base[sizeof...(axes)] = {axes...};
	static constexpr int compute_num_counters() {
		int length = 0;
		for (int arg : {axes...}) {
			length += arg;
		}
		return length - sizeof...(axes) + 1;
	}
	static constexpr int num_counters = compute_num_counters();
	static constexpr int num_axes = (int) sizeof...(axes);
	static constexpr int compute_num_fundament_nodes() {
		int length = 1;
		for (int arg : {axes...}) {
			length *= arg;
		}
		return length;
	}
	static constexpr int num_fundament_nodes = compute_num_fundament_nodes();
};

#endif /* HyperplaneStructure_hpp */
