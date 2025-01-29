//
//  HyperplaneSensor.hpp
//  BaInDaSt
//
//  Created by Francesco Da Dalt on 25.10.23.
//

#ifndef HyperplaneSensor_hpp
#define HyperplaneSensor_hpp

#include "../../Utility/Utility.hpp"
#include <vector>
#include <array>
#include "distinct_counter.h"
#include "MurmurHash3.h"

template <typename key_type, typename support_type, typename hypergrid_type>
class HyperplaneSensor {
	using key_t = key_type;
	using hypergrid_t = hypergrid_type;
	
//	How many axes the LDPC hypergris has. If num_axes=1 then H is the identity matrix
	static constexpr size_t num_axes = hypergrid_type::num_axes;
	const std::array<size_t, num_axes> hypergrid_base;
	const size_t num_counters;
	const size_t num_nodes;
//	Basic hash function to go from key-type to size_t
	static constexpr auto basicHash = std::hash<key_type>{};
//	Array of counters used to store frequency information
	std::vector<support_type> counters;
//	Custom counter to store the sum over all frequencies
	support_type total = 0;
//	Auxiliary sketch to estimate the number of distinct keys in the stream
	hyperloglog_hip::distinct_counter<size_t> distinct_counter;
//	Auxuliary array for the LDPC code
	const std::array<size_t, num_axes> pseudohash_offset;
	
public:
	
	void clamp_min (const support_type min) {
		for (int i = 0; i < num_counters; i++) {
			total += std::max(0.0, min - counters[i]);
			counters[i] = std::max(counters[i], min);
		}
	}
	
	//	Auxuliary function for the LDPC code.
	const std::array<size_t, num_axes> compute_pseudohash_offset() const {
		std::array<size_t, num_axes> ans = {0};
		for (int i = 1; i < num_axes; i++) {
			ans[i] = ans[i - 1] + hypergrid_base[i - 1];
		}
		return ans;
	}
	
	HyperplaneSensor (const hypergrid_type &hg,
					  const size_t num_memorycells_secondary):
	hypergrid_base(hg.axes),
	num_counters(hg.num_counters),
	num_nodes(hg.num_nodes),
	pseudohash_offset(compute_pseudohash_offset()) {
		counters = std::vector<support_type>(num_counters, 0);
		new (&distinct_counter) hyperloglog_hip::distinct_counter<size_t>(num_memorycells_secondary);
	}
	
//	Map keys to their level-1 bins
	void topLevelKeyProjection (const std::vector<key_type>* const keys,
								std::vector<size_t>* const dest,
								const bool bypassTopLevel) const {
		const size_t num_keys = keys->size();
		dest->resize(num_keys);
		for (int i = 0; i < num_keys; i++) {
			const size_t intermediate = basicHash(keys->at(i));
//			128 bit struct for 128 bit hash
			size_t extradest[2];
			if (!bypassTopLevel) {
//				We hash (with seed 42)
				MurmurHash3_x64_128((const void*) &intermediate, sizeof(size_t), 42, (void *) extradest);
//				and then we do modulo
				dest->at(i) = extradest[0] % num_nodes;
			} else {
//				We dont hash and only do modulo operation
				dest->at(i) = intermediate % num_nodes;
			}
		}
	}
	
//	Insert key-value into object
	void insert(const key_type key_raw,
				const support_type value,
				const bool bypassTopLevel) {
//		Hash from key-type to size_t
		const size_t key_aux = basicHash(key_raw);
//		Hash from 64-bit to the intermediate d bins
		size_t key_countdistinct;
		size_t extradest[2];
		MurmurHash3_x64_128(&key_aux, sizeof(size_t), 42, &extradest);
		key_countdistinct = extradest[0];
//		Insert into auxiliar sketch to measure stream cardinality
		distinct_counter.insert(key_countdistinct);
		
		size_t key = key_countdistinct;
		if (bypassTopLevel) {key = (size_t) key_aux;}
		key = key % num_nodes;
		
//		Hash from the d bins to the m counters. If num_axes = 1 this is an identity mapping.
//		This piece of code represents an LDPC code where the d intermediate bins have been arranged on a num_axes-dimensional hypergrid and one of the m counters corresponds to one axis-aligned hyperplane in this hypergrid.
//		Linear dependencies fo hyperplanes have been taken into account and eliminated
		const size_t lower = 1;
		const size_t upper = num_axes;
		const size_t loopsize = upper - lower;
		size_t running_key = key;
		size_t offset = 0;
		const size_t axis = hypergrid_base[0];
		const size_t hash = running_key % axis;
		counters[hash + offset] += value;
		offset += axis;
		running_key /= axis;
		for (size_t idx = 0; idx < loopsize; idx++) {
			const size_t i = idx + lower;
			const size_t axis = hypergrid_base[i];
			const size_t hash = running_key % axis;
			running_key /= axis;
			if (hash > 0) {
				counters[hash - 1 + offset] += value;
			}
			offset += axis - 1;
		}
//		Increase auxiliary counter taht measures total stream volume
		total += value;
	}

	
//	Auxuliary function for the LDPC code.
	const std::array<size_t, num_axes> grid_coordinates(const key_type key) const {
		std::array<size_t, num_axes> ans;
		key_type running_key = key;
		for (int i = 0; i < num_axes; i++) {
			const size_t axis = hypergrid_base[i];
			ans[i] = running_key % axis;
			running_key /= axis;
		}
		return ans;
	}
	
//	Auxuliary function for the LDPC code that for a key k from 0 to d-1 returns the list of indices of counters where this key falls into.
//	Its not really optimized well because it is only used for level-1 sampling and if we care about speed we omit level-1 sampling altogether
	const std::vector<size_t> pseudohash(const key_type key) const {
		auto coordinates = grid_coordinates(key);
		std::vector<size_t> ans;
		ans.push_back(coordinates[0]);
		for (int i = 1; i < num_axes; i++) {
			const size_t coord = coordinates[i];
			if (coord > 0) {
				ans.push_back(coord + pseudohash_offset[i] - i);
			}
		}
		return ans;
	}
	
//	Get pointer to counters
	inline const support_type* const counter_ptr() const {
		return counters.data();
	}
	
//	Get total volume in stream
	inline const support_type get_total() const {
		return total;
	}
	
//	Get estimate fo number of items in stream
	inline const support_type get_distinct_estimate() const {
		return distinct_counter.count();
	}
	
//	Clear the object
	inline void clear () {
//		Clear counters
		std::fill(counters.cbegin(), counters.cend(), (support_type) 0);
//		Clear cardinality sketch
		new (&distinct_counter) hyperloglog_hip::distinct_counter<size_t>();
	}
	
	auto getNummemcellsFast () const {
		return num_counters + distinct_counter.size() / 32;
	}
	
	auto getNummemcells () const {
		return num_counters + distinct_counter.size() / 32;
	}
};

#endif /* HyperplaneSensor_hpp */
