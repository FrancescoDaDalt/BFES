/* 
 * MIT License
 * 
 * Copyright (c) 2014 Takuya Akiba
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



#ifndef HYPERLOGLOG_HIP_DISTINCT_COUNTER_H_
#define HYPERLOGLOG_HIP_DISTINCT_COUNTER_H_

#include <algorithm>
#include <cstdint>
#include <cmath>
#include "dense_array.h"

namespace hyperloglog_hip {
template<typename Key, typename Hash = std::hash<Key>, int NumRegisterBits = 5>
class distinct_counter {
 public:
  typedef Key key_type;
  typedef Hash hash_type;

	
//	i = (int) ((1.4427 * std::log(6.320987654320988 * m)) + 0.5)
	distinct_counter(size_t num_memorycells = 4)
	: num_bucket_bits_(std::max((int) ((1.4427 * std::log(6.320987654320988 * num_memorycells)) + 0.5), 4)), M_(1 << num_bucket_bits_),
	c_(0), s_(1 << num_bucket_bits_) {}
	
//  distinct_counter(size_t num_bucket_bits = 12)
//    : num_bucket_bits_(num_bucket_bits), M_(1 << num_bucket_bits),
//      c_(0), s_(1 << num_bucket_bits) {}
	
	size_t size() const { // Returns size of datastructure in bits
		return 2 * 32 + M_.size();
	}

  void insert(const key_type &v) {
    static constexpr uint64_t num_register_bits = NumRegisterBits;
    static constexpr uint64_t register_limit = (uint64_t(1) << num_register_bits) - 1;

    const uint64_t h = hash_(v) * magic1() + magic2();
    const uint64_t h0 = h & ((uint64_t(1) << num_bucket_bits_) - 1);
    const uint64_t h1 = h >> num_bucket_bits_;

    const uint64_t b_old = M_.get(h0);
    const uint64_t b_new = h1 == 0 ? register_limit :
        std::min(register_limit, uint64_t(1 + __builtin_ctzl(h1)));

    if (b_new > b_old) {
      M_.set(h0, b_new);
      c_ += 1.0 / (s_ / (uint64_t(1) << num_bucket_bits_));
      s_ -= 1.0 / (uint64_t(1) << b_old);
      if (b_new < register_limit) {
        s_ += 1.0 / (uint64_t(1) << b_new);
      }
    }
  }

  size_t count() const {
    return round(c_);
  }

 private:
  const size_t num_bucket_bits_;
  dense_array<NumRegisterBits> M_;
  double c_, s_;
  hash_type hash_;

  static constexpr uint64_t magic1() {
    return 9223372036854775837ULL;
  }

  static constexpr uint64_t magic2() {
    return 1234567890123456789ULL;
  }
};
}  // namespace hyperloglog_hip

#endif  // HYPERLOGLOG_HIP_DISTINCT_COUNTER_H_
