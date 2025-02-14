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

#ifndef HYPERLOGLOG_HIP_DENSE_ARRAY_H_
#define HYPERLOGLOG_HIP_DENSE_ARRAY_H_

#include <algorithm>
#include <climits>
#include <cstdint>
#include <memory>
#include <type_traits>

namespace hyperloglog_hip {
template<size_t NumRegisterBits, typename Value = uint8_t>
class dense_array {
 public:
  typedef Value value_type;
	
  dense_array(size_t num_registers)
    : data_(new value_type[data_length(num_registers)]()),
	num_registers(num_registers) {}
	virtual ~dense_array() {}

  value_type get(size_t pos) {
    const size_t b = pos * num_register_bits();
    const size_t i1 = b / num_value_bits();
    const size_t o1 = b - i1 * num_value_bits();
    const size_t n1 = num_value_bits() - o1;
    value_type v = data_[i1] >> o1;

    if (n1 > num_register_bits()) {
      v &= (value_type(1) << num_register_bits()) - 1;
    }
    else if (n1 < num_register_bits()) {
      const size_t i2 = i1 + 1;
      const size_t n2 = num_register_bits() - n1;
      v |= (data_[i2] & ((value_type(1) << n2) - 1)) << n1;
    }
    return v;
  }

  void set(size_t pos, value_type val) {
    const size_t b = pos * num_register_bits();

    const size_t i1 = b / num_value_bits();
    const size_t o1 = b - i1 * num_value_bits();
    const size_t n1 = std::min(num_value_bits() - o1, num_register_bits());
    data_[i1] &= value_type(-1) ^ (((value_type(1) << n1) - 1) << o1);
    data_[i1] |= val << o1;

    if (n1 < num_register_bits()) {
      const size_t i2 = i1 + 1;
      const size_t n2 = num_register_bits() - n1;
      data_[i2] &= value_type(-1) ^ ((value_type(1) << n2) - 1);
      data_[i2] |= val >> n1;
    }
  }
	
	size_t size() const {
		auto a = data_length(num_registers) * sizeof(value_type) * 8;
		return a;
	}

 private:
  std::shared_ptr<value_type[]> data_;
	const size_t num_registers;

  static constexpr size_t num_register_bits() {
    return NumRegisterBits;
  }

  static constexpr size_t num_value_bits() {
    return sizeof(Value) * CHAR_BIT;
  }

  static constexpr size_t data_length(size_t num_registers) {
    return (num_registers * num_register_bits() + num_value_bits() - 1) / num_value_bits();
  }

  static_assert(std::is_unsigned<value_type>::value,
                "Value should be an unsigned integral type.");

  static_assert(sizeof(value_type) * CHAR_BIT >= NumRegisterBits,
                "Value should have at least NumRegisterBits bits.");
};

template<typename Value>
class dense_array_primitive {
 public:
  typedef Value value_type;

  dense_array_primitive(size_t size) : data_(new value_type[size]()) {}
  virtual ~dense_array_primitive() {}

  value_type get(size_t pos) const {
    return data_[pos];
  }

  void set(size_t pos, value_type val) {
    data_[pos] = val;
  }

 private:
  std::shared_ptr<value_type[]> data_;
};

template<>
class dense_array<8, uint8_t> : public dense_array_primitive<uint8_t> {
 public:
  dense_array(size_t size) : dense_array_primitive(size) {}
};

template<>
class dense_array<16, uint16_t> : public dense_array_primitive<uint16_t> {
 public:
  dense_array(size_t size) : dense_array_primitive(size) {}
};

template<>
class dense_array<32, uint32_t> : public dense_array_primitive<uint32_t> {
 public:
  dense_array(size_t size) : dense_array_primitive(size) {}
};

template<>
class dense_array<64, uint64_t> : public dense_array_primitive<uint64_t> {
 public:
  dense_array(size_t size) : dense_array_primitive(size) {}
};
}  // namespace hyperloglog_hip

#endif  // HYEPRLOGLOG_HIP_DENSE_ARRAY_H_
