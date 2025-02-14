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
 * 
 * Check the references for the algorithm publication and eventual copyrights that apply to it.
 */

#ifndef Stopwatch_hpp
#define Stopwatch_hpp

#include <iostream>
#include <chrono>
#include <string>

class Stopwatch {
	std::chrono::high_resolution_clock::time_point startTime;
	std::chrono::high_resolution_clock::time_point endTime;
	
	std::vector<std::tuple<std::string, double>> data;
	
	std::string name;
	
public:
	
	Stopwatch (std::string name): name(name) {
		
	}
	
	void start() {
		startTime = std::chrono::high_resolution_clock::now();
	}
	
	void end() {
		endTime = std::chrono::high_resolution_clock::now();
	}
	
	auto get_time_data (std::string qualifier) {
		std::cout << "\n-> " +  name + " " + qualifier + " Time" << std::endl;
		data.clear();
		double microsecondstime = ((double) std::chrono::duration<double, std::nano>(endTime - startTime).count()) * 1e-3;
		std::cout << "--> Microseconds: " << microsecondstime << std::endl;
		data.push_back(std::make_tuple("Microseconds_Time_" + name + "_" + qualifier, microsecondstime));
		return &data;
	}
};

#endif /* Stopwatch_hpp */
