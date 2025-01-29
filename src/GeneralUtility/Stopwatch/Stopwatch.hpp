//
//  Stopwatch.hpp
//  braindust
//
//  Created by Francesco Da Dalt on 05.06.24.
//

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
