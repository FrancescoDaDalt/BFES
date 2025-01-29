//
//  DumpContainer.hpp
//  MCMC Sketch
//
//  Created by Francesco Da Dalt on 27.03.22.
//

#ifndef DumpContainer_hpp
#define DumpContainer_hpp

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <vector>
#include <tuple>
#include <string>
#include <ctime>

using namespace std;

// Container for dumping experiment statistics to a csv file
class DumpContainer {
	string dumpname;
	
	int precision;
	
	vector<string> column_names;
	vector<vector<string>> tabulated_data;
	
public:
	DumpContainer(string name, int prec) {
		time_t rawtime;
		struct tm * timeinfo;
		char buffer[80];
		time (&rawtime);
		timeinfo = localtime(&rawtime);
		strftime(buffer, 80, "%Y-%m-%d-%I-%M-%S", timeinfo);
		string a(buffer);
		dumpname = name + "_" + a;
		precision = prec;
		
	};
	
	void insert_row(vector<tuple<string, double>> double_row, vector<tuple<string, string>> string_row) {
		
		if (column_names.size() == 0) {
			for (tuple<string, double> tup : double_row) {
				column_names.push_back(get<0>(tup));
			}
			for (tuple<string, string> tup : string_row) {
				column_names.push_back(get<0>(tup));
			}
		}
		
		int idx = 0;
		
		vector<string> tmp;
		for (tuple<string, double> tup : double_row) {
			assert(get<0>(tup).compare(column_names.at(idx)) == 0);
			stringstream stream;
			stream<<fixed<<setprecision(precision)<<get<1>(tup);
			string s = stream.str();
			tmp.push_back(s);
			idx++;
		}
		for (tuple<string, string> tup : string_row) {
			assert(get<0>(tup).compare(column_names.at(idx)) == 0);
			tmp.push_back(get<1>(tup));
			idx++;
		}
		tabulated_data.push_back(tmp);
	}
	
	void dump_to_csv() {
		ofstream dump_file;
		dump_file.open("Data/Results/" + dumpname + ".csv");
		
		int num_cols = (int) column_names.size();
		
		int idx = 0;
		for (string s: column_names) {
			dump_file<<s;
			if (idx != num_cols - 1) {
				dump_file<<",";
			}
			idx++;
		}
		
		dump_file<<endl;
		
		for (vector<string> strvec : tabulated_data) {
			idx = 0;
			for (string s : strvec) {
				dump_file<<s;
				if (idx != num_cols - 1) {
					dump_file<<",";
				}
				idx++;
			}
			dump_file<<endl;
		}
		dump_file.close();
		
	}
};

#endif /* DumpContainer_hpp */
