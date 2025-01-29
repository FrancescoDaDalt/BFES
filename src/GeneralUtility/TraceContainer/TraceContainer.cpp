//
//  TraceContainer.cpp
//  Bayes Sketch
//
//  Created by Francesco Da Dalt on 17.03.22.
//

#include "TraceContainer.hpp"

using namespace std;

/*
 First constructor, which takes in a trace in form of csv file. The csv file must have three colums: Time, Key, Value
 */
TraceContainer::TraceContainer(string path) {
	tracetype = "File";
	trace = path;
	set<int> id_set;
	fstream fin;
	fin.open(path, ios::in);
	string line, word, temp;
	while (fin >> temp) {
		stringstream s(temp);
		tuple<int, int> tmp;
		getline(s, word, ',');
		getline(s, word, ',');
		id_set.insert(stoi(word));
		get<0>(tmp) = stoi(word);
		getline(s, word, ',');
		get<1>(tmp) = stoi(word);
		
		// The Trace class computes a ground truth when reading in the file
		ground_truth[get<0>(tmp)] += get<1>(tmp);
	}
	
	// For easier handling, store groudn truth in vector form
	get<0>(ground_truth_vectorized).reserve(id_set.size());
	get<1>(ground_truth_vectorized).reserve(id_set.size());
	for (int id : id_set) {
		get<0>(ground_truth_vectorized).push_back(id);
		get<1>(ground_truth_vectorized).push_back(ground_truth[id]);
	}
}

// The csv file must have three colums: Time, Key, Value
TraceContainer::TraceContainer(string path, int num_samples, int seed) {
	tracetype = "File_Subsampled";
	trace = path + "_" + to_string(num_samples) + "_RanSamples";
	
	set<int> id_set;
	fstream fin;
	fin.open(path, ios::in);
	string line, word, temp;
	while (fin >> temp) {
		stringstream s(temp);
		tuple<int, int> tmp;
		getline(s, word, ',');
		getline(s, word, ',');
		id_set.insert(stoi(word));
		get<0>(tmp) = stoi(word);
		getline(s, word, ',');
		get<1>(tmp) = stoi(word);
		
		// The Trace class computes a ground truth when reading in the file
		ground_truth[get<0>(tmp)] += get<1>(tmp);
	}
	
	int num_items_in_stream = id_set.size();
	// For easier handling, store groudn truth in vector form
	get<0>(ground_truth_vectorized).reserve(num_items_in_stream);
	get<1>(ground_truth_vectorized).reserve(num_items_in_stream);
	for (int id : id_set) {
		get<0>(ground_truth_vectorized).push_back(id);
		get<1>(ground_truth_vectorized).push_back(ground_truth[id]);
	}
	
//	Shuffle and discard random flows for subsampling
	std::mt19937 eng1(seed);
	auto eng2 = eng1;
	std::shuffle(begin(get<0>(ground_truth_vectorized)), end(get<0>(ground_truth_vectorized)), eng1);
	std::shuffle(begin(get<1>(ground_truth_vectorized)), end(get<1>(ground_truth_vectorized)), eng2);
	
	int num_to_discard = std::max((int) num_items_in_stream - num_samples, 0);
	get<0>(ground_truth_vectorized).resize(num_items_in_stream - num_to_discard);
	get<1>(ground_truth_vectorized).resize(num_items_in_stream - num_to_discard);
	
//	Update ground truth map
	ground_truth.clear();
	for (int i = 0; i < num_items_in_stream - num_to_discard; i++) {
		ground_truth[get<0>(ground_truth_vectorized).at(i)] = get<1>(ground_truth_vectorized).at(i);
	}
	
}

/*
 Second constructor which is used for syntehtic distributions.
 The arguments specifies how many items we want in the stream.
 */
TraceContainer::TraceContainer(int num_ids) {
	get<0>(ground_truth_vectorized).reserve(num_ids);
	for (int i = 0; i < num_ids; i++) {
		get<0>(ground_truth_vectorized).push_back(i);
	}
}

/*
 Synthesizes an exponential distributed trace.
 Items are i.i.d.
 */
void TraceContainer::synthesizeExponential(double scale, double offset, int seed) {
	tracetype = "Synthetic";
	trace = "Exponential_" + to_string(scale) + "_" + to_string(offset);
	ground_truth.clear();
	mt19937 gen(seed);
	exponential_distribution<> d(1);
	get<1>(ground_truth_vectorized).clear();
	for (int id : get<0>(ground_truth_vectorized)) {
		double sample = d(gen) * scale + offset;
		ground_truth[id] = sample;
		get<1>(ground_truth_vectorized).push_back(sample);
	}
}

/*
 Synthesizes a Pareto distributed trace.
 Items are i.i.d.
 */
void TraceContainer::synthesizePareto(double shape, double scale, double offset, int seed) {
	tracetype = "Synthetic";
	trace = "Pareto_" + to_string(shape) + "_" + to_string(scale) + "_" + to_string(offset);
	
	auto rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, seed);
	
	ground_truth.clear();
	get<1>(ground_truth_vectorized).clear();
	for (int id : get<0>(ground_truth_vectorized)) {
		double sample = (gsl_ran_pareto(rng, shape, 1) - 1) * scale + offset;
		ground_truth[id] = sample;
		get<1>(ground_truth_vectorized).push_back(sample);
	}
	
	gsl_rng_free(rng);
}

/*
 Synthesizes a Normal distributed trace of teh kind trunc(N(scale, scale^2) + low, [low, high]) if scale > 0.
 Items are i.i.d.
 */
void TraceContainer::synthesizeNormal(double scale, double low, double high, int seed) {
	assert(scale > 0);
	tracetype = "Synthetic";
	trace = "Normal_" + to_string(scale) + "_" + to_string(low) + "_" + to_string(high);
	
	auto rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, seed);
	
	ground_truth.clear();
	get<1>(ground_truth_vectorized).clear();
	
	double a = (0 - scale) / scale;
	double b = (high - low - scale) / scale;
	for (int id : get<0>(ground_truth_vectorized)) {
		double sample;
		if (a > 3 and b > 3) {
			sample = std::fmod(gsl_ran_ugaussian_tail(rng, a), b - a) * scale + scale + low;
		} else if (a < -3 and b < -3) {
			sample = std::fmod(gsl_ran_ugaussian_tail(rng, -b), b - a) * -scale + scale + low;
		} else {
			double random_unif_sample = gsl_rng_uniform(rng);
			double truncated_gaussian_sample = scale + scale * gsl_cdf_ugaussian_Pinv(gsl_cdf_ugaussian_P(a) + random_unif_sample * (gsl_cdf_ugaussian_P(b) - gsl_cdf_ugaussian_P(a)));
			sample = truncated_gaussian_sample + low;
		}
		
		ground_truth[id] = sample;
		get<1>(ground_truth_vectorized).push_back(sample);
	}
	
	gsl_rng_free(rng);
}


const std::map<int, double>* const TraceContainer::get_ground_truth_mapping() {return &ground_truth;}

const std::tuple<std::vector<int>, std::vector<double>>* const TraceContainer::get_ground_truth_vector() {return &ground_truth_vectorized;}

void TraceContainer::override_item_ids() {
	
	int num_items_in_stream = get<0>(ground_truth_vectorized).size();
	for (int i = 0; i < num_items_in_stream; i++) {
		get<0>(ground_truth_vectorized).at(i) = i;
	}
	
	ground_truth.clear();
	
	if (get<1>(ground_truth_vectorized).size() == num_items_in_stream) {
		for (int i = 0; i < num_items_in_stream; i++) {
			ground_truth[i] = get<1>(ground_truth_vectorized).at(i);
		}
	}
	
	tracetype = tracetype + "_IdOverride";
}

double TraceContainer::normalize() {
	double mean = std::accumulate(get<1>(ground_truth_vectorized).begin(), get<1>(ground_truth_vectorized).end(), 0.0) / (get<1>(ground_truth_vectorized).size() * 10);
	const int numitems = get<1>(ground_truth_vectorized).size();
	for (int i = 0; i < numitems; i++) {
		get<1>(ground_truth_vectorized).at(i) /= mean;
	}
	for (auto it = ground_truth.begin(); it != ground_truth.end(); it++) {
		it->second /= mean;
	}
	tracetype = tracetype + "_normalized10";
	return mean;
}
