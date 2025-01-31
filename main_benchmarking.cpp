//
//  main.cpp
//  BaInDaSt
//
//  Created by Francesco Da Dalt on 25.10.23.
//

#include "ParetoPrior.hpp"
#include "ExponentialPrior.hpp"
#include "NormalPrior.hpp"
#include "BFESketch.hpp"
#include "CBSketch.hpp"
#include "CMSketch.hpp"
#include "DPSketch.hpp"
#include "NIGPSketch.hpp"
#include "CSketch.hpp"
#include "CUSketch.hpp"
#include "PRSketch.hpp"
#ifndef NO_MOSEK
#include "SeqSketch.hpp"
#endif
#include "CCBSketch.hpp"

#include "EstimationAnalysis.hpp"
#include "DumpContainer.hpp"
#include "TraceContainer.hpp"
#include "Stopwatch.hpp"

#include <chrono>
#include <thread>
#include <iostream>
#include <omp.h>

using namespace std;

constexpr int numthreads = 8;

int main(int argc, const char * argv[]) {
	
	auto rng = gsl_rng_alloc(gsl_rng_mt19937);
	std::random_device rd;
	gsl_rng_set(rng, 42);
	omp_set_num_threads(numthreads);
	
	size_t num_memory_cells = 10000;
	size_t num_items = 1000000;
	constexpr size_t w1 = 10;
	constexpr size_t level_1_period = 2;
	constexpr size_t level_1_burn_in = 10;
	constexpr size_t w2 = 10;
	constexpr size_t level_2_period = 2;
	constexpr size_t level_2_burn_in = 10;
	constexpr size_t fPlusOne = 5;
	constexpr size_t numsamples_max = 1000;
	constexpr bool bypassTopLevel = false;
	using hypergrid_type = hypergrid<1>;
	#ifdef NO_MOSEK
	static_assert((hypergrid_type::num_axes == 1) and "If MOSEK is not provided, need to change the number of axes in the hypergrid to 1: using hypergrid_type = hypergrid<1>;");
	#endif
	
	const int numruns = 10;
	const bool use_synthetic_data = true;
	
//	string commondescription = "Pareto";
//	string commondescription = "Exponential";
	string commondescription = "Normal";
//	string commondescription = "Kosarak";
//	string commondescription = "Retail";
//	string commondescription = "MAWI";
//	string commondescription = "MACCDC";
//	string commondescription = "CAIDA";
//	string commondescription = "UNIV";

	
	DumpContainer bfes_par_dump("BFESketch_ParetoPrior_" + hypergrid_type::description() + "_" + commondescription, 10);
	DumpContainer bfes_exp_dump("BFESketch_ExponentialPrior_" + hypergrid_type::description() + "_" + commondescription, 10);
	DumpContainer bfes_norm_dump("BFESketch_NormalPrior_" + hypergrid_type::description() + "_" + commondescription, 10);
	DumpContainer cbs_dump("CBSketch_" + commondescription, 10);
	DumpContainer cms_dump("CMSketch_" + commondescription, 10);
	DumpContainer ds_dump("DPSketch_" + commondescription, 10);
	DumpContainer nigps_dump("NIGPSketch_" + commondescription, 10);
	DumpContainer cs_dump("CSketch_" + commondescription, 10);
	DumpContainer cus_dump("CUSketch_" + commondescription, 10);
	DumpContainer prs_dump("PRSketch_" + commondescription, 10);
	DumpContainer prsn_dump("PRSketch_Native_" + commondescription, 10);
	#ifndef NO_MOSEK
	DumpContainer seqs_dump("SeqSketch_" + commondescription, 10);
	DumpContainer seqsn_dump("SeqSketch_Native_" + commondescription, 10);
	#endif
	DumpContainer ccbs_dump("CCBSketch_" + commondescription, 10);
	DumpContainer ccbsn_dump("CCBSketch_Native_" + commondescription, 10);
	
	for (int exprun = 0; exprun < numruns; exprun++) {
		
		TraceContainer trace(num_items);
		if (!use_synthetic_data) {
//			Real data
//			new (&trace) TraceContainer ("Data/kosarak.csv");
//			new (&trace) TraceContainer ("Data/retail.csv");
//			new (&trace) TraceContainer ("Data/mawi/mawi_200612301400.csv");
//			new (&trace) TraceContainer ("Data/maccdc/maccdc2012_00000.csv");
//			new (&trace) TraceContainer ("Data/caida/caida_448000-1.csv");
//			new (&trace) TraceContainer ("Data/univ/univ1_1.csv");
			
			num_items = trace.get_ground_truth_mapping()->size();
			num_memory_cells = (int) (0.1 * num_items);
		} else {
//			Synthesize data
//			trace.synthesizePareto(3, 10.0, 1e-3, exprun);
//			trace.synthesizeExponential(10.0, 1e-3, exprun);
			trace.synthesizeNormal(10, 1e-3, 1e+6, exprun);
		}
		const double scale_multiplier = trace.normalize();
		
		const auto* const gt_id_freq = trace.get_ground_truth_vector();
		const auto* const gt_ids = &std::get<0>(*gt_id_freq);
		const auto* const gt_freqs = &std::get<1>(*gt_id_freq);
		auto* const trace_stats = trace.getStats();
		
		trace.writeGroundTruth(commondescription);
		
		
		
//		using level1_priortype = ParetoPrior<double>;
//		level1_priortype level1_prior(1e-3, 1e10, 10.0, 5);
		
//		using level1_priortype = ExponentialPrior<double>;
//		level1_priortype level1_prior(1e-3, 1e10, 1.0);
		
		using level1_priortype = NormalPrior<double>;
		level1_priortype level1_prior(1e-6, 1e6, 1.0, 1.0);
		
		using level2_priortype_pareto = ParetoPrior<double>;
		level2_priortype_pareto level2_prior_pareto(1e-6, 1e6, 10.0, 3);
		
		using level2_priortype_exp = ExponentialPrior<double>;
		level2_priortype_exp  level2_prior_pareto_exp(1e-6, 1e6, 1.0);
		
		using level2_priortype_norm = NormalPrior<double>;
		level2_priortype_norm level2_prior_pareto_norm(1e-6, 1e6, 1.0, 1.0);

//		BFES (Pareto)
//		Initialization
		BFESketch<hypergrid_type, int, double, numthreads, fPlusOne> bfes_par(num_memory_cells, numsamples_max, numsamples_max, bypassTopLevel);
		Stopwatch bfes_par_stopwatch("BFES (Pareto)");
		vector<tuple<string, double>> bfes_par_dump_numvec;
		vector<tuple<string, string>> bfes_par_dump_strvec;
		bfes_par_dump_numvec.push_back(std::make_tuple("RunNr", (double) exprun));
		bfes_par_dump_strvec.insert(bfes_par_dump_strvec.end(), trace_stats->begin(), trace_stats->end());
//		Insert
		bfes_par_stopwatch.start();
		for (int i = 0; i < num_items; i++) {bfes_par.insert(std::make_tuple(gt_ids->at(i), gt_freqs->at(i)));}
		bfes_par_stopwatch.end();
		auto* const bfes_par_insert_stats = bfes_par_stopwatch.get_time_data("Insert");
		bfes_par_dump_numvec.insert(bfes_par_dump_numvec.end(), bfes_par_insert_stats->begin(), bfes_par_insert_stats->end());
//		Query
		bfes_par_stopwatch.start();
		const auto bfes_par_query_ans = bfes_par.batch_query(gt_ids, &level1_prior, &level2_prior_pareto, w1, w2, level_1_period, level_1_burn_in, level_2_period, level_2_burn_in);
		bfes_par_stopwatch.end();
		auto* const bfes_par_query_stats = bfes_par_stopwatch.get_time_data("BatchQuery");
		bfes_par_dump_numvec.insert(bfes_par_dump_numvec.end(), bfes_par_query_stats->begin(), bfes_par_query_stats->end());
//		Analyze
		EstimationAnalysis bfes_par_analysis("BFES (Pareto)", bfes_par_query_ans, gt_freqs);
		auto* const bfes_par_analysis_stats = bfes_par_analysis.getStats();
		auto* const bfes_par_memory_stats = bfes_par.getMemStats();
		bfes_par_dump_numvec.insert(bfes_par_dump_numvec.end(), bfes_par_analysis_stats->begin(), bfes_par_analysis_stats->end());
		bfes_par_dump_numvec.insert(bfes_par_dump_numvec.end(), bfes_par_memory_stats->begin(), bfes_par_memory_stats->end());
//		Cleanup
		bfes_par_dump.insert_row(bfes_par_dump_numvec, bfes_par_dump_strvec);
		delete get<0>(bfes_par_query_ans);
		delete get<1>(bfes_par_query_ans);
		if (get<1>(bfes_par_query_ans) != get<2>(bfes_par_query_ans)) {delete get<2>(bfes_par_query_ans);}

//		BFES (Exp)
//		Initialization
		BFESketch<hypergrid_type, int, double, numthreads, fPlusOne> bfes_exp(num_memory_cells, numsamples_max, numsamples_max, bypassTopLevel);
		Stopwatch bfes_exp_stopwatch("BFES (Exponential)");
		vector<tuple<string, double>> bfes_exp_dump_numvec;
		vector<tuple<string, string>> bfes_exp_dump_strvec;
		bfes_exp_dump_numvec.push_back(std::make_tuple("RunNr", (double) exprun));
		bfes_exp_dump_strvec.insert(bfes_exp_dump_strvec.end(), trace_stats->begin(), trace_stats->end());
//		Insert
		bfes_exp_stopwatch.start();
		for (int i = 0; i < num_items; i++) {bfes_exp.insert(std::make_tuple(gt_ids->at(i), gt_freqs->at(i)));}
		bfes_exp_stopwatch.end();
		auto* const bfes_exp_insert_stats = bfes_exp_stopwatch.get_time_data("Insert");
		bfes_exp_dump_numvec.insert(bfes_exp_dump_numvec.end(), bfes_exp_insert_stats->begin(), bfes_exp_insert_stats->end());
//		Query
		bfes_exp_stopwatch.start();
		const auto bfes_exp_query_ans = bfes_exp.batch_query(gt_ids, &level1_prior, &level2_prior_pareto_exp, w1, w2, level_1_period, level_1_burn_in, level_2_period, level_2_burn_in);
		bfes_exp_stopwatch.end();
		auto* const bfes_exp_query_stats = bfes_exp_stopwatch.get_time_data("BatchQuery");
		bfes_exp_dump_numvec.insert(bfes_exp_dump_numvec.end(), bfes_exp_query_stats->begin(), bfes_exp_query_stats->end());
//		Analyze
		EstimationAnalysis bfes_exp_analysis("BFES (Exponential)", bfes_exp_query_ans, gt_freqs);
		auto* const bfes_exp_analysis_stats = bfes_exp_analysis.getStats();
		auto* const bfes_exp_memory_stats = bfes_exp.getMemStats();
		bfes_exp_dump_numvec.insert(bfes_exp_dump_numvec.end(), bfes_exp_analysis_stats->begin(), bfes_exp_analysis_stats->end());
		bfes_exp_dump_numvec.insert(bfes_exp_dump_numvec.end(), bfes_exp_memory_stats->begin(), bfes_exp_memory_stats->end());
//		Cleanup
		bfes_exp_dump.insert_row(bfes_exp_dump_numvec, bfes_exp_dump_strvec);
		delete get<0>(bfes_exp_query_ans);
		delete get<1>(bfes_exp_query_ans);
		if (get<1>(bfes_exp_query_ans) != get<2>(bfes_exp_query_ans)) {delete get<2>(bfes_exp_query_ans);}


//		BFES (Norm)
//		Initialization
		BFESketch<hypergrid_type, int, double, numthreads, fPlusOne> bfes_norm(num_memory_cells, numsamples_max, numsamples_max, bypassTopLevel);
		Stopwatch bfes_norm_stopwatch("BFES (Normal)");
		vector<tuple<string, double>> bfes_norm_dump_numvec;
		vector<tuple<string, string>> bfes_norm_dump_strvec;
		bfes_norm_dump_numvec.push_back(std::make_tuple("RunNr", (double) exprun));
		bfes_norm_dump_strvec.insert(bfes_norm_dump_strvec.end(), trace_stats->begin(), trace_stats->end());
//		Insert
		bfes_norm_stopwatch.start();
		for (int i = 0; i < num_items; i++) {bfes_norm.insert(std::make_tuple(gt_ids->at(i), gt_freqs->at(i)));}
		bfes_norm_stopwatch.end();
		auto* const bfes_norm_insert_stats = bfes_norm_stopwatch.get_time_data("Insert");
		bfes_norm_dump_numvec.insert(bfes_norm_dump_numvec.end(), bfes_norm_insert_stats->begin(), bfes_norm_insert_stats->end());
//		Query
		bfes_norm_stopwatch.start();
		const auto bfes_norm_query_ans = bfes_norm.batch_query(gt_ids, &level1_prior, &level2_prior_pareto_norm, w1, w2, level_1_period, level_1_burn_in, level_2_period, level_2_burn_in);
		bfes_norm_stopwatch.end();
		auto* const bfes_norm_query_stats = bfes_norm_stopwatch.get_time_data("BatchQuery");
		bfes_norm_dump_numvec.insert(bfes_norm_dump_numvec.end(), bfes_norm_query_stats->begin(), bfes_norm_query_stats->end());
//		Analyze
		EstimationAnalysis bfes_norm_analysis("BFES (Normal)", bfes_norm_query_ans, gt_freqs);
		auto* const bfes_norm_analysis_stats = bfes_norm_analysis.getStats();
		auto* const bfes_norm_memory_stats = bfes_norm.getMemStats();
		bfes_norm_dump_numvec.insert(bfes_norm_dump_numvec.end(), bfes_norm_analysis_stats->begin(), bfes_norm_analysis_stats->end());
		bfes_norm_dump_numvec.insert(bfes_norm_dump_numvec.end(), bfes_norm_memory_stats->begin(), bfes_norm_memory_stats->end());
//		Cleanup
		bfes_norm_dump.insert_row(bfes_norm_dump_numvec, bfes_norm_dump_strvec);
		delete get<0>(bfes_norm_query_ans);
		delete get<1>(bfes_norm_query_ans);
		if (get<1>(bfes_norm_query_ans) != get<2>(bfes_norm_query_ans)) {delete get<2>(bfes_norm_query_ans);}
		
//		CBS
//		Initialization
		CBSketch<int, double> cbs(num_memory_cells);
		Stopwatch cbs_stopwatch("CBS");
		vector<tuple<string, double>> cbs_dump_numvec;
		vector<tuple<string, string>> cbs_dump_strvec;
		cbs_dump_numvec.push_back(std::make_tuple("RunNr", (double) exprun));
		cbs_dump_strvec.insert(cbs_dump_strvec.end(), trace_stats->begin(), trace_stats->end());
//		Insert
		cbs_stopwatch.start();
		for (int i = 0; i < num_items; i++) {cbs.insert(std::make_tuple(gt_ids->at(i), gt_freqs->at(i)));}
		cbs_stopwatch.end();
		auto* const cbs_insert_stats = cbs_stopwatch.get_time_data("Insert");
		cbs_dump_numvec.insert(cbs_dump_numvec.end(), cbs_insert_stats->begin(), cbs_insert_stats->end());
//		Query
		cbs_stopwatch.start();
		const auto cbs_query_ans = cbs.batch_query(gt_ids);
		cbs_stopwatch.end();
		auto* const bcs_query_stats = cbs_stopwatch.get_time_data("BatchQuery");
		cbs_dump_numvec.insert(cbs_dump_numvec.end(), bcs_query_stats->begin(), bcs_query_stats->end());
//		Analyze
		EstimationAnalysis cbs_analysis("CBS", cbs_query_ans, gt_freqs);
		auto* const cbs_analysis_stats = cbs_analysis.getStats();
		auto* const cbs_memory_stats = cbs.getMemStats();
		cbs_dump_numvec.insert(cbs_dump_numvec.end(), cbs_analysis_stats->begin(), cbs_analysis_stats->end());
		cbs_dump_numvec.insert(cbs_dump_numvec.end(), cbs_memory_stats->begin(), cbs_memory_stats->end());
//		Cleanup
		cbs_dump.insert_row(cbs_dump_numvec, cbs_dump_strvec);
		delete get<0>(cbs_query_ans);
		delete get<1>(cbs_query_ans);
		if (get<1>(cbs_query_ans) != get<2>(cbs_query_ans)) {delete get<2>(cbs_query_ans);}
		
//		CMS
//		Initialization
		CMSketch<int, double> cms(num_memory_cells);
		Stopwatch cms_stopwatch("CMS");
		vector<tuple<string, double>> cms_dump_numvec;
		vector<tuple<string, string>> cms_dump_strvec;
		cms_dump_numvec.push_back(std::make_tuple("RunNr", (double) exprun));
		cms_dump_strvec.insert(cms_dump_strvec.end(), trace_stats->begin(), trace_stats->end());
//		Insert
		cms_stopwatch.start();
		for (int i = 0; i < num_items; i++) {cms.insert(std::make_tuple(gt_ids->at(i), gt_freqs->at(i)));}
		cms_stopwatch.end();
		auto* const cms_insert_stats = cms_stopwatch.get_time_data("Insert");
		cms_dump_numvec.insert(cms_dump_numvec.end(), cms_insert_stats->begin(), cms_insert_stats->end());
//		Query
		cms_stopwatch.start();
		const auto cms_query_ans = cms.batch_query(gt_ids);
		cms_stopwatch.end();
		auto* const cms_query_stats = cms_stopwatch.get_time_data("BatchQuery");
		cms_dump_numvec.insert(cms_dump_numvec.end(), cms_query_stats->begin(), cms_query_stats->end());
//		Analyze
		EstimationAnalysis cms_analysis("CMS", cms_query_ans, gt_freqs);
		auto* const cms_analysis_stats = cms_analysis.getStats();
		auto* const cms_memory_stats = cms.getMemStats();
		cms_dump_numvec.insert(cms_dump_numvec.end(), cms_analysis_stats->begin(), cms_analysis_stats->end());
		cms_dump_numvec.insert(cms_dump_numvec.end(), cms_memory_stats->begin(), cms_memory_stats->end());
//		Cleanup
		cms_dump.insert_row(cms_dump_numvec, cms_dump_strvec);
		delete get<0>(cms_query_ans);
		delete get<1>(cms_query_ans);
		if (get<1>(cms_query_ans) != get<2>(cms_query_ans)) {delete get<2>(cms_query_ans);}
		
//		DPS
//		Initialization
		DPSketch<int, double> ds(num_memory_cells, w1 * w2, scale_multiplier);
		Stopwatch ds_stopwatch("DPS");
		vector<tuple<string, double>> ds_dump_numvec;
		vector<tuple<string, string>> ds_dump_strvec;
		ds_dump_numvec.push_back(std::make_tuple("RunNr", (double) exprun));
		ds_dump_strvec.insert(ds_dump_strvec.end(), trace_stats->begin(), trace_stats->end());
//		Insert
		ds_stopwatch.start();
		for (int i = 0; i < num_items; i++) {ds.insert(std::make_tuple(gt_ids->at(i), gt_freqs->at(i)));}
		ds_stopwatch.end();
		auto* const ds_insert_stats = ds_stopwatch.get_time_data("Insert");
		ds_dump_numvec.insert(ds_dump_numvec.end(), ds_insert_stats->begin(), ds_insert_stats->end());
//		Query
		ds_stopwatch.start();
		const auto ds_query_ans = ds.batch_query(gt_ids);
		ds_stopwatch.end();
		auto* const ds_query_stats = ds_stopwatch.get_time_data("BatchQuery");
		ds_dump_numvec.insert(ds_dump_numvec.end(), ds_query_stats->begin(), ds_query_stats->end());
//		Analyze
		EstimationAnalysis ds_analysis("DPS", ds_query_ans, gt_freqs);
		auto* const ds_analysis_stats = ds_analysis.getStats();
		auto* const ds_memory_stats = ds.getMemStats();
		ds_dump_numvec.insert(ds_dump_numvec.end(), ds_analysis_stats->begin(), ds_analysis_stats->end());
		ds_dump_numvec.insert(ds_dump_numvec.end(), ds_memory_stats->begin(), ds_memory_stats->end());
//		Cleanup
		ds_dump.insert_row(ds_dump_numvec, ds_dump_strvec);
		delete get<0>(ds_query_ans);
		delete get<1>(ds_query_ans);
		if (get<1>(ds_query_ans) != get<2>(ds_query_ans)) {delete get<2>(ds_query_ans);}
		
		
//		There is no hope for this working on big streams. Evaluating the modified bessel function a couple times takes longer than running all queries for other sketches. The modified bessel function has to be evaluated many thousand times in order to find alpha. Therefore the query operation on this sketch is so slow that we regard it as a time-out.
////		NIGPS
////		Initialization
//		NIGPSketch<int, double> nigps(num_memory_cells, w1 * w2, scale_multiplier);
//		Stopwatch nigps_stopwatch("NIGPS");
//		vector<tuple<string, double>> nigps_dump_numvec;
//		vector<tuple<string, string>> nigps_dump_strvec;
//		nigps_dump_numvec.push_back(std::make_tuple("RunNr", (double) exprun));
//		nigps_dump_strvec.insert(nigps_dump_strvec.end(), trace_stats->begin(), trace_stats->end());
////		Insert
//		nigps_stopwatch.start();
//		for (int i = 0; i < num_items; i++) {nigps.insert(std::make_tuple(gt_ids->at(i), gt_freqs->at(i)));}
//		nigps_stopwatch.end();
//		auto* const nigps_insert_stats = nigps_stopwatch.get_time_data("Insert");
//		nigps_dump_numvec.insert(nigps_dump_numvec.end(), nigps_insert_stats->begin(), nigps_insert_stats->end());
////		Query
//		nigps_stopwatch.start();
//		const auto nigps_query_ans = nigps.batch_query(gt_ids);
//		nigps_stopwatch.end();
//		auto* const nigps_query_stats = nigps_stopwatch.get_time_data("BatchQuery");
//		nigps_dump_numvec.insert(nigps_dump_numvec.end(), nigps_query_stats->begin(), nigps_query_stats->end());
////		Analyze
//		EstimationAnalysis nigps_analysis("NIGPS", nigps_query_ans, gt_freqs);
//		auto* const nigps_analysis_stats = nigps_analysis.getStats();
//		auto* const nigps_memory_stats = nigps.getMemStats();
//		nigps_dump_numvec.insert(nigps_dump_numvec.end(), nigps_analysis_stats->begin(), nigps_analysis_stats->end());
//		nigps_dump_numvec.insert(nigps_dump_numvec.end(), nigps_memory_stats->begin(), nigps_memory_stats->end());
////		Cleanup
//		nigps_dump.insert_row(nigps_dump_numvec, nigps_dump_strvec);
//		delete get<0>(nigps_query_ans);
//		delete get<1>(nigps_query_ans);
//		if (get<1>(nigps_query_ans) != get<2>(nigps_query_ans)) {delete get<2>(nigps_query_ans);}
		
//		CS
//		Initialization
		CSketch<int, double> cs(num_memory_cells);
		Stopwatch cs_stopwatch("CS");
		vector<tuple<string, double>> cs_dump_numvec;
		vector<tuple<string, string>> cs_dump_strvec;
		cs_dump_numvec.push_back(std::make_tuple("RunNr", (double) exprun));
		cs_dump_strvec.insert(cs_dump_strvec.end(), trace_stats->begin(), trace_stats->end());
//		Insert
		cs_stopwatch.start();
		for (int i = 0; i < num_items; i++) {cs.insert(std::make_tuple(gt_ids->at(i), gt_freqs->at(i)));}
		cs_stopwatch.end();
		auto* const cs_insert_stats = cs_stopwatch.get_time_data("Insert");
		cs_dump_numvec.insert(cs_dump_numvec.end(), cs_insert_stats->begin(), cs_insert_stats->end());
//		Query
		cs_stopwatch.start();
		const auto cs_query_ans = cs.batch_query(gt_ids);
		cs_stopwatch.end();
		auto* const cs_query_stats = cs_stopwatch.get_time_data("BatchQuery");
		cs_dump_numvec.insert(cs_dump_numvec.end(), cs_query_stats->begin(), cs_query_stats->end());
//		Analyze
		EstimationAnalysis cs_analysis("CS", cs_query_ans, gt_freqs);
		auto* const cs_analysis_stats = cs_analysis.getStats();
		auto* const cs_memory_stats = cs.getMemStats();
		cs_dump_numvec.insert(cs_dump_numvec.end(), cs_analysis_stats->begin(), cs_analysis_stats->end());
		cs_dump_numvec.insert(cs_dump_numvec.end(), cs_memory_stats->begin(), cs_memory_stats->end());
//		Cleanup
		cs_dump.insert_row(cs_dump_numvec, cs_dump_strvec);
		delete get<0>(cs_query_ans);
		delete get<1>(cs_query_ans);
		if (get<1>(cs_query_ans) != get<2>(cs_query_ans)) {delete get<2>(cs_query_ans);}
		
//		CUS
//		Initialization
		CUSketch<int, double> cus(num_memory_cells);
		Stopwatch cus_stopwatch("CUS");
		vector<tuple<string, double>> cus_dump_numvec;
		vector<tuple<string, string>> cus_dump_strvec;
		cus_dump_numvec.push_back(std::make_tuple("RunNr", (double) exprun));
		cus_dump_strvec.insert(cus_dump_strvec.end(), trace_stats->begin(), trace_stats->end());
//		Insert
		cus_stopwatch.start();
		for (int i = 0; i < num_items; i++) {cus.insert(std::make_tuple(gt_ids->at(i), gt_freqs->at(i)));}
		cus_stopwatch.end();
		auto* const cus_insert_stats = cus_stopwatch.get_time_data("Insert");
		cus_dump_numvec.insert(cus_dump_numvec.end(), cus_insert_stats->begin(), cus_insert_stats->end());
//		Query
		cus_stopwatch.start();
		const auto cus_query_ans = cus.batch_query(gt_ids);
		cus_stopwatch.end();
		auto* const cus_query_stats = cus_stopwatch.get_time_data("BatchQuery");
		cus_dump_numvec.insert(cus_dump_numvec.end(), cus_query_stats->begin(), cus_query_stats->end());
//		Analyze
		EstimationAnalysis cus_analysis("CUS", cus_query_ans, gt_freqs);
		auto* const cus_analysis_stats = cus_analysis.getStats();
		auto* const cus_memory_stats = cus.getMemStats();
		cus_dump_numvec.insert(cus_dump_numvec.end(), cus_analysis_stats->begin(), cus_analysis_stats->end());
		cus_dump_numvec.insert(cus_dump_numvec.end(), cus_memory_stats->begin(), cus_memory_stats->end());
//		Cleanup
		cus_dump.insert_row(cus_dump_numvec, cus_dump_strvec);
		delete get<0>(cus_query_ans);
		delete get<1>(cus_query_ans);
		if (get<1>(cus_query_ans) != get<2>(cus_query_ans)) {delete get<2>(cus_query_ans);}
		
//		PRS
//		Initialization
		PRSketch<int, double> prs(num_memory_cells);
		Stopwatch prs_stopwatch("PRS");
		vector<tuple<string, double>> prs_dump_numvec;
		vector<tuple<string, string>> prs_dump_strvec;
		prs_dump_numvec.push_back(std::make_tuple("RunNr", (double) exprun));
		prs_dump_strvec.insert(prs_dump_strvec.end(), trace_stats->begin(), trace_stats->end());
//		Insert
		prs_stopwatch.start();
		for (int i = 0; i < num_items; i++) {prs.insert(std::make_tuple(gt_ids->at(i), gt_freqs->at(i)));}
		prs_stopwatch.end();
		auto* const prs_insert_stats = prs_stopwatch.get_time_data("Insert");
		prs_dump_numvec.insert(prs_dump_numvec.end(), prs_insert_stats->begin(), prs_insert_stats->end());
//		Query
		prs_stopwatch.start();
		const auto prs_query_ans = prs.batch_query(gt_ids);
		prs_stopwatch.end();
		auto* const prs_query_stats = prs_stopwatch.get_time_data("BatchQuery");
		prs_dump_numvec.insert(prs_dump_numvec.end(), prs_query_stats->begin(), prs_query_stats->end());
//		Analyze
		EstimationAnalysis prs_analysis("PRS", prs_query_ans, gt_freqs);
		auto* const prs_analysis_stats = prs_analysis.getStats();
		auto* const prs_memory_stats = prs.getMemStats();
		prs_dump_numvec.insert(prs_dump_numvec.end(), prs_analysis_stats->begin(), prs_analysis_stats->end());
		prs_dump_numvec.insert(prs_dump_numvec.end(), prs_memory_stats->begin(), prs_memory_stats->end());
//		Cleanup
		prs_dump.insert_row(prs_dump_numvec, prs_dump_strvec);
		delete get<0>(prs_query_ans);
		delete get<1>(prs_query_ans);
		if (get<1>(prs_query_ans) != get<2>(prs_query_ans)) {delete get<2>(prs_query_ans);}

//		PRSN
//		Initialization
		PRSketch_Native<int, double> prsn(num_memory_cells);
		Stopwatch prsn_stopwatch("PRSN");
		vector<tuple<string, double>> prsn_dump_numvec;
		vector<tuple<string, string>> prsn_dump_strvec;
		prsn_dump_numvec.push_back(std::make_tuple("RunNr", (double) exprun));
		prsn_dump_strvec.insert(prsn_dump_strvec.end(), trace_stats->begin(), trace_stats->end());
//		Insert
		prsn_stopwatch.start();
		for (int i = 0; i < num_items; i++) {prsn.insert(std::make_tuple(gt_ids->at(i), gt_freqs->at(i)));}
		prsn_stopwatch.end();
		auto* const prsn_insert_stats = prsn_stopwatch.get_time_data("Insert");
		prsn_dump_numvec.insert(prsn_dump_numvec.end(), prsn_insert_stats->begin(), prsn_insert_stats->end());
//		Query
		prsn_stopwatch.start();
		const auto prsn_query_ans = prsn.batch_query(gt_ids);
		prsn_stopwatch.end();
		auto* const prsn_query_stats = prsn_stopwatch.get_time_data("BatchQuery");
		prsn_dump_numvec.insert(prsn_dump_numvec.end(), prsn_query_stats->begin(), prsn_query_stats->end());
//		Analyze
		EstimationAnalysis prsn_analysis("PRSN", prsn_query_ans, gt_freqs);
		auto* const prsn_analysis_stats = prsn_analysis.getStats();
		auto* const prsn_memory_stats = prsn.getMemStats();
		prsn_dump_numvec.insert(prsn_dump_numvec.end(), prsn_analysis_stats->begin(), prsn_analysis_stats->end());
		prsn_dump_numvec.insert(prsn_dump_numvec.end(), prsn_memory_stats->begin(), prsn_memory_stats->end());
//		Cleanup
		prsn_dump.insert_row(prsn_dump_numvec, prsn_dump_strvec);
		delete get<0>(prsn_query_ans);
		delete get<1>(prsn_query_ans);
		if (get<1>(prsn_query_ans) != get<2>(prsn_query_ans)) {delete get<2>(prsn_query_ans);}

		#ifndef NO_MOSEK
//		SeqS
//		Initialization
		SeqSketch<int, double> seqs(num_memory_cells);
		Stopwatch seqs_stopwatch("SeqS");
		vector<tuple<string, double>> seqs_dump_numvec;
		vector<tuple<string, string>> seqs_dump_strvec;
		seqs_dump_numvec.push_back(std::make_tuple("RunNr", (double) exprun));
		seqs_dump_strvec.insert(seqs_dump_strvec.end(), trace_stats->begin(), trace_stats->end());
//		Insert
		seqs_stopwatch.start();
		for (int i = 0; i < num_items; i++) {seqs.insert(std::make_tuple(gt_ids->at(i), gt_freqs->at(i)));}
		seqs_stopwatch.end();
		auto* const seqs_insert_stats = seqs_stopwatch.get_time_data("Insert");
		seqs_dump_numvec.insert(seqs_dump_numvec.end(), seqs_insert_stats->begin(), seqs_insert_stats->end());
//		Query
		seqs_stopwatch.start();
		const auto seqs_query_ans = seqs.batch_query(gt_ids);
		seqs_stopwatch.end();
		auto* const seqs_query_stats = seqs_stopwatch.get_time_data("BatchQuery");
		seqs_dump_numvec.insert(seqs_dump_numvec.end(), seqs_query_stats->begin(), seqs_query_stats->end());
//		Analyze
		EstimationAnalysis seqs_analysis("SeqS", seqs_query_ans, gt_freqs);
		auto* const seqs_analysis_stats = seqs_analysis.getStats();
		auto* const seqs_memory_stats = seqs.getMemStats();
		seqs_dump_numvec.insert(seqs_dump_numvec.end(), seqs_analysis_stats->begin(), seqs_analysis_stats->end());
		seqs_dump_numvec.insert(seqs_dump_numvec.end(), seqs_memory_stats->begin(), seqs_memory_stats->end());
//		Cleanup
		seqs_dump.insert_row(seqs_dump_numvec, seqs_dump_strvec);
		delete get<0>(seqs_query_ans);
		delete get<1>(seqs_query_ans);
		if (get<1>(seqs_query_ans) != get<2>(seqs_query_ans)) {delete get<2>(seqs_query_ans);}

//		SeqSN
//		Initialization
		SeqSketch_Native<int, double> seqsn(num_memory_cells);
		Stopwatch seqsn_stopwatch("SeqSN");
		vector<tuple<string, double>> seqsn_dump_numvec;
		vector<tuple<string, string>> seqsn_dump_strvec;
		seqsn_dump_numvec.push_back(std::make_tuple("RunNr", (double) exprun));
		seqsn_dump_strvec.insert(seqsn_dump_strvec.end(), trace_stats->begin(), trace_stats->end());
//		Insert
		seqsn_stopwatch.start();
		for (int i = 0; i < num_items; i++) {seqsn.insert(std::make_tuple(gt_ids->at(i), gt_freqs->at(i)));}
		seqsn_stopwatch.end();
		auto* const seqsn_insert_stats = seqsn_stopwatch.get_time_data("Insert");
		seqsn_dump_numvec.insert(seqsn_dump_numvec.end(), seqsn_insert_stats->begin(), seqsn_insert_stats->end());
//		Query
		seqsn_stopwatch.start();
		const auto seqsn_query_ans = seqsn.batch_query(gt_ids);
		seqsn_stopwatch.end();
		auto* const seqsn_query_stats = seqsn_stopwatch.get_time_data("BatchQuery");
		seqsn_dump_numvec.insert(seqsn_dump_numvec.end(), seqsn_query_stats->begin(), seqsn_query_stats->end());
//		Analyze
		EstimationAnalysis seqsn_analysis("SeqSN", seqsn_query_ans, gt_freqs);
		auto* const seqsn_analysis_stats = seqsn_analysis.getStats();
		auto* const seqsn_memory_stats = seqsn.getMemStats();
		seqsn_dump_numvec.insert(seqsn_dump_numvec.end(), seqsn_analysis_stats->begin(), seqsn_analysis_stats->end());
		seqsn_dump_numvec.insert(seqsn_dump_numvec.end(), seqsn_memory_stats->begin(), seqsn_memory_stats->end());
//		Cleanup
		seqsn_dump.insert_row(seqsn_dump_numvec, seqsn_dump_strvec);
		delete get<0>(seqsn_query_ans);
		delete get<1>(seqsn_query_ans);
		if (get<1>(seqsn_query_ans) != get<2>(seqsn_query_ans)) {delete get<2>(seqsn_query_ans);}
		#endif
//		CCBS
//		Initialization
		CCBSketch<int, double> ccbs(num_memory_cells);
		Stopwatch ccbs_stopwatch("CCBS");
		vector<tuple<string, double>> ccbs_dump_numvec;
		vector<tuple<string, string>> ccbs_dump_strvec;
		ccbs_dump_numvec.push_back(std::make_tuple("RunNr", (double) exprun));
		ccbs_dump_strvec.insert(ccbs_dump_strvec.end(), trace_stats->begin(), trace_stats->end());
//		Insert
		ccbs_stopwatch.start();
		for (int i = 0; i < num_items; i++) {ccbs.insert(std::make_tuple(gt_ids->at(i), gt_freqs->at(i)));}
		ccbs_stopwatch.end();
		auto* const ccbs_insert_stats = ccbs_stopwatch.get_time_data("Insert");
		ccbs_dump_numvec.insert(ccbs_dump_numvec.end(), ccbs_insert_stats->begin(), ccbs_insert_stats->end());
//		Query
		ccbs_stopwatch.start();
		const auto ccbs_query_ans = ccbs.batch_query(gt_ids);
		ccbs_stopwatch.end();
		auto* const ccbs_query_stats = ccbs_stopwatch.get_time_data("BatchQuery");
		ccbs_dump_numvec.insert(ccbs_dump_numvec.end(), ccbs_query_stats->begin(), ccbs_query_stats->end());
//		Analyze
		EstimationAnalysis ccbs_analysis("CCBS", ccbs_query_ans, gt_freqs);
		auto* const ccbs_analysis_stats = ccbs_analysis.getStats();
		auto* const ccbs_memory_stats = ccbs.getMemStats();
		ccbs_dump_numvec.insert(ccbs_dump_numvec.end(), ccbs_analysis_stats->begin(), ccbs_analysis_stats->end());
		ccbs_dump_numvec.insert(ccbs_dump_numvec.end(), ccbs_memory_stats->begin(), ccbs_memory_stats->end());
//		Cleanup
		ccbs_dump.insert_row(ccbs_dump_numvec, ccbs_dump_strvec);
		delete get<0>(ccbs_query_ans);
		delete get<1>(ccbs_query_ans);
		if (get<1>(ccbs_query_ans) != get<2>(ccbs_query_ans)) {delete get<2>(ccbs_query_ans);}
		
//		CCBSN
//		Initialization
		CCBSketch_Native<int, double> ccbsn(num_memory_cells);
		Stopwatch ccbsn_stopwatch("CCBSN");
		vector<tuple<string, double>> ccbsn_dump_numvec;
		vector<tuple<string, string>> ccbsn_dump_strvec;
		ccbsn_dump_numvec.push_back(std::make_tuple("RunNr", (double) exprun));
		ccbsn_dump_strvec.insert(ccbsn_dump_strvec.end(), trace_stats->begin(), trace_stats->end());
//		Insert
		ccbsn_stopwatch.start();
		for (int i = 0; i < num_items; i++) {ccbsn.insert(std::make_tuple(gt_ids->at(i), gt_freqs->at(i)));}
		ccbsn_stopwatch.end();
		auto* const ccbsn_insert_stats = ccbsn_stopwatch.get_time_data("Insert");
		ccbsn_dump_numvec.insert(ccbsn_dump_numvec.end(), ccbsn_insert_stats->begin(), ccbsn_insert_stats->end());
//		Query
		ccbsn_stopwatch.start();
		const auto ccbsn_query_ans = ccbsn.batch_query(gt_ids);
		ccbsn_stopwatch.end();
		auto* const ccbsn_query_stats = ccbsn_stopwatch.get_time_data("BatchQuery");
		ccbsn_dump_numvec.insert(ccbsn_dump_numvec.end(), ccbsn_query_stats->begin(), ccbsn_query_stats->end());
//		Analyze
		EstimationAnalysis ccbsn_analysis("CCBSN", ccbsn_query_ans, gt_freqs);
		auto* const ccbsn_analysis_stats = ccbsn_analysis.getStats();
		auto* const ccbsn_memory_stats = ccbsn.getMemStats();
		ccbsn_dump_numvec.insert(ccbsn_dump_numvec.end(), ccbsn_analysis_stats->begin(), ccbsn_analysis_stats->end());
		ccbsn_dump_numvec.insert(ccbsn_dump_numvec.end(), ccbsn_memory_stats->begin(), ccbsn_memory_stats->end());
//		Cleanup
		ccbsn_dump.insert_row(ccbsn_dump_numvec, ccbsn_dump_strvec);
		delete get<0>(ccbsn_query_ans);
		delete get<1>(ccbsn_query_ans);
		if (get<1>(ccbsn_query_ans) != get<2>(ccbsn_query_ans)) {delete get<2>(ccbsn_query_ans);}
		
	}
	
	
	bfes_par_dump.dump_to_csv();
	bfes_exp_dump.dump_to_csv();
	bfes_norm_dump.dump_to_csv();
//	cbs_dump.dump_to_csv();
//	cms_dump.dump_to_csv();
//	ds_dump.dump_to_csv();
//	cs_dump.dump_to_csv();
//	cus_dump.dump_to_csv();
//	prs_dump.dump_to_csv();
//	prsn_dump.dump_to_csv();
#ifndef NO_MOSEK
//	seqs_dump.dump_to_csv();
//	seqsn_dump.dump_to_csv();
#endif
//	ccbs_dump.dump_to_csv();
//	ccbsn_dump.dump_to_csv();

	return 0;
}
