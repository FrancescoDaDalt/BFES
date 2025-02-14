/* 
 * MIT License
 * 
 * Copyright (c) 2023 Francesco Da Dalt
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

#ifndef InitialPoint_hpp
#define InitialPoint_hpp

#include <Eigen/Dense>
#ifndef NO_MOSEK
#include "fusion.h"

namespace msk = mosek::fusion;
namespace mty = monty;

std::shared_ptr<mty::ndarray<int,1>> nint(const std::vector<int> &X);
std::shared_ptr<mty::ndarray<double,1>> ndou(const std::vector<double> &X);
#endif

// Used to compute initial points for level-2 gibbs sampling. This solver is implemented manually for high performance since the problem to be solved is very simple
template<size_t targetSpace_size, class prior_type, typename counter_type, typename rng_type, typename anchor_type>
void simple_initial_point (const counter_type counter,
								  const typename prior_type::supp_t slack,
								  const prior_type* const prior,
								  rng_type* const rng_local,
								  anchor_type* const anchor) {
	using supp_t = typename prior_type::supp_t;
	int num_counters = 1;
	anchor->resize(targetSpace_size);

	Eigen::Array<double, targetSpace_size, 1> target;
	Eigen::Array<double, targetSpace_size, 1> target_std;

	prior->normPlusBatchMean(counter,
							 &target);
	prior->normPlusBatchVar(counter,
							&target_std);

	target_std = target_std.sqrt();

	Eigen::Array<double, targetSpace_size, 1> weightsAdd = target_std;
	weightsAdd(0) = 0;
	weightsAdd(1) = 0;
	const bool invalidWeightsAdd = std::abs(weightsAdd.sum()) <= 1e-6;
	weightsAdd(0) = (double) invalidWeightsAdd;
	weightsAdd = weightsAdd / weightsAdd.sum();

	Eigen::Array<double, targetSpace_size, 1> weightsSub = target_std;
	weightsSub(0) = 0;
	weightsSub.segment(2, targetSpace_size - 2) = 0;
	const bool invalidWeightsSub = std::abs(weightsSub.sum()) <= 1e-6;
	weightsSub(0) = (double) invalidWeightsSub;
	weightsSub = weightsSub/ weightsSub.sum();

	const double missing_to_constraint = counter - target.sum();
	const bool needToAdd = missing_to_constraint > 0;
	if (needToAdd) {
		Eigen::Array<double, targetSpace_size, 1> delta = missing_to_constraint * weightsAdd;
		*anchor = target + delta;
	} else {
		Eigen::Array<double, targetSpace_size, 1> delta = missing_to_constraint * weightsSub;
		*anchor = target + delta;
	}
}

#ifndef NO_MOSEK
// Used to compute level-1 initial points for gibbs sampling. Solving this problem is not trivial in general so we use MOSEK to solve it
template<class sensor_type, class prior_type, size_t num_anchors, size_t num_threads, typename counter_type, typename rng_type>
void multiple_random_initial_points(const sensor_type* const sensor,
									const size_t targetSpace_size,
									const size_t counters_size,
									const counter_type* const counters,
									Eigen::Array<typename prior_type::supp_t, Eigen::Dynamic, Eigen::Dynamic>* const ans,
									const typename prior_type::supp_t slack,
									const prior_type* const prior,
									rng_type* const rng_local) {
	using supp_t = typename prior_type::supp_t;
	int num_counters = counters_size;
	int num_items = targetSpace_size;
	
	std::vector<int> row_coord;
	std::vector<int> col_coord;
	std::vector<double> values;
	
	std::vector<std::shared_ptr<mty::ndarray<double,1>>> random_anchors_target_mosek;
	
	for (int i = 0; i < num_anchors; i++) {
		auto random_anchor_target_std = prior->randomVectorSample(targetSpace_size, rng_local);
		std::vector<double> random_anchor_target_std_correcttype(random_anchor_target_std->begin(), random_anchor_target_std->end());
		random_anchors_target_mosek.push_back(mty::new_array_ptr<double>(random_anchor_target_std_correcttype));
		delete random_anchor_target_std;
	}
	
	for (int i = 0; i < num_items; i++) {
		const std::vector<size_t> hashes = sensor->pseudohash(i);
		
		row_coord.insert(row_coord.end(), hashes.begin(), hashes.end());
		col_coord.insert(col_coord.end(), hashes.size(), i);
	}
	
	values.resize(row_coord.size());
	std::fill(values.begin(), values.end(), 1.0);
	
	msk::Model::t mdl = new msk::Model("initialpointComputation");
	
	auto _mdl = mty::finally([&](){mdl->dispose();});
	msk::Variable::t aStar_mosek;
	
	if (!std::is_integral<supp_t>::value) {
		aStar_mosek = mdl->variable("aStar", num_items, msk::Domain::unbounded());
	} else {
		aStar_mosek = mdl->variable("aStar", num_items, msk::Domain::integral(msk::Domain::unbounded()));
	}
	
	msk::Variable::t solution_norm = mdl->variable("solution_norm", 1, msk::Domain::unbounded());
	std::vector<double> counters_std(counters_size);
	for (int i = 0; i < counters_size; i++) {
		counters_std[i] = (double) counters->operator()(i);
	}
	auto counters_mosak = mty::new_array_ptr<double>(counters_std);
	
	auto prior_boundaries = prior->get_bounds();
	prior_boundaries.first -= slack;
	prior_boundaries.second += slack;
	if (std::is_integral<supp_t>::value) {
		prior_boundaries.first -= 0.5;
		prior_boundaries.second += 0.5;
	}
	
	auto sens_mat = msk::Matrix::sparse(num_counters, num_items, nint(row_coord), nint(col_coord), ndou(values));
	mdl->constraint(msk::Expr::mul(sens_mat, aStar_mosek), msk::Domain::equalsTo(counters_mosak));
	if (std::isfinite((double) prior_boundaries.first)) {
		mdl->constraint(aStar_mosek, msk::Domain::greaterThan((double) prior_boundaries.first));
	}
	if (std::isfinite((double) prior_boundaries.second)) {
		mdl->constraint(aStar_mosek, msk::Domain::lessThan((double) prior_boundaries.second));
	}
	
	mdl->objective( msk::ObjectiveSense::Minimize, solution_norm);
	mdl->setSolverParam("mioTolRelGap", 1e-6);
	mdl->setSolverParam("intpntCoTolRelGap", 1e-6);
	
	auto multi_models = std::make_shared<mty::ndarray<msk::Model::t, 1>>(mty::shape(num_anchors));
	for (int i = 0; i < num_anchors; i++) {
		(*multi_models)[i] = mdl->clone();
		(*multi_models)[i]->constraint( msk::Expr::vstack(solution_norm, msk::Expr::sub(aStar_mosek, random_anchors_target_mosek[i])), msk::Domain::inQCone());
	}
	
	msk::Model::solveBatch(false, -1.0, num_threads, multi_models);
	Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> aStar_eigen(targetSpace_size, num_anchors);
	
	ans->resize(targetSpace_size, num_anchors);
	for (int i = 0; i < num_anchors; i++) {
		auto sol = (*multi_models)[i]->getVariable("aStar")->level();
		memcpy((void*) (aStar_eigen.data() + i * targetSpace_size), (void*) sol->begin(), num_items * sizeof(double));
		(*multi_models)[i]->dispose();
	}
	if (!std::is_integral<supp_t>::value) {
		*ans = aStar_eigen.template cast<supp_t>();
	} else {
		*ans = aStar_eigen.round().template cast<supp_t>();
	}
}
#endif

#endif /* InitialPoint_hpp */
