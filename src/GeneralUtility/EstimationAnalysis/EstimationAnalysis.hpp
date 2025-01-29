//
//  EstimationAnalysis.hpp
//  MCMC_Sketch
//
//  Created by Francesco Da Dalt on 11.11.22.
//

#ifndef EstimationAnalysis_hpp
#define EstimationAnalysis_hpp

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

#include <Eigen/Dense>

template<typename samples_type, typename index_type, typename gt_type>
class EstimationAnalysis {
	samples_type* const level1_samples;
	samples_type* const level2_samples;
	index_type* const indices;
	const gt_type* gt;
	int num_calibration_levels = 100;
	
	std::vector<double> calibration_levels;
	std::vector<double> calibration_values;
	
	std::vector<double> sharpness_levels;
	std::vector<std::vector<double>> sharpness_values;
	
	std::vector<double> interval_lengths;
	
	std::vector<std::tuple<std::string, double>> stats;
public:
	
	EstimationAnalysis(const std::string name, const std::tuple<index_type*, samples_type*, samples_type*> query_data, const gt_type* const gt):
	level2_samples(std::get<2>(query_data)),
	level1_samples(std::get<1>(query_data)),
	indices(std::get<0>(query_data)),
	gt(gt) {
		
		const int num_items = (int) gt->size();
		const int samples_dim = (int) level2_samples->size();
		std::cout << "\n-> " + name + " Estimation Analysis" << std::endl;
		
		for(int i = 0; i < samples_dim; i++) {
			level2_samples->at(i).erase(std::remove_if(level2_samples->at(i).begin(), level2_samples->at(i).end(), [](double x) {
				return !std::isfinite(x);  // Remove if not finite (NaN, inf, -inf)
			}), level2_samples->at(i).end());
			level1_samples->at(i).erase(std::remove_if(level1_samples->at(i).begin(), level1_samples->at(i).end(), [](double x) {
				return !std::isfinite(x);  // Remove if not finite (NaN, inf, -inf)
			}), level1_samples->at(i).end());
			for(auto x: level2_samples->at(i)) {
				assert(std::isfinite(x));
			}
			for(auto x: level1_samples->at(i)) {
				assert(std::isfinite(x));
			}
			sort(level2_samples->at(i).begin(), level2_samples->at(i).end());
			sort(level1_samples->at(i).begin(), level1_samples->at(i).end());
		}
		
		std::cout << "--> Error Statistics:\n";
		
		auto eigen_gt = Eigen::Map<const Eigen::Array<double, Eigen::Dynamic, 1>>(gt->data(), num_items);
		Eigen::Array<double, Eigen::Dynamic, 1> eigen_gt_sorted = eigen_gt;
		std::sort(eigen_gt_sorted.data(), eigen_gt_sorted.data() + num_items);
		
		Eigen::Array<double, Eigen::Dynamic, 1> corr_maximizer_vec(num_items);
		Eigen::Array<double, Eigen::Dynamic, 1> mse_minimizer_vec(num_items);
		Eigen::Array<double, Eigen::Dynamic, 1> mse_estimation_vec(num_items);
		Eigen::Array<double, Eigen::Dynamic, 1> mae_minimizer_vec(num_items);
		Eigen::Array<double, Eigen::Dynamic, 1> mae_estimation_vec(num_items);
		Eigen::Array<double, Eigen::Dynamic, 1> mundere_minimizer_vec(num_items);
		Eigen::Array<double, Eigen::Dynamic, 1> mundere_estimation_vec(num_items);
		Eigen::Array<double, Eigen::Dynamic, 1> movere_minimizer_vec(num_items);
		Eigen::Array<double, Eigen::Dynamic, 1> movere_estimation_vec(num_items);
		Eigen::Array<double, Eigen::Dynamic, 1> are_minimizer_vec(num_items);
		Eigen::Array<double, Eigen::Dynamic, 1> are_estimation_vec(num_items);
		
		for (int i = 0; i < num_items; i++) {
			double idx = indices->at(i);
			int numsamples_level2 = level2_samples->at(idx).size();
			int numsamples_level1 = level1_samples->at(idx).size();
			auto eigen_estimates_level2 = Eigen::Map<Eigen::Array<double, Eigen::Dynamic, 1>>(level2_samples->at(idx).data(), numsamples_level2);
			auto eigen_estimates_level1 = Eigen::Map<Eigen::Array<double, Eigen::Dynamic, 1>>(level1_samples->at(idx).data(), numsamples_level1);
			
//			Just for visualization
//			for (auto v: level2_samples->at(idx)) {
//				std::cout << v << ", ";
//			}
//			std::cout << std::endl;
			
			
			double corr_maximizer = eigen_estimates_level1.mean();
			
			double mse_minimizer = eigen_estimates_level2.mean();
			double mse_est = (eigen_estimates_level2 - mse_minimizer).square().mean();
			
			double mae_minimizer = eigen_estimates_level2(numsamples_level2 / 2);
			double mae_est = (eigen_estimates_level2 - mae_minimizer).abs().mean();
			
			double mundere_minimizer = eigen_estimates_level2(std::min((int) std::round(numsamples_level2 * 0.909091), numsamples_level2 - 1));
			double mundere_est = ((mundere_minimizer - eigen_estimates_level2) > 0).select((mundere_minimizer - eigen_estimates_level2) * 1, (mundere_minimizer - eigen_estimates_level2) * -10).mean();
			
			double movere_minimizer = eigen_estimates_level2(std::min((int) std::round(numsamples_level2 * 0.0909091), numsamples_level2 - 1));
			double movere_est = ((movere_minimizer - eigen_estimates_level2) > 0).select((mundere_minimizer - eigen_estimates_level2) * 10, (mundere_minimizer - eigen_estimates_level2) * -1).mean();
			
			Eigen::Array<double, Eigen::Dynamic, 1> aux = 1.0 / eigen_estimates_level2;
			std::partial_sum(aux.data(), aux.data() + numsamples_level2, aux.data());
			double thresh = aux(numsamples_level2 - 1) / 2;
			int are_idx = std::lower_bound(aux.data(), aux.data(), thresh) - aux.data();
			double are_minimizer = eigen_estimates_level2(are_idx);
			double are_est = ((eigen_estimates_level2 - are_minimizer) / eigen_estimates_level2).abs().mean();
			
			corr_maximizer_vec(i) = corr_maximizer;
			mse_minimizer_vec(i) = mse_minimizer;
			mse_estimation_vec(i) = mse_est;
			mae_minimizer_vec(i) = mae_minimizer;
			mae_estimation_vec(i) = mae_est;
			mundere_minimizer_vec(i) = mundere_minimizer;
			mundere_estimation_vec(i) = mundere_est;
			movere_minimizer_vec(i) = movere_minimizer;
			movere_estimation_vec(i) = movere_est;
			are_minimizer_vec(i) = are_minimizer;
			are_estimation_vec(i) = are_est;
			
//			Just for visualization
//			std::cout << "mse_minimizer: " << mse_minimizer << std::endl;
//			std::cout << "mse_est: " << mse_est << std::endl;
//			std::cout << "mae_minimizer: " << mae_minimizer << std::endl;
//			std::cout << "mae_est: " << mae_est << std::endl;
//			std::cout << "mundere_minimizer: " << mundere_minimizer << std::endl;
//			std::cout << "mundere_est: " << mundere_est << std::endl;
//			std::cout << "movere_minimizer: " << movere_minimizer << std::endl;
//			std::cout << "movere_est: " << movere_est << std::endl;
//			std::cout << "are_minimizer: " << are_minimizer << std::endl;
//			std::cout << "are_est: " << are_est << std::endl;
//			assert(false);
		}
		
		corr_maximizer_vec = (corr_maximizer_vec - corr_maximizer_vec.mean());
		corr_maximizer_vec /= std::sqrt(corr_maximizer_vec.square().sum());
		
		Eigen::ArrayXd square_errors = (mse_minimizer_vec - eigen_gt).square();
		Eigen::ArrayXd absolute_errors = (mae_minimizer_vec - eigen_gt).abs();
		Eigen::ArrayXd absolute_relative_errors = ((are_minimizer_vec - eigen_gt) / eigen_gt).abs();
		Eigen::ArrayXd underest_penal_errors = ((mundere_minimizer_vec - eigen_gt) > 0).select((mundere_minimizer_vec - eigen_gt) * 1, (mundere_minimizer_vec - eigen_gt) * -10);
		Eigen::ArrayXd overest_penal_errors = ((movere_minimizer_vec - eigen_gt) > 0).select((movere_minimizer_vec - eigen_gt) * 10, (movere_minimizer_vec - eigen_gt) * -1);
		
		double mean_se = (mse_minimizer_vec - eigen_gt).square().mean();
		double neg_r2 = (mse_minimizer_vec - eigen_gt).square().sum() / (eigen_gt.mean() - eigen_gt).square().sum();
		double mean_ae = (mae_minimizer_vec - eigen_gt).abs().mean();
		double neg_a2 = (mae_minimizer_vec - eigen_gt).abs().sum() / (eigen_gt_sorted(num_items / 2) - eigen_gt).abs().sum();
		double mean_undere = ((mundere_minimizer_vec - eigen_gt) > 0).select((mundere_minimizer_vec - eigen_gt) * 1, (mundere_minimizer_vec - eigen_gt) * -10).mean();
		double mean_overe = ((movere_minimizer_vec - eigen_gt) > 0).select((movere_minimizer_vec - eigen_gt) * 10, (movere_minimizer_vec - eigen_gt) * -1).mean();
		double mean_are = ((are_minimizer_vec - eigen_gt) / eigen_gt).abs().mean();
		
		double mean_se_pred = mse_estimation_vec.mean();
		double mean_ae_pred = mae_estimation_vec.mean();
		double mean_undere_pred = mundere_estimation_vec.mean();
		double mean_overe_pred = movere_estimation_vec.mean();
		double mean_are_pred = are_estimation_vec.mean();
		
		std::nth_element(square_errors.data(), square_errors.data() + num_items / 2, square_errors.data() + num_items);
		std::nth_element(absolute_errors.data(), absolute_errors.data() + num_items / 2, absolute_errors.data() + num_items);
		std::nth_element(absolute_relative_errors.data(), absolute_relative_errors.data() + num_items / 2, absolute_relative_errors.data() + num_items);
		std::nth_element(underest_penal_errors.data(), underest_penal_errors.data() + num_items / 2, underest_penal_errors.data() + num_items);
		std::nth_element(overest_penal_errors.data(), overest_penal_errors.data() + num_items / 2, overest_penal_errors.data() + num_items);
		
		double median_se = *(square_errors.data() + num_items / 2);
		double median_ae = *(absolute_errors.data() + num_items / 2);
		double median_are = *(absolute_relative_errors.data() + num_items / 2);
		double median_undere = *(underest_penal_errors.data() + num_items / 2);
		double median_overe = *(overest_penal_errors.data() + num_items / 2);
		
		std::cout << "---> Underestimation-Penalized:\n";
		std::cout << "----> Mean Underest-Pen Error: " << mean_undere << std::endl;
		stats.push_back(std::tuple<std::string, double>("munder_err", mean_undere));
		std::cout << "----> Predicted Mean Underest-Pen Error: " << mean_undere_pred << std::endl;
		stats.push_back(std::tuple<std::string, double>("munder_pred_err", mean_undere_pred));
		std::cout << "----> Median Underest-Pen Error: " << median_undere << std::endl;
		stats.push_back(std::tuple<std::string, double>("med_under_err", median_undere));
		std::cout << "---> Overestimation-Penalized:\n";
		std::cout << "----> Mean Overest-Pen Error: " << mean_overe << std::endl;
		stats.push_back(std::tuple<std::string, double>("mover_err", mean_overe));
		std::cout << "----> Predicted Mean Overest-Pen Error: " << mean_overe_pred << std::endl;
		stats.push_back(std::tuple<std::string, double>("mover_pred_err", mean_overe_pred));
		std::cout << "----> Median Overest-Pen Error: " << median_overe << std::endl;
		stats.push_back(std::tuple<std::string, double>("med_over_err", median_overe));
		std::cout << "---> Mean Square:\n";
		std::cout << "----> Mean Square Error: " << mean_se << std::endl;
		stats.push_back(std::tuple<std::string, double>("ms_err", mean_se));
		std::cout << "----> Predicted Mean Square Error: " << mean_se_pred << std::endl;
		stats.push_back(std::tuple<std::string, double>("ms_pred_err", mean_se_pred));
		std::cout << "----> Negative R2: " << neg_r2 << std::endl;
		stats.push_back(std::tuple<std::string, double>("neg_r2", neg_r2));
		std::cout << "----> Median Square Error: " << median_se << std::endl;
		stats.push_back(std::tuple<std::string, double>("med_s_err", median_se));
		std::cout << "---> Absolute Relative:\n";
		std::cout << "----> Mean AbsRel Error: " << mean_are << std::endl;
		stats.push_back(std::tuple<std::string, double>("mare_err", mean_are));
		std::cout << "----> Predicted Mean AbsRel Error: " << mean_are_pred << std::endl;
		stats.push_back(std::tuple<std::string, double>("mare_pred_err", mean_are_pred));
		std::cout << "----> Median AbsRel Error: " << median_are << std::endl;
		stats.push_back(std::tuple<std::string, double>("med_are_err", median_are));
		std::cout << "---> Mean Absolute:\n";
		std::cout << "----> Mean Abs Error: " << mean_ae << std::endl;
		stats.push_back(std::tuple<std::string, double>("mae_err", mean_ae));
		std::cout << "----> Predicted Mean Abs Error: " << mean_ae_pred << std::endl;
		stats.push_back(std::tuple<std::string, double>("mae_pred_err", mean_ae_pred));
		std::cout << "----> Negative A2: " << neg_a2 << std::endl;
		stats.push_back(std::tuple<std::string, double>("neg_a2", neg_a2));
		std::cout << "----> Median Abs Error: " << median_ae << std::endl;
		stats.push_back(std::tuple<std::string, double>("med_ae_err", median_ae));
		
		Eigen::ArrayXd aux11 = mae_estimation_vec - mae_estimation_vec.mean();
		Eigen::ArrayXd aux12 = (mae_minimizer_vec - eigen_gt).abs();
		aux12 = aux12 - aux12.mean();
		double mae_real_pred_corr = (aux11 * aux12).sum() / sqrt((aux11.square().sum() * aux12.square().sum()));
		
		aux11 = mse_estimation_vec - mse_estimation_vec.mean();
		aux12 = (mse_minimizer_vec - eigen_gt).square();
		aux12 = aux12 - aux12.mean();
		double mse_real_pred_corr = (aux11 * aux12).sum() / sqrt((aux11.square().sum() * aux12.square().sum()));
		
		aux11 = mundere_estimation_vec - mundere_estimation_vec.mean();
		aux12 = ((mundere_minimizer_vec - eigen_gt) > 0).select((mundere_minimizer_vec - eigen_gt) * 1, (mundere_minimizer_vec - eigen_gt) * -10);
		aux12 = aux12 - aux12.mean();
		double mundere_real_pred_corr = (aux11 * aux12).sum() / sqrt((aux11.square().sum() * aux12.square().sum()));
		
		aux11 = movere_estimation_vec - movere_estimation_vec.mean();
		aux12 = ((movere_minimizer_vec - eigen_gt) > 0).select((movere_minimizer_vec - eigen_gt) * 10, (movere_minimizer_vec - eigen_gt) * -1);
		aux12 = aux12 - aux12.mean();
		double movere_real_pred_corr = (aux11 * aux12).sum() / sqrt((aux11.square().sum() * aux12.square().sum()));
		
		
		Eigen::ArrayXd aux1 = corr_maximizer_vec - corr_maximizer_vec.mean();
		Eigen::ArrayXd aux2 = eigen_gt - eigen_gt.mean();
		double linear_corr = (aux1 * aux2).sum() / sqrt((aux2.square().sum() * aux1.square().sum()));
		
		Eigen::ArrayXi aux4 = Eigen::ArrayXi::LinSpaced(num_items, 0, num_items - 1);
		Eigen::ArrayXi aux5 = Eigen::ArrayXi::LinSpaced(num_items, 0, num_items - 1);
		
		std::vector<std::pair<int, double>> aux4_aux;
		aux4_aux.reserve(num_items);
		std::vector<std::pair<int, double>> aux5_aux;
		aux5_aux.reserve(num_items);
		
		for (int i = 0; i < num_items; i++) {
			aux4_aux.push_back(std::make_pair(i, mse_minimizer_vec(i)));
			aux5_aux.push_back(std::make_pair(i, eigen_gt(i)));
		}
		
		std::sort(aux5_aux.begin(), aux5_aux.end(), [](std::pair<int, double> fst, std::pair<int, double> snd){return fst.second < snd.second;});
		std::sort(aux4_aux.begin(), aux4_aux.end(), [](std::pair<int, double> fst, std::pair<int, double> snd){return fst.second < snd.second;});
		
		for (int i = 0; i < gt->size(); i++) {
			aux4(aux4_aux[i].first) = i;
			aux5(aux5_aux[i].first) = i;
		}
		
		Eigen::ArrayXd pt1 = (aux4 - aux5).cast<double>().square() * 6.0;
		double pt2 = (double) num_items  * (((double) num_items) * num_items - 1.0);
		
		double rank_corr = 1.0 - (pt1 / pt2).sum();
		
		std::cout << "--> Correlation Statistics:\n";
		std::cout << "---> Mean Absolute Error Correlation: " << mae_real_pred_corr << std::endl;
		stats.push_back(std::tuple<std::string, double>("maecorr", mae_real_pred_corr));
		std::cout << "---> Mean Square Error Correlation: " << mse_real_pred_corr << std::endl;
		stats.push_back(std::tuple<std::string, double>("msecorr", mse_real_pred_corr));
		std::cout << "---> Mean Overest-Pen Error Correlation: " << movere_real_pred_corr << std::endl;
		stats.push_back(std::tuple<std::string, double>("movercorr", movere_real_pred_corr));
		std::cout << "---> Mean Underest-Pen Error Correlation: " << mundere_real_pred_corr << std::endl;
		stats.push_back(std::tuple<std::string, double>("mundercorr", mundere_real_pred_corr));
		std::cout << "---> Linear Correlation: " << linear_corr << std::endl;
		stats.push_back(std::tuple<std::string, double>("lincorr", linear_corr));
		std::cout << "---> Rank Correlation: " << rank_corr << std::endl;
		stats.push_back(std::tuple<std::string, double>("rankcorr", rank_corr));
		
		
		calibration_levels.reserve(num_calibration_levels);
		calibration_values.reserve(num_calibration_levels);
		for (int i = 0; i < num_calibration_levels; i++) {
			double calib_level = ((double) i) / (num_calibration_levels - 1);
			calibration_levels.push_back(calib_level);
			double calib_value = 0;
			for (int j = 0; j < num_items; j++) {
				double idx = indices->at(j);
				int numsamples_level2 = (int) level2_samples->at(idx).size();
				auto eigen_estimates_level2 = Eigen::Map<Eigen::Array<double, Eigen::Dynamic, 1>>(level2_samples->at(idx).data(), numsamples_level2);
				
				int rank_of_quantile = (int) ((numsamples_level2 - 1) * calib_level);
				double thresh = eigen_estimates_level2(rank_of_quantile);
				calib_value += thresh >= gt->at(j);
			}
			calib_value /= num_items;
			calibration_values.push_back(calib_value);
			std::string aux = "calib_" + std::to_string(calib_level);
			stats.push_back(std::make_tuple(aux, calib_value));
		}
		
		auto calib_diff = Eigen::Map<Eigen::ArrayXd>(calibration_levels.data(), num_calibration_levels) - Eigen::Map<Eigen::ArrayXd>(calibration_values.data(), num_calibration_levels);
		
		std::cout << "--> Calibration Statistics:\n";
		std::cout << "---> Mean Absolute Calibration Bias: " << abs(calib_diff.mean()) << std::endl;
		stats.push_back(std::make_tuple("abs_calib_bias", abs(calib_diff.mean())));
		std::cout << "---> Mean Calibration Bias: " << calib_diff.mean() << std::endl;
		stats.push_back(std::make_tuple("calib_bias", calib_diff.mean()));
		std::cout << "---> Mean Calibration Error: " << calib_diff.abs().mean() << std::endl;
		stats.push_back(std::make_tuple("avg_calib_err", calib_diff.abs().mean()));
	};
	
	auto getStats() const {
		return &stats;
	}
};

#endif /* EstimationAnalysis_hpp */
