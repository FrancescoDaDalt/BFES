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

#ifndef GibbsSamplerEngine_hpp
#define GibbsSamplerEngine_hpp

//Takes a latent base and defines a gibbs sampler that solves A | HA=c
template <typename latentbase_type, typename supp_type>
class GibbsSamplerEngine {
public:
	using supp_t = supp_type;
	using lat_t = typename latentbase_type::template LatentType<supp_t, supp_t, supp_t>::t;
	using delta_t = lat_t;
	
//	static constexpr auto target_dim = latentbase_type::target_dim;
//	static constexpr auto latent_dim = latentbase_type::latent_dim;
	const size_t target_dim;
	const size_t latent_dim;
	
	const latentbase_type* const latentbase;
	
	GibbsSamplerEngine (const latentbase_type* const latentbase):
	latentbase(latentbase),
	target_dim(latentbase->target_dim),
	latent_dim(latentbase->latent_dim) {
		
	};
	
	GibbsSamplerEngine (): latentbase(nullptr) {};
	
//	Takes a current state, a prior and performs one round of gibbs sampling to generate a new sample
	template <typename priorsampler_type, typename targetstate_type, typename rng_type>
	auto sampleNextState_external (const priorsampler_type* const priorsampler,
								   const targetstate_type* const targetstate_source,
								   targetstate_type* const targetstate_target,
								   rng_type* const rng) const {
//		Copy source to target
		*targetstate_target = *targetstate_source;
//		Iterate over latent base
		for (int i = 0; i < latent_dim; i++) {
//			For every latent base vector, perform one sampling step
//			Get the i-th column of N
			const auto* const plane = latentbase->get_col(i);
//			Compute delta based on state, prior and column
			const delta_t d = priorsampler->template sampleHyperplane<delta_t>(plane,
																			   targetstate_target,
																			   rng);
//			Update state by delta times column
			plane->updateState(targetstate_target,
							   d);
		}
	}
};

//Same as GibbsSamplerEngine except that we randomize the choice of base N in order to improve the mixing of the markov chain
template <typename latentbase_type, typename supp_type>
class GibbsSamplerEngineSingleItem {
public:
	using supp_t = supp_type;
	using lat_t = typename latentbase_type::template LatentType<supp_t, supp_t, supp_t>::t;
	using delta_t = lat_t;
	
	static constexpr auto target_dim = latentbase_type::target_dim;
	static constexpr auto latent_dim = latentbase_type::latent_dim;
	
	const latentbase_type* const latentbase;
	
	GibbsSamplerEngineSingleItem (const GibbsSamplerEngineSingleItem& other): latentbase(other.latentbase) {};
	
	GibbsSamplerEngineSingleItem (): latentbase(nullptr) {};
	
	GibbsSamplerEngineSingleItem (const latentbase_type* const latentbase): latentbase(latentbase) {};
	
	template <typename priorsampler_type, typename targetstate_type, typename rng_type>
	auto sampleNextState_external (const priorsampler_type* const priorsampler,
								   const targetstate_type* const targetstate_source,
								   targetstate_type* const targetstate_target,
								   const supp_t counter,
								   rng_type* const rng) const {
		*targetstate_target = *targetstate_source;
		for (int i = 0; i < latent_dim; i++) {
//			Here we pick a random basis vector to improve mixing
			const auto* const plane = latentbase->random_col(rng);
			const lat_t d = priorsampler->template sampleHyperplaneBalanced<delta_t>(plane,
																					 targetstate_target,
																					 counter,
																					 rng);
			plane->updateState(targetstate_target,
							   d);
		}
	}
};

#endif /* GibbsSamplerEngine_hpp */
