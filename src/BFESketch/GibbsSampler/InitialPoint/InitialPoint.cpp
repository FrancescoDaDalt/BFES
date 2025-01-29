//
//  InitialPoint.cpp
//  BaInDaSt
//
//  Created by Francesco Da Dalt on 10.11.23.
//

#include "InitialPoint.hpp"

std::shared_ptr<mty::ndarray<int,1>> nint(const std::vector<int> &X) { return mty::new_array_ptr<int>(X);}
std::shared_ptr<mty::ndarray<double,1>> ndou(const std::vector<double> &X) { return mty::new_array_ptr<double>(X);}
