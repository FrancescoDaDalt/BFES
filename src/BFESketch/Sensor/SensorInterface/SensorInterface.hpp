//
//  SensorInterface.hpp
//  BaInDaSt
//
//  Created by Francesco Da Dalt on 25.10.23.
//

#ifndef SensorInterface_hpp
#define SensorInterface_hpp

#include <stdio.h>
template <typename key_type, typename support_type>
class SensorInterface {
public:
	virtual void insert(const key_type key, const support_type value, const bool nohash) = 0;
	virtual const support_type* const counter_ptr() const = 0;
};

#endif /* SensorInterface_hpp */
