#ifndef __Distribution_h__
#define __Distribution_h__

#include <vector>

class Distribution{

	public:
		virtual double eval(std::vector<double>) = 0;
};

#endif
