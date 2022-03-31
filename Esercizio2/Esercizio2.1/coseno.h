#ifndef __Coseno_h__
#define __Coseno_h__

#include "funzionebase.h"

class Coseno : public FunzioneBase {

	public:
		Coseno();
		~Coseno();

		double eval(double x);
};

#endif
