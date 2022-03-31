#ifndef __F_Function_h__
#define __F_Function_h__

#include "funzionebase.h"

class F_Function : public FunzioneBase {

	public:

		F_Function();
		~F_Function();

		double eval(double x);
};

#endif
