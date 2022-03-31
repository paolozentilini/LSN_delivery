#ifndef __G_Function_h__
#define __G_Function_h__

#include "funzionebase.h"

class G_Function : public FunzioneBase {

	public:

		G_Function();
		~G_Function();

		double eval(double x);
};

#endif
