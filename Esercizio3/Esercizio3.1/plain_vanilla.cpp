#include <cmath>
#include "plain_vanilla.h"

using namespace std;

PlainVanilla::PlainVanilla(double strike, double r, WienerProcess* wiener_process, string option_type){
  _strike=strike;
  _option_type=option_type;
  _wiener_process = wiener_process;
  _r=r;
}

PlainVanilla::~PlainVanilla(){}

double PlainVanilla::pay_off(){
  double asset_price = _wiener_process->asset_price();
  double discount_factor = exp(-1*_r*_wiener_process->get_delivery_time());
  if(_option_type == "Call"){
    return discount_factor*fmax(asset_price -_strike, 0.);
  }else{
    return discount_factor*fmax(_strike - asset_price, 0.);
  }
}
