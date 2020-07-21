#ifndef _CONTR_GRID_H
#define _CONTR_GRID_H
#include <hls_stream.h>
#include <ap_fixed.h>
#include "spline.hpp"
#include "my_types.h"


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


typedef ap_fixed<32, 12> mdrf_value_t;
typedef ap_fixed<16, 6> thresh_value_t;
typedef ap_ufixed<32, 12> gain_value_t;
typedef ap_fixed<24, 6> pos_value_t;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


extern "C" void contr_grid_kernel(estim_event_t *estim_event_dev, PMT_data_t *PMT_data_dev, unsigned int num_events, float mdrf_spline_coefs_dev[(MY + KY) * (MX + KX) * NUM_PMTS], float thresh_spline_coefs_dev[(MY + KY) * (MX + KX)], float gain_values_dev[NUM_PMTS]);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#endif // _CONTR_GRID_H
