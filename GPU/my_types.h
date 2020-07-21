#ifndef _MY_TYPES_H
#define _MY_TYPES_H

#include "spline.hpp"
#include "my_defines.h"


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<int _N> struct _P {
  enum {
    val = _P<_N / 2>::val * 2
  };
};


template<> struct _P<0> {
  enum {
    val = 1
  };
};


template<class _T> struct padding_size {
  enum {
    val = _P<sizeof(_T) - 1>::val - sizeof(_T)
  }; 
};


template<int _N> struct add_padding {
  char padding[_N]; 
};


template<> struct add_padding<0> {
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


struct PMT_data_no_padding_t {
  int16_t val[NUM_PMTS];
};


struct PMT_data_t: PMT_data_no_padding_t, add_padding<padding_size<PMT_data_no_padding_t>::val> {
};


struct calibr_data_t {
  float mdrf[NUM_PMTS][NUM_SAMPL][NUM_SAMPL];
  float thresh[NUM_SAMPL][NUM_SAMPL];
  float gain[NUM_PMTS];
};


struct estim_event_nopadding_t {
  float x_pos;
  float y_pos;
  uint32_t valid;
  float log_like;
};


struct estim_event_t: estim_event_nopadding_t, add_padding<padding_size<estim_event_nopadding_t>::val> {
};


typedef spline_2D<float, float, MX, MY, KX, KY> mdrf_spline_t;
typedef spline_2D<float, float, MX, MY, KX, KY> thresh_spline_t;


struct calibr_funct_t {
  mdrf_spline_t mdrf[NUM_PMTS];
  thresh_spline_t thresh;
  float gain[NUM_PMTS];
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif // _MY_TYPES_H
