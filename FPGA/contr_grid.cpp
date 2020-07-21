#include <hls_stream.h>
#include <hls_math.h>
#include <ap_int.h>
#include "my_defines.h"
#include "my_types.h"
#include "contr_grid.h"


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


static void read_input(hls::stream<PMT_data_t> & PMT_data_stream, PMT_data_t *PMT_data, unsigned int num_events);
static void write_result(estim_event_t *estim_event, hls::stream<estim_event_t> & estim_event_stream, unsigned int num_events);
static void compute(hls::stream<estim_event_t> & estim_event_stream, hls::stream<PMT_data_t> & PMT_data_stream, unsigned int num_events, float mdrf_spline_coefs_dev[(MY + KY) * (MX + KX) * NUM_PMTS], float thresh_spline_coefs_dev[(MY + KY) * (MX + KX)], float gain_values_dev[NUM_PMTS]);
inline double make_double(uint32_t msw, uint32_t lsw = 0x00000000);
template <class _T> _T my_lgamma(_T x);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


extern "C" void contr_grid_kernel(estim_event_t *estim_event_dev, PMT_data_t *PMT_data_dev, unsigned int num_events, float mdrf_spline_coefs_dev[(MY + KY) * (MX + KX) * NUM_PMTS], float thresh_spline_coefs_dev[(MY + KY) * (MX + KX)], float gain_values_dev[NUM_PMTS]) {
  #pragma HLS INTERFACE m_axi port=estim_event_dev offset=slave bundle=gmem0
  #pragma HLS INTERFACE m_axi port=PMT_data_dev offset=slave bundle=gmem1
  #pragma HLS INTERFACE m_axi port=mdrf_spline_coefs_dev offset=slave bundle=gmem2
  #pragma HLS INTERFACE m_axi port=thresh_spline_coefs_dev offset=slave bundle=gmem3
  #pragma HLS INTERFACE m_axi port=gain_values_dev offset=slave bundle=gmem4
  #pragma HLS INTERFACE s_axilite port=estim_event_dev bundle=control
  #pragma HLS INTERFACE s_axilite port=PMT_data_dev bundle=control
  #pragma HLS INTERFACE s_axilite port=num_events bundle=control
  #pragma HLS INTERFACE s_axilite port=mdrf_spline_coefs_dev bundle=control
  #pragma HLS INTERFACE s_axilite port=thresh_spline_coefs_dev bundle=control
  #pragma HLS INTERFACE s_axilite port=gain_values_dev bundle=control
  #pragma HLS INTERFACE s_axilite port=return bundle=control
  #pragma HLS DATA_PACK variable=estim_event_dev field_level
  #pragma HLS DATA_PACK variable=PMT_data_dev field_level
  #pragma HLS STREAM variable=PMT_data_stream depth=32
  #pragma HLS STREAM variable=estim_event_stream depth=32
  #pragma HLS DATAFLOW
  static hls::stream<PMT_data_t> PMT_data_stream("PMT_data_stream");
  static hls::stream<estim_event_t> estim_event_stream("estim_event_stream");
  
  read_input(PMT_data_stream, PMT_data_dev, num_events);
  compute(estim_event_stream, PMT_data_stream, num_events, mdrf_spline_coefs_dev, thresh_spline_coefs_dev, gain_values_dev);
  write_result(estim_event_dev, estim_event_stream, num_events);
  return;
}


static void read_input(hls::stream<PMT_data_t> & PMT_data_stream, PMT_data_t *PMT_data, unsigned int num_events) {
  unsigned int event_index;
  
  for(event_index = 0; event_index < num_events; ++event_index) {
    #pragma HLS PIPELINE II=1
    PMT_data_stream << PMT_data[event_index];
  }
  return;
}


static void write_result(estim_event_t *estim_event, hls::stream<estim_event_t> & estim_event_stream, unsigned int num_events) {
  unsigned int event_index;
  
  for(event_index = 0; event_index < num_events; ++event_index) {
    #pragma HLS PIPELINE II=1
    estim_event_stream >> estim_event[event_index];
  }
  return;
}


static void compute(hls::stream<estim_event_t> & estim_event_stream, hls::stream<PMT_data_t> & PMT_data_stream, unsigned int num_events, float mdrf_spline_coefs_dev[(MY + KY) * (MX + KX) * NUM_PMTS], float thresh_spline_coefs_dev[(MY + KY) * (MX + KX)], float gain_values_dev[NUM_PMTS]) {
  spline_2D<mdrf_value_t, pos_value_t, MX, MY, KX, KY> calibr_funct_mdrf[NUM_PMTS];
  spline_2D<thresh_value_t, pos_value_t, MX, MY, KX, KY> calibr_funct_thresh;
  mdrf_value_t log_like_values[SIZE_CONTR_GRID][SIZE_CONTR_GRID];
  thresh_value_t thresh_spline_coefs[MY + KY][MX + KX];
  mdrf_value_t mdrf_spline_coefs[MY + KY][MX + KX];
  mdrf_value_t log_like, max_log_like;
  gain_value_t gain_values[NUM_PMTS];
  pos_value_t current_x, current_y;
  mdrf_value_t tmp_data[NUM_PMTS];
  int max_index_x, max_index_y;
  pos_value_t test_x, test_y;
  estim_event_t estim_event;
  mdrf_value_t camera_MDRF;
  unsigned int event_index;
  bool inside_x, inside_y;
  int index_x, index_y;
  PMT_data_t PMT_data;
  pos_value_t step;
  int i_x, i_y, n;
  int pmt, iter;
  float coeff;
  
  for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
    for(i_x = 0; i_x < (MX + KX); ++i_x) {
      for(i_y = 0; i_y < (MY + KY); ++i_y) {
        n = MAP_3D(MY + KY, MX + KX, NUM_PMTS, i_y, i_x, pmt);
        coeff = mdrf_spline_coefs_dev[n];
        mdrf_spline_coefs[i_y][i_x] = coeff;
      }
    }
    calibr_funct_mdrf[pmt] = spline_2D<mdrf_value_t, pos_value_t, MX, MY, KX, KY>(mdrf_spline_coefs);    
  }
  for(i_x = 0; i_x < (MX + KX); ++i_x) {
    for(i_y = 0; i_y < (MY + KY); ++i_y) {
      n = MAP_2D(MY + KY, MX + KX, i_y, i_x);
      coeff = thresh_spline_coefs_dev[n];
      thresh_spline_coefs[i_y][i_x] = coeff;
    }
  }
  calibr_funct_thresh = spline_2D<thresh_value_t, pos_value_t, MX, MY, KX, KY>(thresh_spline_coefs);
  for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
    gain_values[pmt] = gain_values_dev[pmt];
  }
  for(event_index = 0; event_index < num_events; ++event_index) {
    PMT_data_stream >> PMT_data;
    for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
      tmp_data[pmt] = mdrf_value_t(PMT_data.val[pmt]) / mdrf_value_t(gain_values_dev[pmt]);
    }
    current_x = current_y = pos_value_t(1) / pos_value_t(2);
    step = (pos_value_t(1) - pos_value_t(0)) / pos_value_t(SIZE_CONTR_GRID);
    for(iter = 0; iter < NUM_CONTR_GRID_ITER; ++iter) {
      for(index_x = 0; index_x < SIZE_CONTR_GRID; ++index_x) {
        test_x = current_x + (pos_value_t(index_x) - pos_value_t(pos_value_t(SIZE_CONTR_GRID - 1) / pos_value_t(2))) * step;
        for(index_y = 0; index_y < SIZE_CONTR_GRID; ++index_y) {
          test_y = current_y + (pos_value_t(index_y) - pos_value_t(pos_value_t(SIZE_CONTR_GRID - 1) / pos_value_t(2))) * step;
          inside_x = (pos_value_t(0) < test_x) && (test_x < pos_value_t(1));
          inside_y = (pos_value_t(0) < test_y) && (test_y < pos_value_t(1));
          if(inside_x && inside_y) {
            log_like = mdrf_value_t(0);
            for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
              camera_MDRF = calibr_funct_mdrf[pmt](test_x, test_y);
              if((tmp_data[pmt] != mdrf_value_t(0)) || (camera_MDRF != mdrf_value_t(0))) {
                log_like += tmp_data[pmt] * hls::log(camera_MDRF) - camera_MDRF;
              }
            }
            log_like_values[index_x][index_y] = log_like;
          } else {
            log_like_values[index_x][index_y] = mdrf_value_t(-1000);
          }
        }
      }
      max_log_like = log_like_values[0][0];
      max_index_x = max_index_y = 0;
      for(index_x = 0; index_x < SIZE_CONTR_GRID; ++index_x) {
        for(index_y = 0; index_y < SIZE_CONTR_GRID; ++index_y) {
          if(max_log_like < log_like_values[index_x][index_y]) {
            max_log_like = log_like_values[index_x][index_y];
            max_index_x = index_x;
            max_index_y = index_y;
          }
        }
      }
      current_x = current_x + (pos_value_t(max_index_x) - (pos_value_t(SIZE_CONTR_GRID - 1) / pos_value_t(2))) * step;
      current_y = current_y + (pos_value_t(max_index_y) - (pos_value_t(SIZE_CONTR_GRID - 1) / pos_value_t(2))) * step;
      step /= pos_value_t(CONTR_FACTOR);
    }
    inside_x = (pos_value_t(0) < current_x) && (current_x < pos_value_t(1));
    inside_y = (pos_value_t(0) < current_y) && (current_y < pos_value_t(1));
    if(inside_x && inside_y) {
      log_like = max_log_like;
      for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
        if(tmp_data[pmt] > mdrf_value_t(0)) {
          log_like -= my_lgamma<mdrf_value_t>(tmp_data[pmt] + mdrf_value_t(1));
        }
      }
      estim_event.valid = log_like > mdrf_value_t(calibr_funct_thresh(current_x, current_y));
      estim_event.log_like = float(log_like);
    } else {
      estim_event.valid = 0;
    }
    estim_event.x_pos = CAMERA_MIN_POS + float(current_x) * (CAMERA_MAX_POS - CAMERA_MIN_POS);
    estim_event.y_pos = CAMERA_MIN_POS + float(current_y) * (CAMERA_MAX_POS - CAMERA_MIN_POS);
    estim_event_stream << estim_event;
  }
  return;
}


double make_double(uint32_t msw, uint32_t lsw) {
  double tmp;
  
  *(((unsigned char *) (& tmp)) + 0) = *(((unsigned char *) (& lsw)) + 0);
  *(((unsigned char *) (& tmp)) + 1) = *(((unsigned char *) (& lsw)) + 1);
  *(((unsigned char *) (& tmp)) + 2) = *(((unsigned char *) (& lsw)) + 2);
  *(((unsigned char *) (& tmp)) + 3) = *(((unsigned char *) (& lsw)) + 3);
  *(((unsigned char *) (& tmp)) + 4) = *(((unsigned char *) (& msw)) + 0);
  *(((unsigned char *) (& tmp)) + 5) = *(((unsigned char *) (& msw)) + 1);
  *(((unsigned char *) (& tmp)) + 6) = *(((unsigned char *) (& msw)) + 2);
  *(((unsigned char *) (& tmp)) + 7) = *(((unsigned char *) (& msw)) + 3);
  return(tmp);
}


// Adapted from https://code.woboq.org/userspace/glibc/sysdeps/ieee754/dbl-64/e_lgamma_r.c.html
template <class _T> _T my_lgamma(_T x) {
  static const _T tc = _T(1.46163214496836224576e+00);
  static const _T tf = _T(-1.21486290535849611461e-01);
  static const _T tt = _T(-3.63867699703950536541e-18);
  static const _T a0 = _T(7.72156649015328655494e-02);
  static const _T a1 = _T(3.22467033424113591611e-01);
  static const _T a2 = _T(6.73523010531292681824e-02);
  static const _T a3 = _T(2.05808084325167332806e-02);
  static const _T a4 = _T(7.38555086081402883957e-03);
  static const _T a5 = _T(2.89051383673415629091e-03);
  static const _T a6 = _T(1.19270763183362067845e-03);
  static const _T a7 = _T(5.10069792153511336608e-04);
  static const _T a8 = _T(2.20862790713908385557e-04);
  static const _T a9 = _T(1.08011567247583939954e-04);
  static const _T a10 = _T(2.52144565451257326939e-05);
  static const _T a11 = _T(4.48640949618915160150e-05);
  static const _T t0 = _T(4.83836122723810047042e-01);
  static const _T t1 = _T(-1.47587722994593911752e-01);
  static const _T t2 = _T(6.46249402391333854778e-02);
  static const _T t3 = _T(-3.27885410759859649565e-02);
  static const _T t4 = _T(1.79706750811820387126e-02);
  static const _T t5 = _T(-1.03142241298341437450e-02);
  static const _T t6 = _T(6.10053870246291332635e-03);
  static const _T t7 = _T(-3.68452016781138256760e-03);
  static const _T t8 = _T(2.25964780900612472250e-03);
  static const _T t9 = _T(-1.40346469989232843813e-03);
  static const _T t10 = _T(8.81081882437654011382e-04);
  static const _T t11 = _T(-5.38595305356740546715e-04);
  static const _T t12 = _T(3.15632070903625950361e-04);
  static const _T t13 = _T(-3.12754168375120860518e-04);
  static const _T t14 = _T(3.35529192635519073543e-04);
  static const _T u0 = _T(-7.72156649015328655494e-02);
  static const _T u1 = _T(6.32827064025093366517e-01);
  static const _T u2 = _T(1.45492250137234768737e+00);
  static const _T u3 = _T(9.77717527963372745603e-01);
  static const _T u4 = _T(2.28963728064692451092e-01);
  static const _T u5 = _T(1.33810918536787660377e-02);
  static const _T v1 = _T(2.45597793713041134822e+00);
  static const _T v2 = _T(2.12848976379893395361e+00);
  static const _T v3 = _T(7.69285150456672783825e-01);
  static const _T v4 = _T(1.04222645593369134254e-01);
  static const _T v5 = _T(3.21709242282423911810e-03);
  static const _T s0 = _T(-7.72156649015328655494e-02);
  static const _T s1 = _T(2.14982415960608852501e-01);
  static const _T s2 = _T(3.25778796408930981787e-01);
  static const _T s3 = _T(1.46350472652464452805e-01);
  static const _T s4 = _T(2.66422703033638609560e-02);
  static const _T s5 = _T(1.84028451407337715652e-03);
  static const _T s6 = _T(3.19475326584100867617e-05);
  static const _T r1 = _T(1.39200533467621045958e+00);
  static const _T r2 = _T(7.21935547567138069525e-01);
  static const _T r3 = _T(1.71933865632803078993e-01);
  static const _T r4 = _T(1.86459191715652901344e-02);
  static const _T r5 = _T(7.77942496381893596434e-04);
  static const _T r6 = _T(7.32668430744625636189e-06);
  static const _T w0 = _T(4.18938533204672725052e-01);
  static const _T w1 = _T(8.33333333333329678849e-02);
  static const _T w2 = _T(-2.77777777728775536470e-03);
  static const _T w3 = _T(7.93650558643019558500e-04);
  static const _T w4 = _T(-5.95187557450339963135e-04);
  static const _T w5 = _T(8.36339918996282139126e-04);
  static const _T w6 = _T(-1.63092934096575273989e-03);
  _T r, y, z, p, p1, p2, p3, q, w;
  int i;
  
  if((x == _T(1)) || (x == _T(2))) {
    // Take care of 1 and 2
    r = _T(0);
  } else {
    if(x < _T(2)) {
      if(x <= _T(0.9)) {
        r = -hls::log(x);
        if(x >= _T(make_double(0x3FE76944))) {
          y = _T(1) - x;
          i = 0;
        } else {
          if(x >= _T(make_double(0x3FCDA661))) {
            y = x - (tc - _T(1));
            i = 1;
          } else {
            y = x;
            i = 2;
          }
        }
      } else {
        r = _T(0);
        if(x >= _T(make_double(0x3FFBB4C3))) {
          y = _T(2) - x;
          i = 0;
        } else {
          if(x >= _T(make_double(0x3FF3B4C4))) {
            y = x - tc;
            i = 1;
          } else {
            y = x - _T(1);
            i = 2;
          }
        }
      }
      switch(i) {
        case 0: {
          z = y * y;
          p1 = a0 + z * (a2 + z * (a4 + z * (a6 + z * (a8 + z * a10))));
          p2 = z * (a1 + z * (a3 + z * (a5 + z * (a7 + z * (a9 + z * a11)))));
          p = y * p1 + p2;
          r += (p - (_T(1) / _T(2)) * y);
          break;
        }
        case 1: {
          z = y * y;
          w = z * y;
          p1 = t0 + w * (t3 + w * (t6 + w * (t9 + w * t12)));
          p2 = t1 + w * (t4 + w * (t7 + w * (t10 + w * t13)));
          p3 = t2 + w * (t5 + w * (t8 + w * (t11 + w * t14)));
          p = z * p1 - (tt - w * (p2 + y * p3));
          r += (tf + p); break;
        }
        case 2: {
          p1 = y * (u0 + y * (u1 + y * (u2 + y * (u3 + y * (u4 + y * u5)))));
          p2 = _T(1) + y * (v1 + y * (v2 + y * (v3 + y * (v4 + y * v5))));
          r += (-(_T(1) / _T(2)) * y + p1 / p2);
        }
      }
    } else {
      // Case x >= 2
      if(x < _T(8)) {
        // Case 2 <= x < 8
        i = int(x);
        y = x - _T(i);
        p = y * (s0 + y * (s1 + y * (s2 + y * (s3 + y * (s4 + y *(s5 + y * s6))))));
        q = _T(1) + y * (r1 + y * (r2 + y * (r3 + y * (r4 + y * (r5 + y * r6)))));
        r = (_T(1) / _T(2)) * y + p / q;
        z = _T(1);
        switch(i) {
          case 7: z *= (y + _T(6));
          case 6: z *= (y + _T(5));
          case 5: z *= (y + _T(4));
          case 4: z *= (y + _T(3));
          case 3: z *= (y + _T(2));
          r += hls::log(z);
        }
      } else {
        // Case x >= 8
        z = _T(1) / x;
        y = z * z;
        w = w0 + z * (w1 + y * (w2 + y * (w3 + y * (w4 + y * (w5 + y * w6)))));
        r = w + (x - _T(1) / _T(2)) * (hls::log(x) - _T(1));
      }
    }
  }
  return(r);
}
