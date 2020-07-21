#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdint>
#include <vector>
#include <chrono>
#include <array>
#include <cmath>
#include "spline.hpp"
#include "my_defines.h"
#include "my_types.h"
#include "my_utils.h"
#include "contr_grid.h"


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void sample_calibr_funct(const struct calibr_funct_t & calibr_funct);
std::vector<estim_event_t, aligned_allocator<estim_event_t>> contr_grid(const std::vector<PMT_data_t, aligned_allocator<PMT_data_t>> & PMT_data, calibr_funct_t calibr_funct);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char **argv) {
  std::vector<estim_event_t, aligned_allocator<estim_event_t>> estim_event;
  std::vector<PMT_data_t, aligned_allocator<PMT_data_t>> PMT_data;
  calibr_funct_t calibr_funct;
  calibr_data_t calibr_data;
  
  calibr_data = get_calibration_data("camera0_79x79_1.5mm_tc99m_mean", "camera0_thresh.dat", "camera0_79x79_1.5mm_tc99m_gains");
  calibr_funct = get_calibration_funct(calibr_data);
  sample_calibr_funct(calibr_funct);
  PMT_data = get_PMT_data("ResPhantom022516-0mm_00.dat");
  PMT_data.resize(10);
  estim_event = contr_grid(PMT_data, calibr_funct);
  write_estim_events(estim_event, "estim_events_FPGA.dat");
  
  unsigned int i;
  for(i = 0; i < estim_event.size(); ++i) {
    std::cout << estim_event[i].x_pos << "  " << estim_event[i].y_pos << "  " << estim_event[i].valid << "  " << estim_event[i].log_like << std::endl;
  }
  return(0);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void sample_calibr_funct(const struct calibr_funct_t & calibr_funct) {
  const int num_sampl_x = 128;
  const int num_sampl_y = 128;
  std::array<std::array<float, num_sampl_y>, num_sampl_x> data;
  pos_value_t pos_x[num_sampl_x];
  pos_value_t pos_y[num_sampl_y];
  std::string filename;
  int nx, ny;
  int i, pmt;
  
  for(i = 0; i < num_sampl_x; ++i) {
    pos_x[i] = pos_value_t(float(i) / float(num_sampl_x - 1));
  }
  for(i = 0; i < num_sampl_y; ++i) {
    pos_y[i] = pos_value_t(float(i) / float(num_sampl_y - 1));
  }
  for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
    for(nx = 0; nx < num_sampl_x; ++nx) {
      for(ny = 0; ny < num_sampl_y; ++ny) {
        data[nx][ny] = calibr_funct.mdrf[pmt](pos_x[nx], pos_y[ny]);
      }
    }
    std::ostringstream ss;
    ss << std::setw(3) << std::setfill('0') << pmt;
    filename = "mdrf_samples_FPGA_" + ss.str() + ".dat";
    write_dat_2d<float, num_sampl_x, num_sampl_y>(data, filename.c_str());
  }
  for(nx = 0; nx < num_sampl_x; ++nx) {
    for(ny = 0; ny < num_sampl_y; ++ny) {
      data[nx][ny] = calibr_funct.thresh(pos_x[nx], pos_y[ny]);
    }
  }
  write_dat_2d<float, num_sampl_x, num_sampl_y>(data, "thresh_samples_FPGA.dat");
  return;
}


std::vector<estim_event_t, aligned_allocator<estim_event_t>> contr_grid(const std::vector<PMT_data_t, aligned_allocator<PMT_data_t>> & PMT_data, calibr_funct_t calibr_funct) {
  std::vector<estim_event_t, aligned_allocator<estim_event_t>> estim_event(PMT_data.size());
  float mdrf_spline_coefs[(MY + KY) * (MX + KX) * NUM_PMTS];
  float thresh_spline_coefs[(MY + KY) * (MX + KX)];
  float tmp_spline_coefs[MY + KY][MX + KX];
  estim_event_t *estim_event_dev;
  unsigned int num_events, index;
  float gain_values[NUM_PMTS];
  PMT_data_t *PMT_data_dev;
  int i_x, i_y, n;
  float coeff;
  int pmt;
  
  num_events = (unsigned int) PMT_data.size();
  std::cout << "Number of events: " << num_events << "." << std::endl;
  auto start = std::chrono::steady_clock::now();
  if(num_events > 0) {
  
    for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
      calibr_funct.mdrf[pmt].get_coefs(tmp_spline_coefs);
      for(i_x = 0; i_x < (MX + KX); ++i_x) {
        for(i_y = 0; i_y < (MY + KY); ++i_y) {
          coeff = tmp_spline_coefs[i_y][i_x];
          n = MAP_3D(MY + KY, MX + KX, NUM_PMTS, i_y, i_x, pmt);
          mdrf_spline_coefs[n] = coeff;
        }
      }
    }
    
    calibr_funct.thresh.get_coefs(tmp_spline_coefs);
    for(i_x = 0; i_x < (MX + KX); ++i_x) {
      for(i_y = 0; i_y < (MY + KY); ++i_y) {
        coeff = tmp_spline_coefs[i_y][i_x];
        n = MAP_2D(MY + KY, MX + KX, i_y, i_x);
        thresh_spline_coefs[n] = coeff;
      }
    }
    
    for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
      gain_values[pmt] = calibr_funct.gain[pmt];
    }
 
    PMT_data_dev = new PMT_data_t[num_events];
    estim_event_dev = new estim_event_t[num_events];
    for(index = 0; index < num_events; ++index) {
      PMT_data_dev[index] = PMT_data[index];
    }
    contr_grid_kernel(estim_event_dev, PMT_data_dev, num_events, mdrf_spline_coefs, thresh_spline_coefs, gain_values);
    
    for(index = 0; index < num_events; ++index) {
      estim_event[index] = estim_event_dev[index];
    }
    delete [] PMT_data_dev;
    delete [] estim_event_dev;
  }
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "Elapsed time: " << diff.count() << " s (" << double(num_events) / diff.count() << " events/s)." << std::endl;
  return(estim_event);
}
