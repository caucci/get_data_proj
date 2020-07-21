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

// To compile: g++ -std=c++11 -pedantic-errors -O2 -Wall -Wdouble-promotion -Wparentheses -Wconversion main.cpp -o main

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void sample_calibr_funct(const calibr_funct_t & calibr_funct);
std::vector<estim_event_t, aligned_allocator<estim_event_t>> contr_grid(const std::vector<PMT_data_t, aligned_allocator<PMT_data_t>> & PMT_data, calibr_funct_t calibr_funct);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char **argv) {
  std::vector<estim_event_t, aligned_allocator<estim_event_t>> estim_event;
  std::vector<PMT_data_t, aligned_allocator<PMT_data_t>> PMT_data;
  calibr_funct_t calibr_funct;
  calibr_data_t calibr_data;
  
  calibr_data = get_calibration_data("../data/camera0_79x79_1.5mm_tc99m_mean", "../data/camera0_thresh.dat", "../data/camera0_79x79_1.5mm_tc99m_gains");
  calibr_funct = get_calibration_funct(calibr_data);
  sample_calibr_funct(calibr_funct);
  PMT_data = get_PMT_data("../data/ResPhantom022516-0mm_00.dat");
  estim_event = contr_grid(PMT_data, calibr_funct);
  write_estim_events(estim_event, "../data/estim_events_CPU.dat");
  return(0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void sample_calibr_funct(const calibr_funct_t & calibr_funct) {
  const int num_sampl_x = 128;
  const int num_sampl_y = 128;
  std::array<std::array<float, num_sampl_y>, num_sampl_x> data;
  float pos_x[num_sampl_x];
  float pos_y[num_sampl_y];
  std::string filename;
  int nx, ny;
  int i, pmt;
  
  for(i = 0; i < num_sampl_x; ++i) {
    pos_x[i] = float(i) / float(num_sampl_x - 1);
  }
  for(i = 0; i < num_sampl_y; ++i) {
    pos_y[i] = float(i) / float(num_sampl_y - 1);
  }
  for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
    for(nx = 0; nx < num_sampl_x; ++nx) {
      for(ny = 0; ny < num_sampl_y; ++ny) {
        data[nx][ny] = calibr_funct.mdrf[pmt](pos_x[nx], pos_y[ny]);
      }
    }
    std::ostringstream ss;
    ss << std::setw(3) << std::setfill('0') << pmt;
    filename = "../data/mdrf_samples_CPU_" + ss.str() + ".dat";
    write_dat_2d<float, num_sampl_x, num_sampl_y>(data, filename.c_str());
  }
  for(nx = 0; nx < num_sampl_x; ++nx) {
    for(ny = 0; ny < num_sampl_y; ++ny) {
      data[nx][ny] = calibr_funct.thresh(pos_x[nx], pos_y[ny]);
    }
  }
  write_dat_2d<float, num_sampl_x, num_sampl_y>(data, "../data/thresh_samples_CPU.dat");
  return;
}


std::vector<estim_event_t, aligned_allocator<estim_event_t>> contr_grid(const std::vector<PMT_data_t, aligned_allocator<PMT_data_t>> & PMT_data, calibr_funct_t calibr_funct) {
  std::vector<estim_event_t, aligned_allocator<estim_event_t>> estim_event(PMT_data.size());
  std::chrono::time_point<std::chrono::steady_clock> start, end;
  float log_like_values[SIZE_CONTR_GRID][SIZE_CONTR_GRID];
  float log_like, max_log_like;
  int max_index_x, max_index_y;
  float current_x, current_y;
  float tmp_data[NUM_PMTS];
  unsigned int event_index;
  unsigned int num_events;
  bool inside_x, inside_y;
  int index_x, index_y;
  float test_x, test_y;
  float camera_MDRF;
  int pmt, iter;
  float step;
  
  num_events = (unsigned int) PMT_data.size();
  std::cout << "Number of events: " << num_events << "." << std::endl;
  start = std::chrono::steady_clock::now();
  for(event_index = 0; event_index < num_events; ++event_index) {
    for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
      tmp_data[pmt] = PMT_data[event_index].val[pmt] / calibr_funct.gain[pmt];
    }
    current_x = current_y = float(1) / float(2);
    step = (float(1) - float(0)) / float(SIZE_CONTR_GRID);
    for(iter = 0; iter < NUM_CONTR_GRID_ITER; ++iter) {
      for(index_x = 0; index_x < SIZE_CONTR_GRID; ++index_x) {
        test_x = current_x + (float(index_x) - (float(SIZE_CONTR_GRID - 1) / 2.00f)) * step;
        for(index_y = 0; index_y < SIZE_CONTR_GRID; ++index_y) {
          test_y = current_y + (float(index_y) - (float(SIZE_CONTR_GRID - 1) / 2.00f)) * step;
          inside_x = (float(0) < test_x) && (test_x < float(1));
          inside_y = (float(0) < test_y) && (test_y < float(1));
          if(inside_x && inside_y) {
            log_like = float(0);
            for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
              camera_MDRF = calibr_funct.mdrf[pmt](test_x, test_y);
              if((tmp_data[pmt] != float(0)) || (camera_MDRF != float(0))) {
                log_like += tmp_data[pmt] * std::log(camera_MDRF) - camera_MDRF;
              }
            }
            log_like_values[index_x][index_y] = log_like;
          } else {
            log_like_values[index_x][index_y] = -HUGE_VALF;
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
      current_x = current_x + (float(max_index_x) - (float(SIZE_CONTR_GRID - 1) / 2.00f)) * step;
      current_y = current_y + (float(max_index_y) - (float(SIZE_CONTR_GRID - 1) / 2.00f)) * step;
      step /= CONTR_FACTOR;
    }
    inside_x = (float(0) < current_x) && (current_x < float(1));
    inside_y = (float(0) < current_y) && (current_y < float(1));
    if(inside_x && inside_y) {
      log_like = max_log_like;
      for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
        if(tmp_data[pmt] > float(0)) {
          log_like -= std::lgamma(tmp_data[pmt] + float(1));
        }
      }
      estim_event[event_index].valid = log_like > calibr_funct.thresh(current_x, current_y);
      estim_event[event_index].log_like = log_like;
    } else {
      estim_event[event_index].valid = 0;
    }
    estim_event[event_index].x_pos = CAMERA_MIN_POS + current_x * (CAMERA_MAX_POS - CAMERA_MIN_POS);
    estim_event[event_index].y_pos = CAMERA_MIN_POS + current_y * (CAMERA_MAX_POS - CAMERA_MIN_POS);
  }
  end = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "Elapsed time: " << diff.count() << " s (" << double(num_events) / diff.count() << " events/s)." << std::endl;
  return(estim_event);
}
