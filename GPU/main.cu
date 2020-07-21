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

// To compile: nvcc -std=c++11 -Xcompiler -O2 -Xcompiler "-Wall -Wdouble-promotion -Wparentheses -Wconversion" main.cu -o main

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#define CUDA_SAFE_CALL(__call) do {											\
  cudaError_t __err = __call;												\
  if(__err != cudaSuccess) {												\
    std::cerr << "CUDA driver error " << cudaGetErrorString(__err) << " while calling " #__call << std::endl;		\
    std::cerr << "File: " << __FILE__ << ", function: " << __FUNCTION__ << ", line: " << __LINE__  << std::endl;	\
    throw std::runtime_error("CUDA driver error");									\
  }															\
} while(0)


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void sample_calibr_funct(const calibr_funct_t & calibr_funct);
std::vector<estim_event_t, aligned_allocator<estim_event_t>> contr_grid(const std::vector<PMT_data_t, aligned_allocator<PMT_data_t>> & PMT_data, calibr_funct_t calibr_funct);
__global__ void contr_grid_kernel(estim_event_t *estim_event_dev, PMT_data_t *PMT_data_dev, mdrf_spline_t *mdrf_dev, thresh_spline_t *thresh_dev, float *gain_dev);


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
  write_estim_events(estim_event, "../data/estim_events_GPU.dat");
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
    filename = "../data/mdrf_samples_GPU_" + ss.str() + ".dat";
    write_dat_2d<float, num_sampl_x, num_sampl_y>(data, filename.c_str());
  }
  for(nx = 0; nx < num_sampl_x; ++nx) {
    for(ny = 0; ny < num_sampl_y; ++ny) {
      data[nx][ny] = calibr_funct.thresh(pos_x[nx], pos_y[ny]);
    }
  }
  write_dat_2d<float, num_sampl_x, num_sampl_y>(data, "../data/thresh_samples_GPU.dat");
  return;
}


std::vector<estim_event_t, aligned_allocator<estim_event_t>> contr_grid(const std::vector<PMT_data_t, aligned_allocator<PMT_data_t>> & PMT_data, calibr_funct_t calibr_funct) {
  std::vector<estim_event_t, aligned_allocator<estim_event_t>> estim_event(PMT_data.size());
  std::chrono::time_point<std::chrono::steady_clock> start, end;
  estim_event_t *estim_event_dev;
  thresh_spline_t *thresh_dev;
  PMT_data_t *PMT_data_dev;
  mdrf_spline_t *mdrf_dev;
  unsigned int num_events;
  float *gain_dev;
  
  num_events = (unsigned int) PMT_data.size();
  std::cout << "Number of events: " << num_events << "." << std::endl;
  start = std::chrono::steady_clock::now();
  if(num_events > 0) {
    CUDA_SAFE_CALL(cudaMalloc(& PMT_data_dev, num_events * sizeof(PMT_data[0])));
    CUDA_SAFE_CALL(cudaMalloc(& estim_event_dev, num_events * sizeof(estim_event[0])));
    CUDA_SAFE_CALL(cudaMalloc(& mdrf_dev, NUM_PMTS * sizeof(*mdrf_dev)));
    CUDA_SAFE_CALL(cudaMalloc(& thresh_dev, sizeof(*thresh_dev)));
    CUDA_SAFE_CALL(cudaMalloc(& gain_dev, NUM_PMTS * sizeof(*gain_dev)));
    CUDA_SAFE_CALL(cudaMemcpy(PMT_data_dev, PMT_data.data(), num_events * sizeof(PMT_data[0]), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(mdrf_dev, calibr_funct.mdrf, NUM_PMTS * sizeof(*mdrf_dev), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(thresh_dev, & calibr_funct.thresh, sizeof(*thresh_dev), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(gain_dev, calibr_funct.gain, NUM_PMTS * sizeof(*gain_dev), cudaMemcpyHostToDevice));
    contr_grid_kernel<<<num_events, dim3(SIZE_CONTR_GRID, SIZE_CONTR_GRID)>>>(estim_event_dev, PMT_data_dev, mdrf_dev, thresh_dev, gain_dev);
    CUDA_SAFE_CALL(cudaMemcpy(estim_event.data(), estim_event_dev, num_events * sizeof(estim_event[0]), cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaFree(PMT_data_dev));
    CUDA_SAFE_CALL(cudaFree(estim_event_dev));
    CUDA_SAFE_CALL(cudaFree(mdrf_dev));
    CUDA_SAFE_CALL(cudaFree(thresh_dev));
    CUDA_SAFE_CALL(cudaFree(gain_dev));
  }
  end = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "Elapsed time: " << diff.count() << " s (" << double(num_events) / diff.count() << " events/s)." << std::endl;
  return(estim_event);
}


__global__ void contr_grid_kernel(estim_event_t *estim_event_dev, PMT_data_t *PMT_data_dev, mdrf_spline_t *mdrf_dev, thresh_spline_t *thresh_dev, float *gain_dev) {
  const unsigned int thread_index = threadIdx.y * SIZE_CONTR_GRID + threadIdx.x;
  __shared__ float log_like_values[SIZE_CONTR_GRID][SIZE_CONTR_GRID];
  __shared__ float tmp_data[NUM_PMTS];
  float log_like, max_log_like;
  int max_index_x, max_index_y;
  float current_x, current_y;
  unsigned int event_index;
  bool inside_x, inside_y;
  int index_x, index_y;
  float test_x, test_y;
  float camera_MDRF;
  int pmt, iter;
  float step;
  
  event_index = blockIdx.x;
  if(thread_index < NUM_PMTS) {
    tmp_data[thread_index] = PMT_data_dev[event_index].val[thread_index] / gain_dev[thread_index];
  }
  __syncthreads();
  current_x = current_y = float(1) / float(2);
  step = (float(1) - float(0)) / float(SIZE_CONTR_GRID);
  for(iter = 0; iter < NUM_CONTR_GRID_ITER; ++iter) {
    test_x = current_x + (float(threadIdx.x) - (float(SIZE_CONTR_GRID - 1) / 2.00f)) * step;
    test_y = current_y + (float(threadIdx.y) - (float(SIZE_CONTR_GRID - 1) / 2.00f)) * step;
    inside_x = (float(0) < test_x) && (test_x < float(1));
    inside_y = (float(0) < test_y) && (test_y < float(1));
    if(inside_x && inside_y) {
      log_like = float(0);
      for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
        camera_MDRF = mdrf_dev[pmt](test_x, test_y);
        if((tmp_data[pmt] != float(0)) || (camera_MDRF != float(0))) {
          log_like += tmp_data[pmt] * logf(camera_MDRF) - camera_MDRF;
        }
      }
      log_like_values[threadIdx.x][threadIdx.y] = log_like;
    } else {
      log_like_values[threadIdx.x][threadIdx.y] = -HUGE_VALF;
    }
    __syncthreads();
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
  if(thread_index == 0) {
    inside_x = (float(0) < current_x) && (current_x < float(1));
    inside_y = (float(0) < current_y) && (current_y < float(1));
    if(inside_x && inside_y) {
      log_like = max_log_like;
      for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
        if(tmp_data[pmt] > float(0)) {
          log_like -= lgammaf(tmp_data[pmt] + float(1));
        }
      }
      estim_event_dev[event_index].valid = log_like > (*thresh_dev)(current_x, current_y);
      estim_event_dev[event_index].log_like = log_like;
    } else {
      estim_event_dev[event_index].valid = 0;
    }
    estim_event_dev[event_index].x_pos = CAMERA_MIN_POS + current_x * (CAMERA_MAX_POS - CAMERA_MIN_POS);
    estim_event_dev[event_index].y_pos = CAMERA_MIN_POS + current_y * (CAMERA_MAX_POS - CAMERA_MIN_POS);
  }
  return;
}
