#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdint>
#include <vector>
#include <chrono>
#include <array>
#include <cmath>
#include <CL/cl.hpp>
#include <unistd.h>
#include "spline.hpp"
#include "my_defines.h"
#include "my_types.h"
#include "my_utils.h"

// To compile: g++ main.cpp -O2 -Wall -std=c++11 -lOpenCL -o main -Wno-unknown-pragmas

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#define OCL_CHECK(error, call)												\
  call;															\
  if(error != CL_SUCCESS) {												\
    std::cout << __FILE__ << ":" << __LINE__ << " Error calling " #call ", error code is: " << error << std::endl;	\
    throw std::runtime_error("OpenCL runtime error code: " + get_error_code_string(error));				\
  }															\


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::string get_error_code_string(cl_int error);
std::vector<unsigned char> read_binary_file(const std::string & xclbin_filename);
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
  write_estim_events(estim_event, "../data/estim_events_FPGA.dat");
  return(0);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::string get_error_code_string(cl_int error) {
  switch(error) {
    case CL_SUCCESS:					return(std::string("CL_SUCCESS"));
    case CL_DEVICE_NOT_FOUND:				return(std::string("CL_DEVICE_NOT_FOUND"));
    case CL_DEVICE_NOT_AVAILABLE:			return(std::string("CL_DEVICE_NOT_AVAILABLE"));
    case CL_COMPILER_NOT_AVAILABLE:			return(std::string("CL_COMPILER_NOT_AVAILABLE"));
    case CL_MEM_OBJECT_ALLOCATION_FAILURE:		return(std::string("CL_MEM_OBJECT_ALLOCATION_FAILURE"));
    case CL_OUT_OF_RESOURCES:				return(std::string("CL_OUT_OF_RESOURCES"));
    case CL_OUT_OF_HOST_MEMORY:				return(std::string("CL_OUT_OF_HOST_MEMORY"));
    case CL_PROFILING_INFO_NOT_AVAILABLE:		return(std::string("CL_PROFILING_INFO_NOT_AVAILABLE"));
    case CL_MEM_COPY_OVERLAP:				return(std::string("CL_MEM_COPY_OVERLAP"));
    case CL_IMAGE_FORMAT_MISMATCH:			return(std::string("CL_IMAGE_FORMAT_MISMATCH"));
    case CL_IMAGE_FORMAT_NOT_SUPPORTED:			return(std::string("CL_IMAGE_FORMAT_NOT_SUPPORTED"));
    case CL_BUILD_PROGRAM_FAILURE:			return(std::string("CL_BUILD_PROGRAM_FAILURE"));
    case CL_MAP_FAILURE:				return(std::string("CL_MAP_FAILURE"));
    case CL_MISALIGNED_SUB_BUFFER_OFFSET:		return(std::string("CL_MISALIGNED_SUB_BUFFER_OFFSET"));
    case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:	return(std::string("CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST"));
    case CL_COMPILE_PROGRAM_FAILURE:			return(std::string("CL_COMPILE_PROGRAM_FAILURE"));
    case CL_LINKER_NOT_AVAILABLE:			return(std::string("CL_LINKER_NOT_AVAILABLE"));
    case CL_LINK_PROGRAM_FAILURE:			return(std::string("CL_LINK_PROGRAM_FAILURE"));
    case CL_DEVICE_PARTITION_FAILED:			return(std::string("CL_DEVICE_PARTITION_FAILED"));
    case CL_KERNEL_ARG_INFO_NOT_AVAILABLE:		return(std::string("CL_KERNEL_ARG_INFO_NOT_AVAILABLE"));
    case CL_INVALID_VALUE:				return(std::string("CL_INVALID_VALUE"));
    case CL_INVALID_DEVICE_TYPE:			return(std::string("CL_INVALID_DEVICE_TYPE"));
    case CL_INVALID_PLATFORM:				return(std::string("CL_INVALID_PLATFORM"));
    case CL_INVALID_DEVICE:				return(std::string("CL_INVALID_DEVICE"));
    case CL_INVALID_CONTEXT:				return(std::string("CL_INVALID_CONTEXT"));
    case CL_INVALID_QUEUE_PROPERTIES:			return(std::string("CL_INVALID_QUEUE_PROPERTIES"));
    case CL_INVALID_COMMAND_QUEUE:			return(std::string("CL_INVALID_COMMAND_QUEUE"));
    case CL_INVALID_HOST_PTR:				return(std::string("CL_INVALID_HOST_PTR"));
    case CL_INVALID_MEM_OBJECT:				return(std::string("CL_INVALID_MEM_OBJECT"));
    case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:		return(std::string("CL_INVALID_IMAGE_FORMAT_DESCRIPTOR"));
    case CL_INVALID_IMAGE_SIZE:				return(std::string("CL_INVALID_IMAGE_SIZE"));
    case CL_INVALID_SAMPLER:				return(std::string("CL_INVALID_SAMPLER"));
    case CL_INVALID_BINARY:				return(std::string("CL_INVALID_BINARY"));
    case CL_INVALID_BUILD_OPTIONS:			return(std::string("CL_INVALID_BUILD_OPTIONS"));
    case CL_INVALID_PROGRAM:				return(std::string("CL_INVALID_PROGRAM"));
    case CL_INVALID_PROGRAM_EXECUTABLE:			return(std::string("CL_INVALID_PROGRAM_EXECUTABLE"));
    case CL_INVALID_KERNEL_NAME:			return(std::string("CL_INVALID_KERNEL_NAME"));
    case CL_INVALID_KERNEL_DEFINITION:			return(std::string("CL_INVALID_KERNEL_DEFINITION"));
    case CL_INVALID_KERNEL:				return(std::string("CL_INVALID_KERNEL"));
    case CL_INVALID_ARG_INDEX:				return(std::string("CL_INVALID_ARG_INDEX"));
    case CL_INVALID_ARG_VALUE:				return(std::string("CL_INVALID_ARG_VALUE"));
    case CL_INVALID_ARG_SIZE:				return(std::string("CL_INVALID_ARG_SIZE"));
    case CL_INVALID_KERNEL_ARGS:			return(std::string("CL_INVALID_KERNEL_ARGS"));
    case CL_INVALID_WORK_DIMENSION:			return(std::string("CL_INVALID_WORK_DIMENSION"));
    case CL_INVALID_WORK_GROUP_SIZE:			return(std::string("CL_INVALID_WORK_GROUP_SIZE"));
    case CL_INVALID_WORK_ITEM_SIZE:			return(std::string("CL_INVALID_WORK_ITEM_SIZE"));
    case CL_INVALID_GLOBAL_OFFSET:			return(std::string("CL_INVALID_GLOBAL_OFFSET"));
    case CL_INVALID_EVENT_WAIT_LIST:			return(std::string("CL_INVALID_EVENT_WAIT_LIST"));
    case CL_INVALID_EVENT:				return(std::string("CL_INVALID_EVENT"));
    case CL_INVALID_OPERATION:				return(std::string("CL_INVALID_OPERATION"));
    case CL_INVALID_GL_OBJECT:				return(std::string("CL_INVALID_GL_OBJECT"));
    case CL_INVALID_BUFFER_SIZE:			return(std::string("CL_INVALID_BUFFER_SIZE"));
    case CL_INVALID_MIP_LEVEL:				return(std::string("CL_INVALID_MIP_LEVEL"));
    case CL_INVALID_GLOBAL_WORK_SIZE:			return(std::string("CL_INVALID_GLOBAL_WORK_SIZE"));
    case CL_INVALID_PROPERTY:				return(std::string("CL_INVALID_PROPERTY"));
    case CL_INVALID_IMAGE_DESCRIPTOR:			return(std::string("CL_INVALID_IMAGE_DESCRIPTOR"));
    case CL_INVALID_COMPILER_OPTIONS:			return(std::string("CL_INVALID_COMPILER_OPTIONS"));
    case CL_INVALID_LINKER_OPTIONS:			return(std::string("CL_INVALID_LINKER_OPTIONS"));
    case CL_INVALID_DEVICE_PARTITION_COUNT:		return(std::string("CL_INVALID_DEVICE_PARTITION_COUNT"));
    case CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR:	return(std::string("CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR"));
    case CL_PLATFORM_NOT_FOUND_KHR:			return(std::string("CL_PLATFORM_NOT_FOUND_KHR"));
    default:						return(std::string("Unknown OpenCL error!"));
  }
}


std::vector<unsigned char> read_binary_file(const std::string & xclbin_filename) {
  std::vector<unsigned char> buff;
  std::ifstream f;
  int num_bytes;
  
  std::cout << "INFO: Reading file " << xclbin_filename << std::endl;
  if(access(xclbin_filename.c_str(), R_OK) != 0) {
    throw std::runtime_error("ERROR: " + xclbin_filename + "xclbin not available; please build!");
  }
  f.open(xclbin_filename.c_str(), std::ifstream::binary);
  f.seekg(0, f.end);
  num_bytes = f.tellg();
  f.seekg(0, f.beg);
  buff.resize(num_bytes);
  f.read(reinterpret_cast<char *>(buff.data()), num_bytes);
  return(buff);
}


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
    filename = "../data/mdrf_samples_FPGA_" + ss.str() + ".dat";
    write_dat_2d<float, num_sampl_x, num_sampl_y>(data, filename.c_str());
  }
  for(nx = 0; nx < num_sampl_x; ++nx) {
    for(ny = 0; ny < num_sampl_y; ++ny) {
      data[nx][ny] = calibr_funct.thresh(pos_x[nx], pos_y[ny]);
    }
  }
  write_dat_2d<float, num_sampl_x, num_sampl_y>(data, "../data/thresh_samples_FPGA.dat");
  return;
}


std::vector<estim_event_t, aligned_allocator<estim_event_t>> contr_grid(const std::vector<PMT_data_t, aligned_allocator<PMT_data_t>> & PMT_data, calibr_funct_t calibr_funct) {
  std::vector<float, aligned_allocator<float>> mdrf_spline_coefs((MY + KY) * (MX + KX) * NUM_PMTS);
  std::vector<estim_event_t, aligned_allocator<estim_event_t>> estim_event(PMT_data.size());
  std::vector<float, aligned_allocator<float>> thresh_spline_coefs((MY + KY) * (MX + KX));
  std::vector<float, aligned_allocator<float>> gain_values(NUM_PMTS);
  std::chrono::time_point<std::chrono::steady_clock> start, end;
  float tmp_spline_coefs[MY + KY][MX + KX];
  std::vector<unsigned char> file_buff;
  std::vector<cl::Platform> platforms;
  std::vector<cl::Device> devices;
  std::string platform_name;
  unsigned int num_events;
  cl::Platform platform;
  cl::Device device;
  cl::Event event;
  int i_x, i_y, n;
  unsigned int i;
  float coeff;
  cl_int err;
  bool found;
  int pmt;
  
  found = false;
  OCL_CHECK(err, err = cl::Platform::get(& platforms));
  for(i = 0; (i < platforms.size()) && !found; i++) {
    platform = platforms[i];
    OCL_CHECK(err, platform_name = platforms[i].getInfo<CL_PLATFORM_NAME>(& err));
    std::cout << "Platform Name: " << platform_name << std::endl;
    std::cout << "Platform Vendor: " << platforms[i].getInfo<CL_PLATFORM_VENDOR>() << std::endl;
    std::cout << "Platform Version: " << platforms[i].getInfo<CL_PLATFORM_VERSION>() << std::endl;
    std::cout << "Platform Profile: " << platforms[i].getInfo<CL_PLATFORM_PROFILE>() << std::endl << std::endl;
    found = (platform_name.c_str() == std::string("Xilinx"));
  }
  if(!found) {
    std::cout << "Error: Failed to find Xilinx platform!" << std::endl;
    return(std::vector<estim_event_t, aligned_allocator<estim_event_t>>());
  }
  OCL_CHECK(err, err = platform.getDevices(CL_DEVICE_TYPE_ACCELERATOR, & devices));
  if(devices.size() < 1) {
    std::cout << "Error: No device found!" << std::endl;
    return(std::vector<estim_event_t, aligned_allocator<estim_event_t>>());
  }
  for(i = 0; i < devices.size(); ++i) {
    device = devices[i];
    std::cout << "\tDevice " << i + 1 << std::endl;
    std::cout << "\tName: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
    std::cout << "\tVendor: " << device.getInfo<CL_DEVICE_VENDOR>() << std::endl;
    std::cout << "\tOpenCL C Version: " << device.getInfo<CL_DEVICE_OPENCL_C_VERSION>() << std::endl;
    std::cout << std::endl;
  }
  device = devices[0];
  OCL_CHECK(err, cl::Context context(device, NULL, NULL, NULL, & err));
  OCL_CHECK(err, cl::CommandQueue q(context, device, CL_QUEUE_PROFILING_ENABLE, & err));
  file_buff = read_binary_file("contr_grid.xclbin");
  cl::Program::Binaries binaries = {{file_buff.data(), file_buff.size()}};
  devices.resize(1);
  OCL_CHECK(err, cl::Program program(context, devices, binaries, NULL, & err));
  OCL_CHECK(err, cl::Kernel kernel(program, "contr_grid_kernel", & err));
  num_events = PMT_data.size();
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
      n = MAP_2D(MX + KX, MY + KY, i_x, i_y);
      thresh_spline_coefs[n] = coeff;
    }
  }
  for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
    gain_values[pmt] = calibr_funct.gain[pmt];
  }
  start = std::chrono::steady_clock::now();
  if(num_events > 0) {
    OCL_CHECK(err, cl::Buffer estim_event_dev(context, CL_MEM_WRITE_ONLY, num_events * sizeof(estim_event[0]), NULL, & err));
    OCL_CHECK(err, cl::Buffer PMT_data_dev(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, num_events * sizeof(PMT_data[0]), const_cast<void *>(reinterpret_cast<void const *>(PMT_data.data())), & err));
    OCL_CHECK(err, cl::Buffer mdrf_spline_coefs_dev(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, (MY + KY) * (MX + KX) * NUM_PMTS * sizeof(mdrf_spline_coefs[0]), mdrf_spline_coefs.data(), & err));
    OCL_CHECK(err, cl::Buffer thresh_spline_coefs_dev(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, (MY + KY) * (MX + KX) * sizeof(thresh_spline_coefs[0]), thresh_spline_coefs.data(), & err));
    OCL_CHECK(err, cl::Buffer gain_values_dev(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, NUM_PMTS * sizeof(gain_values[0]), gain_values.data(), & err));
    OCL_CHECK(err, err = kernel.setArg(0, estim_event_dev));
    OCL_CHECK(err, err = kernel.setArg(1, PMT_data_dev));
    OCL_CHECK(err, err = kernel.setArg(2, num_events));
    OCL_CHECK(err, err = kernel.setArg(3, mdrf_spline_coefs_dev));
    OCL_CHECK(err, err = kernel.setArg(4, thresh_spline_coefs_dev));
    OCL_CHECK(err, err = kernel.setArg(5, gain_values_dev));
    OCL_CHECK(err, err = q.enqueueTask(kernel, NULL, & event));
    OCL_CHECK(err, err = cl::copy(q, estim_event_dev, estim_event.begin(), estim_event.end()));
    OCL_CHECK(err, err = q.finish());
  }
  end = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "Elapsed time: " << diff.count() << " s (" << double(num_events) / diff.count() << " events/s)." << std::endl;
  return(estim_event);
}
