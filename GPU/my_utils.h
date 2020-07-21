#ifndef _MY_UTILS_H
#define _MY_UTILS_H

#include "my_types.h"


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<class _T> class aligned_allocator {
  public:
    typedef _T value_type;
    typedef _T * pointer;
    typedef _T const * const_pointer;
    typedef _T & reference;
    typedef _T const & const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    static pointer allocate(size_type n, std::allocator<void>::const_pointer hint = 0);
    static void deallocate(pointer p, size_type n);
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<class _T> typename aligned_allocator<_T>::pointer aligned_allocator<_T>::allocate(size_type n, std::allocator<void>::const_pointer hint) {
  void *ptr;
  
  if(posix_memalign(& ptr, 4096, n * sizeof(_T))) {
    throw std::bad_alloc();
  }
  return(reinterpret_cast<pointer>(ptr));
}

  
template<class _T> void aligned_allocator<_T>::deallocate(pointer p, size_type n) {
  free(p);
  return;
}


template<class _T1, class _T2> bool operator==(const aligned_allocator<_T1> & lhs, const aligned_allocator<_T2> & rhs) {
  return(true);
}


template<class _T1, class _T2> bool operator!=(const aligned_allocator<_T1> & lhs, const aligned_allocator<_T2> & rhs) {
  return(!(lhs == rhs));
}


template<class _T, int _NX, int _NY> void write_dat_2d(std::array<std::array<_T, _NY>, _NX> data, const char *filename) {
  const uint32_t num_x = _NX;
  const uint32_t num_y = _NY;
  std::ofstream ofs;
  int nx, ny;
  float tmp;
  
  ofs.open(filename, std::ofstream::out | std::ofstream::binary | std::ofstream::trunc);
  if(!ofs) {
    throw std::runtime_error(std::string("Cannot create file ") + std::string(filename));
  }
  ofs.write(reinterpret_cast<const char *>(& num_x), sizeof(num_x));
  ofs.write(reinterpret_cast<const char *>(& num_y), sizeof(num_y));
  for(nx = 0; nx < _NX; ++nx) {
    for(ny = 0; ny < _NY; ++ny) {
      tmp = (float) data[nx][ny];
      ofs.write(reinterpret_cast<const char *>(& tmp), sizeof(tmp));
    }
  }
  ofs.close();
  return;
}


void write_estim_events(const std::vector<estim_event_t, aligned_allocator<estim_event_t>> & estim_events, const char *filename) {
  uint32_t event_index, num_events;
  float x_pos, y_pos;
  std::ofstream ofs;
  uint32_t valid;
  float log_like;
  
  ofs.open(filename, std::ofstream::out | std::ofstream::binary | std::ofstream::trunc);
  if(!ofs) {
    throw std::runtime_error("Cannot create ML estimates file!");
  }
  num_events = (uint32_t) estim_events.size();
  ofs.write(reinterpret_cast<const char *>(& num_events), sizeof(num_events));
  for(event_index = 0; event_index < num_events; ++event_index) {
    x_pos = estim_events[event_index].x_pos;
    y_pos = estim_events[event_index].y_pos;
    valid = estim_events[event_index].valid ? 1 : 0;
    log_like = estim_events[event_index].log_like;
    ofs.write(reinterpret_cast<const char *>(& valid), sizeof(valid));
    ofs.write(reinterpret_cast<const char *>(& x_pos), sizeof(x_pos));
    ofs.write(reinterpret_cast<const char *>(& y_pos), sizeof(y_pos));
    ofs.write(reinterpret_cast<const char *>(& log_like), sizeof(log_like));
  }
  ofs.close();
  return;
}


std::vector<PMT_data_t, aligned_allocator<PMT_data_t>> get_PMT_data(const char *filename) {
  std::vector<PMT_data_t, aligned_allocator<PMT_data_t>> PMT_data;
  int event_index, num_events;
  PMT_data_t tmp_PMT_data;
  std::ifstream LM_file;
  int16_t LM_header[9];
  int i, pmt;
  
  LM_file.open(filename, std::ifstream::in | std::ifstream::binary);
  if(!LM_file) {
    throw std::runtime_error("Cannot open input LM file!");
  }
  LM_file.read(reinterpret_cast<char *>(LM_header), 9 * sizeof(LM_header[0]));
  for(i = 0; i < 9; ++i) {
    LM_header[i] = __builtin_bswap16(LM_header[i]);
  }
  num_events = LM_header[3] * 1000 + LM_header[4];
  PMT_data.resize(num_events);
  for(event_index = 0; event_index < num_events; ++event_index) {
    LM_file.read(reinterpret_cast<char *>(tmp_PMT_data.val), NUM_PMTS * sizeof(tmp_PMT_data.val[0]));
    for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
      tmp_PMT_data.val[pmt] = __builtin_bswap16(tmp_PMT_data.val[pmt]);
      PMT_data[event_index].val[pmt] = std::max(int16_t(0), tmp_PMT_data.val[pmt]);
    }
  }
  LM_file.close();
  return(PMT_data);
}


calibr_data_t get_calibration_data(const char *mdrf_filename, const char *thresh_filename, const char *gain_filename) {
  calibr_data_t calibr_data;
  std::ifstream ifs;
  int row, col, pmt;
  float tmp;
  
  ifs.open(mdrf_filename, std::ifstream::in | std::ifstream::binary);
  if(!ifs) {
    throw std::runtime_error("Cannot open MDRF calibration file!");
  }
  for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
    for(row = 0; row < NUM_SAMPL; ++row) {
      for(col = 0; col < NUM_SAMPL; ++col) {
        ifs.read(reinterpret_cast<char *>(& tmp), sizeof(tmp));
        calibr_data.mdrf[pmt][row][col] = tmp;
      }
    }
  }
  ifs.close();
  ifs.open(thresh_filename, std::ifstream::in | std::ifstream::binary);
  if(!ifs) {
    throw std::runtime_error("Cannot open threshold file!");
  }
  for(row = 0; row < NUM_SAMPL; ++row) {
    for(col = 0; col < NUM_SAMPL; ++col) {
      ifs.read(reinterpret_cast<char *>(& tmp), sizeof(tmp));
      calibr_data.thresh[row][col] = tmp;
    }
  }
  ifs.close();
  ifs.open(gain_filename, std::ifstream::in | std::ifstream::binary);
  if(!ifs) {
    throw std::runtime_error("Cannot open gain file!");
  }
  for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
    ifs.read(reinterpret_cast<char *>(& tmp), sizeof(tmp));
    calibr_data.gain[pmt] = tmp;
  }
  ifs.close();
  for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
    for(row = 0; row < NUM_SAMPL; ++row) {
      for(col = 0; col < NUM_SAMPL; ++col) {
        calibr_data.mdrf[pmt][row][col] /= calibr_data.gain[pmt];
      }
    }
  }
  return(calibr_data);
}


calibr_funct_t get_calibration_funct(const calibr_data_t & calibr_data) {
  calibr_funct_t calibr_funct;
  float pos[NUM_SAMPL];
  int i, pmt;
  
  for(i = 0; i < NUM_SAMPL; ++i) {
    pos[i] = float(i) / float(NUM_SAMPL - 1);
  }
  for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
    calibr_funct.mdrf[pmt] = spap2<float, float, MX, MY, KX, KY, NUM_SAMPL, NUM_SAMPL>(pos, pos, calibr_data.mdrf[pmt]);
  }
  calibr_funct.thresh = spap2<float, float, MX, MY, KX, KY, NUM_SAMPL, NUM_SAMPL>(pos, pos, calibr_data.thresh);
  for(pmt = 0; pmt < NUM_PMTS; ++pmt) {
    calibr_funct.gain[pmt] = calibr_data.gain[pmt];
  }
  return(calibr_funct);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#endif // _MY_UTILS_H
