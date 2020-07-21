#ifndef _SPLINE_HPP
#define _SPLINE_HPP

#include <iostream>


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<class _V, class _C, int _M, int _K> class spline_1D {
  public:
    spline_1D();
    spline_1D(const _V my_coefs[_M + _K]);
    _V operator()(const _C & x) const;
    void get_coefs(_V output[_M + _K]) const;
    
  private:
    _V coefs[_M + _K];
};


template<class _V, class _C, int _MX, int _MY, int _KX, int _KY> class spline_2D {
  public:
    spline_2D();
    spline_2D(const _V my_coefs[_MY + _KY][_MX + _KX]);
    _V operator()(const _C & x, const _C & y) const;
    void get_coefs(_V output[_MY + _KY][_MX + _KX]) const;
    
  private:
    _V coefs[_MY + _KY][_MX + _KX];
};


template<class _V, class _C, int _MX, int _MY, int _MZ, int _KX, int _KY, int _KZ> class spline_3D {
  public:
    spline_3D();
    spline_3D(const _V my_coefs[_MZ + _KZ][_MY + _KY][_MX + _KX]);
    _V operator()(const _C & x, const _C & y, const _C & z) const;
    void get_coefs(_V output[_MZ + _KZ][_MY + _KY][_MX + _KX]) const;
    
  private:
    _V coefs[_MZ + _KZ][_MY + _KY][_MX + _KX];
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<class _V, class _C, int _M, int _K> spline_1D<_V, _C, _M, _K> spapi(const _C x[_M + _K], const _V v[_M + _K]);
template<class _V, class _C, int _M, int _K, int _L> spline_1D<_V, _C, _M, _K> spap2(const _C x[_L], const _V v[_L]);
template<class _V, class _C, int _MX, int _MY, int _KX, int _KY> spline_2D<_V, _C, _MX, _MY, _KX, _KY> spapi(const _C x[_MX + _KX], const _C y[_MY + _KY], const _V v[_MY + _KY][_MX + _KX]);
template<class _V, class _C, int _MX, int _MY, int _KX, int _KY, int _LX, int _LY> spline_2D<_V, _C, _MX, _MY, _KX, _KY> spap2(const _C x[_LX], const _C y[_LY], const _V v[_LY][_LX]);
template<class _V, class _C, int _MX, int _MY, int _MZ, int _KX, int _KY, int _KZ> spline_3D<_V, _C, _MX, _MY, _MZ, _KX, _KY, _KZ> spapi(const _C x[_MX + _KX], _C y[_MY + _KY], const _C z[_MZ + _KZ], const _V v[_MZ + _KZ][_MY + _KY][_MX + _KX]);
template<class _V, class _C, int _MX, int _MY, int _MZ, int _KX, int _KY, int _KZ, int _LX, int _LY, int _LZ> spline_3D<_V, _C, _MX, _MY, _MZ, _KX, _KY, _KZ> spap2(const _C x[_LX], const _C y[_LY], const _C z[_LZ], const _V v[_LZ][_LY][_LX]);
template<class _C, int _M, int _K> inline int find_span(const _C & x);
template<class _C, int _M, int _K> void evaluate_basis(_C basis[_M], const _C & x, int ell);
template<class _C, int N> void get_inv(_C inv[N][N], const _C matr[N][N]);
template<class _C, int _M, int _K> void get_inv_interp_matr(_C inv[_M + _K][_M + _K], const _C x[_M + _K]);
template<class _C, int _M, int _K, int _L> void get_inv_approx_matr(_C inv[_M + _K][_M + _K], const _C colmat[_L][_M], const int t[_L]);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<class _V, class _C, int _M, int _K> spline_1D<_V, _C, _M, _K>::spline_1D() {
  int i;
  
  for(i = 0; i < (_M + _K); ++i) {
    coefs[i] = _V(_C(0));
  }
}


template<class _V, class _C, int _M, int _K> spline_1D<_V, _C, _M, _K>::spline_1D(const _V my_coefs[_M + _K]) {
  int i;
  
  for(i = 0; i < (_M + _K); ++i) {
    #pragma HLS UNROLL
    #pragma HLS PIPELINE II=1
    coefs[i] = my_coefs[i];
  }
}


template<class _V, class _C, int _M, int _K> _V spline_1D<_V, _C, _M, _K>::operator()(const _C & x) const {
  #pragma HLS ARRAY_PARTITION variable=coefs complete dim=0
  _C basis[_M];
  int ell;
  int i;
  _V s;
  
  s = _V(_C(0));
  ell = find_span<_C, _M, _K>(x);
  if(ell >= 0) {
    evaluate_basis<_C, _M, _K>(basis, x, ell);
    for(i = 0; i < _M; ++i) {
      #pragma HLS UNROLL
      #pragma HLS PIPELINE II=1
      s += coefs[i + ell] * _V(basis[i]);
    }
  }
  return(s);
}


template<class _V, class _C, int _M, int _K> void spline_1D<_V, _C, _M, _K>::get_coefs(_V output[_M + _K]) const {
  int i;
  
  for(i = 0; i < (_M + _K); ++i) {
    output[i] = coefs[i];
  }
  return;
}


template<class _V, class _C, int _MX, int _MY, int _KX, int _KY> spline_2D<_V, _C, _MX, _MY, _KX, _KY>::spline_2D() {
  int i_x, i_y;
  
  for(i_y = 0; i_y < (_MY + _KY); ++i_y) {
    for(i_x = 0; i_x < (_MX + _KX); ++i_x) {
      coefs[i_y][i_x] = _V(_C(0));
    }
  }
}


template<class _V, class _C, int _MX, int _MY, int _KX, int _KY> spline_2D<_V, _C, _MX, _MY, _KX, _KY>::spline_2D(const _V my_coefs[_MY + _KY][_MX + _KX]) {
  int i_x, i_y;
  
  for(i_y = 0; i_y < (_MY + _KY); ++i_y) {
    #pragma HLS UNROLL
    #pragma HLS PIPELINE II=1
    for(i_x = 0; i_x < (_MX + _KX); ++i_x) {
      #pragma HLS UNROLL
      #pragma HLS PIPELINE II=1
      coefs[i_y][i_x] = my_coefs[i_y][i_x];
    }
  }
}


template<class _V, class _C, int _MX, int _MY, int _KX, int _KY> _V spline_2D<_V, _C, _MX, _MY, _KX, _KY>::operator()(const _C & x, const _C & y) const {
  #pragma HLS ARRAY_PARTITION variable=coefs complete dim=0
  _C basis_x[_MX];
  _C basis_y[_MY];
  int ell_x, ell_y;
  int i_x, i_y;
  _V s;
  
  s = _V(_C(0));
  ell_x = find_span<_C, _MX, _KX>(x);
  ell_y = find_span<_C, _MY, _KY>(y);
  if((ell_x >= 0) && (ell_y >= 0)) {
    evaluate_basis<_C, _MX, _KX>(basis_x, x, ell_x);
    evaluate_basis<_C, _MY, _KY>(basis_y, y, ell_y);
    for(i_x = 0; i_x < _MX; ++i_x) {
      #pragma HLS UNROLL
      #pragma HLS PIPELINE II=1
      for(i_y = 0; i_y < _MY; ++i_y) {
        #pragma HLS UNROLL
        #pragma HLS PIPELINE II=1
        s += coefs[i_y + ell_y][i_x + ell_x] * _V(basis_x[i_x] * basis_y[i_y]);
      }
    }
  }
  return(s);
}


template<class _V, class _C, int _MX, int _MY, int _KX, int _KY> void spline_2D<_V, _C, _MX, _MY, _KX, _KY>::get_coefs(_V output[_MY + _KY][_MX + _KX]) const {
  int i_x, i_y;
  
  for(i_y = 0; i_y < (_MY + _KY); ++i_y) {
    for(i_x = 0; i_x < (_MX + _KX); ++i_x) {
      output[i_y][i_x] = coefs[i_y][i_x];
    }
  }
  return;
}


template<class _V, class _C, int _MX, int _MY, int _MZ, int _KX, int _KY, int _KZ> spline_3D<_V, _C, _MX, _MY, _MZ, _KX, _KY, _KZ>::spline_3D() {
  int i_x, i_y, i_z;
  
  for(i_z = 0; i_z < (_MZ + _KZ); ++i_z) {
    for(i_y = 0; i_y < (_MY + _KY); ++i_y) {
      for(i_x = 0; i_x < (_MX + _KX); ++i_x) {
        coefs[i_z][i_y][i_x] = _V(_C(0));
      }
    }
  }
}


template<class _V, class _C, int _MX, int _MY, int _MZ, int _KX, int _KY, int _KZ> spline_3D<_V, _C, _MX, _MY, _MZ, _KX, _KY, _KZ>::spline_3D(const _V my_coefs[_MZ + _KZ][_MY + _KY][_MX + _KX]) {
  int i_x, i_y, i_z;
  
  for(i_z = 0; i_z < (_MZ + _KZ); ++i_z) {
    #pragma HLS UNROLL
    #pragma HLS PIPELINE II=1
    for(i_y = 0; i_y < (_MY + _KY); ++i_y) {
      #pragma HLS UNROLL
      #pragma HLS PIPELINE II=1
      for(i_x = 0; i_x < (_MX + _KX); ++i_x) {
        #pragma HLS UNROLL
        #pragma HLS PIPELINE II=1
        coefs[i_z][i_y][i_x] = my_coefs[i_z][i_y][i_x];
      }
    }
  }
}


template<class _V, class _C, int _MX, int _MY, int _MZ, int _KX, int _KY, int _KZ> _V spline_3D<_V, _C, _MX, _MY, _MZ, _KX, _KY, _KZ>::operator()(const _C & x, const _C & y, const _C & z) const {
  #pragma HLS ARRAY_PARTITION variable=coefs complete dim=0
  _C basis_x[_MX];
  _C basis_y[_MY];
  _C basis_z[_MZ];
  int ell_x, ell_y, ell_z;
  int i_x, i_y, i_z;
  _V s;
  
  s = _V(_C(0));
  ell_x = find_span<_C, _MX, _KX>(x);
  ell_y = find_span<_C, _MY, _KY>(y);
  ell_z = find_span<_C, _MZ, _KZ>(z);
  if((ell_x >= 0) && (ell_y >= 0) && (ell_z >= 0)) {
    evaluate_basis<_C, _MX, _KX>(basis_x, x, ell_x);
    evaluate_basis<_C, _MY, _KY>(basis_y, y, ell_y);
    evaluate_basis<_C, _MZ, _KZ>(basis_z, z, ell_z);
    for(i_x = 0; i_x < _MX; ++i_x) {
      #pragma HLS UNROLL
      #pragma HLS PIPELINE II=1
      for(i_y = 0; i_y < _MY; ++i_y) {
        #pragma HLS UNROLL
        #pragma HLS PIPELINE II=1
        for(i_z = 0; i_z < _MZ; ++i_z) {
          #pragma HLS UNROLL
          #pragma HLS PIPELINE II=1
          s += coefs[i_z + ell_z][i_y + ell_y][i_x + ell_x] * _V(basis_x[i_x] * basis_y[i_y] * basis_z[i_z]);
        }
      }
    }
  }
  return(s);
}


template<class _V, class _C, int _MX, int _MY, int _MZ, int _KX, int _KY, int _KZ> void spline_3D<_V, _C, _MX, _MY, _MZ, _KX, _KY, _KZ>::get_coefs(_V output[_MZ + _KZ][_MY + _KY][_MX + _KX]) const {
  int i_x, i_y, i_z;
  
  for(i_z = 0; i_z < (_MZ + _KZ); ++i_z) {
    for(i_y = 0; i_y < (_MY + _KY); ++i_y) {
      for(i_x = 0; i_x < (_MX + _KX); ++i_x) {
        output[i_z][i_y][i_x] = coefs[i_z][i_y][i_x];
      }
    }
  }
  return;
}


template<class _V, class _C, int _M, int _K> spline_1D<_V, _C, _M, _K> spapi(const _C x[_M + _K], const _V v[_M + _K]) {
  _C inv[_M + _K][_M + _K];
  _V coefs[_M + _K];
  int i, j;
  _V sum;
  
  get_inv_interp_matr<_C, _M, _K>(inv, x);
  for(i = 0; i < (_M + _K); ++i) {
    sum = _V(_C(0));
    for(j = 0; j < (_M + _K); ++j) {
      sum += _V(inv[i][j]) * v[j];
    }
    coefs[i] = sum;
  }
  return(spline_1D<_V, _C, _M, _K>(coefs));
}


template<class _V, class _C, int _M, int _K, int _L> spline_1D<_V, _C, _M, _K> spap2(const _C x[_L], const _V v[_L]) {
  _C inv[_M + _K][_M + _K];
  _C colmat[_L][_M];
  _V vectB[_M + _K];
  _V coefs[_M + _K];
  int t[_L];
  int ell;
  int i;
  int j;
  _V sum;

  for(i = 0; i < (_M + _K); ++i) {
    vectB[i] = _V(_C(0));
  }
  for(ell = 0; ell < _L; ++ell) {
    t[ell] = find_span<_C, _M, _K>(x[ell]);
    evaluate_basis<_C, _M, _K>(colmat[ell], x[ell], t[ell]);
  }
  for(ell = 0; ell < _L; ++ell) {
    for(i = 0; i < _M; ++i) {
      vectB[i + t[ell]] += _V(colmat[ell][i]) * v[ell];
    }
  }
  get_inv_approx_matr<_C, _M, _K, _L>(inv, colmat, t);
  for(j = 0; j < (_M + _K); ++j) {
    sum = _V(_C(0));
    for(i = 0; i < (_M + _K); ++i) {
      sum += _V(inv[j][i]) * vectB[i];
    }
    coefs[j] = sum;
  }
  return(spline_1D<_V, _C, _M, _K>(coefs));
}


template<class _V, class _C, int _MX, int _MY, int _KX, int _KY> spline_2D<_V, _C, _MX, _MY, _KX, _KY> spapi(const _C x[_MX + _KX], const _C y[_MY + _KY], const _V v[_MY + _KY][_MX + _KX]) {
  _C inv_x[_MX + _KX][_MX + _KX];
  _C inv_y[_MY + _KY][_MY + _KY];
  _V new_c[_MY + _KY][_MX + _KX];
  _V coefs[_MY + _KY][_MX + _KX];
  int i_x, i_y;
  int j_x, j_y;
  _V sum;

  get_inv_interp_matr<_C, _MX, _KX>(inv_x, x);
  for(i_y = 0; i_y < (_MY + _KY); ++i_y) {
    for(j_x = 0; j_x < (_MX + _KX); ++j_x) {
      sum = _V(_C(0));
      for(i_x = 0; i_x < (_MX + _KX); ++i_x) {
        sum += _V(inv_x[j_x][i_x]) * v[i_y][i_x];
      }
      new_c[i_y][j_x] = sum;
    }
  }
  get_inv_interp_matr<_C, _MY, _KY>(inv_y, y);
  for(i_x = 0; i_x < (_MX + _KX); ++i_x) {
    for(j_y = 0; j_y < (_MY + _KY); ++j_y) {
      sum = _V(_C(0));
      for(i_y = 0; i_y < (_MY + _KY); ++i_y) {
        sum += _V(inv_y[j_y][i_y]) * new_c[i_y][i_x];
      }
      coefs[j_y][i_x] = sum;
    }
  }
  return(spline_2D<_V, _C, _MX, _MY, _KX, _KY>(coefs));
}


template<class _V, class _C, int _MX, int _MY, int _KX, int _KY, int _LX, int _LY> spline_2D<_V, _C, _MX, _MY, _KX, _KY> spap2(const _C x[_LX], const _C y[_LY], const _V v[_LY][_LX]) {
  _C inv_x[_MX + _KX][_MX + _KX];
  _C inv_y[_MY + _KY][_MY + _KY];
  _V new_c[_MY + _KY][_MX + _KX];
  _V coefs[_MY + _KY][_MX + _KX];
  _V new_v[_MY + _KY][_LX];
  _C colmat_x[_LX][_MX];
  _C colmat_y[_LY][_MY];
  _V vectB[_MX + _KX];
  int t_x[_LX];
  int t_y[_LY];
  int ell_x, ell_y;
  int i_x, i_y;
  int j_x, j_y;
  _V sum;
  
  for(i_y = 0; i_y < (_MY + _KY); ++i_y) {
    for(ell_x = 0; ell_x < _LX; ++ell_x) {
      new_v[i_y][ell_x] = _V(_C(0));
    }
  }
  for(ell_y = 0; ell_y < _LY; ++ell_y) {
    t_y[ell_y] = find_span<_C, _MY, _KY>(y[ell_y]);
    evaluate_basis<_C, _MY, _KY>(colmat_y[ell_y], y[ell_y], t_y[ell_y]);
  }
  for(ell_y = 0; ell_y < _LY; ++ell_y) {
    for(i_y = 0; i_y < _MY; ++i_y) {
      for(ell_x = 0; ell_x < _LX; ++ell_x) {
        new_v[i_y + t_y[ell_y]][ell_x] += _V(colmat_y[ell_y][i_y]) * v[ell_y][ell_x];
      }
    }
  }
  for(ell_x = 0; ell_x < _LX; ++ell_x) {
    t_x[ell_x] = find_span<_C, _MX, _KX>(x[ell_x]);
    evaluate_basis<_C, _MX, _KX>(colmat_x[ell_x], x[ell_x], t_x[ell_x]);
  }
  get_inv_approx_matr<_C, _MX, _KX, _LX>(inv_x, colmat_x, t_x);
  for(i_y = 0; i_y < (_MY + _KY); ++i_y) {
    for(i_x = 0; i_x < (_MX + _KX); ++i_x) {
      vectB[i_x] = _V(_C(0));
    }
    for(ell_x = 0; ell_x < _LX; ++ell_x) {
      for(i_x = 0; i_x < _MX; ++i_x) {
        vectB[i_x + t_x[ell_x]] += _V(colmat_x[ell_x][i_x]) * new_v[i_y][ell_x];
      }
    }
    for(j_x = 0; j_x < (_MX + _KX); ++j_x) {
      sum = _V(_C(0));
      for(i_x = 0; i_x < (_MX + _KX); ++i_x) {
        sum += _V(inv_x[j_x][i_x]) * vectB[i_x];
      }
      new_c[i_y][j_x] = sum;
    }
  }
  get_inv_approx_matr<_C, _MY, _KY, _LY>(inv_y, colmat_y, t_y);
  for(j_x = 0; j_x < (_MX + _KX); ++j_x) {
    for(j_y = 0; j_y < (_MY + _KY); ++j_y) {
      sum = _V(_C(0));
      for(i_y = 0; i_y < (_MY + _KY); ++i_y) {
        sum += _V(inv_y[j_y][i_y]) * new_c[i_y][j_x];
      }
      coefs[j_y][j_x] = sum;
    }
  }
  return(spline_2D<_V, _C, _MX, _MY, _KX, _KY>(coefs));
}


template<class _V, class _C, int _MX, int _MY, int _MZ, int _KX, int _KY, int _KZ> spline_3D<_V, _C, _MX, _MY, _MZ, _KX, _KY, _KZ> spapi(const _C x[_MX + _KX], _C y[_MY + _KY], const _C z[_MZ + _KZ], const _V v[_MZ + _KZ][_MY + _KY][_MX + _KX]) {
  _V new_c[_MZ + _KZ][_MY + _KY][_MX + _KX];
  _V new_d[_MZ + _KZ][_MY + _KY][_MX + _KX];
  _V coefs[_MZ + _KZ][_MY + _KY][_MX + _KX];
  _C inv_x[_MX + _KX][_MX + _KX];
  _C inv_y[_MY + _KY][_MY + _KY];
  _C inv_z[_MZ + _KZ][_MZ + _KZ];
  int i_x, i_y, i_z;
  int j_x, j_y, j_z;
  _V sum;

  get_inv_interp_matr<_C, _MX, _KX>(inv_x, x);
  for(i_y = 0; i_y < (_MY + _KY); ++i_y) {
    for(i_z = 0; i_z < (_MZ + _KZ); ++i_z) {
      for(j_x = 0; j_x < (_MX + _KX); ++j_x) {
        sum = _V(_C(0));
        for(i_x = 0; i_x < (_MX + _KX); ++i_x) {
          sum += _V(inv_x[j_x][i_x]) * v[i_z][i_y][i_x];
        }
        new_c[i_z][i_y][j_x] = sum;
      }
    }
  }
  get_inv_interp_matr<_C, _MY, _KY>(inv_y, y);
  for(i_x = 0; i_x < (_MX + _KX); ++i_x) {
    for(i_z = 0; i_z < (_MZ + _KZ); ++i_z) {
      for(j_y = 0; j_y < (_MY + _KY); ++j_y) {
        sum = _V(_C(0));
        for(i_y = 0; i_y < (_MY + _KY); ++i_y) {
          sum += _V(inv_y[j_y][i_y]) * new_c[i_z][i_y][i_x];
        }
        new_d[i_z][j_y][i_x] = sum;
      }
    }
  }
  get_inv_interp_matr<_C, _MZ, _KZ>(inv_z, z);
  for(i_x = 0; i_x < (_MX + _KX); ++i_x) {
    for(i_y = 0; i_y < (_MY + _KY); ++i_y) {
      for(j_z = 0; j_z < (_MZ + _KZ); ++j_z) {
        sum = _V(_C(0));
        for(i_z = 0; i_z < (_MZ + _KZ); ++i_z) {
          sum += _V(inv_z[j_z][i_z]) * new_d[i_z][i_y][i_x];
        }
        coefs[j_z][i_y][i_x] = sum;
      }
    }
  }
  return(spline_3D<_V, _C, _MX, _MY, _MZ, _KX, _KY, _KZ>(coefs));
}


template<class _V, class _C, int _MX, int _MY, int _MZ, int _KX, int _KY, int _KZ, int _LX, int _LY, int _LZ> spline_3D<_V, _C, _MX, _MY, _MZ, _KX, _KY, _KZ> spap2(const _C x[_LX], const _C y[_LY], const _C z[_LZ], const _V v[_LZ][_LY][_LX]) {
  _V coefs[_MZ + _KZ][_MY + _KY][_MX + _KX];
  _V new_c[_MZ + _KZ][_MY + _KY][_MX + _KX];
  _V new_d[_MZ + _KZ][_MY + _KY][_MX + _KX];
  _V new_v[_MZ + _KZ][_MY + _KY][_LX];
  _C inv_x[_MX + _KX][_MX + _KX];
  _C inv_y[_MY + _KY][_MY + _KY];
  _C inv_z[_MZ + _KZ][_MZ + _KZ];
  _C colmat_x[_LX][_MX];
  _C colmat_y[_LY][_MY];
  _C colmat_z[_LZ][_MZ];
  _V vectB[_MX + _KX];
  int t_x[_LX];
  int t_y[_LY];
  int t_z[_LZ];
  int ell_x, ell_y, ell_z;
  int i_x, i_y, i_z;
  int j_x, j_y, j_z;
  _V sum, tmp;
  
  for(i_z = 0; i_z < (_MZ + _KZ); ++i_z) {
    for(i_y = 0; i_y < (_MY + _KY); ++i_y) {
      for(ell_x = 0; ell_x < _LX; ++ell_x) {
        new_v[i_z][i_y][ell_x] = _V(_C(0));
      }
    }
  }
  for(ell_z = 0; ell_z < _LZ; ++ell_z) {
    t_z[ell_z] = find_span<_C, _MZ, _KZ>(z[ell_z]);
    evaluate_basis<_C, _MZ, _KZ>(colmat_z[ell_z], z[ell_z], t_z[ell_z]);
  }
  for(ell_y = 0; ell_y < _LY; ++ell_y) {
    t_y[ell_y] = find_span<_C, _MY, _KY>(y[ell_y]);
    evaluate_basis<_C, _MY, _KY>(colmat_y[ell_y], y[ell_y], t_y[ell_y]);
  }
  for(ell_y = 0; ell_y < _LY; ++ell_y) {
    for(ell_z = 0; ell_z < _LZ; ++ell_z) {
      for(i_y = 0; i_y < _MY; ++i_y) {
        for(i_z = 0; i_z < _MZ; ++i_z) {
          tmp = _V(colmat_z[ell_z][i_z] * colmat_y[ell_y][i_y]);
          for(ell_x = 0; ell_x < _LX; ++ell_x) {
            new_v[i_z + t_z[ell_z]][i_y + t_y[ell_y]][ell_x] += tmp * v[ell_z][ell_y][ell_x];
          }
        }
      }
    }
  }
  for(ell_x = 0; ell_x < _LX; ++ell_x) {
    t_x[ell_x] = find_span<_C, _MX, _KX>(x[ell_x]);
    evaluate_basis<_C, _MX, _KX>(colmat_x[ell_x], x[ell_x], t_x[ell_x]);
  }
  get_inv_approx_matr<_C, _MX, _KX, _LX>(inv_x, colmat_x, t_x);
  for(i_y = 0; i_y < (_MY + _KY); ++i_y) {
    for(i_z = 0; i_z < (_MZ + _KZ); ++i_z) {
      for(i_x = 0; i_x < (_MX + _KX); ++i_x) {
        vectB[i_x] = _V(_C(0));
      }
      for(ell_x = 0; ell_x < _LX; ++ell_x) {
        for(i_x = 0; i_x < _MX; ++i_x) {
          vectB[i_x + t_x[ell_x]] += _V(colmat_x[ell_x][i_x]) * new_v[i_z][i_y][ell_x];
        }
      }
      for(j_x = 0; j_x < (_MX + _KX); ++j_x) {
        sum = _V(_C(0));
        for(i_x = 0; i_x < (_MX + _KX); ++i_x) {
          sum += _V(inv_x[j_x][i_x]) * vectB[i_x];
        }
        new_c[i_z][i_y][j_x] = sum;
      }
    }
  }
  get_inv_approx_matr<_C, _MY, _KY, _LY>(inv_y, colmat_y, t_y);
  for(j_x = 0; j_x < (_MX + _KX); ++j_x) {
    for(i_z = 0; i_z < (_MZ + _KZ); ++i_z) {
      for(j_y = 0; j_y < (_MY + _KY); ++j_y) {
        sum = _V(_C(0));
        for(i_y = 0; i_y < (_MY + _KY); ++i_y) {
          sum += _V(inv_y[j_y][i_y]) * new_c[i_z][i_y][j_x];
        }
        new_d[i_z][j_y][j_x] = sum;
      }
    }
  }
  get_inv_approx_matr<_C, _MZ, _KZ, _LZ>(inv_z, colmat_z, t_z);
  for(j_x = 0; j_x < (_MX + _KX); ++j_x) {
    for(j_y = 0; j_y < (_MY + _KY); ++j_y) {
      for(j_z = 0; j_z < (_MZ + _KZ); ++j_z) {
        sum = _V(_C(0));
        for(i_z = 0; i_z < (_MZ + _KZ); ++i_z) {
          sum += _V(inv_z[j_z][i_z]) * new_d[i_z][j_y][j_x];
        }
        coefs[j_z][j_y][j_x] = sum;
      }
    }
  }
  return(spline_3D<_V, _C, _MX, _MY, _MZ, _KX, _KY, _KZ>(coefs));
}


template<class _C, int _M, int _K> int find_span(const _C & x) {
  // Do not force inlining of this function as it might crash Vivado HLS!
  return(((x < _C(0)) || (x > _C(1))) ? -1 : ((x != _C(1)) ? int(x * _C(_K + 1)) : _K));
}


template<class _C, int _M, int _K> void evaluate_basis(_C basis[_M], const _C & x, int ell) {
  #pragma HLS ARRAY_PARTITION variable=basis complete
  _C saved, tmp, t;
  _C buff[2][_M];
  int m, j;
  
  if(ell >= 0) {
    buff[0][0] = _C(1);
    for(m = 1; m < (_M - 1); ++m) {
      #pragma HLS PIPELINE II=1
      saved = _C(0);
      for(j = 0; j < m; ++j) {
        tmp = buff[(~m) & 1][j] * _C(_K + 1) / _C(m);
        t = _C(ell + j + 1) / _C(_K + 1);
        buff[m & 1][j] = saved + (t - x) * tmp;
        t = _C(ell + j - m + 1) / _C(_K + 1);
        saved = (x - t) * tmp;
      }
      buff[m & 1][m] = saved;
    }
    saved = _C(0);
    for(j = 0; j < (_M - 1); ++j) {
      #pragma HLS PIPELINE II=1
      tmp = buff[_M & 1][j] * _C(_K + 1) / _C(_M - 1);
      t = _C(ell + j + 1) / _C(_K + 1);
      basis[j] = saved + (t - x) * tmp;
      t = _C(ell + j - _M + 2) / _C(_K + 1);
      saved = (x - t) * tmp;
    }
    basis[_M - 1] = saved;
  }
  return;
}


template<class _C, int N> void get_inv(_C inv[N][N], const _C matr[N][N]) {
  _C tmp_matr[N][N];
  _C tmp_col[N];
  int perm[N];
  _C pivot, coeff;
  int i, j, k;
  
  for(i = 0; i < N; ++i) {
    for(j = 0; j < N; ++j) {
      tmp_matr[i][j] = matr[i][j];
      inv[i][j] = _C(0);
    }
  }
  for(i = 0; i < N; ++i) {
    perm[i] = i;
  }
  for(j = 0; j < (N - 1); ++j) {
    pivot = _C(0);
    k = j;
    for(i = j; i < N; ++i) {
      coeff = tmp_matr[i][j];
      if(coeff < _C(0)) {
        coeff = -coeff;
      }
      if(coeff > pivot) {
        pivot = coeff;
        k = i;
      }
    }
    if(k != j) {
      i = perm[k];
      perm[k] = perm[j];
      perm[j] = i;
      for(i = 0; i < N; ++i) {
        coeff = tmp_matr[k][i];
        tmp_matr[k][i] = tmp_matr[j][i];
        tmp_matr[j][i] = coeff;
      }
    }
    for(i = (j + 1); i < N; ++i) {
      tmp_matr[i][j] /= tmp_matr[j][j];
      for(k = (j + 1); k < N; ++k) {
        tmp_matr[i][k] -= tmp_matr[j][k] * tmp_matr[i][j];
      }
    }
  }
  for(k = 0; k < N; ++k) {
    for(i = 0; i < N; ++i) {
      tmp_col[i] = (k == perm[i]) ? _C(1) : _C(0);
      for(j = 0; j < i; ++j) {
        tmp_col[i] -= tmp_col[j] * tmp_matr[i][j];
      }
    }
    for(i = (N - 1); i >= 0; --i) {
      for(j = (i + 1); j < N; ++j) {
        tmp_col[i] -= inv[j][k] * tmp_matr[i][j];
      }
      inv[i][k] = tmp_col[i] / tmp_matr[i][i];
    }
  }
  return;
}


template<class _C, int _M, int _K> void get_inv_interp_matr(_C inv[_M + _K][_M + _K], const _C x[_M + _K]) {
  _C mat[_M + _K][_M + _K];
  _C basis[_M];
  int i, j, t;
  
  for(i = 0; i < (_M + _K); ++i) {
    for(j = 0; j < (_M + _K); ++j) {
      mat[i][j] = _C(0);
    }
    t = find_span<_C, _M, _K>(x[i]);
    evaluate_basis<_C, _M, _K>(basis, x[i], t);
    for(j = 0; j < _M; ++j) {
      mat[i][j + t] = basis[j];
    }
  }
  get_inv<_C, _M + _K>(inv, mat);
  return;
}


template<class _C, int _M, int _K, int _L> void get_inv_approx_matr(_C inv[_M + _K][_M + _K], const _C colmat[_L][_M], const int t[_L]) {
  _C mat[_M + _K][_M + _K];
  int ell, i, j;
  _C tmp;
  
  for(i = 0; i < (_M + _K); ++i) {
    for(j = 0; j < (_M + _K); ++j) {
      mat[i][j] = _C(0);
    }
  }
  for(ell = 0; ell < _L; ++ell) {
    for(i = 0; i < _M; ++i) {
      for(j = 0; j <= i; ++j) {
        tmp = colmat[ell][i] * colmat[ell][j];
        mat[i + t[ell]][j + t[ell]] += tmp;
      }
    }
  }
  for(i = 0; i < (_M + _K); ++i) {
    for(j = (i + 1); j < (_M + _K); ++j) {
      mat[i][j] = mat[j][i];
    }
  }
  get_inv<_C, _M + _K>(inv, mat);
  return;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#endif // _SPLINE_HPP
