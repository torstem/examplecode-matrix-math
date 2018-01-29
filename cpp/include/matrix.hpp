#ifndef __MATRIX_HPP__
#define __MATRIX_HPP__

typedef double FP;

namespace vector_math {

template <typename M>
struct matrix_base;

template <typename M1, typename M2>
struct mat_mat_add : public matrix_base<mat_mat_add<M1, M2> >
{
  static const unsigned int NUM_ROWS = M1::NUM_ROWS;
  static const unsigned int NUM_COLS = M1::NUM_COLS;

  typedef mat_mat_add<M1, M2> self_type;

  const M1 &m1;
  const M2 &m2;

  mat_mat_add(const M1 &m1_, const M2 &m2_) : m1(m1_), m2(m2_)  {
    static_assert(M1::NUM_ROWS == M2::NUM_ROWS, "M1::NUM_ROWS != M2::NUM_ROWS");
    static_assert(M1::NUM_COLS == M2::NUM_COLS, "M1::NUM_COLS != M2::NUM_COLS");
  }

  FP operator() (unsigned int row, unsigned int col) const
  {
    return m1(row, col) + m2(row, col);
  }
};

template <typename M1, typename M2>
struct mat_mat_sub : public matrix_base<mat_mat_sub<M1, M2> >
{
  static const unsigned int NUM_ROWS = M1::NUM_ROWS;
  static const unsigned int NUM_COLS = M1::NUM_COLS;

  typedef mat_mat_sub<M1, M2> self_type;

  const M1 &m1;
  const M2 &m2;

  mat_mat_sub(const M1 &m1_, const M2 &m2_) : m1(m1_), m2(m2_) {
    static_assert(M1::NUM_ROWS == M2::NUM_ROWS, "M1::NUM_ROWS != M2::NUM_ROWS");
    static_assert(M1::NUM_COLS == M2::NUM_COLS, "M1::NUM_COLS != M2::NUM_COLS");
  }

  FP operator() (unsigned int row, unsigned int col) const
  {
    return m1(row, col) - m2(row, col);
  }

  FP operator[] (unsigned int row) const
  {
    return (*this)(row, 0);
  }
};
  
template <typename M1, typename M2>
struct mat_mat_mul : public matrix_base<mat_mat_mul<M1, M2> >
{
  static const unsigned int NUM_ROWS = M1::NUM_ROWS;
  static const unsigned int NUM_COLS = M2::NUM_COLS;

  typedef mat_mat_mul<M1, M2> self_type;
  
  const M1 &m1;
  const M2 &m2;

  mat_mat_mul(const M1 &m1_, const M2 &m2_) : m1(m1_), m2(m2_)  {
    static_assert(M1::NUM_COLS == M2::NUM_ROWS, "M1::NUM_COLS != M2::NUM_ROWS");
  }

  FP operator() (unsigned int row, unsigned int col) const
  {
    FP sum = 0.0;
    for (int i=0; i<M2::NUM_ROWS; i++) {
      sum += m1(row, i) * m2(i, col);
    }
    return sum;
  }

  FP operator[] (unsigned int row) const
  {
    return (*this)(row, 0);
  }
};

template <typename M>
struct matrix_scalar_mul : public matrix_base<matrix_scalar_mul<M> >
{
  static const unsigned int NUM_ROWS = M::NUM_ROWS;
  static const unsigned int NUM_COLS = M::NUM_COLS;
  
  const M &m;
  FP val;
  
  matrix_scalar_mul(const M &m_, FP val_) : m(m_), val(val_) {}
    
  FP operator() (unsigned int row, unsigned int col) const {
    return val * m(row, col);
  }

  FP operator[] (unsigned int row) const
  {
    return (*this)(row, 0);
  }
};
  
template <typename MC>
struct matrix_base
{
  typedef MC CHILD_TYPE;

  template <typename OTHER_TYPE>
  mat_mat_add<CHILD_TYPE, OTHER_TYPE> operator+(const OTHER_TYPE &m)
  {
    return mat_mat_add<CHILD_TYPE, OTHER_TYPE>(static_cast<CHILD_TYPE&>(*this), m);
  }

  template <typename OTHER_TYPE>
  mat_mat_sub<CHILD_TYPE, OTHER_TYPE> operator-(const OTHER_TYPE &m)
  {
    return mat_mat_sub<CHILD_TYPE, OTHER_TYPE>(static_cast<CHILD_TYPE&>(*this), m);
  }
  
  template <typename OTHER_TYPE>
  mat_mat_mul<CHILD_TYPE, OTHER_TYPE> operator*(const OTHER_TYPE &m)
  {
    return mat_mat_mul<CHILD_TYPE, OTHER_TYPE>(static_cast<CHILD_TYPE&>(*this), m);
  }

  matrix_scalar_mul<CHILD_TYPE> operator*(FP value)
  {
    return matrix_scalar_mul<CHILD_TYPE>(static_cast<CHILD_TYPE&>(*this), value);
  }
};
  
  
template <unsigned int NR, unsigned int NC, unsigned int RS = NC>
struct matrix : public matrix_base<matrix<NR,NC,RS> >
{
  static const unsigned int NUM_ROWS = NR;  
  static const unsigned int NUM_COLS = NC;
  static const unsigned int ROW_STRIDE = RS;

  typedef matrix<NR, NC, RS> self_type;
  
  FP *data;
  matrix(FP *data_) : data(data_) {  }
  
  FP operator()(unsigned int row, unsigned int col) const
  {
    return data[row * ROW_STRIDE + col];
  }
  
  FP& operator()(unsigned int row, unsigned int col)
  {
    return data[row * ROW_STRIDE + col];
  }
  
  FP operator[](unsigned int row) const
  {
    return (*this)(row, 0);
  }

  FP& operator[](unsigned int row)
  {
    return (*this)(row, 0);
  }

  matrix<NR, 1, NR> col(unsigned int i)
  {
    return matrix<NR, 1, NR>(&data[i]);
  }
};

template <unsigned int NR, unsigned int NC>
struct matws : public matrix_base<matws<NR,NC> >
{
  static const unsigned int NUM_ROWS = NR;
  static const unsigned int NUM_COLS = NC;

  typedef matws<NR, NC> self_type;

  FP data[NR * NC];
  matws() { }

  FP operator()(unsigned int row, unsigned int col) const
  {
    return data[row * ROW_STRIDE + col];
  }

  FP& operator()(unsigned int row, unsigned int col)
  {
    return data[row * ROW_STRIDE + col];
  }

  FP operator[](unsigned int row) const
  {
    return (*this)(row, 0);
  }

  FP& operator[](unsigned int row)
  {
    return (*this)(row, 0);
  }

  matrix<NR, 1, NR> col(unsigned int i)
  {
    return matrix<NR, 1, NR>(&data[i]);
  }
};

template <unsigned int NR>
struct e
{
  static const unsigned int NUM_ROWS = NR;
  static const unsigned int NUM_COLS = 1;

  unsigned int N;
  
  e(unsigned int N_) : N(N_) {}
  
  FP operator[] (unsigned int i) const
  {
    if (i == N) {
      return 1;
    }
    return 0;
  }

  FP operator() (unsigned int row, unsigned int col) const
  {
    return (*this)[row];
  }  
};

template <unsigned int N>
struct diag_t : public matrix_base<diag_t<N> >
{
  static const unsigned int NUM_ROWS = N;
  static const unsigned int NUM_COLS = N;

  FP *data;
  diag_t(FP *data_) : data(data_)
  {}

  FP operator() (unsigned int row, unsigned int col) const
  {
    if (row == col) {
      return data[row];
    }

    return 0.0;
  }
};
  
template <typename M>
struct col_t
{
  static const unsigned int NUM_ROWS = M::NUM_ROWS;
  static const unsigned int NUM_COLS = 1;

  const M &m;
  unsigned int col_num;
   col_t(const M &m_, unsigned int col_num_) : m(m_), col_num(col_num_)
  {}

  FP operator() (unsigned int row_num, unsigned int not_used) const
  {
    return m(row_num, col_num);
  }

  FP operator[] (unsigned int row_num) const
  {
    return m(row_num, col_num);
  }
};

template <typename M>
struct transpose_t : public matrix_base<transpose_t<M> >
{
  static const unsigned int NUM_ROWS = M::NUM_COLS;
  static const unsigned int NUM_COLS = M::NUM_ROWS;

  typedef transpose_t<M> self_type;
  
  const M &m;
  transpose_t(const M &m_) : m(m_)
  {
  }

  FP operator() (unsigned int row, unsigned int col) const
  {
    return m(col, row);
  }

  FP & operator() (unsigned int row, unsigned int col)
  {
    return m(col, row);
  }
};


template <typename M>
col_t<M> col(const M &m, unsigned int col_num) {
  return col_t<M>(m, col_num);
}

template <typename M1, typename M2>
void copy_assign(M1 &m1, const M2 &m2)
{
  static_assert(M1::NUM_ROWS == M2::NUM_ROWS, "copy_assign. M1::NUM_ROWS == M2::NUM_ROWS failed.");
  static_assert(M1::NUM_COLS == M2::NUM_COLS, "copy_assign. M1::NUM_COLS == M2::NUM_COLS failed.");

  FP tmp[M1::NUM_ROWS * M1::NUM_COLS];
  for (unsigned int i=0; i<M2::NUM_ROWS; i++) {
    for (unsigned int j=0; j<M2::NUM_COLS; j++) {
      tmp[i * M2::NUM_COLS + j] = m2(i, j);
    }
  }
  for (unsigned int i=0; i<M1::NUM_ROWS; i++) {
    for (unsigned int j=0; j<M1::NUM_COLS; j++) {
      m1(i, j) = tmp[i * M1::NUM_COLS + j];
    }
  }  
}

template <typename M1, typename M2>
void assign(M1 &m1, const M2 &m2)
{
  static_assert(M1::NUM_ROWS == M2::NUM_ROWS, "assign. M1::NUM_ROWS == M2::NUM_ROWS failed.");
  static_assert(M1::NUM_COLS == M2::NUM_COLS, "assign. M1::NUM_COLS == M2::NUM_COLS failed.");

  for (unsigned int i=0; i<M2::NUM_ROWS; i++) {
    for (unsigned int j=0; j<M2::NUM_COLS; j++) {
      m1(i, j) = m2(i, j);
    }
  }
}


inline  FP identity(unsigned int row, unsigned int col)
{
  if (row == col) {
    return 1.0;
  }
  return 0.0;
}

template <typename M>
transpose_t<M> transpose(const M &m)
{
  return transpose_t<M>(m);
}

template <typename M>
transpose_t<M> T(const M &m)
{
  return transpose_t<M>(m);
}
  
template <unsigned int N>
diag_t<N> diag(FP *data) {
  return diag_t<N>(data);
}

template <typename V1, typename V2>
FP dot(const V1 &v1, const V2 &v2)
{
  FP sum = 0.0;
  for (unsigned int i=0; i<V1::NUM_ROWS; i++) {
    sum += v1[i] * v2[i];
  }
  return sum;
}

template <typename M>
FP to_scalar(const M &m) {
  static_assert(M::NUM_COLS == 1, "M::NUM_COLS != 1");
  static_assert(M::NUM_ROWS == 1, "M::NUM_ROWS != 1");
  return m(0, 0);
}

}


#endif
