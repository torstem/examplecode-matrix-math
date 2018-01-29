#ifndef __MATRIX_HPP__
#define __MATRIX_HPP__

template <typename M>
struct base;

template <typename M1, typename M2>
struct mat_mat_add : public base<mat_mat_add<M1, M2> >
{
  const M1 &m1;
  const M2 &m2;

  mat_mat_add(const M1 &m1_, const M2 &m2_) : m1(m1_), m2(m2_)
  {  }

  float operator()(unsigned int row, unsigned int col) const
  {
    return m1(row, col) + m2(row, col);
  }
};

template <typename M1, typename M2>
struct mat_mat_sub : public base<mat_mat_sub<M1, M2> >
{
  const M1 &m1;
  const M2 &m2;

  mat_mat_sub(const M1 &m1_, const M2 &m2_) : m1(m1_), m2(m2_)
  {  }

  float operator()(unsigned int row, unsigned int col) const
  {
    return m1(row, col) - m2(row, col);
  }
};

template <typename M>
struct base
{
  typedef M CHILD_TYPE;

  template <typename OTHER_TYPE>
  mat_mat_add<CHILD_TYPE, OTHER_TYPE> operator+(const OTHER_TYPE &other)
  {
    return mat_mat_add<CHILD_TYPE, OTHER_TYPE>(static_cast<CHILD_TYPE&>(*this), other);
  }

  template <typename OTHER_TYPE>
  mat_mat_sub<CHILD_TYPE, OTHER_TYPE> operator-(const OTHER_TYPE &other)
  {
    return mat_mat_sub<CHILD_TYPE, OTHER_TYPE>(static_cast<CHILD_TYPE&>(*this), other);
  }  
};


struct matrix : public base<matrix>
{
  float *data;
  matrix(float *data_) : data(data_)
  {
  }

  float operator()(unsigned int row, unsigned int col) const
  {
    return data[row * 3 + col];
  }
  
};

#endif
