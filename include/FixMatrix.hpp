#ifndef __FIXMATRIX_HPP__
#define __FIXMATRIX_HPP__
#include <stdexcept>

//start of container
template<class T, std::size_t ... RestD>
struct FixMatrix;

template<class T, std::size_t PrimaryD >
struct FixMatrix<T, PrimaryD> {
  typedef T type[PrimaryD];
  type data;
  T& at(std::size_t index) {
    if(index >= PrimaryD) {
      throw std::out_of_range("FixMatrix");
    }
    return data[index];
  }

  static size_t dimension() {
    return 1;
  }
};

template<class T, std::size_t PrimaryD, std::size_t ... RestD >
struct FixMatrix<T, PrimaryD, RestD...> {
  typedef FixMatrix<T, RestD...> sub_dimension_Matrix;
  typedef sub_dimension_Matrix type[PrimaryD];
  type data;

  template<typename... RestI>
  T& at(std::size_t index, RestI... rest) {
    if(index >= PrimaryD) {
      throw std::out_of_range("FixMatrix");
    }
    return data[index].at(rest...);
  }

  static size_t dimension() {
    return 1 + sizeof...(RestD);
  }
  
};
#endif //__FIXMATRIX_HPP__
