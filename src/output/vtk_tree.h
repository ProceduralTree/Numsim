#ifndef VTK_TREE_H_
#define VTK_TREE_H_

#include "grid/densetree.h"
#include "grid/sparsegrid.h"
#include "pde/system.h"
#include "utils/index.h"
#include <cstddef>
#include <cstdint>
#include <vector>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkXMLImageDataWriter.h>
constexpr void set_index(vtkSmartPointer<vtkDoubleArray> array, size_t idx, Index I, const DenseTree::DenseTree& tree)
{
  size_t index = DenseTree::get_dense_index(tree, I);
  array->SetValue(idx, static_cast<double>(index));
};

template <typename T>
constexpr void set_value(vtkSmartPointer<vtkDoubleArray> array, size_t idx, Index I, const SparseGrid2D<T>& data)
{
  array->SetValue(idx, static_cast<double>(data[I]));
};

constexpr void set_pressure(vtkSmartPointer<vtkDoubleArray> array, size_t idx, Index I, const PDESystem& sys)
{
  array->SetValue(idx, interpolate_p(sys, sys.p, I));
};

constexpr void set_velocity(vtkSmartPointer<vtkDoubleArray> array, size_t idx, Index I, const PDESystem& system)
{
  std::array<double, 3> velocityVector;
  velocityVector[0] = interpolate_u(system, system.u, I);
  velocityVector[1] = interpolate_v(system, system.v, I);
  velocityVector[2] = 0.0;
  array->SetTuple(idx, velocityVector.data());
};

template <typename Operator, typename... Args>
constexpr void write_data(Operator&& O, std::string name, int num_components, const DenseTree::DenseTree& tree, vtkSmartPointer<vtkImageData> dataSet, bool interpolate, Args&&... args)
{
  vtkSmartPointer<vtkDoubleArray> array = vtkDoubleArray::New();
  array->SetNumberOfTuples(dataSet->GetNumberOfPoints());
  array->SetName(name.c_str());
  array->SetNumberOfComponents(num_components);
  array->SetNumberOfTuples(dataSet->GetNumberOfPoints());
  size_t idx = 0;
  for (uint16_t j = interpolate ? 1 : 0; j < (1ULL << tree.maxDepth) - interpolate ? 1 : 0; j++)
    for (uint16_t i = interpolate ? 1 : 0; i < (1ULL << tree.maxDepth) - interpolate ? 1 : 0; i++)
    {
      {

        Index I = { i, j, tree.maxDepth };
        std::forward<Operator>(O)(array, idx++, I, std::forward<Args>(args)...);
      }
    }
  dataSet->GetPointData()->AddArray(array);
};

constexpr void write_depth(std::string name, const DenseTree::DenseTree& tree, vtkSmartPointer<vtkImageData> dataSet)
{
  ProfileScope("Write Field");
  vtkSmartPointer<vtkIntArray> array = vtkIntArray::New();
  array->SetNumberOfTuples(dataSet->GetNumberOfPoints());
  array->SetName(name.c_str());
  array->SetNumberOfComponents(1);
  array->SetNumberOfTuples(dataSet->GetNumberOfPoints());

  // double* ptr = static_cast<double*>(array->GetVoidPointer(0));
  size_t idx = 0;
  for (uint16_t j = 0; j < (1ULL << tree.maxDepth); j++)
    for (uint16_t i = 0; i < (1ULL << tree.maxDepth); i++)
    {
      {

        Index I = { i, j, tree.maxDepth };
        Zindex Z = Zindex(I);
        size_t index = DenseTree::get_dense_index(tree, Z);

        array->SetValue(idx++, tree._index_cache[index].depth);
      }
    }
  dataSet->GetPointData()->AddArray(array);
};

template <typename T>
void write_field(std::string name, const DenseTree::DenseTree& tree, const std::vector<T>& data, uint16_t level, vtkSmartPointer<vtkImageData> dataSet)
{
  ProfileScope("Write Field");
  vtkSmartPointer<vtkDoubleArray> array = vtkDoubleArray::New();
  array->SetNumberOfTuples(dataSet->GetNumberOfPoints());
  array->SetName(name.c_str());
  array->SetNumberOfComponents(1);
  array->SetNumberOfTuples(dataSet->GetNumberOfPoints());

  // double* ptr = static_cast<double*>(array->GetVoidPointer(0));
  size_t idx = 0;
  for (uint16_t j = 0; j < (1ULL << tree.maxDepth); j++)
    for (uint16_t i = 0; i < (1ULL << tree.maxDepth); i++)
    {
      {

        Index I = { i, j, level };
        Zindex Z = Zindex(I);
        size_t index = DenseTree::get_dense_index(tree, Z);

        array->SetValue(idx++, static_cast<double>(data[index]));
      }
    }
  dataSet->GetPointData()->AddArray(array);
};
constexpr void write(std::string name, const PDESystem& system, vtkSmartPointer<vtkImageData> dataSet)
{
  write_data(set_pressure, "Pressure", 1, system.boundary.tree, dataSet, true, system);
  write_data(set_velocity, "Velocity", 2, system.boundary.tree, dataSet, true, system);
};
template <typename T>
void write_field(std::string name, const DenseTree::DenseTree& tree, const std::vector<T>& data, vtkSmartPointer<vtkImageData> dataSet)
{
  write_field<T>(name, tree, data, tree.maxDepth, dataSet);
};
void save_dataset(vtkSmartPointer<vtkImageData> dataSet);
vtkSmartPointer<vtkImageData> init(const DenseTree::DenseTree& tree, bool interpolate);
#endif // VTK_TREE_H_
