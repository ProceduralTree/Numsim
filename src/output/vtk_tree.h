#ifndef VTK_TREE_H_
#define VTK_TREE_H_

#include "grid/densetree.h"
#include "grid/sparsegrid.h"
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkXMLImageDataWriter.h>

template <typename T>
void write_field(std::string name, const SparseGrid2D<T>& data, vtkSmartPointer<vtkImageData> dataSet)
{
  ProfileScope("Write Field");
  vtkSmartPointer<vtkIntArray> array = vtkIntArray::New();
  array->SetNumberOfTuples(dataSet->GetNumberOfPoints());
  array->SetName(name.c_str());
  array->SetNumberOfComponents(1);
  array->SetNumberOfTuples(dataSet->GetNumberOfPoints());

  // double* ptr = static_cast<double*>(array->GetVoidPointer(0));
  size_t idx = 0;
  for (uint16_t j = 0; j < (1ULL << data.tree.maxDepth); j++)
    for (uint16_t i = 0; i < (1ULL << data.tree.maxDepth); i++)
    {
      {

        Index I = { i, j, data.tree.maxDepth };
        Zindex Z = Zindex(I);
        size_t index = DenseTree::get_dense_index(data.tree, Z);

        array->SetValue(idx++, data.tree._index_cache[index].depth);
      }
    }
  dataSet->GetPointData()->AddArray(array);
};
void save_dataset(vtkSmartPointer<vtkImageData> dataSet);
vtkSmartPointer<vtkImageData> init(const DenseTree::DenseTree& tree);
#endif // VTK_TREE_H_
