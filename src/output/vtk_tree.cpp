#include "vtk_tree.h"
#include "grid/zindex.h"
#include <cstddef>
#include <cstdint>
#include <output/vtk_util.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkXMLImageDataWriter.h>

vtkSmartPointer<vtkImageData> init(DenseTree::DenseTree tree)
{
  vtkSmartPointer<vtkImageData> dataSet = vtkSmartPointer<vtkImageData>::New();
  dataSet->SetOrigin(0, 0, 0);

  size_t resolution = (1 << tree.maxDepth);
  // set spacing of mesh
  const double dx
    = 1. / resolution;
  const double dy = 1. / resolution;
  const double dz = 1;
  dataSet->SetSpacing(dx, dy, dz);

  // set number of points in each dimension, 1 cell in z direction
  dataSet->SetDimensions(
    resolution, resolution, 1);

  return dataSet;
};

void save_dataset(vtkSmartPointer<vtkImageData> dataSet)
{
  static auto vtkWriter_ = vtkSmartPointer<vtkXMLImageDataWriter>::New();
  static int fileNumber = 0;
  set_filename(vtkWriter_, fileNumber);
  dataSet->Squeeze();
  vtkWriter_->SetInputData(dataSet);
  vtkWriter_->SetDataModeToBinary(); // set file mode to binary files:
  vtkWriter_->Write();
};

void write_tree(DenseTree::DenseTree tree)
{
  auto dataSet = init(tree);
  vtkSmartPointer<vtkDoubleArray> array = vtkDoubleArray::New();
  array->SetNumberOfTuples(dataSet->GetNumberOfPoints());
  array->SetName("Indices");
  array->SetNumberOfComponents(1);
  array->SetNumberOfTuples(dataSet->GetNumberOfPoints());

  size_t idx = 0;
  for (uint16_t j = 0; j < (1ULL << tree.maxDepth); j++)
    for (uint16_t i = 0; i < (1ULL << tree.maxDepth); i++)
    {
      {

        Index I = { i, j, tree.maxDepth };
        Zindex Z = Zindex(I);
        size_t index = DenseTree::get_dense_index(tree, Z);
        array->SetValue(idx++, static_cast<double>(index));
      }
    }
  dataSet->GetPointData()->AddArray(array);
  save_dataset(dataSet);
}
