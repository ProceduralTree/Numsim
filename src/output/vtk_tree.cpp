#include "vtk_tree.h"
#include "grid/densetree.h"
#include "grid/sparsegrid.h"
#include "grid/zindex.h"
#include "utils/profiler.h"
#include <cstddef>
#include <cstdint>
#include <output/vtk_util.h>
#include <string>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkXMLImageDataWriter.h>

vtkSmartPointer<vtkImageData> init(const DenseTree::DenseTree& tree, bool interpolate)
{
  vtkSmartPointer<vtkImageData> dataSet = vtkSmartPointer<vtkImageData>::New();
  dataSet->SetOrigin(0, 0, 0);

  size_t resolution = (1 << tree.maxDepth);
  // set spacing of mesh
  const double dx
    = 1. / resolution;
  const double dy = 1. / resolution;
  const double dz = 1;
  if (interpolate)
    resolution -= 2;

  dataSet->SetSpacing(dx, dy, dz);

  // set number of points in each dimension, 1 cell in z direction
  dataSet->SetDimensions(
    resolution, resolution, 1);

  return dataSet;
};

void save_dataset(vtkSmartPointer<vtkImageData> dataSet)
{
  ProfileScope("Save to file");
  static auto vtkWriter_ = vtkSmartPointer<vtkXMLImageDataWriter>::New();
  static int fileNumber = 0;
  set_filename(vtkWriter_, fileNumber);
  dataSet->Squeeze();
  vtkWriter_->SetInputData(dataSet);
  vtkWriter_->SetDataModeToBinary(); // set file mode to binary files:
  vtkWriter_->Write();
};
