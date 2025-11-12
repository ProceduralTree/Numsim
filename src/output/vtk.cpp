#include "vtk.h"
#include "utils/index.h"
#include <stdio.h>
#include <vtkImageData.h>

void set_filename(const vtkSmartPointer<vtkXMLImageDataWriter> writer, int& fileNumber)
{

  // Assemble the filename
  std::stringstream fileName;
  fileName << "out/output_" << std::setw(4) << setfill('0') << fileNumber << "."
           << writer->GetDefaultFileExtension();
  // increment file no.
  // assign the new file name to the output vtkWriter_
  writer->SetFileName(fileName.str().c_str());
  fileNumber++;
};

vtkSmartPointer<vtkImageData> initialize_dataset(const PDESystem& system)
{
  vtkSmartPointer<vtkImageData> dataSet = vtkSmartPointer<vtkImageData>::New();
  dataSet->SetOrigin(0, 0, 0);

  // set spacing of mesh
  const double dx = system.h.x;
  const double dy = system.h.y;
  const double dz = 1;
  dataSet->SetSpacing(dx, dy, dz);

  // set number of points in each dimension, 1 cell in z direction
  dataSet->SetDimensions(
    system.size_x, system.size_y, 1); // we want to have points at each corner of each cell

  return dataSet;
};

double write_pressure(vtkSmartPointer<vtkDoubleArray> arrayPressure, const PDESystem& system)
{

  double index = 0; // index for the vtk data structure, will be incremented
                    // in the inner loop
  for (int j = system.begin.y; j <= system.end.y; j++)
  {
    for (int i = system.begin.x; i <= system.end.x; i++, index++)
    {
      arrayPressure->SetValue(index, system.p[i, j]);
    }
  }
  return index;
};

void write_vtk(const PDESystem& system, double time)
{
  static int fileNumber = 0;
  static auto vtkWriter_ = vtkSmartPointer<vtkXMLImageDataWriter>::New();
  set_filename(vtkWriter_, fileNumber);
  auto dataSet = initialize_dataset(system);

  // initialize data set that will be output to the file
  dataSet->SetOrigin(0, 0, 0);

  // set spacing of mesh
  //
  vtkSmartPointer<vtkDoubleArray> arrayPressure = vtkDoubleArray::New();
  // the pressure is a scalar which means the number of components is 1
  arrayPressure->SetNumberOfComponents(1);
  // Set the number of pressure values and allocate memory for it. We
  // already know the number, it has to be the same as there are nodes in
  // the mesh.
  arrayPressure->SetNumberOfTuples(dataSet->GetNumberOfPoints());
  arrayPressure->SetName("pressure");
  double index = write_pressure(arrayPressure, system);
  dataSet->GetPointData()->AddArray(arrayPressure);
  assert(index == dataSet->GetNumberOfPoints());

  // loop over the nodes of the mesh and assign the interpolated p values
  // in the vtk data structure we only consider the cells that are the
  // actual computational domain, not the helper values in the "halo"

  // now, we should have added as many values as there are points in the
  // vtk data structure

  // add the field variable to the data set

  // add velocity field variable
  // ---------------------------
  vtkSmartPointer<vtkDoubleArray> arrayVelocity = vtkDoubleArray::New();

  // here we have two components (u,v), but ParaView will only allow
  // vector glyphs if we have an ℝ^3 vector, therefore we use a
  // 3-dimensional vector and set the 3rd component to zero
  arrayVelocity->SetNumberOfComponents(3);

  // set the number of values
  arrayVelocity->SetNumberOfTuples(dataSet->GetNumberOfPoints());

  arrayVelocity->SetName("velocity");

  // loop over the mesh where p is defined and assign the values in the
  // vtk data structure
  index = 0; // index for the vtk data structure
  Index I = { 0, 0 };
  Offset o = { 0, 0 };
  for (int j = system.begin.y; j <= system.end.y; j++)
  {
    for (int i = system.begin.x; i <= system.end.x; i++, index++)
    {
      I = { static_cast<uint16_t>(i), static_cast<uint16_t>(j) };
      o = { 1, 0 };
      std::array<double, 3> velocityVector;
      velocityVector[0] = interpolate_at(system, system.u, I, o);
      o = { 0, 1 };
      velocityVector[1] = interpolate_at(system, system.v, I, o);
      velocityVector[2] = 0.0;

      arrayVelocity->SetTuple(index, velocityVector.data());
    }
  }
  // now, we should have added as many values as there are points in the
  // vtk data structure
  assert(index == dataSet->GetNumberOfPoints());

  // add the field variable to the data set
  dataSet->GetPointData()->AddArray(arrayVelocity);

  // add current time
  vtkSmartPointer<vtkDoubleArray> arrayTime = vtkDoubleArray::New();
  arrayTime->SetName("TIME");
  arrayTime->SetNumberOfTuples(1);
  arrayTime->SetTuple1(0, time);
  dataSet->GetFieldData()->AddArray(arrayTime);

  // Remove unused memory
  dataSet->Squeeze();

  // Write the data
  vtkWriter_->SetInputData(dataSet);

  // vtkWriter_->SetDataModeToAscii();     // comment this in to get ascii
  // text files: those can be checked in an editor
  vtkWriter_->SetDataModeToBinary(); // set file mode to binary files:
                                     // smaller file sizes

  // finally write out the data
  vtkWriter_->Write();
}
