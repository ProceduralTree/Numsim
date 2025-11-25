#include "vtk.h"
#include "pde/system.h"
#include "utils/index.h"
#include "utils/profiler.h"
#include <stdio.h>
#include <vtkImageData.h>

static auto vtkWriter_ = vtkSmartPointer<vtkXMLImageDataWriter>::New();
void checkForDir(const std::string& name)
{
  if (std::filesystem::is_directory(name))
    return;
  std::filesystem::create_directory(name);
}
void set_filename(const vtkSmartPointer<vtkXMLImageDataWriter> writer, int& fileNumber)
{

  // Assemble the filename
  std::stringstream fileName;
  fileName << "out/output_" << std::setw(4) << setfill('0') << fileNumber << "."
           << writer->GetDefaultFileExtension();
  checkForDir("out");
  // increment file no.
  // assign the new file name to the output vtkWriter_
  writer->SetFileName(fileName.str().c_str());
  fileNumber++;
}

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
    system.settings.nCells[0] + 1, system.settings.nCells[1] + 1, 1); // we want to have points at each corner of each cell

  return dataSet;
}

double write_pressure(vtkSmartPointer<vtkDoubleArray> arrayPressure, const PDESystem& system)
{

  double index = 0; // index for the vtk data structure, will be incremented
                    // in the inner loop
  Index I;
  for (int j = system.begin.y - 1; j <= system.end.y; j++)
  {
    for (int i = system.begin.x - 1; i <= system.end.x; i++, index++)
    {
      I = { static_cast<uint16_t>(i), static_cast<uint16_t>(j) };
      arrayPressure->SetValue(index, interpolate_p(system, system.p, I));
    }
  }
  return index;
}

void write_vtk(const PDESystem& system, double time)
{
  Scope scope("Output VTK");
  static int fileNumber = 0;
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
  Index I;
  for (int j = system.begin.y - 1; j <= system.end.y; j++)
  {
    for (int i = system.begin.x - 1; i <= system.end.x; i++, index++)
    {
      I = { static_cast<uint16_t>(i), static_cast<uint16_t>(j) };
      std::array<double, 3> velocityVector;
      velocityVector[0] = interpolate_u(system, system.u, I);
      velocityVector[1] = interpolate_v(system, system.v, I);
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
