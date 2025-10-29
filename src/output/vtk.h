#ifndef VTK_H_
#define VTK_H_

#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkXMLImageDataWriter.h>

#include <pde/system.h>

#include <grid/grid.h>

struct ParaviewOutput
{
  vtkSmartPointer<vtkXMLImageDataWriter> vtkWriter_; //< vtk writer to write ImageData
  int fileNo_;
};

void write_vtk(const PDESystem& system, double time);

#endif // VTK_H_
