#ifndef VTK_H_
#define VTK_H_

#include "vtkDoubleArray.h"
#include "vtkImageData.h"
#include "vtkPointData.h"

#include "system.h"

#include "grid.h"

void write_vtk(const PDESystem &system, double time);

#endif // VTK_H_
