#ifndef VTK_UTIL_H_
#define VTK_UTIL_H_

#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkXMLImageDataWriter.h>
void checkForDir(const std::string& name);

void set_filename(const vtkSmartPointer<vtkXMLImageDataWriter> writer, int& fileNumber);
#endif // VTK_UTIL_H_
