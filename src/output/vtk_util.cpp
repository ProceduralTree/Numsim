#include "vtk_util.h"
#include <filesystem>

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
