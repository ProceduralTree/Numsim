// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause

#ifndef vtkCxxABIConfigure_h
#define vtkCxxABIConfigure_h

#define VTK_HAS_CXXABI_DEMANGLE

#ifdef VTK_HAS_CXXABI_DEMANGLE
#  include <cxxabi.h>
#endif

#endif //vtkCxxABIConfigure_h
