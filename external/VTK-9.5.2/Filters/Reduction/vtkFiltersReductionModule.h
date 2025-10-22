
#ifndef VTKFILTERSREDUCTION_EXPORT_H
#define VTKFILTERSREDUCTION_EXPORT_H

#ifdef VTKFILTERSREDUCTION_STATIC_DEFINE
#  define VTKFILTERSREDUCTION_EXPORT
#  define VTKFILTERSREDUCTION_NO_EXPORT
#else
#  ifndef VTKFILTERSREDUCTION_EXPORT
#    ifdef FiltersReduction_EXPORTS
        /* We are building this library */
#      define VTKFILTERSREDUCTION_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define VTKFILTERSREDUCTION_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef VTKFILTERSREDUCTION_NO_EXPORT
#    define VTKFILTERSREDUCTION_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef VTKFILTERSREDUCTION_DEPRECATED
#  define VTKFILTERSREDUCTION_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef VTKFILTERSREDUCTION_DEPRECATED_EXPORT
#  define VTKFILTERSREDUCTION_DEPRECATED_EXPORT VTKFILTERSREDUCTION_EXPORT VTKFILTERSREDUCTION_DEPRECATED
#endif

#ifndef VTKFILTERSREDUCTION_DEPRECATED_NO_EXPORT
#  define VTKFILTERSREDUCTION_DEPRECATED_NO_EXPORT VTKFILTERSREDUCTION_NO_EXPORT VTKFILTERSREDUCTION_DEPRECATED
#endif

/* NOLINTNEXTLINE(readability-avoid-unconditional-preprocessor-if) */
#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef VTKFILTERSREDUCTION_NO_DEPRECATED
#    define VTKFILTERSREDUCTION_NO_DEPRECATED
#  endif
#endif

/* VTK-HeaderTest-Exclude: vtkFiltersReductionModule.h */

/* Include ABI Namespace */
#include "vtkABINamespace.h"

#endif /* VTKFILTERSREDUCTION_EXPORT_H */
