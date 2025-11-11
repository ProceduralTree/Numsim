#ifndef PRESSURESOLVERS_H_
#define PRESSURESOLVERS_H_

#include <pde/system.h>
#include <utils/index.h>

void sor_step(PDESystem& system, Index I);
void gauss_seidel_step(PDESystem& system, Index I);
#endif // PRESSURESOLVERS_H_
