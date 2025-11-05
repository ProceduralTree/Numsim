#pragma once
struct PDESystem;

namespace Writer {
void writeTextVelocity(const PDESystem& system, double currentTime);

void writeTextPressure(const PDESystem& system);
}
