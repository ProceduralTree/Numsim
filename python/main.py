import pandas as pd
import sys
import os
import subprocess

setting_template = """
# Settings file for numsim program
# Run ./numsim lid_driven_cavity.txt

# Problem description
physicalSizeX = 2.0   # physical size of the domain
physicalSizeY = 2.0
endTime = 10.0        # duration of the simulation
re = 100             # Reynolds number
gX = 0.0              # external forces, set to (gX,gY) = (0,-9.81) to account for gravity
gY = 0.0

# Dirichlet boundary conditions
dirichletBottomX = 0
dirichletBottomY = 0
dirichletTopX    = 1
dirichletTopY    = 0
dirichletLeftX   = 0
dirichletLeftY   = 0
dirichletRightX  = 0
dirichletRightY  = 0

# Discretization parameters
nCellsX = {}          # number of cells in x and y direction
nCellsY = {} 
useDonorCell = true   # if donor cell discretization should be used, possible values: true false
alpha = 0.5           # factor for donor-cell scheme, 0 is equivalent to central differences
tau = 0.5             # safety factor for time step width
maximumDt = 0.1       # maximum values for time step width

# Solver parameters
pressureSolver = SOR  # which pressure solver to use, possible values: GaussSeidel SOR CG
omega = 1.6           # overrelaxation factor, only for SOR solver
epsilon = 1e-5        # tolerance for 2-norm of residual
maximumNumberOfIterations = 1e4    # maximum number of iterations in the solver

"""



if __name__ == "__main__":
    for n in [2**i for i in range(4,8)]:

        with open("/tmp/settings.txt" , "w") as file:
            file.write(setting_template.format(n,n))

        subprocess.run(["../build/numsim" , "/tmp/settings.txt"] , capture_output=True)


        df = pd.read_csv("Profiler.csv")
        df.columns = ["name" , "t-begin" , "t-end" , "calltime[ms]" , "#calls"]
        df["calltime[ms]"] = (df["calltime[ms]"].str.replace("ns", "", regex=False).astype(float) * 1e-6)

        df["#calls"] = df["#calls"] + 1
        group = df.groupby(by="name").agg({"calltime[ms]" : 'mean' , "#calls" : 'sum'})

        print(f"\nGridsize: {n}x{n}\n")
        print(group)
        print("\n\n")
