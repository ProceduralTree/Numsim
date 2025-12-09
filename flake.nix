{
  description = "Numsim C++ project with VTK + MPI built via CMake";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.05";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      ...
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs { inherit system; };
      in
      {

        packages.default = pkgs.stdenv.mkDerivation {
          pname = "numsim";
          version = "0.1";
          src = self;

          # Tools needed at build time
          nativeBuildInputs = [
            pkgs.ninja
            pkgs.cmake
            pkgs.gcc
          ];

          # Libraries available at runtime
          buildInputs = [
            pkgs.vtk
            pkgs.openmpi
          ];

          CXX = "${pkgs.gcc}/bin/g++";

          cmakeFlags = [
            "-DCMAKE_CXX_STANDARD=23"
          ];

          # match your original flags
          NIX_ENFORCE_NO_NATIVE = 0;
          NIX_CFLAGS_COMPILE = [
            "-march=znver4"
            "-mbmi2"
            "-Wall"
            "-pedantic"
            "-O3"
            "-fopenmp"
            "-funroll-loops"
            "-ftree-vectorize"
          ];
        };

        devShells.default = pkgs.mkShell {
          nativeBuildInputs = [
            pkgs.ninja
            pkgs.cmake
            pkgs.gcc
            pkgs.openmpi
          ];

          buildInputs = [
            pkgs.vtk
          ];

          shellHook = ''
            # Unset NIX_ENFORCE_NO_NATIVE so -march=native works
            unset NIX_ENFORCE_NO_NATIVE
            echo "[Numsim dev shell]"
          '';
        };
      }
    );
}
