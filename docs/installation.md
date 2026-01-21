---
parent: CaptureFOAM v1.0 User Guide
nav_order: 1
---

# Installation

---

## Dependencies
CaptureFOAM is built upon [OpenFOAM](https://openfoam.org/), a free open-source computational fluid dynamics software. CaptureFOAM is currently compatible with [OpenFOAM v12](https://github.com/OpenFOAM/OpenFOAM-12), which will need to be compiled and sourced prior to beginning the CaptureFOAM installation process. Installation instructions for OpenFOAM can be found [here](https://openfoam.org/download/source/).

---

## Compiling CaptureFOAM

Once OpenFOAM v12 is built, CaptureFOAM is ready to be downloaded and compiled via the following steps:

1. Clone CaptureFOAM from the repository to your OpenFOAM project directory `WM_PROJECT_USER_DIR`:
   ```bash
   cd $WM_PROJECT_USER_DIR
   git clone https://github.com/ORNL/AdditiveFOAM.git
   ```

2. Build the CaptureFOAM libraries:
   ```bash
   cd $WM_PROJECT_USER_DIR/CaptureFOAM
   ./Allwmake
   ```

If the `Allwmake` script completes without errors, CaptureFOAM is ready to use.
