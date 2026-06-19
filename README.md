
---

# TheMeshProject

**TheMeshProject** is the engineering codebase that follows the educational arc
of *Before the Mesh*. It implements practical geometry, meshing, and solver
components using the `btm-framework` as a foundation.

---

## Repository Structure

LEGAL/              # legal and licensing information  
projects/           # contains individual tutorial project folders  
external/           # Contains the btm-framework submodule  
runtime/            # contains sample models, shaders and other runtime required files  
CMakeLists.txt  

---

## clone the repo

git clone https://github.com/cosfer65/themeshproject.git  
cd themeshproject  
git submodule update --init --recursive

## Framework Dependency

The project depends on the framework via:

external/btm-framework/

The project is currently based on version v0.3.5 of btm-framework  
to checkout this version of the framework:  
cd external/btm-framework/  
git checkout v0.3.5


## Building the project

mkdir build  
cd build  
cmake ..  
cmake --build .  
