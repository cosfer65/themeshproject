
---

# TheMeshProject

**TheMeshProject** is the engineering codebase that follows the educational arc
of *Before the Mesh*. It implements practical geometry, meshing, and solver
components using the `btm-framework` as a foundation.

This repository is intentionally minimal at the start. It grows progressively as
the tutorials introduce new concepts and tools.

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

The project is currently based on version v0.3.3 of btm-framework
