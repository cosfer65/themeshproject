
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

## Framework Dependency

The project depends on the framework via:

external/btm-framework/

The project is currently based on version v0.3.2 of btm-framework

To update the submodule:

## bash
cd external/btm-framework

git fetch

git checkout v0.3.2
