
---

# TheMeshProject

**TheMeshProject** is the engineering codebase that follows the educational arc
of *Before the Mesh*. It implements practical geometry, meshing, and solver
components using the `btm-framework` as a foundation.

This repository is intentionally minimal at the start. It grows progressively as
the tutorials introduce new concepts and tools.

---

## Repository Structure

themeshproject/
include/            # Public headers for meshcore
src/                # Implementation files
apps/               # Future standalone tools
external/           # Contains the btm-framework submodule
CMakeLists.txt

---

## Framework Dependency

The project depends on the framework via:

external/btm-framework/

Initialize or update the submodule:

## bash
git submodule update --init --recursive
