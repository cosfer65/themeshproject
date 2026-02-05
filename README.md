# TheMeshProject — C++ Mesh Viewer

TheMeshProject is a lightweight, open-source C++ application for loading, visualizing, and analyzing 3D triangle meshes. It provides interactive inspection tools along with robust differential‑geometry operators for computing curvature quantities. The project is designed to be easy to download, compile, and extend, making it useful for students, researchers, and developers working with surface geometry.

---

## Features

### 🔹 Mesh Visualization
- Load and display common triangle‑mesh formats  
- Interactive camera controls (rotate, pan, zoom)  
- Smooth shading and normal visualization  

### 🔹 Curvature Computation
The application computes per‑vertex:
- **Gaussian curvature**  
- **Mean curvature**  
- **Principal curvatures** (minimum and maximum)  
- **Principal directions**  

Curvature values can be visualized using color maps or directional fields.

### 🔹 Lightweight & Portable
- Written in standard C++  
- Minimal external dependencies  
- Fully compatible with Visual Studio  
- Clear, modular code structure for easy extension  

---

## Project Structure

Each module is self‑contained, making it straightforward to modify or replace components.

---

## Building the Project

### Requirements
- Visual Studio (2026 Community Edition was used)
- A C++20‑compatible compiler
- OpenGL and GLFW/GLUT (depending on your configuration)

### Steps
1. Clone the repository  
2. Open the provided Visual Studio solution (`.sln`)  
3. Build the project in Debug or Release mode  
4. Run the executable and load a mesh file from the UI  

No additional setup is required.

---

## Sample Meshes
The `runtime/resources/models/` folder includes a few example meshes for testing curvature visualization. You can also load your own `.obj` or other supported formats.

---

## Purpose & Vision
TheMeshProject aims to provide a clear, accessible reference implementation for mesh visualization and curvature analysis in C++. It is both a practical tool and a learning resource for anyone exploring geometry processing, differential operators, or mesh‑based computation.

Contributions, suggestions, and extensions are welcome.

---

## License
This project is released under the MIT License. See `LICENSE` for details.
