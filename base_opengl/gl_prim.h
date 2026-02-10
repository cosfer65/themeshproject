#ifndef __primitives__
#define __primitives__

#include "vector.h"
#include "gl_mesh.h"
#include "gl_material.h"
#include "gl_shaders.h"

using namespace base_math;

namespace base_opengl {
	/**
	 * @brief Base class for OpenGL drawable primitives.
	 *
	 * Encapsulates mesh data, transformation matrices, material, and rendering state.
	 * Provides methods for setting transformation, draw mode/type, and rendering.
	 */
	class gl_prim {
	protected:
		mesh_data m_mesh_data;   ///< Mesh data for the primitive.
		GLuint vao;              ///< Vertex Array Object handle.

		/**
		 * @brief Specifies how polygons will be rasterized.
		 * Accepted values: GL_POINT, GL_LINE, GL_FILL.
		 */
		GLenum draw_mode;

		/**
		 * @brief Specifies the OpenGL primitive type for drawing.
		 * Common values: GL_TRIANGLES, GL_LINES, GL_PATCHES, etc.
		 */
		GLenum draw_type;

		/**
		 * @brief If true, uses glDrawElements; otherwise, uses glDrawArrays.
		 */
		bool draw_elements;

		fvec3 position;           ///< Object position in world space.
		fvec3 scale;              ///< Object scale factors.
		fvec3 rotation;           ///< Object rotation angles (radians).

		/**
		 * @brief If nonzero, enables per-vertex color in the shader.
		 */
		int use_vertex_color;

        fvec3 m_color;           ///< Base color of the primitive.
		cg_material* m_material; ///< Material for shading the primitive.
		GLuint m_texture;        ///< Texture handle for the primitive.

	public:
		fmat4 rmat;  ///< Local rotation matrix.
		fmat4 tmat;  ///< Translation matrix.
		fmat4 smat;  ///< Scaling matrix.
		bool force_black;        ///< If true, forces the primitive to render in black (e.g., for wireframe).
        fmat4 view_matrix; ///< View matrix for the primitive (optional, can be set externally).

		/**
		 * @brief Constructs a gl_prim with default transformation and rendering state.
		 */
		gl_prim() {
			vao = 0;
			scale = fvec3(1);
			position = fvec3(0);
			rotation = fvec3(0);
			draw_mode = GL_FILL;
			draw_type = GL_TRIANGLES;
			draw_elements = true;
			// matrices are row-major!
			rmat = rotation_matrix<float>(rotation.x(), 0, 0, 1) * rotation_matrix<float>(rotation.y(), 0, 1, 0) * rotation_matrix<float>(rotation.z(), 1, 0, 0);
			smat = scaling_matrix(scale.x(), scale.y(), scale.z());
			tmat = translation_matrix(position.x(), position.y(), position.z());
			use_vertex_color = 0;
			m_texture = 0;
            m_material = nullptr;
            m_color = fvec3(0.8f, 0.8f, 0.8f);
            force_black = false;
            view_matrix.loadIdentity();
		}

		int get_use_vertex_color() const {
			return use_vertex_color;
        }

		int set_use_vertex_color(int v) {
            int ret = use_vertex_color;
			use_vertex_color = v;
			return ret;
        }

		/**
		 * @brief Virtual destructor.
		 */
		virtual ~gl_prim() {}

		/**
		 * @brief Initializes the primitive from mesh data.
		 * @param mesh Pointer to mesh data.
		 * @param drmode Polygon rasterization mode (GL_FILL, GL_LINE, etc.).
		 * @param dr_el If true, use glDrawElements; otherwise, use glDrawArrays.
		 */
		virtual void create_from_mesh(mesh_data* mesh, GLenum drmode = GL_FILL, bool dr_el = true);

		/**
		 * @brief Sets the polygon rasterization mode.
		 * @param dm OpenGL polygon mode (GL_FILL, GL_LINE, GL_POINT).
		 */
		void set_draw_mode(GLenum dm) {
			draw_mode = dm;
		}

		/**
		 * @brief Sets the texture for the primitive.
		 * @param tex OpenGL texture handle.
		 */
		void set_texture(GLuint tex) {
			m_texture = tex;
		}

		/**
		 * @brief Sets the material for the primitive.
		 * @param m Pointer to cg_material.
		 */
		void set_material(cg_material* m) {
			m_material = m;
		}

		void set_color(const fvec3& col) {
			m_color = col;
        }

		/**
		 * @brief Gets the current material.
		 * @return Pointer to cg_material.
		 */
		const cg_material* material() {
			return m_material;
		}

		/**
		 * @brief Renders the primitive using the provided shader.
		 * @param _shader Pointer to the shader program.
		 */
		virtual void render(gl_shader* _shader) {
			if (!vao) return;

            _shader->set_int("object_or_vertex_color", 0);  // object color by default
			if (force_black)
				_shader->set_vec4("object_color", fvec4(0, 0, 0, 1)); // black for wireframe
            else
				_shader->set_vec4("object_color", fvec4(m_color, 1));
			if (use_vertex_color) {
                _shader->set_int("object_or_vertex_color", 1); // use per-vertex color
			}
			else if (m_material) {
				m_material->apply(_shader);
				_shader->set_vec4("object_color", fvec4(m_material->get_diffuse(), 1));
			}

			// position object
			fmat4 ob_matrix = tmat * rmat * smat;
			ob_matrix = ob_matrix.transpose();    // convert to column wise for OpenGL!
			ob_matrix = view_matrix * ob_matrix;

			// pass transformation to shader
			_shader->set_mat4("model", ob_matrix);

			glBindVertexArray(vao);
			if (draw_elements)
			{
				// setup drawing
				if (draw_type == GL_TRIANGLES) {
					glFrontFace(GL_CCW);
					glPolygonMode(GL_FRONT_AND_BACK, draw_mode);
				}
				glDrawElements(draw_type, (unsigned int)m_mesh_data.num_indices, GL_UNSIGNED_INT, 0);
			}
			else
			{
				// setup drawing
				glPatchParameteri(GL_PATCH_VERTICES, 4);
				glPolygonMode(GL_FRONT_AND_BACK, draw_mode);
				glDrawArrays(GL_PATCHES, 0, (unsigned int)m_mesh_data.num_indices);
			}
			glBindVertexArray(0);
		}

		/**
		 * @brief Sets the OpenGL primitive type for drawing.
		 * @param dt OpenGL draw type (GL_TRIANGLES, GL_LINES, etc.).
		 */
		void set_draw_type(GLenum dt) {
			draw_type = dt;
		}

		// virtual void step_simulation(float fElapsed) {}

		/**
		 * @brief Sets position, scale, and rotation in one call.
		 * @param _p Position vector.
		 * @param _s Scale vector.
		 * @param _r Rotation vector (radians).
		 */
		void set_all(const fvec3& _p, const fvec3& _s, const fvec3& _r) {
			position = _p;
			scale = _s;
			rotation = _r;
			tmat = translation_matrix(position.x(), position.y(), position.z());
			smat = scaling_matrix(scale.x(), scale.y(), scale.z());
			rmat = rotation_matrix<float>(rotation.x(), 1, 0, 0) * rotation_matrix<float>(rotation.y(), 0, 1, 0) * rotation_matrix<float>(rotation.z(), 0, 0, 1);
		}

		/**
		 * @brief Sets the rotation to the specified vector.
		 * @param _r Rotation vector (radians).
		 */
		void rotate_to(const fvec3& _r) {
			rotation = _r;
			rmat = rotation_matrix<float>(rotation.x(), 1, 0, 0) * rotation_matrix<float>(rotation.y(), 0, 1, 0) * rotation_matrix<float>(rotation.z(), 0, 0, 1);
		}

		/**
		 * @brief Sets the rotation to the specified angles.
		 * @param x Rotation around X axis (radians).
		 * @param y Rotation around Y axis (radians).
		 * @param z Rotation around Z axis (radians).
		 */
		void rotate_to(float x, float y, float z) {
			rotation = fvec3(x, y, z);
			rmat = rotation_matrix<float>(rotation.x(), 1, 0, 0) * rotation_matrix<float>(rotation.y(), 0, 1, 0) * rotation_matrix<float>(rotation.z(), 0, 0, 1);
		}

		/**
		 * @brief Adds the specified vector to the current rotation.
		 * @param _r Rotation vector to add (radians).
		 */
		void rotate_by(const fvec3& _r) {
			rotation += _r;
			rmat = rotation_matrix<float>(rotation.x(), 1, 0, 0) * rotation_matrix<float>(rotation.y(), 0, 1, 0) * rotation_matrix<float>(rotation.z(), 0, 0, 1);
		}

		/**
		 * @brief Adds the specified angles to the current rotation.
		 * @param x Rotation around X axis (radians).
		 * @param y Rotation around Y axis (radians).
		 * @param z Rotation around Z axis (radians).
		 */
		void rotate_by(float x, float y, float z) {
			rotation += fvec3(x, y, z);
			rmat = rotation_matrix<float>(rotation.x(), 1, 0, 0) * rotation_matrix<float>(rotation.y(), 0, 1, 0) * rotation_matrix<float>(rotation.z(), 0, 0, 1);
		}

		/**
		 * @brief Sets the position to the specified vector.
		 * @param _r Position vector.
		 */
		void move_to(const fvec3& _r) {
			position = _r;
			tmat = translation_matrix(position.x(), position.y(), position.z());
		}

		/**
		 * @brief Sets the position to the specified coordinates.
		 * @param x X coordinate.
		 * @param y Y coordinate.
		 * @param z Z coordinate.
		 */
		void move_to(float x, float y, float z) {
			position = fvec3(x, y, z);
			tmat = translation_matrix(x, y, z);
		}

		/**
		 * @brief Adds the specified vector to the current position.
		 * @param _r Position vector to add.
		 */
		void move_by(const fvec3& _r) {
			position += _r;
			tmat = translation_matrix(position.x(), position.y(), position.z());
		}

		/**
		 * @brief Adds the specified values to the current position.
		 * @param x X increment.
		 * @param y Y increment.
		 * @param z Z increment.
		 */
		void move_by(float x, float y, float z) {
			position += fvec3(x, y, z);
			tmat = translation_matrix(position.x(), position.y(), position.z());
		}

		/**
		 * @brief Virtual method for creating the primitive geometry.
		 * @param drmode Polygon rasterization mode.
		 * @param dr_el If true, use glDrawElements; otherwise, use glDrawArrays.
		 */
		virtual void create(GLenum drmode = GL_FILL, bool dr_el = true) {
		}

		/**
		 * @brief Sets the scale factors.
		 * @param x Scale along X axis.
		 * @param y Scale along Y axis.
		 * @param z Scale along Z axis.
		 */
		void set_scale(float x, float y, float z) {
			scale = fvec3(x, y, z);
			smat = scaling_matrix(x, y, z);
		}

		/**
		 * @brief Sets the scale vector.
		 * @param _s Scale vector.
		 */
		void set_scale(const fvec3& _s) {
			scale = _s;
			smat = scaling_matrix(scale.x(), scale.y(), scale.z());
		}

		/**
		 * @brief Sets the scale along the X axis.
		 * @param _s Scale value.
		 */
		void set_xscale(float _s) {
			scale.x() = _s;
			smat = scaling_matrix(scale.x(), scale.y(), scale.z());
		}

		/**
		 * @brief Sets the scale along the Y axis.
		 * @param _s Scale value.
		 */
		void set_yscale(float _s) {
			scale.y() = _s;
			smat = scaling_matrix(scale.x(), scale.y(), scale.z());
		}

		/**
		 * @brief Sets the scale along the Z axis.
		 * @param _s Scale value.
		 */
		void set_zscale(float _s) {
			scale.z() = _s;
			smat = scaling_matrix(scale.x(), scale.y(), scale.z());
		}

		/**
		 * @brief Gets the current position vector.
		 * @return Reference to position vector.
		 */
		fvec3& get_position() {
			return position;
		}

		/**
		 * @brief Gets the current scale vector.
		 * @return Reference to scale vector.
		 */
		fvec3& get_scale() {
			return scale;
		}

		/**
		 * @brief Gets the current rotation vector.
		 * @return Reference to rotation vector.
		 */
		fvec3& get_rotation() {
			return rotation;
		}
	};

	/**
	 * @brief Creates a gl_prim from a half_edge_mesh.
	 * @param ms Pointer to half_edge_mesh.
	 * @param drmode OpenGL draw mode (default: GL_LINE).
	 * @param dr_el Whether to use element drawing (default: true).
	 * @return Pointer to the created gl_prim.
	 */
	gl_prim* create_prim(half_edge_mesh<float>* ms, GLenum drmode=GL_LINE, bool dr_el=true);

	/**
	 * @brief Creates a gl_prim from a gl_mesh.
	 * @param ms Pointer to gl_mesh.
	 * @param drmode OpenGL draw mode (default: GL_LINE).
	 * @param dr_el Whether to use element drawing (default: true).
	 * @return Pointer to the created gl_prim.
	 */
	gl_prim* create_prim(gl_mesh* ms, GLenum drmode=GL_LINE, bool dr_el=true);

	/**
	 * @brief Creates a cone primitive.
	 * @param drmode Polygon rasterization mode.
	 * @param dr_el If true, use glDrawElements; otherwise, use glDrawArrays.
	 * @return Pointer to the created gl_prim.
	 */
	gl_prim* create_cone(GLenum drmode = GL_LINE, bool dr_el = true);

	/**
	 * @brief Creates a cube primitive.
	 * @param drmode Polygon rasterization mode.
	 * @param dr_el If true, use glDrawElements; otherwise, use glDrawArrays.
	 * @return Pointer to the created gl_prim.
	 */
	gl_prim* create_cube(GLenum drmode = GL_LINE, bool dr_el = true);

	/**
	 * @brief Creates a cylinder primitive.
	 * @param drmode Polygon rasterization mode.
	 * @param dr_el If true, use glDrawElements; otherwise, use glDrawArrays.
	 * @return Pointer to the created gl_prim.
	 */
	gl_prim* create_cylinder(GLenum drmode = GL_LINE, bool dr_el = true);

	/**
	 * @brief Creates a dodecahedron primitive.
	 * @param drmode Polygon rasterization mode.
	 * @param dr_el If true, use glDrawElements; otherwise, use glDrawArrays.
	 * @return Pointer to the created gl_prim.
	 */
	gl_prim* create_dodecahedron(GLenum drmode = GL_LINE, bool dr_el = true);

	/**
	 * @brief Creates an icosahedron primitive.
	 * @param drmode Polygon rasterization mode.
	 * @param dr_el If true, use glDrawElements; otherwise, use glDrawArrays.
	 * @return Pointer to the created gl_prim.
	 */
	gl_prim* create_icosahedron(GLenum drmode = GL_LINE, bool dr_el = true);

	/**
	 * @brief Creates an octahedron primitive.
	 * @param drmode Polygon rasterization mode.
	 * @param dr_el If true, use glDrawElements; otherwise, use glDrawArrays.
	 * @return Pointer to the created gl_prim.
	 */
	gl_prim* create_octa(GLenum drmode = GL_LINE, bool dr_el = true);

	/**
	 * @brief Creates a pentagonal primitive.
	 * @param drmode Polygon rasterization mode.
	 * @param dr_el If true, use glDrawElements; otherwise, use glDrawArrays.
	 * @return Pointer to the created gl_prim.
	 */
	gl_prim* create_penta(GLenum drmode = GL_LINE, bool dr_el = true);

	/**
	 * @brief Creates a plane primitive.
	 * @param drmode Polygon rasterization mode.
	 * @param dr_el If true, use glDrawElements; otherwise, use glDrawArrays.
	 * @return Pointer to the created gl_prim.
	 */
	gl_prim* create_plane(GLenum drmode = GL_LINE, bool dr_el = true);

	/**
	 * @brief Creates a sphere primitive.
	 * @param drmode Polygon rasterization mode.
	 * @param dr_el If true, use glDrawElements; otherwise, use glDrawArrays.
	 * @return Pointer to the created gl_prim.
	 */
	gl_prim* create_sphere(GLenum drmode = GL_LINE, bool dr_el = true);

	/**
	 * @brief Creates a tetrahedron primitive.
	 * @param drmode Polygon rasterization mode.
	 * @param dr_el If true, use glDrawElements; otherwise, use glDrawArrays.
	 * @return Pointer to the created gl_prim.
	 */
	gl_prim* create_tetra(GLenum drmode = GL_LINE, bool dr_el = true);

	/**
	 * @brief Creates a torus primitive.
	 * @param drmode Polygon rasterization mode.
	 * @param dr_el If true, use glDrawElements; otherwise, use glDrawArrays.
	 * @return Pointer to the created gl_prim.
	 */
	gl_prim* create_torus(GLenum drmode = GL_LINE, bool dr_el = true);

	/**
	 * @brief Creates a Universal Coordinate System (XYZ axes with arrows).
	 * @param drmode Polygon rasterization mode.
	 * @param dr_el If true, use glDrawElements; otherwise, use glDrawArrays.
	 * @return Pointer to the created gl_prim.
	 */
	gl_prim* create_UCS(GLenum drmode = GL_LINE, bool dr_el = true);
}

#endif // __primitives__
