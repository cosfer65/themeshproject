#define _CRT_SECURE_NO_WARNINGS
#include "gl_graphics.h"
#include "themeshview.h"
#include "themeshmodel.h"

#include "gl_camera.h"
#include "gl_shaders.h"
#include "gl_light.h"
#include "gl_prim.h"
#include "gl_font.h"

#include "visual_objects.h"
#include "arcball.h"

#include "resource.h"

using namespace base_math;
using namespace base_opengl;

struct mesh_view_private {
	std::unique_ptr<gl_camera> m_cam;                   ///< gl_camera used to view the scene and compute view/projection.
	std::unique_ptr<gl_shader> m_shader;                ///< Shader program used for mesh and helper rendering.
	std::unique_ptr<gl_light> m_light;                  ///< Primary scene light affecting shading.

	std::unique_ptr<UCS_view> m_ucs_view;               ///< UCS (user coordinate system) view used to render the axes widget.

	// {{ the loaded model and the acompanying visualization objects
	std::vector<std::unique_ptr<gl_prim>> m_draw_parts; ///< Drawable mesh parts converted to OpenGL primitives.
	// }}

	bool face_normals_created = false;                                  ///< Flag to track whether face normals visualization has been generated.
	std::vector<std::unique_ptr<visual_objects>> m_face_normals;        ///< Visualization of per-face normals as line segments.
	bool vertex_normals_created = false;                                ///< Flag to track whether vertex normals visualization has been generated.
	std::vector<std::unique_ptr<visual_objects>> m_vertex_normals;      ///< Visualization of per-vertex normals as line segments.
	bool curvature_directions_created = false;                          ///< Flag to track whether curvature directions visualization has been generated.
	std::vector<std::unique_ptr<visual_objects>> m_model_curvatures_k1; ///< Visualization of minimum principal curvature directions per vertex.
	std::vector<std::unique_ptr<visual_objects>> m_model_curvatures_k2; ///< Visualization of maximum principal curvature directions per vertex.

	void invalidate_visualizations() {
		face_normals_created = false;
		vertex_normals_created = false;
		curvature_directions_created = false;
	}

	std::unique_ptr<arcball> m_arcball;                                 ///< Arcball for mouse interaction
	gl_font* font2D;
};

meshViewWindow::meshViewWindow() : glViewWindow() {
	m_private = new mesh_view_private();
	m_private->m_ucs_view.reset(new UCS_view());
	m_private->m_arcball.reset(new arcball(800, 600));
	m_private->m_cam.reset(new gl_camera(fvec3(0, 0, 50), fvec3(0, 0, 0), fvec3(0, 1, 0)));
	m_private->m_cam->set_fov(dtr(30.f));
}

meshViewWindow::~meshViewWindow() {
	delete m_private;
}

void meshViewWindow::init() {
	m_private->m_ucs_view->initialize();

	m_private->m_shader.reset(new gl_shader);
	m_private->m_shader->add_file(GL_VERTEX_SHADER, "resources/shaders/mesh_tools_VertexShader.glsl");
	m_private->m_shader->add_file(GL_FRAGMENT_SHADER, "resources/shaders/mesh_tools_FragmentShader.glsl");
	m_private->m_shader->load();

	m_private->m_light.reset(new gl_light(gl_light::SPOTLIGHT));
	m_private->m_light->set_position(fvec3(-20, 20, 20));
	m_private->m_light->set_ambient(fvec3(0.75f));
	m_private->m_light->set_diffuse(fvec3(0.5f));
	m_private->m_light->set_specular(fvec3(0.1f));

	m_private->font2D = get_font_manager().create_font("Consolas", "Consolas", 12);

	// OpenGL initialization
	//glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_MULTISAMPLE);
	init_done = true;
}

void meshViewWindow::reset_view() {
	m_private->m_cam->setup(fvec3(0, 0, 50), fvec3(0, 0, 0), fvec3(0, 1, 0));
	m_private->m_arcball->reset();
	m_private->m_model_curvatures_k1.clear();
	m_private->m_model_curvatures_k2.clear();
	m_view_state.reset();
	m_private->invalidate_visualizations();
	create_model_view();
}

void meshViewWindow::create_model_view() {
	m_private->m_draw_parts.clear();

	if (m_model) {
		const auto& parts = m_model->m_parts;
		m_private->m_draw_parts.reserve(parts.size());

		for (base_math::mesh<double>* part : parts) {
			part->computeFaceProperties();
			mesh_data md;
			collect_mesh_data(part, md);
			auto prim = std::make_unique<gl_prim>();
			prim->create_from_mesh(&md, GL_FILL);
			m_private->m_draw_parts.push_back(std::move(prim));
		}
	}
}

bool meshViewWindow::generate_face_normals_view() {
	if (!m_model) return false; // no model, cannot generate view
	if (m_private->face_normals_created) return true; // already created, no need to regenerate

	m_private->m_face_normals.clear();
	if (m_model) {
		const auto& parts = m_model->m_parts; // get_parts();
		m_private->m_face_normals.reserve(parts.size());

		for (base_math::mesh<double>* part : parts) {
			float vl = float(part->getAverageEdgeLength()) * 0.5f;
			part->computeFaceProperties();
			auto mo = std::make_unique<visual_objects>();
			for (const auto& face : part->faces) {
				fvec3 face_normal = to_fvec3(face.second->normal);
				fvec3 face_center = to_fvec3(face.second->center);
				fvec3 normal_end = face_center + face_normal * vl;
				mo->add_vector(face_center, normal_end);
			}
			mo->create_prim();
			mo->get_prim()->set_color(fvec3(0.f, 0.f, 1.f));
			m_private->m_face_normals.push_back(std::move(mo));
		}
	}
	m_private->face_normals_created = true;
	return true;
}

bool meshViewWindow::generate_vertex_normals_view() {
	if (!m_model) return false; // no model, cannot generate view
	if (m_private->vertex_normals_created) return true; // already created, no need to regenerate

	m_private->m_vertex_normals.clear();
	const auto& parts = m_model->m_parts;
	m_private->m_vertex_normals.reserve(parts.size());
	for (base_math::mesh<double>* part : parts) {
		float vl = float(part->getAverageEdgeLength()) * 0.5f;
		part->computeFaceNormals();
		part->computeFaceProperties();
		part->computeVertexNormals();
		auto mo = std::make_unique<visual_objects>();
		for (const auto& vertex : part->vertices) {
			fvec3 vertex_normal = to_fvec3(vertex.second->normal);
			fvec3 vertex_position = to_fvec3(vertex.second->position);
			fvec3 normal_end = vertex_position + vertex_normal * vl;
			mo->add_vector(vertex_position, normal_end);
		}
		mo->create_prim();
		mo->get_prim()->set_color(fvec3(1.f, 0.f, 0.f));
		m_private->m_vertex_normals.push_back(std::move(mo));
	}
	m_private->vertex_normals_created = true;
	return true;
}

bool meshViewWindow::generate_principal_curvature_view() {
	if (!m_model || !m_model->curvatures_computed()) return false; // no model or curvatures not computed, cannot generate view
	if (m_private->curvature_directions_created) return true; // already created, no need to regenerate

	m_private->m_model_curvatures_k1.clear();
	m_private->m_model_curvatures_k2.clear();

	const auto& parts = m_model->m_parts; // get_parts();
	m_private->m_model_curvatures_k1.reserve(parts.size());
	m_private->m_model_curvatures_k2.reserve(parts.size());

	for (base_math::mesh<double>* part : parts) {
		double vl = part->getAverageEdgeLength() * 0.5f;
		auto mo_min = std::make_unique<visual_objects>();
		auto mo_max = std::make_unique<visual_objects>();
		// for each vertex, add curvature vectors
		for (const auto& v : part->vertices) {
			dvec3 p = v.second->position;
			dvec3 k1_end = p + v.second->curvature_info.k_1_dir * vl;
			dvec3 k2_end = p + v.second->curvature_info.k_2_dir * vl;
			// add min curvature vector in red
			mo_min->add_vector(to_fvec3(p), to_fvec3(k1_end));
			// add max curvature vector in green
			mo_max->add_vector(to_fvec3(p), to_fvec3(k2_end));
		}
		mo_min->create_prim();
		mo_min->get_prim()->rotate_to(m_private->m_draw_parts[0]->get_rotation());
		mo_min->get_prim()->set_color(fvec3(1.f, 0.f, 0.f));
		mo_max->create_prim();
		mo_max->get_prim()->rotate_to(m_private->m_draw_parts[0]->get_rotation());
		mo_max->get_prim()->set_color(fvec3(0.f, 1.f, 0.f));
		m_private->m_model_curvatures_k1.push_back(std::move(mo_min));
		m_private->m_model_curvatures_k2.push_back(std::move(mo_max));
	}
	m_private->curvature_directions_created = true;
	return true;
}

bool meshViewWindow::generate_gaussian_curvature_map_view() {
	if (!m_model || !m_model->curvatures_computed()) return false; // no model or curvatures not computed, cannot generate view
	m_model->prepare_gaussian_curvature_preview();
	return true;
}

bool meshViewWindow::generate_k1_k2_map_view() {
	if (!m_model || !m_model->curvatures_computed()) return false; // no model or curvatures not computed, cannot generate view
	m_model->prepare_k1_k2_preview();
	return true;
}

bool meshViewWindow::generate_mean_curvature_map_view() {
	if (!m_model || !m_model->curvatures_computed()) return false; // no model or curvatures not computed, cannot generate view
	m_model->prepare_mean_curvature_preview();
	return true;
}

bool meshViewWindow::toggle_map_view(int map_type) {
	std::vector<bool> map_view_flags = { m_view_state.show_gaussian_curvature, m_view_state.show_mean_curvature, m_view_state.show_k1_k2_preview };

	if (map_type < 1 || map_type > 3) {
		return false; // invalid map type
	}

	if (map_view_flags[map_type - 1]) {
		// map view is currently active, so disable it
		switch (map_type) {
		case MAP_TYPE_GAUSSIAN_CURVATURE:
			m_view_state.show_gaussian_curvature = false;
			break;
		case MAP_TYPE_MEAN_CURVATURE:
			m_view_state.show_mean_curvature = false;
			break;
		case MAP_TYPE_K1_K2_PREVIEW:
			m_view_state.show_k1_k2_preview = false;
			break;
		}
	}
	else {
		// map view is currently inactive, so enable it
		bool generation_success = false;
		switch (map_type) {
		case MAP_TYPE_GAUSSIAN_CURVATURE:
			generation_success = generate_gaussian_curvature_map_view();
			if (generation_success) {
				m_view_state.show_gaussian_curvature = true;
				m_view_state.show_mean_curvature = false;
				m_view_state.show_k1_k2_preview = false;
			}
			break;
		case MAP_TYPE_MEAN_CURVATURE:
			generation_success = generate_mean_curvature_map_view();
			if (generation_success) {
				m_view_state.show_gaussian_curvature = false;
				m_view_state.show_mean_curvature = true;
				m_view_state.show_k1_k2_preview = false;
			}
			break;
		case MAP_TYPE_K1_K2_PREVIEW:
			generation_success = generate_k1_k2_map_view();
			if (generation_success) {
				m_view_state.show_gaussian_curvature = false;
				m_view_state.show_mean_curvature = false;
				m_view_state.show_k1_k2_preview = true;
				break;
			}
		}
		if (generation_success) {
			create_model_view(); // rebuild model view to reflect the new curvature map visualization
		}
	}
	return true;
}

LRESULT meshViewWindow::OnSize(int width, int height) {
	m_private->m_arcball->resize((float)width, (float)height);
	m_private->m_cam->set_aspect(width, height);
	if (m_private->m_ucs_view.get() != nullptr)
		m_private->m_ucs_view->resize_window(width, height);
	return 0;
}

int meshViewWindow::onCommand(int cmd) {
	if (!init_done) return 0;
	if (m_model == nullptr && ID_FILE_LOADMODEL != cmd) return 0;

	switch (cmd) {
	case ID_FILE_LOADMODEL:
		break;
	case ID_EDIT_FLIPMESH:
		m_model->flip_all_faces();
		m_model->curvatures_computed() = false; // flipping invalidates curvature directions, so we need to recalculate them if the user wants to view them again
		m_view_state.reset(); // reset view state to disable any active visualizations that may not be valid after flipping
		m_private->invalidate_visualizations();
		create_model_view(); // rebuild model view to reflect the flipped geometry
		break;
	case ID_CALCULATE_MESHCURVATURES:
		m_model->compute_curvatures();
		m_view_state.reset(); // reset view state to disable any active visualizations that may not be valid after flipping
		m_private->invalidate_visualizations();
		break;
	case ID_VIEW_MESH:
		m_view_state.show_wireframe = !m_view_state.show_wireframe;
		break;
	case ID_VIEW_FACENORMALS:
		generate_face_normals_view();
		m_view_state.show_face_normals = !m_view_state.show_face_normals;
		break;
	case ID_VIEW_VERTEXNORNALS:
		if (generate_vertex_normals_view()) {
			m_view_state.show_vertex_normals = !m_view_state.show_vertex_normals;
		}
		break;
	case ID_VIEW_PRINCIPALCURVATURES:
		if (generate_principal_curvature_view()) {
			m_view_state.show_principal_curvatures = !m_view_state.show_principal_curvatures;
		}
		break;

	case ID_VIEW_GAUSS_CURVATURE:
		toggle_map_view(MAP_TYPE_GAUSSIAN_CURVATURE);
		break;
	case ID_VIEW_MEAN_CURVATURE:
		toggle_map_view(MAP_TYPE_MEAN_CURVATURE);
		break;
	case ID_VIEW_K1_K2_COMBO:
		toggle_map_view(MAP_TYPE_K1_K2_PREVIEW);
		break;
	}

	return 0;
}

void meshViewWindow::onMouseWheel(int delta, unsigned __int64 extra_btn) {
	m_private->m_cam->zoom(float(delta));
}

void meshViewWindow::onMouseMove(int x, int y, unsigned __int64 extra) {
	m_private->m_arcball->drag(float(x), float(y));
	if (dragging) {
		int deltax = x - last_mouse_x;
		int deltay = y - last_mouse_y;

		m_private->m_cam->pan(float(deltax), float(deltay)); // adjust the camera position based on mouse movement

		last_mouse_x = x;
		last_mouse_y = y;
	}
}
void meshViewWindow::onLMouseDown(int x, int y, unsigned __int64 extra) {
	m_private->m_arcball->beginDrag(float(x), float(y));
}
void meshViewWindow::onLMouseUp(int x, int y, unsigned __int64 extra) {
	m_private->m_arcball->endDrag();
}
void meshViewWindow::onRMouseDown(int x, int y, unsigned __int64 extra) {
	dragging = true;
	last_mouse_x = x;
	last_mouse_y = y;
}
void meshViewWindow::onRMouseUp(int x, int y, unsigned __int64 extra) {
	dragging = false;
}

void meshViewWindow::render_scene(fmat4& cam_matrix, fmat4& rot_mat, gl_shader* shdr)
{
	shdr->use();
	m_private->m_light->apply(shdr);
	// set the combined view matrix
	shdr->set_mat4("camera", cam_matrix);

	// render the objects using color
	shdr->set_vec4("object_color", fvec4(.9f, .9f, .9f, 1.f));

	// render the mesh parts with the current rotation applied
	{
		// render filled polygons first
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glEnable(GL_DEPTH_TEST);
		for (const auto& part : m_private->m_draw_parts) {
			part->set_draw_mode(GL_FILL);
			part->force_black = false;
			part->view_matrix = rot_mat;   // apply the current arcball rotation to the mesh parts
			//part->set_use_vertex_color(0); // ensure vertex color is disabled by default
			part->set_use_vertex_color(m_view_state.show_curvture_map() ? 1 : 0); // ensure vertex color is disabled by default
			part->render(shdr);
		}
		if (m_view_state.show_wireframe) {
			// then render wireframe on top
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glEnable(GL_POLYGON_OFFSET_LINE);
			glPolygonOffset(-1.0f, -1.0f); // pull lines toward camera
			for (const auto& part : m_private->m_draw_parts) {
				part->set_draw_mode(GL_LINE);
				part->force_black = true; // render wireframe in black
				part->set_use_vertex_color(0); // ensure vertex color is disabled for wireframe
				part->render(shdr);
			}
		}
	}

	if (m_view_state.show_face_normals) {
		for (const auto& part : m_private->m_face_normals) {
			part->get_prim()->view_matrix = rot_mat;  // apply the same rotation to the normals visualization
			part->render(m_private->m_shader.get());
		}
	}
	if (m_view_state.show_vertex_normals) {
		for (const auto& part : m_private->m_vertex_normals) {
			part->get_prim()->view_matrix = rot_mat;  // apply the same rotation to the vertex normals visualization
			part->render(m_private->m_shader.get());
		}
	}

	if (m_view_state.show_principal_curvatures) {
		glEnable(GL_POLYGON_OFFSET_LINE);
		glPolygonOffset(-1.0f, -1.0f); // pull lines toward camera
		for (const auto& part : m_private->m_model_curvatures_k1) {
			part->get_prim()->view_matrix = rot_mat;  // apply the same rotation to the curvature visualization
			part->render(m_private->m_shader.get());
		}
		for (const auto& part : m_private->m_model_curvatures_k2) {
			part->get_prim()->view_matrix = rot_mat;  // apply the same rotation to the curvature visualization
			part->render(m_private->m_shader.get());
		}
	}

	m_private->m_shader->end();
	print_info();
}

void meshViewWindow::print_info() {
	m_private->font2D->set_color(fvec4(1));
	m_private->font2D->set_position(5, 5);
	m_private->font2D->render((!m_model || m_model->get_name().empty()) ? "No model loaded" : m_model->get_name().c_str());	
	if (!m_model) return;

	int ypos = 20;
	if (m_view_state.show_principal_curvatures) {
		m_private->font2D->set_position(5, ypos);
		m_private->font2D->render("Red: k1 direction, Green: k2 direction");
		ypos += 15;
	}
	if (m_view_state.show_curvture_map()) {
		m_private->font2D->set_position(5, ypos);
		if (m_view_state.show_gaussian_curvature) {
			m_private->font2D->render("Gaussian curvature map (negative=negative(saddle), blue=near zero, green=positive(sphere-like))");
		}
		else if (m_view_state.show_mean_curvature) {
			m_private->font2D->render("Mean curvature map (red=convex, blue=flat, green=concave)");
		}
		else if (m_view_state.show_k1_k2_preview) {
			m_private->font2D->render("K1/K2 preview (red=concave ell, blue=concave cyl, green=saddle, yellow=flat, magenta=convex cyl, cyan=convex ell)");
		}
	}
}

void meshViewWindow::render() {
	if (!hWnd || !hDC || !hGLRC) return;
	begin_render();

	m_private->m_cam->set_viewport();
	glClearColor(.25f, .25f, .25f, 1.f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (init_done) {
		fmat4 cam_matrix = m_private->m_cam->perspective();
		fmat4 rot_mat(m_private->m_arcball->rotation());

		render_scene(cam_matrix, rot_mat, m_private->m_shader.get());

		// render coordinate system arrows
		if (m_private->m_ucs_view.get() != nullptr) {
			m_private->m_ucs_view->set_user_rotation(rot_mat); // apply the same rotation to the UCS view
			m_private->m_ucs_view->render();
		}
	}

	end_render();
}