#include <application.h>
#include <mesh.h>
#include <camera.h>
#include <material.h>
#include <memory>
#include <iostream>
#include <stack>

#define _USE_MATH_DEFINES
#include <math.h>

#include "atmosphere_model.h"
#include "constants.h"

// Uniform buffer data structure.
struct GlobalUniforms
{
    DW_ALIGNED(16) glm::mat4 view;
    DW_ALIGNED(16) glm::mat4 projection;
	DW_ALIGNED(16) glm::mat4 inv_view;
	DW_ALIGNED(16) glm::mat4 inv_proj;
	DW_ALIGNED(16) glm::vec4 camera_pos;
};

#define CAMERA_FAR_PLANE 10000.0f
#define RELFECTION_MAP_SIZE 512

class AtmosphericScattering : public dw::Application
{
protected:
    
    // -----------------------------------------------------------------------------------------------------------------------------------
    
	bool init(int argc, const char* argv[]) override
	{
		// Create GPU resources.
		if (!create_shaders())
			return false;

		if (!create_uniform_buffer())
			return false;

		// Create camera.
		create_camera();
		create_cube();
		create_meshes();

		initialize_atmosphere();
		initialize_reflection_map();

		return true;
	}

	// -----------------------------------------------------------------------------------------------------------------------------------

	void update(double delta) override
	{
		// Update camera.
        update_camera();

		update_global_uniforms(m_global_uniforms);

		ImGui::SliderAngle("Sun Angle", &m_sun_angle, 0.0f, 180.0f);

		m_debug_draw.sphere(0.1f, glm::vec3(1000.0f), glm::vec3(0.0f, 1.0f, 0.0f));

		render_reflection_map();
		render_scene();
		render_sky();

        // Render debug draw.
        m_debug_draw.render(nullptr, m_width, m_height, m_debug_mode ? m_debug_camera->m_view_projection : m_main_camera->m_view_projection);
	}

	// -----------------------------------------------------------------------------------------------------------------------------------

	void shutdown() override
	{
		dw::Mesh::unload(m_floor_mesh);
		dw::Mesh::unload(m_bunny_mesh);
	}

	// -----------------------------------------------------------------------------------------------------------------------------------

	void window_resized(int width, int height) override
	{
		// Override window resized method to update camera projection.
		m_main_camera->update_projection(60.0f, 0.1f, CAMERA_FAR_PLANE, float(m_width) / float(m_height));
        m_debug_camera->update_projection(60.0f, 0.1f, CAMERA_FAR_PLANE * 2.0f, float(m_width) / float(m_height));
	}

	// -----------------------------------------------------------------------------------------------------------------------------------
    
    void key_pressed(int code) override
    {
        // Handle forward movement.
        if(code == GLFW_KEY_W)
            m_heading_speed = m_camera_speed;
        else if(code == GLFW_KEY_S)
            m_heading_speed = -m_camera_speed;
        
        // Handle sideways movement.
        if(code == GLFW_KEY_A)
            m_sideways_speed = -m_camera_speed;
        else if(code == GLFW_KEY_D)
            m_sideways_speed = m_camera_speed;
    }
    
    // -----------------------------------------------------------------------------------------------------------------------------------
    
    void key_released(int code) override
    {
        // Handle forward movement.
        if(code == GLFW_KEY_W || code == GLFW_KEY_S)
            m_heading_speed = 0.0f;
        
        // Handle sideways movement.
        if(code == GLFW_KEY_A || code == GLFW_KEY_D)
            m_sideways_speed = 0.0f;
    }
    
    // -----------------------------------------------------------------------------------------------------------------------------------

	void mouse_pressed(int code) override
	{
		// Enable mouse look.
		if (code == GLFW_MOUSE_BUTTON_RIGHT)
			m_mouse_look = true;
	}

	// -----------------------------------------------------------------------------------------------------------------------------------

	void mouse_released(int code) override
	{
		// Disable mouse look.
		if (code == GLFW_MOUSE_BUTTON_RIGHT)
			m_mouse_look = false;
	}

	// -----------------------------------------------------------------------------------------------------------------------------------

protected:

	// -----------------------------------------------------------------------------------------------------------------------------------

	dw::AppSettings intial_app_settings() override
	{
		dw::AppSettings settings;

		settings.resizable = true;
		settings.maximized = false;
		settings.refresh_rate = 60;
		settings.major_ver = 4;
		settings.width = 1280;
		settings.height = 720;
		settings.title = "Bruneton Atmospheric Scattering";

		return settings;
	}

	// -----------------------------------------------------------------------------------------------------------------------------------

private:

	// -----------------------------------------------------------------------------------------------------------------------------------

	bool create_shaders()
	{
		std::vector<std::string> defines;

		if (m_use_luminance == LUMINANCE::NONE)
			defines.push_back("RADIANCE_API_ENABLED");

		if (m_use_combined_textures)
			defines.push_back("COMBINED_SCATTERING_TEXTURES");

		{
			m_mesh_vs = std::unique_ptr<dw::Shader>(dw::Shader::create_from_file(GL_VERTEX_SHADER, "shader/mesh_vs.glsl", defines));
			m_mesh_fs = std::unique_ptr<dw::Shader>(dw::Shader::create_from_file(GL_FRAGMENT_SHADER, "shader/mesh_fs.glsl", defines));

			if (!m_mesh_vs || !m_mesh_fs)
			{
				DW_LOG_FATAL("Failed to create Shaders");
				return false;
			}

			// Create general shader program
			dw::Shader* shaders[] = { m_mesh_vs.get(), m_mesh_fs.get() };
			m_mesh_program = std::make_unique<dw::Program>(2, shaders);

			if (!m_mesh_program)
			{
				DW_LOG_FATAL("Failed to create Shader Program");
				return false;
			}

			m_mesh_program->uniform_block_binding("u_GlobalUBO", 0);
		}

		{
			m_reflection_map_vs = std::unique_ptr<dw::Shader>(dw::Shader::create_from_file(GL_VERTEX_SHADER, "shader/reflection_map_vs.glsl", defines));
			m_reflection_map_fs = std::unique_ptr<dw::Shader>(dw::Shader::create_from_file(GL_FRAGMENT_SHADER, "shader/reflection_map_fs.glsl", defines));

			if (!m_reflection_map_vs || !m_reflection_map_fs)
			{
				DW_LOG_FATAL("Failed to create Shaders");
				return false;
			}

			// Create general shader program
			dw::Shader* shaders[] = { m_reflection_map_vs.get(), m_reflection_map_fs.get() };
			m_reflection_map_program = std::make_unique<dw::Program>(2, shaders);

			if (!m_reflection_map_program)
			{
				DW_LOG_FATAL("Failed to create Shader Program");
				return false;
			}

			m_reflection_map_program->uniform_block_binding("u_GlobalUBO", 0);
		}

		{
			m_cubemap_vs = std::unique_ptr<dw::Shader>(dw::Shader::create_from_file(GL_VERTEX_SHADER, "shader/cubemap_vs.glsl"));
			m_cubemap_fs = std::unique_ptr<dw::Shader>(dw::Shader::create_from_file(GL_FRAGMENT_SHADER, "shader/cubemap_fs.glsl"));

			if (!m_cubemap_vs || !m_cubemap_fs)
			{
				DW_LOG_FATAL("Failed to create Shaders");
				return false;
			}

			// Create general shader program
			dw::Shader* shaders[] = { m_cubemap_vs.get(), m_cubemap_fs.get() };
			m_cubemap_program = std::make_unique<dw::Program>(2, shaders);

			if (!m_cubemap_program)
			{
				DW_LOG_FATAL("Failed to create Shader Program");
				return false;
			}
		}

		return true;
	}

	// -----------------------------------------------------------------------------------------------------------------------------------

	bool create_uniform_buffer()
	{
        // Create uniform buffer for global data
        m_global_ubo = std::make_unique<dw::UniformBuffer>(GL_DYNAMIC_DRAW, sizeof(GlobalUniforms));
        
		return true;
	}

	// -----------------------------------------------------------------------------------------------------------------------------------

	void create_camera()
	{
        m_main_camera = std::make_unique<dw::Camera>(60.0f, 0.1f, CAMERA_FAR_PLANE, float(m_width) / float(m_height), glm::vec3(0.0f, 2.0f, 10.0f), glm::vec3(0.0f, 0.0, -1.0f));
        m_debug_camera = std::make_unique<dw::Camera>(60.0f, 0.1f, CAMERA_FAR_PLANE * 2.0f, float(m_width) / float(m_height), glm::vec3(0.0f, 2.0f, 10.0f), glm::vec3(0.0f, 0.0, -1.0f));
	}

	// -----------------------------------------------------------------------------------------------------------------------------------

	void render_mesh(dw::Mesh* mesh, const std::unique_ptr<dw::Program>& program, glm::vec3 color, glm::mat4 model)
	{
		program->set_uniform("model", model);
		program->set_uniform("color", color);

		// Bind vertex array.
		mesh->mesh_vertex_array()->bind();

		dw::SubMesh* submeshes = mesh->sub_meshes();

		for (uint32_t i = 0; i < mesh->sub_mesh_count(); i++)
		{
			dw::SubMesh& submesh = submeshes[i];
			// Issue draw call.
			glDrawElementsBaseVertex(GL_TRIANGLES, submesh.index_count, GL_UNSIGNED_INT, (void*)(sizeof(unsigned int) * submesh.base_index), submesh.base_vertex);
		}
	}

	// -----------------------------------------------------------------------------------------------------------------------------------

	void render_scene()
	{
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_CULL_FACE);

		m_mesh_program->use();

		m_model->bind_rendering_uniforms(m_mesh_program.get());

		m_mesh_program->set_uniform("exposure", m_use_luminance != LUMINANCE::NONE ? m_exposure * 1e-5f : m_exposure);
		m_mesh_program->set_uniform("earth_center", glm::vec3(0.0f, -kBottomRadius / kLengthUnitInMeters, 0.0f));
		m_mesh_program->set_uniform("sun_size", glm::vec2(tan(kSunAngularRadius), cos(kSunAngularRadius)));
		m_mesh_program->set_uniform("sun_direction", glm::normalize(glm::vec3(0.0f, sin(m_sun_angle), cos(m_sun_angle))));
		m_mesh_program->set_uniform("white_point", m_white_point);

		m_global_ubo->bind_base(0);

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glViewport(0, 0, m_width, m_height);

		glClearColor(1.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glm::mat4 m = glm::mat4(1.0f);
		
		glm::vec3 floor_scale = glm::vec3(1000.0f);
		glm::vec3 bunny_scale = glm::vec3(2.0f);

		render_mesh(m_floor_mesh, m_mesh_program, glm::vec3(0.5f, 0.5f, 0.5f), glm::translate(m, m_floor_pos) * glm::scale(glm::mat4(1.0f), floor_scale));
		render_mesh(m_bunny_mesh, m_mesh_program, glm::vec3(0.8f, 0.8f, 0.8f), glm::translate(m, m_bunny_pos) * glm::scale(glm::mat4(1.0f), bunny_scale));
	}

	// -----------------------------------------------------------------------------------------------------------------------------------
	
	void render_sky()
	{
		glDisable(GL_CULL_FACE);

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glViewport(0, 0, m_width, m_height);

		glDepthFunc(GL_LEQUAL);

		m_cubemap_program->use();

		m_cube_vao->bind();

		m_cubemap_program->set_uniform("projection", m_global_uniforms.projection);
		m_cubemap_program->set_uniform("view", glm::mat4(glm::mat3(m_global_uniforms.view)));
		m_cubemap_program->set_uniform("exposure", m_use_luminance != LUMINANCE::NONE ? m_exposure * 1e-5f : m_exposure);
		m_cubemap_program->set_uniform("white_point", m_white_point);

		if (m_cubemap_program->set_uniform("s_Skybox", 0))
			m_reflection_map->bind(0);
			
		glDrawArrays(GL_TRIANGLES, 0, 36);

		glDepthFunc(GL_LESS);
	}
    
    // -----------------------------------------------------------------------------------------------------------------------------------
    
    void update_global_uniforms(const GlobalUniforms& global)
    {
        void* ptr = m_global_ubo->map(GL_WRITE_ONLY);
        memcpy(ptr, &global, sizeof(GlobalUniforms));
        m_global_ubo->unmap();
    }

    // -----------------------------------------------------------------------------------------------------------------------------------
    
    void update_transforms(dw::Camera* camera)
    {
        // Update camera matrices.
        m_global_uniforms.view = camera->m_view;
        m_global_uniforms.projection = camera->m_projection;
		m_global_uniforms.inv_view = glm::inverse(camera->m_view);
		m_global_uniforms.inv_proj = glm::inverse(camera->m_projection);
		m_global_uniforms.camera_pos = glm::vec4(camera->m_position, 0.0f);
    }

    // -----------------------------------------------------------------------------------------------------------------------------------
    
    void update_camera()
    {
        dw::Camera* current = m_main_camera.get();
        
        if (m_debug_mode)
            current = m_debug_camera.get();
        
        float forward_delta = m_heading_speed * m_delta;
        float right_delta = m_sideways_speed * m_delta;
        
        current->set_translation_delta(current->m_forward, forward_delta);
        current->set_translation_delta(current->m_right, right_delta);

		double d = 1 - exp(log(0.5) * m_springness * m_delta_seconds);

		m_camera_x = m_mouse_delta_x * m_camera_sensitivity;
		m_camera_y = m_mouse_delta_y * m_camera_sensitivity;
        
        if (m_mouse_look)
        {
            // Activate Mouse Look
            current->set_rotatation_delta(glm::vec3((float)(m_camera_y),
                                                    (float)(m_camera_x),
                                                    (float)(0.0f)));
        }
        else
        {
            current->set_rotatation_delta(glm::vec3((float)(0),
                                                    (float)(0),
                                                    (float)(0)));
        }
        
        current->update();

		update_transforms(current);
    }

	// -----------------------------------------------------------------------------------------------------------------------------------

	void initialize_atmosphere()
	{
		m_sun_angle = glm::radians(180.0f);

		// Values from "Reference Solar Spectral Irradiance: ASTM G-173", ETR column
		// (see http://rredc.nrel.gov/solar/spectra/am1.5/ASTMG173/ASTMG173.html),
		// summed and averaged in each bin (e.g. the value for 360nm is the average
		// of the ASTM G-173 values for all wavelengths between 360 and 370nm).
		// Values in W.m^-2.
		int lambda_min = 360;
		int lambda_max = 830;

		double kSolarIrradiance[] =
		{
			1.11776, 1.14259, 1.01249, 1.14716, 1.72765, 1.73054, 1.6887, 1.61253,
			1.91198, 2.03474, 2.02042, 2.02212, 1.93377, 1.95809, 1.91686, 1.8298,
			1.8685, 1.8931, 1.85149, 1.8504, 1.8341, 1.8345, 1.8147, 1.78158, 1.7533,
			1.6965, 1.68194, 1.64654, 1.6048, 1.52143, 1.55622, 1.5113, 1.474, 1.4482,
			1.41018, 1.36775, 1.34188, 1.31429, 1.28303, 1.26758, 1.2367, 1.2082,
			1.18737, 1.14683, 1.12362, 1.1058, 1.07124, 1.04992
		};

		// Values from http://www.iup.uni-bremen.de/gruppen/molspec/databases/
		// referencespectra/o3spectra2011/index.html for 233K, summed and averaged in
		// each bin (e.g. the value for 360nm is the average of the original values
		// for all wavelengths between 360 and 370nm). Values in m^2.
		double kOzoneCrossSection[] = 
		{
			1.18e-27, 2.182e-28, 2.818e-28, 6.636e-28, 1.527e-27, 2.763e-27, 5.52e-27,
			8.451e-27, 1.582e-26, 2.316e-26, 3.669e-26, 4.924e-26, 7.752e-26, 9.016e-26,
			1.48e-25, 1.602e-25, 2.139e-25, 2.755e-25, 3.091e-25, 3.5e-25, 4.266e-25,
			4.672e-25, 4.398e-25, 4.701e-25, 5.019e-25, 4.305e-25, 3.74e-25, 3.215e-25,
			2.662e-25, 2.238e-25, 1.852e-25, 1.473e-25, 1.209e-25, 9.423e-26, 7.455e-26,
			6.566e-26, 5.105e-26, 4.15e-26, 4.228e-26, 3.237e-26, 2.451e-26, 2.801e-26,
			2.534e-26, 1.624e-26, 1.465e-26, 2.078e-26, 1.383e-26, 7.105e-27
		};

		// From https://en.wikipedia.org/wiki/Dobson_unit, in molecules.m^-2.
		double kDobsonUnit = 2.687e20;
		// Maximum number density of ozone molecules, in m^-3 (computed so at to get
		// 300 Dobson units of ozone - for this we divide 300 DU by the integral of
		// the ozone density profile defined below, which is equal to 15km).
		double kMaxOzoneNumberDensity = 300.0 * kDobsonUnit / 15000.0;
		// Wavelength independent solar irradiance "spectrum" (not physically
		// realistic, but was used in the original implementation).
		double kConstantSolarIrradiance = 1.5;
		double kTopRadius = 6420000.0;
		double kRayleigh = 1.24062e-6;
		double kRayleighScaleHeight = 8000.0;
		double kMieScaleHeight = 1200.0;
		double kMieAngstromAlpha = 0.0;
		double kMieAngstromBeta = 5.328e-3;
		double kMieSingleScatteringAlbedo = 0.9;
		double kMiePhaseFunctionG = 0.8;
		double kGroundAlbedo = 0.1;
		double max_sun_zenith_angle = (m_use_half_precision ? 102.0 : 120.0) / 180.0 * M_PI;

		DensityProfileLayer* rayleigh_layer = new DensityProfileLayer("rayleigh", 0.0, 1.0, -1.0 / kRayleighScaleHeight, 0.0, 0.0);
		DensityProfileLayer* mie_layer = new DensityProfileLayer("mie", 0.0, 1.0, -1.0 / kMieScaleHeight, 0.0, 0.0);

		// Density profile increasing linearly from 0 to 1 between 10 and 25km, and
		// decreasing linearly from 1 to 0 between 25 and 40km. This is an approximate
		// profile from http://www.kln.ac.lk/science/Chemistry/Teaching_Resources/
		// Documents/Introduction%20to%20atmospheric%20chemistry.pdf (page 10).
		std::vector<DensityProfileLayer*> ozone_density;
		ozone_density.push_back(new DensityProfileLayer("absorption0", 25000.0, 0.0, 0.0, 1.0 / 15000.0, -2.0 / 3.0));
		ozone_density.push_back(new DensityProfileLayer("absorption1", 0.0, 0.0, 0.0, -1.0 / 15000.0, 8.0 / 3.0));

		std::vector<double> wavelengths;
		std::vector<double> solar_irradiance;
		std::vector<double> rayleigh_scattering;
		std::vector<double> mie_scattering;
		std::vector<double> mie_extinction;
		std::vector<double> absorption_extinction;
		std::vector<double> ground_albedo;

		for (int l = lambda_min; l <= lambda_max; l += 10)
		{
			double lambda = l * 1e-3;  // micro-meters
			double mie = kMieAngstromBeta / kMieScaleHeight * pow(lambda, -kMieAngstromAlpha);

			wavelengths.push_back(l);

			if (m_use_constant_solar_spectrum)
				solar_irradiance.push_back(kConstantSolarIrradiance);
			else
				solar_irradiance.push_back(kSolarIrradiance[(l - lambda_min) / 10]);

			rayleigh_scattering.push_back(kRayleigh * pow(lambda, -4));
			mie_scattering.push_back(mie * kMieSingleScatteringAlbedo);
			mie_extinction.push_back(mie);
			absorption_extinction.push_back(m_use_ozone ? kMaxOzoneNumberDensity * kOzoneCrossSection[(l - lambda_min) / 10] : 0.0);
			ground_albedo.push_back(kGroundAlbedo);
		}

		m_model = std::make_unique<AtmosphereModel>();

		m_model->m_half_precision = m_use_half_precision;
		m_model->m_combine_scattering_textures = m_use_combined_textures;
		m_model->m_use_luminance = m_use_luminance;
		m_model->m_wave_lengths = wavelengths;
		m_model->m_solar_irradiance = solar_irradiance;
		m_model->m_sun_angular_radius = kSunAngularRadius;
		m_model->m_bottom_radius = kBottomRadius;
		m_model->m_top_radius = kTopRadius;
		m_model->m_rayleigh_density = rayleigh_layer;
		m_model->m_rayleigh_scattering = rayleigh_scattering;
		m_model->m_mie_density = mie_layer;
		m_model->m_mie_scattering = mie_scattering;
		m_model->m_mie_extinction = mie_extinction;
		m_model->m_mie_phase_function_g = kMiePhaseFunctionG;
		m_model->m_absorption_density = ozone_density;
		m_model->m_absorption_extinction = absorption_extinction;
		m_model->m_ground_albedo = ground_albedo;
		m_model->m_max_sun_zenith_angle = max_sun_zenith_angle;
		m_model->m_length_unit_in_meters = kLengthUnitInMeters;

		int num_scattering_orders = 4;

		m_model->initialize(num_scattering_orders);

		double white_point_r = 1.0;
		double white_point_g = 1.0;
		double white_point_b = 1.0;

		if (m_do_white_balance)
		{
			m_model->convert_spectrum_to_linear_srgb(white_point_r, white_point_g, white_point_b);

			double white_point = (white_point_r + white_point_g + white_point_b) / 3.0;
			white_point_r /= white_point;
			white_point_g /= white_point;
			white_point_b /= white_point;
		}

		m_white_point = glm::vec3(float(white_point_r), float(white_point_g), float(white_point_b));
	}

	// -----------------------------------------------------------------------------------------------------------------------------------

	void create_cube()
	{
		float cube_vertices[] =
		{
			// back face
			-1.0f, -1.0f, -1.0f,  0.0f,  0.0f, -1.0f, 0.0f, 0.0f, // bottom-left
			1.0f,  1.0f, -1.0f,  0.0f,  0.0f, -1.0f, 1.0f, 1.0f, // top-right
			1.0f, -1.0f, -1.0f,  0.0f,  0.0f, -1.0f, 1.0f, 0.0f, // bottom-right         
			1.0f,  1.0f, -1.0f,  0.0f,  0.0f, -1.0f, 1.0f, 1.0f, // top-right
			-1.0f, -1.0f, -1.0f,  0.0f,  0.0f, -1.0f, 0.0f, 0.0f, // bottom-left
			-1.0f,  1.0f, -1.0f,  0.0f,  0.0f, -1.0f, 0.0f, 1.0f, // top-left
			// front face
			-1.0f, -1.0f,  1.0f,  0.0f,  0.0f,  1.0f, 0.0f, 0.0f, // bottom-left
			1.0f, -1.0f,  1.0f,  0.0f,  0.0f,  1.0f, 1.0f, 0.0f, // bottom-right
			1.0f,  1.0f,  1.0f,  0.0f,  0.0f,  1.0f, 1.0f, 1.0f, // top-right
			1.0f,  1.0f,  1.0f,  0.0f,  0.0f,  1.0f, 1.0f, 1.0f, // top-right
			-1.0f,  1.0f,  1.0f,  0.0f,  0.0f,  1.0f, 0.0f, 1.0f, // top-left
			-1.0f, -1.0f,  1.0f,  0.0f,  0.0f,  1.0f, 0.0f, 0.0f, // bottom-left
			// left face
			-1.0f,  1.0f,  1.0f, -1.0f,  0.0f,  0.0f, 1.0f, 0.0f, // top-right
			-1.0f,  1.0f, -1.0f, -1.0f,  0.0f,  0.0f, 1.0f, 1.0f, // top-left
			-1.0f, -1.0f, -1.0f, -1.0f,  0.0f,  0.0f, 0.0f, 1.0f, // bottom-left
			-1.0f, -1.0f, -1.0f, -1.0f,  0.0f,  0.0f, 0.0f, 1.0f, // bottom-left
			-1.0f, -1.0f,  1.0f, -1.0f,  0.0f,  0.0f, 0.0f, 0.0f, // bottom-right
			-1.0f,  1.0f,  1.0f, -1.0f,  0.0f,  0.0f, 1.0f, 0.0f, // top-right
			// right face
			1.0f,  1.0f,  1.0f,  1.0f,  0.0f,  0.0f, 1.0f, 0.0f, // top-left
			1.0f, -1.0f, -1.0f,  1.0f,  0.0f,  0.0f, 0.0f, 1.0f, // bottom-right
			1.0f,  1.0f, -1.0f,  1.0f,  0.0f,  0.0f, 1.0f, 1.0f, // top-right         
			1.0f, -1.0f, -1.0f,  1.0f,  0.0f,  0.0f, 0.0f, 1.0f, // bottom-right
			1.0f,  1.0f,  1.0f,  1.0f,  0.0f,  0.0f, 1.0f, 0.0f, // top-left
			1.0f, -1.0f,  1.0f,  1.0f,  0.0f,  0.0f, 0.0f, 0.0f, // bottom-left     
			// bottom face
			-1.0f, -1.0f, -1.0f,  0.0f, -1.0f,  0.0f, 0.0f, 1.0f, // top-right
			1.0f, -1.0f, -1.0f,  0.0f, -1.0f,  0.0f, 1.0f, 1.0f, // top-left
			1.0f, -1.0f,  1.0f,  0.0f, -1.0f,  0.0f, 1.0f, 0.0f, // bottom-left
			1.0f, -1.0f,  1.0f,  0.0f, -1.0f,  0.0f, 1.0f, 0.0f, // bottom-left
			-1.0f, -1.0f,  1.0f,  0.0f, -1.0f,  0.0f, 0.0f, 0.0f, // bottom-right
			-1.0f, -1.0f, -1.0f,  0.0f, -1.0f,  0.0f, 0.0f, 1.0f, // top-right
			// top face
			-1.0f,  1.0f, -1.0f,  0.0f,  1.0f,  0.0f, 0.0f, 1.0f, // top-left
			1.0f,  1.0f , 1.0f,  0.0f,  1.0f,  0.0f, 1.0f, 0.0f, // bottom-right
			1.0f,  1.0f, -1.0f,  0.0f,  1.0f,  0.0f, 1.0f, 1.0f, // top-right     
			1.0f,  1.0f,  1.0f,  0.0f,  1.0f,  0.0f, 1.0f, 0.0f, // bottom-right
			-1.0f,  1.0f, -1.0f,  0.0f,  1.0f,  0.0f, 0.0f, 1.0f, // top-left
			-1.0f,  1.0f,  1.0f,  0.0f,  1.0f,  0.0f, 0.0f, 0.0f  // bottom-left        
		};

		m_cube_vbo = std::make_unique<dw::VertexBuffer>(GL_STATIC_DRAW, sizeof(cube_vertices), (void*)cube_vertices);

		dw::VertexAttrib attribs[] =
		{
			{ 3,GL_FLOAT, false, 0, },
			{ 3,GL_FLOAT, false, sizeof(float) * 3 },
			{ 2,GL_FLOAT, false, sizeof(float) * 6 }
		};

		m_cube_vao = std::make_unique<dw::VertexArray>(m_cube_vbo.get(), nullptr, sizeof(float) * 8, 3, attribs);

		if (!m_cube_vbo || !m_cube_vao)
		{
			DW_LOG_FATAL("Failed to create Vertex Buffers/Arrays");
		}
	}

	// -----------------------------------------------------------------------------------------------------------------------------------
	
	void create_meshes()
	{
		m_floor_mesh = dw::Mesh::load("mesh/plane.obj", false);
		m_bunny_mesh = dw::Mesh::load("mesh/bunny.obj", false);
	}

	// -----------------------------------------------------------------------------------------------------------------------------------

	void initialize_reflection_map()
	{
		glm::mat4 capture_projection = glm::perspective(glm::radians(90.0f), 1.0f, 0.1f, 10.0f);
		glm::mat4 capture_views[] =
		{
			glm::lookAt(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(1.0f,  0.0f,  0.0f), glm::vec3(0.0f, -1.0f,  0.0f)),
			glm::lookAt(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(-1.0f,  0.0f,  0.0f), glm::vec3(0.0f, -1.0f,  0.0f)),
			glm::lookAt(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f,  1.0f,  0.0f), glm::vec3(0.0f,  0.0f,  1.0f)),
			glm::lookAt(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, -1.0f,  0.0f), glm::vec3(0.0f,  0.0f, -1.0f)),
			glm::lookAt(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f,  0.0f,  1.0f), glm::vec3(0.0f, -1.0f,  0.0f)),
			glm::lookAt(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f,  0.0f, -1.0f), glm::vec3(0.0f, -1.0f,  0.0f))
		};

		for (int i = 0; i < 6; i++)
			m_reflection_views[i] = capture_projection * capture_views[i];

		m_reflection_map = std::make_unique<dw::TextureCube>(RELFECTION_MAP_SIZE, RELFECTION_MAP_SIZE, 1, 1, GL_RGB16F, GL_RGB, GL_HALF_FLOAT);
		m_reflection_map->set_min_filter(GL_LINEAR);
		m_reflection_map->set_mag_filter(GL_LINEAR);

		glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);

		m_reflection_fbos.resize(6);

		for (int i = 0; i < 6; i++)
		{
			m_reflection_fbos[i] = std::make_unique<dw::Framebuffer>();
			m_reflection_fbos[i]->attach_render_target(0, m_reflection_map.get(), i, 0, 0);
		}
	}

	// -----------------------------------------------------------------------------------------------------------------------------------

	void render_reflection_map()
	{
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_CULL_FACE);

		m_reflection_map_program->use();

		m_model->bind_rendering_uniforms(m_reflection_map_program.get());

		m_reflection_map_program->set_uniform("exposure", m_use_luminance != LUMINANCE::NONE ? m_exposure * 1e-5f : m_exposure);
		m_reflection_map_program->set_uniform("earth_center", glm::vec3(0.0f, -kBottomRadius / kLengthUnitInMeters, 0.0f));
		m_reflection_map_program->set_uniform("sun_size", glm::vec2(tan(kSunAngularRadius), cos(kSunAngularRadius)));
		m_reflection_map_program->set_uniform("sun_direction", glm::normalize(glm::vec3(0.0f, sin(m_sun_angle), cos(m_sun_angle))));
		m_reflection_map_program->set_uniform("white_point", m_white_point);

		m_global_ubo->bind_base(0);

		for (int i = 0; i < 6; i++)
		{
			m_reflection_map_program->set_uniform("view_projection", m_reflection_views[i]);

			m_reflection_fbos[i]->bind();
			glViewport(0, 0, RELFECTION_MAP_SIZE, RELFECTION_MAP_SIZE);

			glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			m_cube_vao->bind();

			glDrawArrays(GL_TRIANGLES, 0, 36);
		}
	}
    
    // -----------------------------------------------------------------------------------------------------------------------------------
    
private:
	// General GPU resources.
    std::unique_ptr<dw::UniformBuffer> m_global_ubo;
	std::unique_ptr<dw::Shader>  m_mesh_vs;
	std::unique_ptr<dw::Shader>  m_mesh_fs;
	std::unique_ptr<dw::Program> m_mesh_program;

    // Camera.
    std::unique_ptr<dw::Camera> m_main_camera;
    std::unique_ptr<dw::Camera> m_debug_camera;

	// Reflection Map
	std::unique_ptr<dw::TextureCube> m_reflection_map;
	std::vector<std::unique_ptr<dw::Framebuffer>> m_reflection_fbos;
	glm::mat4 m_reflection_views[6];
	std::unique_ptr<dw::Shader>  m_reflection_map_vs;
	std::unique_ptr<dw::Shader>  m_reflection_map_fs;
	std::unique_ptr<dw::Program> m_reflection_map_program;

	std::unique_ptr<dw::VertexBuffer> m_cube_vbo;
	std::unique_ptr<dw::VertexArray> m_cube_vao;

	std::unique_ptr<dw::Shader>  m_cubemap_vs;
	std::unique_ptr<dw::Shader>  m_cubemap_fs;
	std::unique_ptr<dw::Program> m_cubemap_program;

	dw::Mesh* m_floor_mesh;
	dw::Mesh* m_bunny_mesh;

	glm::vec3 m_floor_pos = glm::vec3(0.0f);
	glm::vec3 m_bunny_pos = glm::vec3(0.0f, -0.1f, 0.0f);
    
	// Uniforms.
    GlobalUniforms m_global_uniforms;

    // Camera controls.
    bool m_mouse_look = false;
    bool m_debug_mode = false;
    float m_heading_speed = 0.0f;
    float m_sideways_speed = 0.0f;
    float m_camera_sensitivity = 0.05f;
    float m_camera_speed = 0.01f;

	// Camera orientation.
	float m_camera_x;
	float m_camera_y;
	float m_springness = 1.0f;

	// Atmosphere settings
	glm::vec3 m_white_point;
	float m_sun_angle = 0.0f;
	float kSunAngularRadius = 0.00935f / 2.0f;
	float kBottomRadius = 6360000.0f;
	float kLengthUnitInMeters = 1000.0f;
	bool m_use_constant_solar_spectrum = false;
	bool m_use_ozone = true;
	bool m_use_combined_textures = true;
	bool m_use_half_precision = false;
	bool m_do_white_balance = false;
	LUMINANCE m_use_luminance = LUMINANCE::NONE;
	float m_exposure = 10.0f;
	std::unique_ptr<AtmosphereModel> m_model;
};

DW_DECLARE_MAIN(AtmosphericScattering)
