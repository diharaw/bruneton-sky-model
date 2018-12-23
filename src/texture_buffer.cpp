#include "texture_buffer.h"
#include "constants.h"
#include <macros.h>

// -----------------------------------------------------------------------------------------------------------------------------------

TextureBuffer::TextureBuffer(bool half_precision)
{
	m_delta_irradiance_texture = nullptr;
	m_delta_rayleigh_scattering_texture = nullptr;
	m_delta_mie_scattering_texture = nullptr;
	m_delta_scattering_density_texture = nullptr;
	m_delta_multiple_scattering_texture = nullptr;
	m_transmittance_array[0] = nullptr;
	m_transmittance_array[1] = nullptr;
	m_irradiance_array[0] = nullptr;
	m_irradiance_array[1] = nullptr;
	m_scattering_array[0] = nullptr;
	m_scattering_array[1] = nullptr;
	m_optional_single_mie_scattering_array[0] = nullptr;
	m_optional_single_mie_scattering_array[1] = nullptr;

    // 16F precision for the transmittance gives artifacts. Always use full.
    //Also using full for irradiance as they original code did.

    new_texture_2d_array(m_transmittance_array, CONSTANTS::TRANSMITTANCE_WIDTH, CONSTANTS::TRANSMITTANCE_HEIGHT, false);
    new_texture_2d_array(m_irradiance_array, CONSTANTS::IRRADIANCE_WIDTH, CONSTANTS::IRRADIANCE_HEIGHT, false);
    new_texture_3d_array(m_scattering_array, CONSTANTS::SCATTERING_WIDTH, CONSTANTS::SCATTERING_HEIGHT, CONSTANTS::SCATTERING_DEPTH, half_precision);
    new_texture_3d_array(m_optional_single_mie_scattering_array, CONSTANTS::SCATTERING_WIDTH, CONSTANTS::SCATTERING_HEIGHT, CONSTANTS::SCATTERING_DEPTH, half_precision);
    
    m_delta_irradiance_texture = new_texture_2d(CONSTANTS::IRRADIANCE_WIDTH, CONSTANTS::IRRADIANCE_HEIGHT, false);
    m_delta_rayleigh_scattering_texture = new_texture_3d(CONSTANTS::SCATTERING_WIDTH, CONSTANTS::SCATTERING_HEIGHT, CONSTANTS::SCATTERING_DEPTH, half_precision);
    m_delta_mie_scattering_texture = new_texture_3d(CONSTANTS::SCATTERING_WIDTH, CONSTANTS::SCATTERING_HEIGHT, CONSTANTS::SCATTERING_DEPTH, half_precision);
    m_delta_scattering_density_texture = new_texture_3d(CONSTANTS::SCATTERING_WIDTH, CONSTANTS::SCATTERING_HEIGHT, CONSTANTS::SCATTERING_DEPTH, half_precision);

    // delta_multiple_scattering_texture is only needed to compute scattering
    // order 3 or more, while delta_rayleigh_scattering_texture and
    // delta_mie_scattering_texture are only needed to compute double scattering.
    // Therefore, to save memory, we can store delta_rayleigh_scattering_texture
    // and delta_multiple_scattering_texture in the same GPU texture.
    m_delta_multiple_scattering_texture = m_delta_rayleigh_scattering_texture;
}

// -----------------------------------------------------------------------------------------------------------------------------------

TextureBuffer::~TextureBuffer()
{
	DW_SAFE_DELETE(m_delta_irradiance_texture);
	DW_SAFE_DELETE(m_delta_rayleigh_scattering_texture);
	DW_SAFE_DELETE(m_delta_mie_scattering_texture);
	DW_SAFE_DELETE(m_delta_scattering_density_texture);
	DW_SAFE_DELETE(m_transmittance_array[0]);
	DW_SAFE_DELETE(m_transmittance_array[1]);
	DW_SAFE_DELETE(m_irradiance_array[0]);
	DW_SAFE_DELETE(m_irradiance_array[1]);
	DW_SAFE_DELETE(m_scattering_array[0]);
	DW_SAFE_DELETE(m_scattering_array[1]);
	DW_SAFE_DELETE(m_optional_single_mie_scattering_array[0]);
	DW_SAFE_DELETE(m_optional_single_mie_scattering_array[1]);
}				   

// -----------------------------------------------------------------------------------------------------------------------------------

void TextureBuffer::clear(dw::Program* program_2d, dw::Program* program_3d)
{
    clear_texture(program_2d, m_delta_irradiance_texture);
    clear_texture(program_3d, m_delta_rayleigh_scattering_texture);
    clear_texture(program_3d, m_delta_mie_scattering_texture);
    clear_texture(program_3d, m_delta_scattering_density_texture);
    clear_array(program_2d, m_transmittance_array);
    clear_array(program_2d, m_irradiance_array);
    clear_array(program_3d, m_scattering_array);
    clear_array(program_3d, m_optional_single_mie_scattering_array);
}

// -----------------------------------------------------------------------------------------------------------------------------------

void TextureBuffer::new_texture_2d_array(dw::Texture** arr, int width, int height, bool half_precision)
{
    arr[0] = new_texture_2d(width, height, half_precision);
    arr[1] = new_texture_2d(width, height, half_precision);
}

// -----------------------------------------------------------------------------------------------------------------------------------

void TextureBuffer::new_texture_3d_array(dw::Texture** arr, int width, int height, int depth, bool half_precision)
{
    arr[0] = new_texture_3d(width, height, depth, half_precision);
    arr[1] = new_texture_3d(width, height, depth, half_precision);
}

// -----------------------------------------------------------------------------------------------------------------------------------

dw::Texture* TextureBuffer::new_texture_2d(int width, int height, bool half_precision)
{
	dw::Texture* texture = new dw::Texture2D(width, height, 1, 1, 1, half_precision ? GL_RGBA16F : GL_RGBA32F, GL_RGBA, half_precision ? GL_HALF_FLOAT : GL_FLOAT);
	texture->set_min_filter(GL_LINEAR);
	texture->set_wrapping(GL_CLAMP_TO_EDGE, GL_CLAMP_TO_EDGE, GL_CLAMP_TO_EDGE);

	return texture;
}

// -----------------------------------------------------------------------------------------------------------------------------------

dw::Texture* TextureBuffer::new_texture_3d(int width, int height, int depth, bool half_precision)
{
	dw::Texture* texture = new dw::Texture3D(width, height, depth, 1, half_precision ? GL_RGBA16F : GL_RGBA32F, GL_RGBA, half_precision ? GL_HALF_FLOAT : GL_FLOAT);
	texture->set_min_filter(GL_LINEAR);
	texture->set_wrapping(GL_CLAMP_TO_EDGE, GL_CLAMP_TO_EDGE, GL_CLAMP_TO_EDGE);

	return texture;
}

// -----------------------------------------------------------------------------------------------------------------------------------

void TextureBuffer::clear_array(dw::Program* program, dw::Texture** arr)
{
	clear_texture(program, arr[0]);
	clear_texture(program, arr[1]);
}

// -----------------------------------------------------------------------------------------------------------------------------------

void TextureBuffer::clear_texture(dw::Program* program, dw::Texture* tex)
{
	tex->bind_image(0, 0, 0, GL_READ_WRITE, tex->internal_format());
	program->use();

	if (tex->target() == GL_TEXTURE_3D)
	{
		dw::Texture3D* tex_3d = (dw::Texture3D*)tex;

		int width = tex_3d->width();
		int height = tex_3d->height();
		int depth = tex_3d->depth();

		glDispatchCompute(width / CONSTANTS::NUM_THREADS, height / CONSTANTS::NUM_THREADS, depth / CONSTANTS::NUM_THREADS);
	}
	else
	{
		dw::Texture2D* tex_2d = (dw::Texture2D*)tex;

		int width = tex_2d->width();
		int height = tex_2d->height();

		glDispatchCompute(width / CONSTANTS::NUM_THREADS, height / CONSTANTS::NUM_THREADS, 1);
	}
}

// -----------------------------------------------------------------------------------------------------------------------------------