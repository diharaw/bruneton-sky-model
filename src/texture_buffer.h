#include "ogl.h"

struct TextureBuffer
{
    dw::Texture* m_delta_irradiance_texture;
    dw::Texture* m_delta_rayleigh_scattering_texture;
    dw::Texture* m_delta_mie_scattering_texture;
    dw::Texture* m_delta_scattering_density_texture;
    dw::Texture* m_delta_multiple_scattering_texture;
    dw::Texture* m_transmittance_array[2];
    dw::Texture* m_irradiance_array[2];
    dw::Texture* m_scattering_array[2];
    dw::Texture* m_optional_single_mie_scattering_array[2];

    TextureBuffer(bool half_precision);
    ~TextureBuffer();

    void clear(dw::Program* program);
	void new_texture_2d_array(dw::Texture** arr, int width, int height, bool half_precision);
    void new_texture_3d_array(dw::Texture** arr, int width, int height, int depth, bool half_precision);

    static dw::Texture* new_render_texture_2d(int width, int height, bool half_precision);
    static dw::Texture* new_render_texture_3d(int width, int height, int depth, bool half_precision);

private:
    void clear_array(dw::Program* program, dw::Texture** arr);
    void clear_texture(dw::Program* program, dw::Texture* arr);
};