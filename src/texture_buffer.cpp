TextureBuffer::TextureBuffer(bool half_precision)
{
    // 16F precision for the transmittance gives artifacts. Always use full.
    //Also using full for irradiance as they original code did.

    new_texture_2d_array(m_transmittance_array, CONSTANTS.TRANSMITTANCE_WIDTH, CONSTANTS.TRANSMITTANCE_HEIGHT, false);
    new_texture_2d_array(m_irradiance_array, CONSTANTS.IRRADIANCE_WIDTH, CONSTANTS.IRRADIANCE_HEIGHT, false);
    new_texture_3d_array(m_scattering_array, CONSTANTS.SCATTERING_WIDTH, CONSTANTS.SCATTERING_HEIGHT, CONSTANTS.SCATTERING_DEPTH, half_precision);
    new_texture_3d_array(m_optional_single_mie_scattering_array, CONSTANTS.SCATTERING_WIDTH, CONSTANTS.SCATTERING_HEIGHT, CONSTANTS.SCATTERING_DEPTH, half_precision);
    
    m_delta_irradiance_texture = new_render_texture_2d(CONSTANTS.IRRADIANCE_WIDTH, CONSTANTS.IRRADIANCE_HEIGHT, false);
    m_delta_rayleigh_scattering_texture = new_render_texture_3d(CONSTANTS.SCATTERING_WIDTH, CONSTANTS.SCATTERING_HEIGHT, CONSTANTS.SCATTERING_DEPTH, half_precision);
    m_delta_mie_scattering_texture = new_render_texture_3d(CONSTANTS.SCATTERING_WIDTH, CONSTANTS.SCATTERING_HEIGHT, CONSTANTS.SCATTERING_DEPTH, half_precision);
    m_delta_scattering_density_texture = new_render_texture_3d(CONSTANTS.SCATTERING_WIDTH, CONSTANTS.SCATTERING_HEIGHT, CONSTANTS.SCATTERING_DEPTH, half_precision);

    // delta_multiple_scattering_texture is only needed to compute scattering
    // order 3 or more, while delta_rayleigh_scattering_texture and
    // delta_mie_scattering_texture are only needed to compute double scattering.
    // Therefore, to save memory, we can store delta_rayleigh_scattering_texture
    // and delta_multiple_scattering_texture in the same GPU texture.
    m_delta_multiple_scattering_texture = m_delta_rayleigh_scattering_texture;
}

TextureBuffer::~TextureBuffer()
{

}

void TextureBuffer::clear(Program* program)
{
    clear_texture(program, m_delta_irradiance_texture);
    clear_texture(program, m_delta_rayleigh_scattering_texture);
    clear_texture(program, m_delta_mie_scattering_texture);
    clear_texture(program, m_delta_scattering_density_texture);
    clear_array(program, m_transmittance_array);
    clear_array(program, m_irradiance_array);
    clear_array(program, m_scattering_array);
    clear_array(program, m_optional_single_mie_scattering_array);
}

void TextureBuffer::new_texture_2d_array(Texture** arr, int width, int height, bool half_precision)
{
    arr[0] = new_render_texture_2d(width, height, half_precision);
    arr[1] = new_render_texture_2d(width, height, half_precision);
}

void TextureBuffer::new_texture_dd_array(Texture** arr, int width, int height, int depth, bool half_precision)
{
    arr[0] = new_render_texture_3d(width, height, depth, half_precision);
    arr[1] = new_render_texture_3d(width, height, depth, half_precision);
}

Texture* TextureBuffer::new_render_texture_2d(int width, int height, bool half_precision)
{

}

Texture* TextureBuffer::new_render_texture_3d(int width, int height, int depth, bool half_precision)
{

}

void TextureBuffer::clear_array(Program* program, Texture** arr)
{

}

void TextureBuffer::clear_texture(Program* program, Texture* arr)
{

}