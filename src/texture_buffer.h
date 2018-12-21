
struct TextureBuffer
{
    Texture* m_delta_irradiance_texture;
    Texture* m_delta_rayleigh_scattering_texture;
    Texture* m_delta_mie_scattering_texture;
    Texture* m_delta_scattering_density_texture;
    Texture* m_delta_multiple_scattering_texture;
    Texture* m_transmittance_array[2];
    Texture* m_irradiance_array[2];
    Texture* m_scattering_array[2];
    Texture* m_optional_single_mie_scattering_array[2];

    TextureBuffer(bool half_precision);
    ~TextureBuffer();

    void clear(Program* program);
    void new_texture_2d_array(Texture** arr, int width, int height, bool half_precision)
    void new_texture_dd_array(Texture** arr, int width, int height, int depth, bool half_precision);

    static Texture* new_render_texture_2d(int width, int height, bool half_precision);
    static Texture* new_render_texture_3d(int width, int height, int depth, bool half_precision);

private:
    void clear_array(Program* program, Texture** arr);
    void clear_texture(Program* program, Texture* arr);
};