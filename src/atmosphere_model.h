
/// <summary>
/// An atmosphere layer of width 'width' (in m), and whose density is defined as
/// 'exp_term' * exp('exp_scale' * h) + 'linear_term' * h + 'constant_term',
/// clamped to [0,1], and where h is the altitude (in m). 'exp_term' and
/// 'constant_term' are unitless, while 'exp_scale' and 'linear_term' are in m^-1.
/// </summary>
struct DensityProfileLayer
{
    std::string name;
    double width;
    double exp_term;
    double exp_scale;
    double linear_term;
    double constant_term;

    DensityProfileLayer(std::string _name, double _width, double _exp_term, double _exp_scale, double _linear_term, double _constant_term)
    {
        name = _name;
        width = _width;
        exp_term = _exp_term;
        exp_scale = _exp_scale;
        linear_term = _linear_term;
        constant_term = _constant_term;
    }
};

class AtmosphereModel
{
private:
    const double kLambdaR = 680.0;
    const double kLambdaG = 550.0;
    const double kLambdaB = 440.0;

    const int kLambdaMin = 360;
    const int kLambdaMax = 830;

public:
    /// <summary>
    /// The wavelength values, in nanometers, and sorted in increasing order, for
    /// which the solar_irradiance, rayleigh_scattering, mie_scattering,
    /// mie_extinction and ground_albedo samples are provided. If your shaders
    /// use luminance values (as opposed to radiance values, see above), use a
    /// large number of wavelengths (e.g. between 15 and 50) to get accurate
    /// results (this number of wavelengths has absolutely no impact on the
    /// shader performance).
    /// </summary>
    std::vector<double> m_wave_lengths;

    /// <summary>
    /// The solar irradiance at the top of the atmosphere, in W/m^2/nm. This
    /// vector must have the same size as the wavelengths parameter.
    /// </summary>
    std::vector<double> m_solar_irradiance;

    /// <summary>
    /// The sun's angular radius, in radians. Warning: the implementation uses
    /// approximations that are valid only if this value is smaller than 0.1.
    /// </summary>
    double m_sun_angular_radius;

    /// <summary>
    /// The distance between the planet center and the bottom of the atmosphere in m.
    /// </summary>
    double m_bottom_radius;

    /// <summary>
    /// The distance between the planet center and the top of the atmosphere in m.
    /// </summary>
    double m_top_radius;

    /// <summary>
    /// The density profile of air molecules, i.e. a function from altitude to
    /// dimensionless values between 0 (null density) and 1 (maximum density).
    /// Layers must be sorted from bottom to top. The width of the last layer is
    /// ignored, i.e. it always extend to the top atmosphere boundary. At most 2
    /// layers can be specified.
    /// </summary>
    DensityProfileLayer m_rayleigh_density;

    /// <summary>
    /// The scattering coefficient of air molecules at the altitude where their
    /// density is maximum (usually the bottom of the atmosphere), as a function
    /// of wavelength, in m^-1. The scattering coefficient at altitude h is equal
    /// to 'rayleigh_scattering' times 'rayleigh_density' at this altitude. This
    /// vector must have the same size as the wavelengths parameter.
    /// </summary>
    std::vector<double> m_rayleigh_scattering;

    /// <summary>
    /// The density profile of aerosols, i.e. a function from altitude to
    /// dimensionless values between 0 (null density) and 1 (maximum density).
    /// Layers must be sorted from bottom to top. The width of the last layer is
    /// ignored, i.e. it always extend to the top atmosphere boundary. At most 2
    /// layers can be specified.
    /// </summary>
    DensityProfileLayer m_mie_density;

    /// <summary>
    /// The scattering coefficient of aerosols at the altitude where their
    /// density is maximum (usually the bottom of the atmosphere), as a function
    /// of wavelength, in m^-1. The scattering coefficient at altitude h is equal
    /// to 'mie_scattering' times 'mie_density' at this altitude. This vector
    /// must have the same size as the wavelengths parameter.
    /// </summary>
    std::vector<double> m_mie_scattering;

    /// <summary>
    /// The extinction coefficient of aerosols at the altitude where their
    /// density is maximum (usually the bottom of the atmosphere), as a function
    /// of wavelength, in m^-1. The extinction coefficient at altitude h is equal
    /// to 'mie_extinction' times 'mie_density' at this altitude. This vector
    /// must have the same size as the wavelengths parameter.
    /// </summary>
    std::vector<double> m_mie_extinction;

    /// <summary>
    /// The asymetry parameter for the Cornette-Shanks phase function for the aerosols.
    /// </summary>
    double m_mie_phase_function_g;

    /// <summary>
    /// The density profile of air molecules that absorb light (e.g. ozone), i.e.
    /// a function from altitude to dimensionless values between 0 (null density)
    /// and 1 (maximum density). Layers must be sorted from bottom to top. The
    /// width of the last layer is ignored, i.e. it always extend to the top
    /// atmosphere boundary. At most 2 layers can be specified.
    /// </summary>
    std::vector<DensityProfileLayer> m_absorption_density;

    /// <summary>
    /// The extinction coefficient of molecules that absorb light (e.g. ozone) at
    /// the altitude where their density is maximum, as a function of wavelength,
    /// in m^-1. The extinction coefficient at altitude h is equal to
    /// 'absorption_extinction' times 'absorption_density' at this altitude. This
    /// vector must have the same size as the wavelengths parameter.
    /// </summary>
    std::vector<double> m_absorption_extinction;

    /// <summary>
    /// The average albedo of the ground, as a function of wavelength. This
    /// vector must have the same size as the wavelengths parameter.
    /// </summary>
    std::vector<double> m_ground_albedo;

    /// <summary>
    /// The maximum Sun zenith angle for which atmospheric scattering must be
    /// precomputed, in radians (for maximum precision, use the smallest Sun
    /// zenith angle yielding negligible sky light radiance values. For instance,
    /// for the Earth case, 102 degrees is a good choice for most cases (120
    /// degrees is necessary for very high exposure values).
    /// </summary>
    double m_max_sun_zenith_angle;

    /// <summary>
    /// The length unit used in your shaders and meshes. This is the length unit
    /// which must be used when calling the atmosphere model shader functions.
    /// </summary>
    double m_length_unit_in_meters;

    /// <summary>
    /// Use radiance or luminance mode.
    /// </summary>
    LUMINANCE m_use_luminance;

    /// <summary>
    /// The number of wavelengths for which atmospheric scattering must be
    /// precomputed (the temporary GPU memory used during precomputations, and
    /// the GPU memory used by the precomputed results, is independent of this
    /// number, but the precomputation time is directly proportional to this number):
    /// - if this number is less than or equal to 3, scattering is precomputed
    /// for 3 wavelengths, and stored as irradiance values. Then both the
    /// radiance-based and the luminance-based API functions are provided (see
    /// the above note).
    /// - otherwise, scattering is precomputed for this number of wavelengths
    /// (rounded up to a multiple of 3), integrated with the CIE color matching
    /// functions, and stored as illuminance values. Then only the
    /// luminance-based API functions are provided (see the above note).
    /// </summary>
    int m_num_precomputed_wavelengths;

    inline int num_precomputed_wavelengths() { return m_use_luminance == LUMINANCE.PRECOMPUTED ? 15 : 3; }

    /// <summary>
    /// Whether to pack the (red component of the) single Mie scattering with the
    /// Rayleigh and multiple scattering in a single texture, or to store the
    /// (3 components of the) single Mie scattering in a separate texture.
    /// </summary>
    bool m_combine_scattering_textures;

    /// <summary>
    /// Whether to use half precision floats (16 bits) or single precision floats
    /// (32 bits) for the precomputed textures. Half precision is sufficient for
    /// most cases, except for very high exposure values.
    /// </summary>
    bool m_half_precision;

    Texture* m_transmittance_texture;

    Texture* m_scattering_texture;

    Texture* m_irradiance_texture;

    Texture* m_optional_single_mie_scattering_texture;

    Program* m_compute_program;

public:
    AtmosphereModel();
    ~AtmosphereModel();

    void initialize(int num_scattering_orders);
    void bind_rendering_uniforms(Program* program);
    void convert_spectrum_to_linear_srgb(double& r, double& g, double& b);
   
private:
    double coeff(double lambda, int component);
    void bind_compute_uniforms(Program* program, double lambdas[], double luminance_from_radiance[]);
    void bind_density_layer(Program* program, DensityProfileLayer layer);
    glm::vec3 to_vector(std::vector<double> wavelengths, std::vector<double> v, std::vector<double> lambdas, double scale);
    void sky_sun_radiance_to_luminance(glm::vec3& sky_spectral_radiance_to_luminance, glm::vec3& sun_spectral_radiance_to_luminance);
    void precompute(TextureBuffer buffer, double lambdas[], double luminance_from_radiance[], bool blend, int num_scattering_orders);
    void swap(RenderTexture arr[]);
    float[] to_matrix(double arr[]);

    static double cie_color_matching_function_table_value(double wavelength, int column);
    static double interpolate(std::vector<double> wavelengths, std::vector<double> wavelength_function, double wavelength);
    static void compute_spectral_radiance_to_luminance_factors(std::vector<double> wavelengths, std::vector<double> solar_irradiance, double lambda_power, double& k_r, double& k_g, double& k_b);
};