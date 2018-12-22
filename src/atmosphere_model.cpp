#include "atmosphere_model.h"

// -----------------------------------------------------------------------------------------------------------------------------------

AtmosphereModel::AtmosphereModel()
{

}

// -----------------------------------------------------------------------------------------------------------------------------------

AtmosphereModel::~AtmosphereModel()
{

}

// -----------------------------------------------------------------------------------------------------------------------------------

void AtmosphereModel::initialize(int num_scattering_orders)
{

}

// -----------------------------------------------------------------------------------------------------------------------------------

void AtmosphereModel::bind_rendering_uniforms(dw::Program* program)
{

}

// -----------------------------------------------------------------------------------------------------------------------------------

void AtmosphereModel::convert_spectrum_to_linear_srgb(double& r, double& g, double& b)
{

}

// -----------------------------------------------------------------------------------------------------------------------------------

double AtmosphereModel::coeff(double lambda, int component)
{

}

// -----------------------------------------------------------------------------------------------------------------------------------

void AtmosphereModel::bind_compute_uniforms(dw::Program* program, double lambdas[], double luminance_from_radiance[])
{

}

// -----------------------------------------------------------------------------------------------------------------------------------

void AtmosphereModel::bind_density_layer(dw::Program* program, DensityProfileLayer* layer)
{

}

// -----------------------------------------------------------------------------------------------------------------------------------

void AtmosphereModel::sky_sun_radiance_to_luminance(glm::vec3& sky_spectral_radiance_to_luminance, glm::vec3& sun_spectral_radiance_to_luminance)
{

}

// -----------------------------------------------------------------------------------------------------------------------------------

void AtmosphereModel::precompute(TextureBuffer buffer, double lambdas[], double luminance_from_radiance[], bool blend, int num_scattering_orders)
{

}

// -----------------------------------------------------------------------------------------------------------------------------------

void AtmosphereModel::swap(dw::Texture** arr)
{

}

// -----------------------------------------------------------------------------------------------------------------------------------

glm::vec3 AtmosphereModel::to_vector(const std::vector<double>& wavelengths, const std::vector<double>& v, const std::vector<double>& lambdas, double scale)
{

}

// -----------------------------------------------------------------------------------------------------------------------------------

glm::mat4 AtmosphereModel::to_matrix(double arr[])
{

}

// -----------------------------------------------------------------------------------------------------------------------------------

double AtmosphereModel::cie_color_matching_function_table_value(double wavelength, int column)
{

}

// -----------------------------------------------------------------------------------------------------------------------------------

double AtmosphereModel::interpolate(const std::vector<double>& wavelengths, const std::vector<double>& wavelength_function, double wavelength)
{

}

// -----------------------------------------------------------------------------------------------------------------------------------

void AtmosphereModel::compute_spectral_radiance_to_luminance_factors(const std::vector<double>& wavelengths, const std::vector<double>& solar_irradiance, double lambda_power, double& k_r, double& k_g, double& k_b)
{

}

// -----------------------------------------------------------------------------------------------------------------------------------