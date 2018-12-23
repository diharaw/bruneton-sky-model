#include <constants.glsl>
#include <uniforms.glsl>
#include <utility.glsl>
#include <transmittance_functions.glsl>
#include <scattering_functions.glsl>
#include <irradiance_functions.glsl>

// ------------------------------------------------------------------
// INPUTS -----------------------------------------------------------
// ------------------------------------------------------------------

layout (local_size_x = LOCAL_SIZE, local_size_y = LOCAL_SIZE, local_size_z = 1) in;

// ------------------------------------------------------------------
// IMAGES -----------------------------------------------------------
// ------------------------------------------------------------------

layout (binding = 0, rgba32f) uniform image3D delta_multiple_scattering;
layout (binding = 1, rgba32f) uniform image3D scattering_read;
layout (binding = 2, rgba32f) uniform image3D scattering_write;

// ------------------------------------------------------------------
// UNIFORMS ---------------------------------------------------------
// ------------------------------------------------------------------

uniform vec4 blend;
uniform int layer;

uniform sampler2D transmittance;
uniform sampler3D delta_scattering_density;

// ------------------------------------------------------------------
// MAIN -------------------------------------------------------------
// ------------------------------------------------------------------

void main()
{
    ivec3 coord = ivec3(gl_GlobalInvocationID);
    coord.z = layer;

    vec3 frag_coord = coord + vec3(0.5, 0.5, 0.5);

    float nu;
    vec3 delta_multiple_scattering_value = ComputeMultipleScatteringTexture(transmittance, delta_scattering_density, frag_coord, nu);

    imageStore(delta_multiple_scattering, coord, vec4(delta_multiple_scattering_value, 1.0));

    imageStore(scattering_write, coord, vec4(RadianceToLuminance(delta_multiple_scattering_value) / RayleighPhaseFunction(nu), 0.0));

    if(blend[1] == 1)
        imageStore(scattering_write, coord, imageLoad(scattering_write, coord) + imageLoad(scattering_read, coord));
}

// ------------------------------------------------------------------