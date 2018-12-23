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

layout (binding = 0, rgba32f) uniform image3D delta_rayleigh_scattering;
layout (binding = 1, rgba32f) uniform image3D delta_mie_scattering;
layout (binding = 2, rgba32f) uniform image3D scattering_read;
layout (binding = 3, rgba32f) uniform image3D scattering_write;
layout (binding = 4, rgba32f) uniform image3D single_mie_scattering_read;
layout (binding = 5, rgba32f) uniform image3D single_mie_scattering_write;

// ------------------------------------------------------------------
// UNIFORMS ---------------------------------------------------------
// ------------------------------------------------------------------

uniform vec4 blend;
uniform int layer;
uniform sampler2D transmittance;

// ------------------------------------------------------------------
// MAIN -------------------------------------------------------------
// ------------------------------------------------------------------

void main()
{
    ivec3 coord = ivec3(gl_GlobalInvocationID);
    coord.z = layer;
    vec3 frag_coord = coord + vec3(0.5,0.5,0.5);

    vec3 delta_rayleigh, delta_mie;
    ComputeSingleScatteringTexture(transmittance, frag_coord, delta_rayleigh, delta_mie);

    imageStore(delta_rayleigh_scattering, coord, vec4(delta_rayleigh, 1));

    imageStore(delta_mie_scattering, coord, vec4(delta_mie, 1));

    imageStore(scattering_write, coord, vec4(RadianceToLuminance(delta_rayleigh), RadianceToLuminance(delta_mie).r));

    imageStore(single_mie_scattering_write, coord, vec4(RadianceToLuminance(delta_mie), 1));

    if(blend[2] == 1)
        imageStore(scattering_write, coord, imageLoad(scattering_write, coord) + imageLoad(scattering_read, coord));

    if(blend[3] == 1)
        imageStore(single_mie_scattering_write, coord, imageLoad(single_mie_scattering_write, coord) + imageLoad(single_mie_scattering_read, coord));
}

// ------------------------------------------------------------------