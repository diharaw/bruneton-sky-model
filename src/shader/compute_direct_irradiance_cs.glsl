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
// INPUT ------------------------------------------------------------
// ------------------------------------------------------------------

layout (binding = 0, rgba32f) uniform image2D irradiance_read;
layout (binding = 1, rgba32f) uniform image2D irradiance_write;
layout (binding = 2, rgba32f) uniform image2D delta_irradiance;

// ------------------------------------------------------------------
// UNIFORMS ---------------------------------------------------------
// ------------------------------------------------------------------

uniform vec4 blend;

uniform sampler2D transmittance;

// ------------------------------------------------------------------
// MAIN -------------------------------------------------------------
// ------------------------------------------------------------------

void main()
{
    ivec2 coord = ivec2(gl_GlobalInvocationID.xy);
    vec2 frag_coord = coord + vec2(0.5, 0.5);
    
    imageStore(delta_irradiance, coord, vec4(ComputeDirectIrradianceTexture(transmittance, frag_coord), 1.0));
    imageStore(irradiance_write, coord, vec4(0.0, 0.0, 0.0, 1.0));

    if(blend[1] == 1)
      imageStore(irradiance_write, coord, imageLoad(irradiance_write, coord) + imageLoad(irradiance_read, coord));
}

// ------------------------------------------------------------------