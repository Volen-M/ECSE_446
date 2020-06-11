/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#version 330 core


#define PI 3.14159265359
#define MAX_NUM_EMITTER_TRIANGLES 40 // Max number of emitter triangles allowed (tuned for A4)
uniform float emitterVertices[MAX_NUM_EMITTER_TRIANGLES*3*3]; // Need to specify a size at compile time (max size is 512 floats)

uniform int nbTriangles; // Use nbTriangles to index the buffer correctly
uniform vec3 lightIrradiance;
uniform vec3 albedo;
uniform vec2 windowSize; // [width, height]

uniform sampler2D cvTerm; // Creates a 2D texture we can sample from (range x,y = [0,1])

in vec3 vNormal;
in vec3 vPos;

out vec3 color;

// Compute edge (v1--v2) contribution
float getEdgeContrib(vec3 v1, vec3 v2, vec3 pos) {
	// Adapt your getEdgeContrib code from the offline part
	float value = 0.f;
// TODO(A4): Implement this
	float dotSolution = dot(((v1 - pos) / length(v1 - pos)), ((v2 - pos) / length(v2 - pos)));
	float theta = acos(dotSolution);
	vec3 gamma = cross(v2-pos,v1-pos) / length(cross(v2-pos,v1-pos));
	value = dot(gamma, vNormal) * theta;
	return value;
}


void main()
{	
	color = vec3(0);

	// 1) Extract vertices of triangles from `emitterVertices` buffer using `nbTriangles`
	// 2) Calculate G term
	// 3) Subtract modification term for G after extracting it from texture (use built-in `texture()` function)
	//	    e.g. `vec3 delta = texture(cvTerm, coords).xyz;`

    // TODO(A4): Implement this
	float E_poly = 0.f;
	for (int i = 0; i < nbTriangles; i++) { 
		vec3 v1 = vec3(emitterVertices[i*9+0], emitterVertices[i*9+1], emitterVertices[i*9+2]);
		vec3 v2 = vec3(emitterVertices[i*9+3], emitterVertices[i*9+4], emitterVertices[i*9+5]);
		vec3 v3 = vec3(emitterVertices[i*9+6], emitterVertices[i*9+7], emitterVertices[i*9+8]);

		E_poly += getEdgeContrib(v1, v2, vPos);
		E_poly += getEdgeContrib(v2, v3, vPos);
		E_poly += getEdgeContrib(v3, v1, vPos);
	}
	vec3 G = lightIrradiance*albedo/PI/(2*PI)*E_poly;
	vec2 texCoords = vec2(gl_FragCoord.x/windowSize.x, 1- gl_FragCoord.y/windowSize.y);
	vec3 delta = texture(cvTerm, texCoords).xyz;
	color = G-delta;

}

