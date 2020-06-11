/*
    This file is part of TinyRender, an educative PBR system.
    Designed for ECSE 446/546 Realistic Image Synthesis, McGill University.

    Copyright (c) 2018 by Derek Nowrouzezahrai and others.
*/
#version 330 core

#define M_PI       3.14159265358979323846f
#define INV_PI     0.31830988618379067154f
#define INV_TWOPI  0.15915494309189533577f

uniform vec3 camPos;
uniform vec3 lightPos;
uniform vec3 lightInt;
uniform vec3 rho_d;
uniform vec3 rho_s;
uniform float exponent;

in vec3 vNormal;
in vec3 vPos;
out vec3 color;

void main() {
	vec3 ptToLight = lightPos - vPos;
	vec3 ptToCam = camPos - vPos;
    vec3 reflection = reflect((ptToLight), (vNormal));
    float specularAngle = dot(normalize(reflection), normalize(ptToCam));

	vec3 phongBSDFEval = rho_d * INV_PI + rho_s * (exponent + 2.0f) * INV_TWOPI * pow(specularAngle,exponent);
	color = phongBSDFEval * dot(normalize(ptToLight), normalize(vNormal)) * lightInt/ pow(length(ptToLight),2);
	
}