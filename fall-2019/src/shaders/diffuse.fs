/*
    This file is part of TinyRender, an educative PBR system.
    Designed for ECSE 446/546 Realistic Image Synthesis, McGill University.

    Copyright (c) 2018 by Derek Nowrouzezahrai and others.
*/
#version 330 core

#define M_PI       3.14159265358979323846f
#define INV_PI     0.31830988618379067154f

uniform vec3 albedo;
uniform vec3 lightInt;
uniform vec3 camPos;
uniform vec3 lightPos;

in vec3 vNormal;
in vec3 vPos;
out vec3 color;

void main() {
	vec3 ptToLight = lightPos - vPos;
	vec3 ptToCam = camPos - vPos;
	color = albedo * INV_PI * lightInt / pow(length(ptToLight), 2) * dot(normalize(vNormal), normalize(ptToLight));
}