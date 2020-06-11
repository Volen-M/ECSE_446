/*
    This file is part of TinyRender, an educative PBR system.
    Designed for ECSE 446/546 Realistic Image Synthesis, McGill University.

    Copyright (c) 2018 by Derek Nowrouzezahrai and others.
*/
#version 330 core

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform mat4 normalMat;

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
out vec3 vPos;
out vec3 vNormal;

void main() {
	vec4 position4 = projection*view*model*vec4(position, 1.0);
	vPos = vec3(model*vec4(position, 1.0));
	gl_Position = position4;
	vNormal =   vec3((normalMat * vec4(normal, 1.0)));
}