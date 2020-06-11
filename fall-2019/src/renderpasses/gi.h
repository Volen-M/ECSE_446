/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <core/core.h>
#include "core/renderpass.h"
#include "tiny_obj_loader.h"
#include "integrators/path.h"

TR_NAMESPACE_BEGIN

/**
 * Global Illumination baking renderpass.
 */
struct GIPass : RenderPass {
    GLuint shader{0};

    GLuint modelMatUniform{0};
    GLuint viewMatUniform{0};
    GLuint projectionMatUniform{0};

    int m_samplePerVertex;

    std::unique_ptr<PathTracerIntegrator> m_ptIntegrator;

    explicit GIPass(const Scene& scene) : RenderPass(scene) {
        m_ptIntegrator = std::unique_ptr<PathTracerIntegrator>(new PathTracerIntegrator(scene));
        m_ptIntegrator->m_maxDepth = scene.config.integratorSettings.gi.maxDepth;
        m_ptIntegrator->m_rrProb = scene.config.integratorSettings.gi.rrProb;
        m_ptIntegrator->m_rrDepth = scene.config.integratorSettings.gi.rrDepth;
        m_samplePerVertex = scene.config.integratorSettings.gi.samplesByVertex;
    }

    virtual void buildVBO(size_t objectIdx) override {
        GLObject& obj = objects[objectIdx];

        // TODO(A5): Implement this
		obj.nVerts = scene.getObjectNbVertices(objectIdx);
		obj.vertices.resize(obj.nVerts * N_ATTR_PER_VERT); // 3 for pos, 3 for normals


		for (int i = 0; i < obj.nVerts; i++) {
			Sampler sampler = Sampler(260746982);

			size_t vertexIdx = i; //???? or +1
			v3f position = scene.getObjectVertexPosition(objectIdx, vertexIdx);
			v3f normal = scene.getObjectVertexNormal(objectIdx, vertexIdx);

			//Create a SurfaceInteraction and manually populate its attributes.
			//set an arbitrary ray direction (wo) and shift the shading point position by e along the normal to avoid self-intersections during integration.
			SurfaceInteraction info = SurfaceInteraction();
			info.p = position + normal * Epsilon;
			info.wo = v3f(0, 0, 1);
			info.shapeID = objectIdx;
			info.frameNg = Frame(normal);
			info.frameNs = Frame(normal);
			info.primID = scene.getPrimitiveID(vertexIdx);
			info.matID = scene.getMaterialID(objectIdx, info.primID);


			Ray ray = Ray(v3f(0.f), v3f(0.f));
			v3f colour(0.f);

			for (int j = 0; j < m_samplePerVertex; j++) {
				colour += m_ptIntegrator->renderExplicit(ray, sampler, info) / (float)m_samplePerVertex;
			}

			//Attributes
			obj.vertices[i * N_ATTR_PER_VERT + 0] = position.x;
			obj.vertices[i * N_ATTR_PER_VERT + 1] = position.y;
			obj.vertices[i * N_ATTR_PER_VERT + 2] = position.z;


			obj.vertices[i * N_ATTR_PER_VERT + 3] = colour.x;
			obj.vertices[i * N_ATTR_PER_VERT + 4] = colour.y;
			obj.vertices[i * N_ATTR_PER_VERT + 5] = colour.z;
			
		}

        // VBO
        glGenVertexArrays(1, &obj.vao);
        glBindVertexArray(obj.vao);

        glGenBuffers(1, &obj.vbo);
        glBindBuffer(GL_ARRAY_BUFFER, obj.vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     sizeof(GLfloat) * obj.nVerts * N_ATTR_PER_VERT,
                     (GLvoid*) (&obj.vertices[0]),
                     GL_STATIC_DRAW);
    }

    bool init(const Config& config) override {
        RenderPass::init(config);

        // Create shader
        GLuint vs = compileShader("gi.vs", GL_VERTEX_SHADER);
        GLuint fs = compileShader("gi.fs", GL_FRAGMENT_SHADER);
        shader = compileProgram(vs, fs);
        glDeleteShader(vs);
        glDeleteShader(fs);

        // Create uniforms
        modelMatUniform = GLuint(glGetUniformLocation(shader, "model"));
        viewMatUniform = GLuint(glGetUniformLocation(shader, "view"));
        projectionMatUniform = GLuint(glGetUniformLocation(shader, "projection"));

        // Create vertex buffers
        objects.resize(scene.worldData.shapes.size());
        for (size_t i = 0; i < objects.size(); i++) {
            buildVBO(i);
            buildVAO(i);
        }

        return true;
    }

    void cleanUp() override {
        // Delete vertex buffers
        for (size_t i = 0; i < objects.size(); i++) {
            glDeleteBuffers(1, &objects[i].vbo);
            glDeleteVertexArrays(1, &objects[i].vao);
        }

        RenderPass::cleanUp();
    }

    void render() override {
        glBindFramebuffer(GL_FRAMEBUFFER, postprocess_fboScreen);
        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);

        // TODO(A5): Implement this
		//A1 ripoff 
		glUseProgram(shader);

		// Update camera
		glm::mat4 model, view, projection;
		camera.Update();
		camera.GetMatricies(projection, view, model);

		// Pass uniforms
		glUniformMatrix4fv(modelMatUniform, 1, GL_FALSE, &(modelMat[0][0]));
		glUniformMatrix4fv(viewMatUniform, 1, GL_FALSE, &(view[0][0]));
		glUniformMatrix4fv(projectionMatUniform, 1, GL_FALSE, &(projection[0][0]));

		// Draw
		for (auto& object : objects) {

			glBindVertexArray(object.vao);
			glDrawArrays(GL_TRIANGLES, 0, object.nVerts);
			glBindVertexArray(0);
		}

		RenderPass::render();
    }

};

TR_NAMESPACE_END
