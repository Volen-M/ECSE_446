/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once
#include <random>
#include "bsdfs/phong.h"

TR_NAMESPACE_BEGIN

/**
 * Reflection occlusion integrator
 */
struct ROIntegrator : Integrator {

    float m_exponent;

    explicit ROIntegrator(const Scene& scene) : Integrator(scene) {
        m_exponent = scene.config.integratorSettings.ro.exponent;
    }

    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }


    v3f render(const Ray& ray, Sampler& sampler) const override {
		v3f Li(0.f);

		// TODO(A3): Implement this

		//Retrieve the light position p and intensity I.
		SurfaceInteraction info;

		//Intersect your ray with the scene geometry
		if (scene.bvh->intersect(ray, info)) {
			SurfaceInteraction shadowInfo;
			v2f sample = sampler.next2D();

			//Math.h stuff			
			v3f wi = Warp::squareToPhongLobe(sample, m_exponent);
			float pdf = Warp::squareToPhongLobePdf(wi, m_exponent);
			float cosAlpha = fmax(wi.z, 0);


			v3f wr = reflect(info.wo);
			wi = glm::toMat4(glm::quat(v3f(0, 0, 1), wr)) * v4f(wi, 1);
			wi = glm::normalize(wi); //Normalized on different line since it would darken it if not as it'd be normalized as v4f
			v3f wiWorld = glm::normalize(info.frameNs.toWorld(wi));

			//CosTheta
			float cosTheta = fmax(wi.z, 0.f);

			//Check for the shadow intersect
			Ray shadowRay = Ray(info.p + Epsilon, wiWorld);
			if (scene.bvh->intersect(shadowRay, shadowInfo)) {
				return Li;
			}
			else {
				Li += cosTheta * (m_exponent + 2) * INV_TWOPI * max(0.0f, pow(cosAlpha, m_exponent)) / pdf;
			}
		}

		return Li;
    }
};

TR_NAMESPACE_END