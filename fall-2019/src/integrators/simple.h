/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Simple direct illumination integrator.
 */
struct SimpleIntegrator : Integrator {
    explicit SimpleIntegrator(const Scene& scene) : Integrator(scene) { }

    v3f render(const Ray& ray, Sampler& sampler) const override {
		v3f Li(0.f);

		// TODO(A2): Implement this
		//Retrieve the light position p and intensity I.
		v3f position = scene.getFirstLightPosition();
		v3f intensity = scene.getFirstLightIntensity();
		SurfaceInteraction info;
		SurfaceInteraction infoShadow;


		//Intersect your ray with the scene geometry.
		if (!(scene.bvh->intersect(ray, info))) {
			return Li;
		}
		v3f shadowRayDirection = glm::normalize(position - info.p);

		Ray shadowRay = Ray(info.p, shadowRayDirection, Epsilon, glm::distance(info.p, position));

		//Check for the shadow intersect
		if (scene.bvh->intersect(shadowRay, infoShadow)) {
			return Li;
		}

		//If an intersection i is found, retrieve the hit surface material using getBSDF(i).
		const BSDF* bsdf = getBSDF(info);
		//Map the incoming direction i.wi to local coordinates by using the hit point frameNs.toLocal() transform.
		info.wi = info.frameNs.toLocal(position - info.p);
		info.wi = glm::normalize(info.wi);

		//Evaluate the BRDF locally and the light, and set the Li term accordingly.
		v3f lightIntensity = bsdf->eval(info);
		v3f multiplier = intensity / glm::length2(info.p - position);

		Li = multiplier * lightIntensity;
		return Li;
    }
};

TR_NAMESPACE_END