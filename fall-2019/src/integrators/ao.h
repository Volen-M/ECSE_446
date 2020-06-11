/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Ambient occlusion integrator
 */
struct AOIntegrator : Integrator {

	// Use this in your switch statement to select the sampling type 
	ESamplingType m_samplingStrategy;

    explicit AOIntegrator(const Scene& scene) : Integrator(scene) { 
		m_samplingStrategy = scene.config.integratorSettings.ao.sampling_type;
	}

    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f Li(0.f);
		
		/*
		Use the m_sampling_type variable to set wi and the corresponding pdf 
		appropriately for sphere, hemisphere, or cosine sampling.

		You can use a switch statement or an if/else block.

		The m_sampling_type variable is an enum. The different values of the enum 
		can be accessed through:
		ESamplingType::ESpherical
		ESamplingType::EHemispherical
		ESamplingType::ECosineHemispherical
		*/
		
        // TODO(A3): Implement this
	

		//Retrieve the light position p and intensity I.
		SurfaceInteraction info;
		SurfaceInteraction infoShadow;
		v3f shadowRayDirection;
		float divisor;
		p2f sample = sampler.next2D();

		//Intersect your ray with the scene geometry.
		if (!(scene.bvh->intersect(ray, info))) {
			return Li;
		}

		if (m_samplingStrategy == ESamplingType::ESpherical) {
			shadowRayDirection = Warp::squareToUniformSphere(sample);
			divisor = Warp::squareToUniformSpherePdf();
		}
		else if (m_samplingStrategy == ESamplingType::EHemispherical) {
			shadowRayDirection = Warp::squareToUniformHemisphere(sample);
			divisor = Warp::squareToUniformHemispherePdf(shadowRayDirection);
		}
		else if (m_samplingStrategy == ESamplingType::ECosineHemispherical) {
			shadowRayDirection = Warp::squareToCosineHemisphere(sample);
			divisor = Warp::squareToCosineHemispherePdf(shadowRayDirection);
		}

		Ray shadowRay = Ray(info.p + Epsilon, info.frameNs.toWorld(shadowRayDirection), Epsilon, scene.aabb.getBSphere().radius / 2.0f);
		//Check for the shadow intersect
		if (scene.bvh->intersect(shadowRay, infoShadow)) {
			return Li;
		}
		info.wi = shadowRayDirection;
		float cosTheta = Frame::cosTheta(shadowRayDirection);
		if (cosTheta > 0.0f) {
			Li = Li + 1.0f * cosTheta / divisor * INV_PI;
		}
		return Li;
    }
};

TR_NAMESPACE_END