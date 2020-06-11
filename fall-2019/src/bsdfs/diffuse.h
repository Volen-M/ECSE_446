/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Perfectly diffuse, Lambertian reflectance model
 */
struct DiffuseBSDF : BSDF {
    std::unique_ptr<Texture < v3f>> albedo;

    DiffuseBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
        const tinyobj::material_t& mat = scene.materials[matID];

        if (mat.diffuse_texname.empty())
            albedo = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.diffuse)));
        else
            albedo = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.diffuse_texname));

        components.push_back(EDiffuseReflection);

        combinedType = 0;
        for (size_t i = 0; i < components.size(); ++i)
            combinedType |= components[i];
    }

    inline float getExponent(const SurfaceInteraction& i) const override { return 1.f; }

    v3f eval(const SurfaceInteraction& i) const override {
        v3f val(0.f);

		if (Frame::cosTheta(i.wo) < 0.0f || Frame::cosTheta(i.wi) < 0.0f) {
			v3f val(0.f);
			return val;
		}
		else {

			v3f rho = albedo->eval(worldData, i);
			v3f val = rho * INV_PI * Frame::cosTheta(i.wi);
			return val;
		}

        return val;
    }

    float pdf(const SurfaceInteraction& i) const override {
        float pdf = 0.f;

        // TODO(A3): Implement this
		pdf = Warp::squareToCosineHemispherePdf(i.wi);
		return pdf;
    }

    v3f sample(SurfaceInteraction& i, Sampler& sampler, float* pdf) const override {
        v3f val(0.f);

        // TODO(A3): Implement this
		v2f sample = sampler.next2D();
		val = Warp::squareToCosineHemisphere(sample);
        i.wi = val;
		*pdf = DiffuseBSDF::pdf(i);
		val = eval(i) / *pdf;
        return val;
    }

    std::string toString() const override { return "Diffuse"; }
};

TR_NAMESPACE_END