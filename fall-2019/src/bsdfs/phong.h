/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Modified Phong reflectance model
 */
struct PhongBSDF : BSDF {

    std::unique_ptr<Texture < v3f>> specularReflectance;
    std::unique_ptr<Texture < v3f>> diffuseReflectance;
    std::unique_ptr<Texture < float>> exponent;
    float specularSamplingWeight;
    float scale;

    PhongBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
        const tinyobj::material_t& mat = scene.materials[matID];

        if (mat.specular_texname.empty())
            specularReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.specular)));
        else
            specularReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.specular_texname));

        if (mat.diffuse_texname.empty())
            diffuseReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.diffuse)));
        else
            diffuseReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.diffuse_texname));

        exponent = std::unique_ptr<Texture<float>>(new ConstantTexture1f(mat.shininess));

        //get scale value to ensure energy conservation
        v3f maxValue = specularReflectance->getMax() + diffuseReflectance->getMax();
        float actualMax = max(max(maxValue.x, maxValue.y), maxValue.z);
        scale = actualMax > 1.0f ? 0.99f * (1.0f / actualMax) : 1.0f;

        float dAvg = getLuminance(diffuseReflectance->getAverage() * scale);
        float sAvg = getLuminance(specularReflectance->getAverage() * scale);
        specularSamplingWeight = sAvg / (dAvg + sAvg);

        components.push_back(EGlossyReflection);
        components.push_back(EDiffuseReflection);

        combinedType = 0;
        for (unsigned int component : components)
            combinedType |= component;
    }

    inline float getExponent(const SurfaceInteraction& i) const override {
        return exponent->eval(worldData, i);
    }

    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }

    v3f eval(const SurfaceInteraction& i) const override {
        v3f val(0.f);

        // TODO(A2): Implement this
		if (Frame::cosTheta(i.wo) < 0.0f || Frame::cosTheta(i.wi) < 0.0f) {
			return val;
		}
		else {
			float angle = glm::dot(glm::normalize(PhongBSDF::reflect(i.wo)), i.wi);
			v3f rhoD = diffuseReflectance->eval(worldData, i);
			v3f rhoS = specularReflectance->eval(worldData, i);
			float exp = exponent->eval(worldData, i);
			float cosAlpha = fmax(angle, 0.0f);
			cosAlpha = fmin(cosAlpha, 1.f);
			val = rhoD * INV_PI + rhoS * (exp + 2.0f) * INV_TWOPI * pow(cosAlpha, exp);

			return val * glm::max(0.f, Frame::cosTheta(i.wi)) * scale;
		}
    }

    float pdf(const SurfaceInteraction& i) const override {
        float pdf = 0.f;

        // TODO(A3): Implement this
		v3f phongWi = glm::toMat4(glm::quat(reflect(i.wo), v3f(0, 0, 1))) * v4f(i.wi, 1);
		float pdfS = Warp::squareToPhongLobePdf(phongWi, exponent->eval(worldData, i));
		float pdfD = Warp::squareToCosineHemispherePdf(glm::normalize(i.wi));

		//PDF should implement a linear combination of glossy and diffuse BRDFs
		pdfS = fmax(pdfS, 0);
		pdfD = fmax(pdfD, 0);
		pdf = pdfS * specularSamplingWeight + pdfD * (1.f - specularSamplingWeight);
		return pdf;
    }

    v3f sample(SurfaceInteraction& i, Sampler& sampler, float* pdf) const override {
        v3f val(0.f);

        // TODO(A3): Implement this
		v2f sample = sampler.next2D();
		if (sample.x < (1.f - specularSamplingWeight)) {

			//Sample in between 0 and 1- specularSamplingWeight
			v2f sample2 = v2f(sample.x / (1.f - specularSamplingWeight), sample.y);

			//Diffuse Code 
			i.wi = Warp::squareToCosineHemisphere(sample2);

			*pdf = PhongBSDF::pdf(i);
			
			if (*pdf > 0.f) {
				val = eval(i) / *pdf;
			}
			else {
				val = v3f(0.f);
			}
			//End Diffuse
		}
		else {
			//Sample in between 1-specularSamplingWeight and 1
			v2f sample2 = v2f((sample.x - 1.f + specularSamplingWeight) / specularSamplingWeight, sample.y);
			float exp = exponent->eval(worldData, i);

			//Phong Code
			v3f tempWi = Warp::squareToPhongLobe(sample2, exp);
			i.wi = glm::toMat4(glm::quat(v3f(0, 0, 1), reflect(i.wo))) * v4f(tempWi, 1);
			*pdf = PhongBSDF::pdf(i);
			if (*pdf > 0.f) {
				val = eval(i) / *pdf;
			}
			else {
				val = v3f(0.f);
			}
			//End Phong
		}
		return val;
    }

    std::string toString() const override { return "Phong"; }
};

TR_NAMESPACE_END