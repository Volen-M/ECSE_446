/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <tiny_obj_loader.h>
#define RAY_EPS_CV 1e-5 // Use when setting min and max dist for ray in control variates code
TR_NAMESPACE_BEGIN

/**
 * Direct illumination integrator for polygonal light sources
 * Follows Arvo '94.
 */
struct PolygonalIntegrator : Integrator {

	float m_alpha;             // Control variates "strength"
	size_t m_visSamples;       // # of samples to estimate h - alpha*g
	bool m_traceShadows;       // Trace shadows or not
	EPolygonalMethod m_method; // Method to use (Arvo, or control variates)

	std::vector<std::vector<v3f>> m_triangles; // Data structure to store triangles

    explicit PolygonalIntegrator(const Scene& scene) : Integrator(scene) {
        m_alpha = scene.config.integratorSettings.poly.alpha;
        m_visSamples = scene.config.integratorSettings.poly.visSamples;
        m_traceShadows = scene.config.integratorSettings.poly.traceShadows;
        m_method = scene.config.integratorSettings.poly.method;

		/**
		 * 1) Get # of triangles on emitter
		 * 2) Store vertices in m_triangles
		 */
// TODO(A4): Implement this
		size_t shapeID = scene.getFirstLight();
		auto shape = scene.worldData.shapes[shapeID];
		int triangleAmt = shape.mesh.indices.size()/3;
		//printf("Indices: %d", shape.mesh.indices.size());

		
		for (int i = 0; i < triangleAmt; i++) {
			std::vector<v3f> tri(3); // Create empty triangle (3 vertices / triangle)
			tri[0] = scene.getObjectVertexPosition(shapeID, 0+i*3); //add vertices to triangle
			tri[1] = scene.getObjectVertexPosition(shapeID, 1+i*3); //add vertices to triangle
			tri[2] = scene.getObjectVertexPosition(shapeID, 2+i*3); //add vertices to triangle
			m_triangles.push_back(tri); // Append triangle vector to vector
		}
		

    }

    /// Reflect
    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }

    /**
     * === PHONG BONUS ONLY ===
     * Compute the following integral:
     *    T(a, b, n, x) = \int_0^x [a \cos(\theta) + b \sin(\theta)]ˆn d\theta
     * Uses a recurrent relation (see Snyder's note, 1996)
     *
     * Series function:
     *    T_sum(a, b, n, x) = \sum_{i=0}ˆ{(n-1)/2} T(a, b, 2i+1, x)
     * assuming n is _odd_
     */
    float cosineSinePowerIntegralSum(float a, float b, int exp, float theta) const {
        if (exp % 2 == 0) exp += 1; // Make exponent odd
        float Tsum = 0.f;

		// Implementing this function may be useful if you attempt the bonus

        // TODO(A4): Implement this

        return Tsum;
    }

    /**
     * Compute edge (v1--v2) contribution
	 * The exp term is only needed if you attempt the bonus, otherwise, you can ignore it
     */
    float getEdgeContrib(const v3f& v1, const v3f& v2, const SurfaceInteraction& i, int exp = 0) const {
        float contrib = 0.f;

        // TODO(A4): Implement this
		
		float dotSolution = glm::dot(((v1 - i.p) / glm::length(v1 - i.p)), ((v2 - i.p) / glm::length(v2 - i.p)));
		float theta = acosf(dotSolution);
		v3f gamma = glm::cross(v2-i.p,v1-i.p) / glm::length(glm::cross(v2-i.p,v1-i.p)); //if there's an error it's here, not sure about the magniture thingy 
		contrib = glm::dot(gamma, i.frameNs.n) * theta;
		return contrib;
    }
	   

    /// Direct illumination using Arvo '94 analytic solution for polygonal lights
    v3f renderAnalytic(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A4): Implement this
		//	1. intersect your eye rays with the scene geometry
		SurfaceInteraction info;
		float E_poly = 0.f;
		if (scene.bvh->intersect(ray, info)) {
			if (getEmission(info) != v3f(0.f, 0.f, 0.f)) {

				Lr += getEmitterByID(getEmitterIDByShapeID(info.shapeID)).getRadiance();
				
			}
			else {
				//	2. if there's an intersection i, loop through the vector of triangles and compute Equation by summing the
				//	contribution of each subtended edge(implement the getEdgeContrib() function)
				for (int i = 0; i < m_triangles.size(); i++) { 
					v3f v1 = m_triangles.at(i)[0];
					v3f v2 = m_triangles.at(i)[1];
					v3f v3 = m_triangles.at(i)[2];
					E_poly += getEdgeContrib(v1, v2, info);
					E_poly += getEdgeContrib(v2, v3, info);
					E_poly += getEdgeContrib(v3, v1, info);
				}
				//Intesity divided by Area
				//	4. compute and return the diffuse outgoing radiance at .
				float emPDF;
				size_t id = scene.getFirstLight();
				size_t emId = getEmitterIDByShapeID(id);
				const Emitter& em = getEmitterByID(emId);
				v3f E = em.getPower() / em.area;

				//	3. retrieve the albedo at the surface point and compute the light's incident radiance
				info.wi = v3f(0, 0, 1);
				v3f albedoOverPi = getBSDF(info)->eval(info);
				//albedoOverPi = v3f(1.f);
				//E= v3f(1.f);

				Lr = E * albedoOverPi * E_poly * INV_TWOPI ;

			}
		}
		
        return clampBelow(Lr, 0.f);
    }

    /**
     * Stand-alone estimator for h - alpha*g (with primary ray)
     * Trace a primary ray, check for emitter hit, and then call `estimateVisDiff()`
     * Used by polygonal render pass
     */
    v3f estimateVisDiffRealTime(const Ray& ray, Sampler& sampler, const Emitter& em) {
        v3f D(0.f);

        SurfaceInteraction hit;
        if (!scene.bvh->intersect(ray, hit)) return D;

        const BSDF* bsdf = getBSDF(hit);
        if (bsdf->isEmissive()) return D;

        hit.wi = v3f(0, 0, 1); // Trick to get 1/pi * albedo without cosine term
        D = estimateVisDiff(sampler, hit, em);

        return D;
    }

    /// Stand-alone estimator for h - alpha*g (without primary ray)
	/// Use RAY_EPS_CV when setting min and max dist for shadow ray
    v3f estimateVisDiff(Sampler& sampler, SurfaceInteraction& i, const Emitter& em) const {
        v3f sum(0.f);
		v3f h(0.f);
		v3f g(0.f);
        // TODO(A4): Implement this
		for (int j = 0; j < m_visSamples; j++) {

			//Begin handout code
			float emPdf;
			size_t id = selectEmitter(sampler.next(), emPdf);
			const Emitter& em = getEmitterByID(id);
			//end handout code

			//argument initialization
			p2f sample = sampler.next2D();
			p3f pShading = i.p;
			v3f emitterCenter = scene.getShapeCenter(em.shapeID);
			float emitterRadius = scene.getShapeRadius(em.shapeID);
			v3f wiW;//not set here

			//sample
			v3f pe; // Point on emitter
			v3f ne; // Surface normal at point
			float pdf; // PDF of choosing point
			sampleEmitterPosition(sampler, em, ne, pe, pdf); // Sample mesh uniformly
			wiW = glm::normalize(pe - pShading);

			//update Wi for the rest
			i.wi = glm::normalize(i.frameNs.toLocal(wiW));
			//info.wi = wiW;

			//v3f multiplier = getBSDF(info)->sample(info, sampler.next2D(), &pdf);
			Ray shadowRay = Ray(i.p + Epsilon, i.frameNs.toWorld(i.wi), Epsilon);
			SurfaceInteraction infoShadow;


			v3f fr(0.f);
			float G = 0.f;
			v3f Le(0.f);

			//-------- THIS IS H------------------
			if (scene.bvh->intersect(shadowRay, infoShadow)) {



				//Li can be replaced by Le
				Le = getEmission(infoShadow);
				if (Le != v3f(0.f)) {

					//printf("It reaches here");
					G = -1.f * glm::normalizeDot(wiW, ne);
					//float G = glm::dot(info.wi, ne);

					//Full G with visibility
					G = fmax(G, 0.f) / glm::distance2(pe, i.p);
					//fr (brdf)
					fr = getBSDF(i)->eval(i);

					h += Le * G * fr / pdf / emPdf;
				}
			}
			//-------------------------------------

			//-------- THIS IS G ---------------- but written as "g" here
			G = -1.f * glm::normalizeDot(wiW, ne);
			//float G = glm::dot(info.wi, ne);

			//Full G with visibility
			G = fmax(G, 0.f) / glm::distance2(pe, i.p);

			fr = getBSDF(i)->eval(i);
			Le = em.getRadiance();
			g += m_alpha*Le * G * fr / pdf / emPdf;
			//-------------------------------------
		}

        return sum = (h - g)/m_visSamples;
    }

    /// Control variates using Arvo '94 for direct illumination; ray trace shadows
	
    v3f renderControlVariates(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A4): Implement this
		v3f aG = m_alpha * renderAnalytic(ray, sampler);
		v3f h_minus_aG = v3f(0.f);
		SurfaceInteraction info;
		size_t id = scene.getFirstLight();
		size_t emId = getEmitterIDByShapeID(id);
		const Emitter& em = getEmitterByID(emId);
		if (scene.bvh->intersect(ray, info)) {
			if (getEmission(info) != v3f(0.f, 0.f, 0.f)) {

				Lr += getEmitterByID(getEmitterIDByShapeID(info.shapeID)).getRadiance();
			}
			else {
				h_minus_aG = estimateVisDiff(sampler, info, em);
			}
		}
        return Lr = clampBelow(aG + h_minus_aG, 0.f);
    }

    /// Direct illumination using surface area sampling
    v3f renderArea(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A4): Implement this
		SurfaceInteraction info;
		SurfaceInteraction infoShadow;
		if (scene.bvh->intersect(ray, info)) {


			if (getEmission(info) != v3f(0.f, 0.f, 0.f)) {

				Lr += getEmitterByID(getEmitterIDByShapeID(info.shapeID)).getRadiance();
				/**/
			}
			else {


				//Begin handout code
				float emPdf;
				size_t id = selectEmitter(sampler.next(), emPdf);
				const Emitter& em = getEmitterByID(id);
				//end handout code

				//argument initialization
				p2f sample = sampler.next2D();
				p3f pShading = info.p;
				v3f emitterCenter = scene.getShapeCenter(em.shapeID);
				float emitterRadius = scene.getShapeRadius(em.shapeID);
				v3f wiW;//not set here

				//sample
				v3f pe; // Point on emitter
				v3f ne; // Surface normal at point
				float pdf; // PDF of choosing point
				sampleEmitterPosition(sampler, em, ne, pe, pdf); // Sample mesh uniformly
				wiW = glm::normalize(pe - pShading);

				//update Wi for the rest
				info.wi = glm::normalize(info.frameNs.toLocal(wiW));
				//info.wi = wiW;

				//v3f multiplier = getBSDF(info)->sample(info, sampler.next2D(), &pdf);
				Ray shadowRay = Ray(info.p + Epsilon, info.frameNs.toWorld(info.wi), Epsilon);
				SurfaceInteraction infoShadow;

				if (m_traceShadows) {
					if (scene.bvh->intersect(shadowRay, infoShadow)) {


						//Li can be replaced by Le
						v3f Le = getEmission(infoShadow);
						if (Le != v3f(0.f)) {

							float G = -1.f * glm::normalizeDot(wiW, ne);
							//float G = glm::dot(info.wi, ne);

							//Full G with visibility
							G = fmax(G, 0.f) / glm::distance2(pe, info.p);
							//fr (brdf)
							v3f fr = getBSDF(info)->eval(info);

							Lr += Le * G * fr / pdf / emPdf;
						}
					}
				}
				else {
					float G = -1.f * glm::normalizeDot(wiW, ne);
					//float G = glm::dot(info.wi, ne);

					//Full G with visibility
					G = fmax(G, 0.f) / glm::distance2(pe, info.p);
					
					v3f fr = getBSDF(info)->eval(info);
					v3f Le = em.getRadiance();
					Lr +=  Le* G * fr / pdf / emPdf;
				}
			}

		}
		return Lr;
    }

    /// Branch to corresponding method
    v3f render(const Ray& ray, Sampler& sampler) const override {
        switch (m_method) {
            case EPolygonalMethod::ESurfaceArea:
                return PolygonalIntegrator::renderArea(ray, sampler);
                break;
            case EPolygonalMethod::EControlVariates:
                return PolygonalIntegrator::renderControlVariates(ray, sampler);
                break;
            default:
                return PolygonalIntegrator::renderAnalytic(ray, sampler);
                break;
        }
    }

};

TR_NAMESPACE_END