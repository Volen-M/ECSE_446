/*
	This file is part of TinyRender, an educative rendering system.

	Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
	Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Path tracer integrator
 */
	struct PathTracerIntegrator : Integrator {
	explicit PathTracerIntegrator(const Scene& scene) : Integrator(scene) {
		m_isExplicit = scene.config.integratorSettings.pt.isExplicit;
		m_maxDepth = scene.config.integratorSettings.pt.maxDepth;
		m_rrDepth = scene.config.integratorSettings.pt.rrDepth;
		m_rrProb = scene.config.integratorSettings.pt.rrProb;
	}


	v3f renderImplicit(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit) const {
		v3f Li(0.f);

		// TODO(A5): Implement this
		float pdf;
		bool lightHit = false;

		Li = v3f(1.0f); //Multiplication is applied to this
		for (int depthCounter = 0; depthCounter <= m_maxDepth; depthCounter++) {
			if (getEmission(hit) == v3f(0.f)) {
				v3f colour = getBSDF(hit)->sample(hit, sampler, &pdf);

				Li = Li * colour; //Colour Reflection
				if (depthCounter + 1 <= m_maxDepth) {
					v3f direction = glm::normalize(hit.frameNs.toWorld(hit.wi));
					Ray newRay = Ray(hit.p + Epsilon, direction, Epsilon);

					if (scene.bvh->intersect(newRay, hit)) {
						continue;
					}
					else {
						break;
					}
				}
			}
			else {
				if (hit.wo.z > 0) {
					Li = Li * getEmission(hit);
					lightHit = true;
					break;
				}
				else {
					Li = v3f(0.f);
					break;
				}
			}
		}

		if (!lightHit) {
			Li = v3f(0.f);
		}

		return Li;

	}

	v3f renderExplicit(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit) const {
		v3f Li(0.f);
		v3f Ldir(0.f);
		v3f Lind(0.f);

		// TODO(A5): Implement this
		if (getEmission(hit) != v3f(0.f)) {
			Ldir = v3f(1.f) * getEmission(hit);
			Lind = v3f(0.f);
		}
		else if (m_maxDepth == 0) {
			Ldir = v3f(0.f);
			Lind = v3f(0.f);
		}
		else if (m_maxDepth == 1) {
			//Do only direct illumination
			Ldir = getLdir(ray, sampler, hit);
			Lind = v3f(0.f);
		}
		else {
			//Do both direct and indirect illumination
			Ldir = getLdir(ray, sampler, hit);
			Lind = getLind(ray, sampler, hit, 1);
		}

		/*if (std::isnan(Lind) == true) {
			Lind = v3f(0.f);
		}
		else if (std::isnan(Ldir) == true ) {
			Ldir = v3f(0.f);
		}*/
		
		Li = Ldir + Lind;
		return Li;
	}

	v3f render(const Ray& ray, Sampler& sampler) const override {
		Ray r = ray;
		SurfaceInteraction hit;

		if (scene.bvh->intersect(r, hit)) {
			if (m_isExplicit)
				return this->renderExplicit(ray, sampler, hit);
			else
				return this->renderImplicit(ray, sampler, hit);
		}
		return v3f(0.0);
	}

	v3f getLdir(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit) const {
		v3f Lr(0.f);
		//-----------------------A4 Code----------------------------

		SurfaceInteraction infoShadow;
		//Begin handout code
		float emPdf;
		size_t id = selectEmitter(sampler.next(), emPdf);
		const Emitter& em = getEmitterByID(id);
		//end handout code

		v3f pos;//not set here
		v3f ne;//not set here
		v3f wiW;//not set here
		float pdf;//not set here
		//--------------------End A4 Code--------------------------

		//Instead of sampleSpherebyArea
		sampleEmitterPosition(sampler, em, ne, pos, pdf);
		//Missing wiW
		wiW = glm::normalize(pos - hit.p);


		//-------------------A4 Code 2---------------------------------
		//update Wi for the rest
		hit.wi = glm::normalize(hit.frameNs.toLocal(wiW));
		Ray shadowRay = Ray(hit.p + Epsilon, wiW + Epsilon);
		if (glm::dot(wiW, ne) < 0.f) {
			if (scene.bvh->intersect(shadowRay, infoShadow)) {

				//Why does this work with -1
				float G = -1 * glm::normalizeDot(wiW, ne);
				//float G = glm::dot(info.wi, ne);

				//Full G with visibility
				G = fmax(G, 0.f) / glm::distance2(pos, hit.p);

				//Li can be replaced by Le
				v3f Le = getEmission(infoShadow);

				//fr (brdf)
				v3f fr = getBSDF(hit)->eval(hit);
				Lr += Le * fr * G / pdf / emPdf;
				//Lr += Le * fr*G / pdf / (float)m_emitterSamples / emPdf;

			}
		}
		//-------------------End A4 Code 2--------------------------
		return Lr;
	}

	v3f getLind(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit, int depth) const {
		
		//-----  Variable init for RR
		float rrPdf = 1.f; //Needs to be default 1.f since m_rrProb is sometimes undefined
		bool isRR = (m_maxDepth == -1);
		//-----

		//----- End Conditions Check and RR Prob check
		if (isRR && depth >= m_rrDepth) {
				rrPdf = m_rrProb;
				if (sampler.next() > m_rrProb) {
					return v3f(0.f);
				}
		}
		else if (!isRR && depth >= m_maxDepth) {
				return v3f(0.f);
		}
		//------

		//----- Variable init
		v3f colour(1.f);
		v3f emission = v3f(0.f);
		SurfaceInteraction newHit;
		Ray newRay = Ray(v3f(0.f),v3f(0.f));
		//-----

		do {
			float pdf;
			colour = getBSDF(hit)->sample(hit, sampler, &pdf);
			v3f direction = hit.frameNs.toWorld(hit.wi);
			direction = glm::normalize(direction);
			newRay = Ray(hit.p, direction, Epsilon);

			if (scene.bvh->intersect(newRay, newHit)) {
				emission = getEmission(newHit);
			}
			else {
				return v3f(0.f);
			}
		} while (emission != v3f(0.f));

		return (1.f / rrPdf) * colour * ( getLdir(newRay, sampler, newHit) + getLind(newRay, sampler, newHit, depth + 1) );
	}

	int m_maxDepth;     // Maximum number of bounces
	int m_rrDepth;      // When to start Russian roulette
	float m_rrProb;     // Russian roulette probability
	bool m_isExplicit;  // Implicit or explicit
};

TR_NAMESPACE_END
