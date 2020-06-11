/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Direct illumination integrator with MIS
 */
struct DirectIntegrator : Integrator {
    explicit DirectIntegrator(const Scene& scene) : Integrator(scene) {
        m_emitterSamples = scene.config.integratorSettings.di.emitterSamples;
        m_bsdfSamples = scene.config.integratorSettings.di.bsdfSamples;
        m_samplingStrategy = scene.config.integratorSettings.di.samplingStrategy;
    }

    static inline float balanceHeuristic(float nf, float fPdf, float ng, float gPdf) {
        float f = nf * fPdf, g = ng * gPdf;
        return f / (f + g);
    }

    void sampleSphereByCosineHemisphere(const p2f& sample,
                                        const v3f& n,
                                        const p3f& pShading,
                                        const v3f& emitterCenter,
                                        float emitterRadius,
                                        v3f& wiW,
                                        float& pdf) const {
        // TODO(A3): Implement this
		//sampling position y on the sphere
		v3f sampleY = Warp::squareToCosineHemisphere(sample) * emitterRadius;

		wiW = sampleY;

		//and the corresponding PDF
		//Why does manual divide by distance squared?
		pdf = Warp::squareToCosineHemispherePdf(sampleY);
    }

    void sampleSphereByArea(const p2f& sample,
                            const p3f& pShading,
                            const v3f& emitterCenter,
                            float emitterRadius,
                            v3f& pos,
                            v3f& ne,
                            v3f& wiW,
                            float& pdf) const {
        // TODO(A3): Implement this
		//sampling position y on the sphere
		v3f sampleY = Warp::squareToUniformSphere(sample) * emitterRadius; 

		//to set the direction wi=x->y
		pos = emitterCenter + sampleY;
		wiW = glm::normalize(pos - pShading);

		//the normal at y,
		ne = glm::normalize(sampleY);

		//and the corresponding PDF
		//Why does manual divide by distance squared?
		pdf = Warp::squareToUniformSpherePdf() / pow(emitterRadius, 2);  
    }

    void sampleSphereBySolidAngle(const p2f& sample,
                                  const p3f& pShading,
                                  const v3f& emitterCenter,
                                  float emitterRadius,
                                  v3f& wiW,
                                  float& pdf) const {
        // TODO(A3): Implement this
			//Manual inspired section 

		//IN manual called: "Computer Coordinate system for sphere sampling"
		v3f pCDir = glm::normalize(emitterCenter - pShading); //P to center norm dir
		Frame frame = Frame(pCDir);

		//Compute theta and phi values for sample in cone>
		float sinThetaMax2 = emitterRadius * emitterRadius / glm::distance2(pShading, emitterCenter);
		float cosThetaMax = std::sqrt(std::max((float)0, 1 - sinThetaMax2));
		float cosTheta = (1 - sample.x) + sample.x * cosThetaMax;
		float sinTheta = std::sqrt(std::max((float)0, 1 - cosTheta * cosTheta));
		float phi = sample.y * 2 * M_PI;

		//Compute angle alpha from center of sphere to sampled point on surface
		float dc = glm::distance(pShading, emitterCenter);
		float ds = dc * cosTheta - std::sqrt(std::max((float)0, emitterRadius * emitterRadius - dc * dc * sinTheta * sinTheta));
		float cosAlpha = (dc * dc + emitterRadius * emitterRadius - ds * ds) / (2 * dc * emitterRadius);
		float sinAlpha = std::sqrt(std::max((float)0, 1 - cosAlpha * cosAlpha));

		//Compute surface normal and sampled point on sphere
		//Vector3f nObj = SphericalDirection(sinAlpha, cosAlpha, phi,-wcX, -wcY, -wc);
		v3f nObj = v3f(sinAlpha * std::cos(phi), sinAlpha * std::sin(phi), cosAlpha);
		v3f pObj = emitterRadius * v3f(nObj.x, nObj.y, nObj.z);

		//Can use instead math.h warp but it works here
		wiW = glm::normalize(frame.toWorld(Warp::squareToUniformCone(sample, cosThetaMax)));
		pdf = Warp::squareToUniformConePdf(cosThetaMax); 
    }

    v3f renderArea(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A3): Implement this
		SurfaceInteraction info;
		SurfaceInteraction infoShadow;
		if (scene.bvh->intersect(ray, info)) {


			if (getEmission(info) != v3f(0.f, 0.f, 0.f)) {

				Lr += getEmitterByID(getEmitterIDByShapeID(info.shapeID)).getRadiance();
				/**/
			}
			else {

				for (int i = 0; i < m_emitterSamples; i++) {

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
					v3f pos;//not set here
					v3f ne;//not set here
					v3f wiW;//not set here
					float pdf;//not set here

					//sample
					sampleSphereByArea(sample, pShading, emitterCenter, emitterRadius, pos, ne, wiW, pdf);

					//update Wi for the rest
					info.wi = glm::normalize(info.frameNs.toLocal(wiW));
					//info.wi = wiW;

					//v3f multiplier = getBSDF(info)->sample(info, sampler.next2D(), &pdf);
					Ray shadowRay = Ray(info.p + Epsilon, info.frameNs.toWorld(info.wi), Epsilon);


					if (scene.bvh->intersect(shadowRay, infoShadow)) {
						//Lr =sum of Li*fr*G*visibility /pdf/samples

						//Jacobian (G) without ||x-y||^2
						float G = glm::normalizeDot(wiW, ne);
						//float G = glm::dot(info.wi, ne);

						//Full G with visibility
						G = fmax(G, 0.f) / glm::distance2(pos, info.p);
									

						//Li can be replaced by Le
						v3f Le = getEmission(infoShadow);

						//fr (brdf)
						v3f fr = getBSDF(info)->eval(info);

						Lr += Le * fr * G / pdf / (float)m_emitterSamples / emPdf;

					}
				}
			}

		}
		return Lr;

    }

    v3f renderCosineHemisphere(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);
		SurfaceInteraction info;
		SurfaceInteraction infoShadow;
		if (scene.bvh->intersect(ray, info)) {


			if (getEmission(info) != v3f(0.f, 0.f, 0.f)) {

				Lr += getEmitterByID(getEmitterIDByShapeID(info.shapeID)).getRadiance();
				/**/
			}
			else {

				for (int i = 0; i < m_bsdfSamples; i++) {

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
					v3f pos;//not set here
					v3f ne;//not set here
					v3f wiW;//not set here
					float pdf;//not set here

					//sample
					sampleSphereByCosineHemisphere(sample, info.frameNs.n, pShading, emitterCenter, emitterRadius, wiW, pdf);

					//update Wi for the rest
					info.wi = glm::normalize(wiW);
					//info.wi = wiW;

					//v3f multiplier = getBSDF(info)->sample(info, sampler.next2D(), &pdf);
					Ray shadowRay = Ray(info.p + Epsilon, info.frameNs.toWorld(wiW), Epsilon);


					if (scene.bvh->intersect(shadowRay, infoShadow)) {
						//Lr =sum of Li*fr*G*visibility /pdf/samples

						//Jacobian (G) without ||x-y||^2
						float G = glm::normalizeDot(wiW, ne);
						//float G = glm::dot(info.wi, ne);

						//Full G with visibility
						G = fmax(G, 0.f) / glm::distance2(pos, info.p);				
						G = 1.f;//G is not taken in consideration in this section (Not sure why)
						//Li can be replaced by Le
						v3f Le = getEmission(infoShadow);

						//fr (brdf)
						v3f fr = getBSDF(info)->eval(info);

						Lr += Le * fr * G / pdf / (float)m_bsdfSamples / emPdf;
						//Lr += Le * fr*G / pdf / (float)m_emitterSamples;

					}
				}
			}

		}
		return Lr;
    }

    v3f renderBSDF(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A3): Implement this
		SurfaceInteraction info;
		SurfaceInteraction infoShadow;
		float pdf;
		if (scene.bvh->intersect(ray, info)) {


			if (getBSDF(info)->isEmissive()) {

				Lr += getEmitterByID(getEmitterIDByShapeID(info.shapeID)).getRadiance();

			}
			else {

				for (int i = 0; i < m_bsdfSamples; i++) {
					p2f sample = sampler.next2D();
					v3f multiplier = getBSDF(info)->sample(info, sampler, &pdf);
					Ray shadowRay = Ray(info.p + Epsilon, info.frameNs.toWorld(info.wi), Epsilon);

					
					if (scene.bvh->intersect(shadowRay, infoShadow)) {
						
						Lr += getEmission(infoShadow) * multiplier / (float)m_bsdfSamples;

					}
				}
			}

		}
		return Lr;
    }

    v3f renderSolidAngle(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A3): Implement this
		SurfaceInteraction info;
		SurfaceInteraction infoShadow;
		if (scene.bvh->intersect(ray, info)) {


			if (getEmission(info) != v3f(0.f, 0.f, 0.f)) {

				Lr += getEmitterByID(getEmitterIDByShapeID(info.shapeID)).getRadiance();
				/**/
			}
			else {

				for (int i = 0; i < m_emitterSamples; i++) {

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
					v3f pos;//not set here
					v3f ne;//not set here
					v3f wiW;//not set here
					float pdf;//not set here

					//sample
					sampleSphereBySolidAngle(sample, pShading, emitterCenter, emitterRadius, wiW, pdf);

					//update Wi for the rest
					info.wi = glm::normalize(info.frameNs.toLocal(wiW));
					//info.wi = wiW;

					//v3f multiplier = getBSDF(info)->sample(info, sampler.next2D(), &pdf);
					Ray shadowRay = Ray(info.p + Epsilon, info.frameNs.toWorld(info.wi), Epsilon);


					if (scene.bvh->intersect(shadowRay, infoShadow)) {
						//Lr =sum of Li*fr*G*visibility /pdf/samples

						//Jacobian (G) without ||x-y||^2
						float G = glm::normalizeDot(wiW, ne);
						//float G = glm::dot(info.wi, ne);

						//Full G with visibility
						G = fmax(G, 0.f) / glm::distance2(pos, info.p);

						//Li can be replaced by Le
						v3f Le = getEmission(infoShadow);

						//fr (brdf)
						v3f fr = getBSDF(info)->eval(info);

						Lr += Le * fr / pdf / (float)m_emitterSamples / emPdf;
						//Lr += Le * fr*G / pdf / (float)m_emitterSamples;

					}
				}
			}

		}
		return Lr;
    }

    v3f renderMIS(const Ray& ray, Sampler& sampler) const {

        v3f Lr(0.f);

        // TODO(A4): Implement this
		SurfaceInteraction info;
		float pdf1;
		float pdf2;

		if (scene.bvh->intersect(ray, info)) {



			if (getEmission(info) != v3f(0.f, 0.f, 0.f)) {

				Lr += getEmitterByID(getEmitterIDByShapeID(info.shapeID)).getRadiance();
				/**/
			}
			else {



				//SAME AS BEFORE BUT WITH WEIGHT
				for (int i = 0; i < m_emitterSamples; i++) {

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
					v3f pos;//not set here
					v3f ne;//not set here
					v3f wiW;//not set here

					//sample
					sampleSphereBySolidAngle(sample, pShading, emitterCenter, emitterRadius, wiW, pdf1);


					//update Wi for the rest
					info.wi = glm::normalize(info.frameNs.toLocal(wiW));
					//info.wi = wiW;

					//v3f multiplier = getBSDF(info)->sample(info, sampler.next2D(), &pdf);
					Ray shadowRay = Ray(info.p + Epsilon, info.frameNs.toWorld(info.wi), Epsilon);
					SurfaceInteraction infoShadow;


					if (scene.bvh->intersect(shadowRay, infoShadow)) {

						

						//Li can be replaced by Le
						v3f Le = getEmission(infoShadow);
						if (Le != v3f(0.f)) {


							//fr (brdf)
							v3f fr = getBSDF(info)->eval(info);


							float fn = m_emitterSamples;
							float fPdf = pdf1 * getEmitterPdf(em);
							float gn = m_bsdfSamples;
							float gPdf = getBSDF(info)->pdf(info);


							float weightEmitter = balanceHeuristic(fn, fPdf, gn, gPdf);

							Lr += weightEmitter * Le * fr / pdf1 / (float)m_emitterSamples / emPdf;
						}
					}
				}
				
				for (int i = 0; i < m_bsdfSamples; i++) {

					if (scene.bvh->intersect(ray, info)) {
						p2f sample = sampler.next2D();
						v3f multiplier = getBSDF(info)->sample(info, sampler, &pdf2);
						Ray shadowRay = Ray(info.p + Epsilon, glm::normalize(info.frameNs.toWorld(info.wi)), Epsilon);
						SurfaceInteraction infoShadow;


						if (scene.bvh->intersect(shadowRay, infoShadow)) {

							v3f emission = getEmission(infoShadow);
							if (emission != v3f(0.f)) {
								p3f pShading = info.p;
								size_t emID = getEmitterIDByShapeID(infoShadow.shapeID);
								const Emitter& em = getEmitterByID(emID);
								float emPdf = getEmitterPdf(em);
								v3f emitterCenter = scene.getShapeCenter(em.shapeID);
								float emitterRadius = scene.getShapeRadius(em.shapeID);

								float sinThetaMax = emitterRadius * emitterRadius / pow(glm::distance(pShading, emitterCenter),2);
								float cosThetaMax = std::sqrt(std::max((float)0, 1 - sinThetaMax));
								

								float fn = m_emitterSamples;
								float fPdf = Warp::squareToUniformConePdf(cosThetaMax) * emPdf;
								float gn = m_bsdfSamples;
								float gPdf = pdf2;

								float weightBSDF = balanceHeuristic(gn, gPdf, fn, fPdf);
								Lr += weightBSDF * getEmission(infoShadow) * multiplier / (float)m_bsdfSamples;
							}
						}
					}
				}
			}


		}
		return Lr;
    }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        if (m_samplingStrategy == ESamplingStrategy::EMIS)
            return this->renderMIS(ray, sampler);
        else if (m_samplingStrategy == ESamplingStrategy::EArea)
            return this->renderArea(ray, sampler);
        else if (m_samplingStrategy == ESamplingStrategy::ESolidAngle)
            return this->renderSolidAngle(ray, sampler);
        else if (m_samplingStrategy == ESamplingStrategy::ECosineHemisphere)
            return this->renderCosineHemisphere(ray, sampler);
        else
            return this->renderBSDF(ray, sampler);
    }

    size_t m_emitterSamples;     // Number of emitter samples
    size_t m_bsdfSamples;        // Number of BSDF samples
    ESamplingStrategy m_samplingStrategy;   // Sampling strategy to use
};

TR_NAMESPACE_END