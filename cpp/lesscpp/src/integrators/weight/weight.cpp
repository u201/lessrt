#include <mitsuba/render/scene.h>
//#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

class Weight : public SamplingIntegrator {
public:
	Weight(const Properties &props) : SamplingIntegrator(props) {
		Spectrum defaultColor;
		defaultColor.fromLinearRGB(0.2f, 0.5f, 0.2f);
		m_color = props.getSpectrum("color", defaultColor);
	}
	
	/// Unserialize from a binary data stream
	Weight(Stream *stream, InstanceManager *manager) : SamplingIntegrator(stream, manager) {
		m_color = Spectrum(stream);
		m_maxDist = stream->readFloat();
	}

	/// Serialize to a binary data stream
	void serialize(Stream *stream, InstanceManager *manager) const {
		SamplingIntegrator::serialize(stream, manager);
		m_color.serialize(stream);
		stream->writeFloat(m_maxDist);
	}

	/// Query for an unbiased estimate of the radiance along <tt>r</tt>
	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
		if (rRec.rayIntersect(r)) {
			//Float distance = rRec.its.t;
			//return Spectrum(1.0f - distance / m_maxDist) * m_color;
			//std::cout << rRec.its.toString() << std::endl;
			const BSDF *bsdf = rRec.its.getBSDF();
			Intersection its_tmp;
			its_tmp.p = rRec.its.p;
			BSDFSamplingRecord bRecref(its_tmp, Vector(0, 0, 1), Vector(0, 0, 1));
			Spectrum r = bsdf->eval(bRecref);

			//std::cout << "samplingRecord: " << r.toString() << std::endl;
			//std::cout << r.toString() << std::endl;
			
			//Float cosTheta = absDot(rRec.its.wi, rRec.its.shFrame.n);
			Float cosTheta = absDot(normalize(cameraPosition - rRec.its.p), rRec.its.shFrame.n);
			//std::cout << "camera position" << std::endl;
			//std::cout << cameraPosition.toString() << std::endl;
			//std::cout << "intersection point" << std::endl;
			//std::cout << rRec.its.toString() << std::endl;

			r = bsdf->getDiffuseReflectance(rRec.its) / M_PI_DBL * cosTheta;
			//std::cout << bsdf->getDiffuseReflectance(rRec.its).toString() << std::endl;
			//std::cout << "cosTheta = " << cosTheta << std::endl;
			//Spectrum r = Spectrum(0.5) / M_PI_DBL * cosTheta;
			//std::cout << bsdf->getDiffuseReflectance(rRec.its).toString() << std::endl;
			//std::cout << r.toString() << std::endl;	
			return r;
		}
		return Spectrum(0.0f);
	}

	/// Preprocess function -- called on the initiating machine
	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job, int sceneResID, int cameraResID, int samplerResID) {
		SamplingIntegrator::preprocess(scene, queue, job, sceneResID, cameraResID, samplerResID);

		std::cout << scene->toString() << std::endl;

		const AABB &sceneAABB = scene->getAABB();

		//Point cameraPosition = scene->getSensor()->getWorldTransform()->eval(0).transformAffine(Point(0.0f));
		cameraPosition = scene->getSensor()->getWorldTransform()->eval(0).transformAffine(Point(0.0f));
		m_maxDist = -std::numeric_limits<Float>::infinity();

		for (int i = 0; i < 8; ++i) {
			m_maxDist = std::max(m_maxDist, (cameraPosition - sceneAABB.getCorner(i)).length());
		}

		return true;
	}

	MTS_DECLARE_CLASS()

private:
	Spectrum m_color;
	Float m_maxDist;
	Point cameraPosition;
};

MTS_IMPLEMENT_CLASS_S(Weight, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(Weight, "first order scattering receive")
MTS_NAMESPACE_END