#include <mitsuba/render/scene.h>
//#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

class MyIntegrator : public SamplingIntegrator {
public:
	MyIntegrator(const Properties &props) : SamplingIntegrator(props) {
		Spectrum defaultColor;
		defaultColor.fromLinearRGB(0.2f, 0.5f, 0.2f);
		m_color = props.getSpectrum("color", defaultColor);
	}
	
	/// Unserialize from a binary data stream
	MyIntegrator(Stream *stream, InstanceManager *manager) : SamplingIntegrator(stream, manager) {
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
			Float distance = rRec.its.t;
			//Float distance = (rRec.its.p - Point(0, 5000, 0)).length();
			//return Spectrum(1.0f - distance / m_maxDist) * m_color;

			//std::cout << rRec.its.toString() << std::endl;
			//std::cout << Spectrum(distance).toString() << std::endl;
			return Spectrum(distance);

		}
		return Spectrum(0.0f);
	}

	/// Preprocess function -- called on the initiating machine
	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job, int sceneResID, int cameraResID, int samplerResID) {
		SamplingIntegrator::preprocess(scene, queue, job, sceneResID, cameraResID, samplerResID);

		const AABB &sceneAABB = scene->getAABB();

		Point cameraPosition = scene->getSensor()->getWorldTransform()->eval(0).transformAffine(Point(0.0f));
		m_maxDist = -std::numeric_limits<Float>::infinity();

		Float minDist = std::numeric_limits<Float>::infinity();
		for (int i = 0; i < 8; ++i) {
			m_maxDist = std::max(m_maxDist, (cameraPosition - sceneAABB.getCorner(i)).length());
			minDist = std::min(m_maxDist, (cameraPosition - sceneAABB.getCorner(i)).length());
		}
		std::cout << "m_maxDist = " << m_maxDist << std::endl;
		std::cout << "minDist   = " << minDist << std::endl;

		return true;
	}

	MTS_DECLARE_CLASS()

private:
	Spectrum m_color;
	Float m_maxDist;
};

MTS_IMPLEMENT_CLASS_S(MyIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(MyIntegrator, "A contrived integrator")
MTS_NAMESPACE_END