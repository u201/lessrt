/*
author: Jianbo QI
date: 2018.4.3
Implementation of a real sun emitter
*/

#include <mitsuba/render/emitter.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

class SunSphereLight : public Emitter {
public:
	SunSphereLight(const Properties &props) : Emitter(props) {
		m_type |= EOnSurface;

		if (props.hasProperty("toWorld"))
			Log(EError, "Found a 'toWorld' transformation -- this is not "
				"allowed -- the area light inherits this transformation from "
				"its parent shape");

		m_radiance = props.getSpectrum("radiance", Spectrum::getD65());
		m_power = Spectrum(0.0f); /// Don't know the power yet
	}

	AreaLight(Stream *stream, InstanceManager *manager)
		: Emitter(stream, manager) {
		m_radiance = Spectrum(stream);
		m_power = Spectrum(stream);
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Emitter::serialize(stream, manager);
		m_radiance.serialize(stream);
		m_power.serialize(stream);
	}

	Spectrum samplePosition(PositionSamplingRecord &pRec,
		const Point2 &sample, const Point2 *extra) const {
		m_shape->samplePosition(pRec, sample);
		return m_power;
	}

	Spectrum evalPosition(const PositionSamplingRecord &pRec) const {
		return m_radiance * M_PI;
	}

	Spectrum eval(const Intersection &its, const Vector &d) const {
		if (dot(its.shFrame.n, d) <= 0)
			return Spectrum(0.0f);
		else
			return m_radiance;
	}

	Float pdfPosition(const PositionSamplingRecord &pRec) const {
		return m_shape->pdfPosition(pRec);
	}

	Spectrum sampleDirection(DirectionSamplingRecord &dRec,
		PositionSamplingRecord &pRec,
		const Point2 &sample, const Point2 *extra) const {
		Vector local = warp::squareToCosineHemisphere(sample);
		dRec.d = Frame(pRec.n).toWorld(local);
		dRec.pdf = warp::squareToCosineHemispherePdf(local);
		dRec.measure = ESolidAngle;
		return Spectrum(1.0f);
	}

	Spectrum evalDirection(const DirectionSamplingRecord &dRec,
		const PositionSamplingRecord &pRec) const {
		Float dp = dot(dRec.d, pRec.n);

		if (dRec.measure != ESolidAngle || dp < 0)
			dp = 0.0f;

		return Spectrum(INV_PI * dp);
	}

	Float pdfDirection(const DirectionSamplingRecord &dRec,
		const PositionSamplingRecord &pRec) const {
		Float dp = dot(dRec.d, pRec.n);

		if (dRec.measure != ESolidAngle || dp < 0)
			dp = 0.0f;

		return INV_PI * dp;
	}

	Spectrum sampleRay(Ray &ray,
		const Point2 &spatialSample,
		const Point2 &directionalSample,
		Float time) const {
		PositionSamplingRecord pRec(time);
		m_shape->samplePosition(pRec, spatialSample);
		Vector local = warp::squareToCosineHemisphere(directionalSample);
		ray.setTime(time);
		ray.setOrigin(pRec.p);
		ray.setDirection(Frame(pRec.n).toWorld(local));
		return m_power;
	}

	Spectrum sampleDirect(DirectSamplingRecord &dRec,
		const Point2 &sample) const {
		m_shape->sampleDirect(dRec, sample);

		/* Check that the emitter and reference position are oriented correctly
		with respect to each other. Note that the >= 0 check
		for 'refN' is intentional -- those sampling requests that specify
		a reference point within a medium or on a transmissive surface
		will set dRec.refN = 0, hence they should always be accepted. */
		if (dot(dRec.d, dRec.refN) >= 0 && dot(dRec.d, dRec.n) < 0 && dRec.pdf != 0) {
			return m_radiance / dRec.pdf;
		}
		else {
			dRec.pdf = 0.0f;
			return Spectrum(0.0f);
		}
	}

	Float pdfDirect(const DirectSamplingRecord &dRec) const {
		/* Check that the emitter and receiver are oriented correctly
		with respect to each other. */
		if (dot(dRec.d, dRec.refN) >= 0 && dot(dRec.d, dRec.n) < 0) {
			return m_shape->pdfDirect(dRec);
		}
		else {
			return 0.0f;
		}
	}

	void setParent(ConfigurableObject *parent) {
		Emitter::setParent(parent);

		if (parent->getClass()->derivesFrom(MTS_CLASS(Shape))) {
			Shape *shape = static_cast<Shape *>(parent);
			if (m_shape == shape || shape->isCompound())
				return;

			if (m_shape != NULL)
				Log(EError, "An area light cannot be parent of multiple shapes");

			m_shape = shape;
			m_shape->configure();
			m_power = m_radiance * M_PI * m_shape->getSurfaceArea();
		}
		else {
			Log(EError, "An area light must be child of a shape instance");
		}
	}

	AABB getAABB() const {
		return m_shape->getAABB();
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "AreaLight[" << endl
			<< "  radiance = " << m_radiance.toString() << "," << endl
			<< "  samplingWeight = " << m_samplingWeight << "," << endl
			<< "  surfaceArea = ";
		if (m_shape)
			oss << m_shape->getSurfaceArea();
		else
			oss << "<no shape attached!>";
		oss << "," << endl
			<< "  medium = " << indent(m_medium.toString()) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
protected:
	Spectrum m_radiance, m_power;
};


MTS_IMPLEMENT_CLASS_S(AreaLight, false, Emitter)
MTS_EXPORT_PLUGIN(AreaLight, "Area light");
MTS_NAMESPACE_END
