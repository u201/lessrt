/*
This file is part of Mitsuba, a physically based rendering system.

Copyright (c) 2007-2014 by Wenzel Jakob and others.

Mitsuba is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License Version 3
as published by the Free Software Foundation.

Mitsuba is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/sensor.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/core/track.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{pushbroom}{Orthographic camera}
* From Orthographic camera
*/
class PushbroomCamera : public ProjectiveCamera {
public:
	PushbroomCamera(const Properties &props)
		: ProjectiveCamera(props) {
		m_type |= EDeltaDirection | EOrthographicCamera | EPositionSampleMapsToPixels;
	}

	PushbroomCamera(Stream *stream, InstanceManager *manager)
		: ProjectiveCamera(stream, manager) {
		configure();
	}

	void configure() {
		ProjectiveCamera::configure();

		const Vector2i &filmSize = m_film->getSize();
		const Vector2i &cropSize = m_film->getCropSize();
		const Point2i  &cropOffset = m_film->getCropOffset();

		Vector2 relSize((Float)cropSize.x / (Float)filmSize.x,
			(Float)cropSize.y / (Float)filmSize.y);
		Point2 relOffset((Float)cropOffset.x / (Float)filmSize.x,
			(Float)cropOffset.y / (Float)filmSize.y);

		/**
		* These do the following (in reverse order):
		*
		* 1. Create transform from camera space to [-1,1]x[-1,1]x[0,1] clip
		*    coordinates (not taking account of the aspect ratio yet)
		*
		* 2+3. Translate and scale to shift the clip coordinates into the
		*    range from zero to one, and take the aspect ratio into account.
		*
		* 4+5. Translate and scale the coordinates once more to account
		*     for a cropping window (if there is any)
		*/
		m_cameraToSample =
			Transform::scale(Vector(1.0f / relSize.x, 1.0f / relSize.y, 1.0f))
			* Transform::translate(Vector(-relOffset.x, -relOffset.y, 0.0f))
			* Transform::scale(Vector(-0.5f, -0.5f*m_aspect, 1.0f))
			* Transform::translate(Vector(-1.0f, -1.0f / m_aspect, 0.0f))
			* Transform::orthographic(m_nearClip, m_farClip);

		m_sampleToCamera = m_cameraToSample.inverse();

		/* Position differentials on the near plane */
		m_dx = m_sampleToCamera(Point(m_invResolution.x, 0.0f, 0.0f))
			- m_sampleToCamera(Point(0.0f));
		m_dy = m_sampleToCamera(Point(0.0f, m_invResolution.y, 0.0f))
			- m_sampleToCamera(Point(0.0f));

		/* Clip-space transformation for OpenGL */
		m_clipTransform = Transform::translate(
			Vector((1 - 2 * relOffset.x) / relSize.x - 1,
			-(1 - 2 * relOffset.y) / relSize.y + 1, 0.0f)) *
			Transform::scale(Vector(1.0f / relSize.x, 1.0f / relSize.y, 1.0f));

		const Transform &trafo = m_worldTransform->eval(0.0f);

		m_invSurfaceArea = 1.0f / (
			trafo(m_sampleToCamera(Vector(1, 0, 0))).length() *
			trafo(m_sampleToCamera(Vector(0, 1, 0))).length());

		m_scale = trafo(Vector(0, 0, 1)).length();
	}

	Spectrum sampleRay(Ray &ray, const Point2 &pixelSample,
		const Point2 &otherSample, Float timeSample) const {
		ray.time = sampleTime(timeSample);
		const Transform &trafo = m_worldTransform->eval(ray.time);

		/* Compute the corresponding position on the
		near plane (in local camera space) */
		/*Point nearP = m_sampleToCamera.transformAffine(Point(
		pixelSample.x * m_invResolution.x,
		pixelSample.y * m_invResolution.y, 0.0f));*/

		const Vector2i &filmSize = m_film->getSize();
		Point nearP = m_sampleToCamera.transformAffine(Point(
			filmSize.x*0.5*m_invResolution.x,
			pixelSample.y * m_invResolution.y, 0.0f));

		ray.setOrigin(trafo.transformAffine(
			Point(nearP.x, nearP.y, 0.0f)));
		/*Vector dir_o;
		Float scale_x = trafo.getMatrix().m[0][0];
		if (pixelSample.x > filmSize.x*0.5){
			dir_o = normalize(Vector(-1 / scale_x, 0, -1 / (1 - 2 * pixelSample.x / Float(filmSize.x)) / m_scanHFOV_tan));
		}
		else if (pixelSample.x == filmSize.x*0.5){
			dir_o = Vector(0, 0, 1);
		}
		else{
			dir_o = normalize(Vector(1 / scale_x, 0, 1 / (1 - 2 * pixelSample.x / Float(filmSize.x)) / m_scanHFOV_tan));
		}
		ray.setDirection(normalize(trafo(dir_o)));*/
		Float height = trafo.getMatrix().m[1][3];
		Point near_sample = m_sampleToCamera.transformAffine(Point(
			pixelSample.x*m_invResolution.x,
			pixelSample.y * m_invResolution.y, 0.0f));
		Point near_sample_real = trafo.transformAffine(near_sample) + trafo(Vector(0, 0, 1))*height;
		Vector dir = near_sample_real - ray.o;

		ray.setDirection(normalize(dir));
		ray.mint = m_nearClip;
		ray.maxt = m_farClip;

		return Spectrum(1.0f);
	}

	Spectrum sampleRayDifferential(RayDifferential &ray, const Point2 &pixelSample,
		const Point2 &otherSample, Float timeSample) const {
		ray.time = sampleTime(timeSample);
		const Transform &trafo = m_worldTransform->eval(ray.time);

		/* Compute the corresponding position on the
		near plane (in local camera space) */
		const Vector2i &filmSize = m_film->getSize();
		Point nearP = m_sampleToCamera.transformAffine(Point(
			filmSize.x*0.5*m_invResolution.x,
			pixelSample.y * m_invResolution.y, 0.0f));
		nearP.z = 0.0f;
		ray.setOrigin(trafo.transformAffine(nearP));
		//Vector dir_o;
		//Float scale_x = trafo.getMatrix().m[0][0];
		//if (pixelSample.x > filmSize.x*0.5){
		//	dir_o = normalize(Vector(-1 / scale_x, 0, -1 / (1 - 2 * pixelSample.x / Float(filmSize.x)) / m_scanHFOV_tan));
		//}
		//else if (pixelSample.x == filmSize.x*0.5){
		//	dir_o = Vector(0, 0, 1);
		//}
		//else{
		//	dir_o = normalize(Vector(1 / scale_x, 0, 1 / (1 - 2 * pixelSample.x / Float(filmSize.x)) / m_scanHFOV_tan));
		//}
		//cout << trafo.toString() << endl;
		Float height = trafo.getMatrix().m[1][3];
		Point near_sample = m_sampleToCamera.transformAffine(Point(
			pixelSample.x*m_invResolution.x,
			pixelSample.y * m_invResolution.y, 0.0f));
		Vector dir = (trafo.transformAffine(near_sample) + trafo(Vector(0, 0, 1))*height) - ray.o;
		ray.setDirection(normalize(dir));

		ray.mint = m_nearClip;
		ray.maxt = m_farClip;
		Vector rxdir = trafo.transformAffine(near_sample + m_dx) + trafo(Vector(0, 0, 1))*height - ray.o;
		ray.rxDirection = normalize(rxdir);
		Vector rydir = trafo.transformAffine(near_sample + m_dy) + trafo(Vector(0, 0, 1))*height - ray.o;
		ray.ryDirection = trafo(rydir);
		ray.rxOrigin = ray.ryOrigin = ray.o;
		ray.hasDifferentials = true;
		return Spectrum(1.0f);
	}

	Spectrum samplePosition(PositionSamplingRecord &pRec,
		const Point2 &sample, const Point2 *extra) const {
		const Transform &trafo = m_worldTransform->eval(pRec.time);

		Point samplePos(sample.x, sample.y, 0.0f);

		if (extra) {
			/* The caller wants to condition on a specific pixel position */
			samplePos.x = (extra->x + sample.x) * m_invResolution.x;
			samplePos.y = (extra->y + sample.y) * m_invResolution.y;
		}

		pRec.uv = Point2(samplePos.x * m_resolution.x,
			samplePos.y * m_resolution.y);

		Point nearP = m_sampleToCamera.transformAffine(samplePos);

		nearP.z = 0.0f;
		pRec.p = trafo.transformAffine(nearP);
		pRec.n = trafo(Vector(0.0f, 0.0f, 1.0f));
		pRec.pdf = m_invSurfaceArea;
		pRec.measure = EArea;
		return Spectrum(1.0f);
	}

	Spectrum evalPosition(const PositionSamplingRecord &pRec) const {
		return Spectrum((pRec.measure == EArea) ? m_invSurfaceArea : 0.0f);
	}

	Float pdfPosition(const PositionSamplingRecord &pRec) const {
		return (pRec.measure == EArea) ? m_invSurfaceArea : 0.0f;
	}

	Spectrum sampleDirection(DirectionSamplingRecord &dRec,
		PositionSamplingRecord &pRec, const Point2 &sample,
		const Point2 *extra) const {

		dRec.d = pRec.n;
		dRec.measure = EDiscrete;
		dRec.pdf = 1.0f;

		return Spectrum(1.0f);
	}

	Float pdfDirection(const DirectionSamplingRecord &dRec,
		const PositionSamplingRecord &pRec) const {
		return (dRec.measure == EDiscrete) ? 1.0f : 0.0f;
	}

	Spectrum evalDirection(const DirectionSamplingRecord &dRec,
		const PositionSamplingRecord &pRec) const {
		return Spectrum((dRec.measure == EDiscrete) ? 1.0f : 0.0f);
	}

	Spectrum sampleDirect(DirectSamplingRecord &dRec, const Point2 &) const {
		const Transform &trafo = m_worldTransform->eval(dRec.time);

		dRec.n = trafo(Vector(0, 0, 1));
		Float scale = dRec.n.length();

		Point localP = trafo.inverse().transformAffine(dRec.ref);
		localP.z *= scale;

		Point sample = m_cameraToSample.transformAffine(localP);

		if (sample.x < 0 || sample.x > 1 || sample.y < 0 ||
			sample.y > 1 || sample.z < 0 || sample.z > 1) {
			dRec.pdf = 0.0f;
			return Spectrum(0.0f);
		}

		dRec.p = trafo.transformAffine(Point(localP.x, localP.y, 0.0f));
		dRec.n /= scale;
		dRec.d = -dRec.n;
		dRec.dist = localP.z;
		dRec.uv = Point2(sample.x * m_resolution.x,
			sample.y * m_resolution.y);
		dRec.pdf = 1.0f;
		dRec.measure = EDiscrete;

		return Spectrum(m_invSurfaceArea);
	}

	Float pdfDirect(const DirectSamplingRecord &dRec) const {
		return (dRec.measure == EDiscrete) ? 1.0f : 0.0f;
	}

	bool getSamplePosition(const PositionSamplingRecord &pRec,
		const DirectionSamplingRecord &dRec, Point2 &samplePosition) const {
		const Transform &trafo = m_worldTransform->eval(pRec.time);

		Point localP = trafo.inverse()(pRec.p);
		Point sample = m_cameraToSample.transformAffine(localP);

		if (sample.x < 0 || sample.x > 1 || sample.y < 0 || sample.y > 1)
			return false;

		samplePosition = Point2(sample.x * m_resolution.x,
			sample.y * m_resolution.y);
		return true;
	}

	Transform getProjectionTransform(const Point2 &apertureSample,
		const Point2 &aaSample) const {
		Point2 offset(
			2.0f * m_invResolution.x * (aaSample.x - 0.5f),
			2.0f * m_invResolution.y * (aaSample.y - 0.5f));

		return m_clipTransform *
			Transform::translate(Vector(offset.x, offset.y, 0.0f)) *
			Transform::scale(Vector(1.0f, m_aspect, 1.0f)) *
			Transform::glOrthographic(m_nearClip, m_farClip)*
			Transform::scale(Vector(1.0f, 1.0f, m_scale));
	}

	AABB getAABB() const {
		AABB bounds;
		bounds.expandBy(m_sampleToCamera(Point(0, 0, 0)));
		bounds.expandBy(m_sampleToCamera(Point(1, 1, 0)));

		return m_worldTransform->getSpatialBounds(bounds);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "PushbroomCamera[" << endl
			<< "  nearClip = " << m_nearClip << "," << endl
			<< "  farClip = " << m_farClip << "," << endl
			<< "  worldTransform = " << indent(m_worldTransform.toString()) << "," << endl
			<< "  sampler = " << indent(m_sampler->toString()) << "," << endl
			<< "  film = " << indent(m_film->toString()) << "," << endl
			<< "  medium = " << indent(m_medium.toString()) << "," << endl
			<< "  shutterOpen = " << m_shutterOpen << "," << endl
			<< "  shutterOpenTime = " << m_shutterOpenTime << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Transform m_cameraToSample;
	Transform m_sampleToCamera;
	Transform m_clipTransform;
	Float m_invSurfaceArea, m_scale;
	Vector m_dx, m_dy;
};

MTS_IMPLEMENT_CLASS_S(PushbroomCamera, false, ProjectiveCamera)
MTS_EXPORT_PLUGIN(PushbroomCamera, "PushbroomCamera camera");
MTS_NAMESPACE_END
