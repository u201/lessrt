#pragma once

#if !defined(_PHOTONRT_H_)
#define _PHOTONRT_H_

#include <mitsuba/render/photonproc.h>
#include <mitsuba/render/range.h>
#include <mitsuba/render/renderjob.h>
#include <mitsuba/core/bitmap.h>
MTS_NAMESPACE_BEGIN

//结果保存到图像中
class CapturePhotonWorkResult :public ImageBlock {
public:
	inline CapturePhotonWorkResult(const Vector2i &res, const ReconstructionFilter *filter)
		: ImageBlock(Bitmap::ESpectrum, res, filter) {
		setOffset(Point2i(0, 0));
		setSize(res);
		m_range = new RangeWorkUnit();

		m_downwellingWorkResult = new ImageBlock(Bitmap::ESpectrum, res, filter);
		m_upwellingWorkResult = new ImageBlock(Bitmap::ESpectrum, res, filter);
		m_PhtonsEachProcess = 0;

	}

	inline const RangeWorkUnit *getRangeWorkUnit() const {
		return m_range.get();
	}

	inline void setRangeWorkUnit(const RangeWorkUnit *range) {
		m_range->set(range);
	}

	/* Work unit implementation */
	void load(Stream *stream);
	void save(Stream *stream) const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~CapturePhotonWorkResult() { }
protected:
	ref<RangeWorkUnit> m_range;
public:
	ref<ImageBlock> m_downwellingWorkResult;
	ref<ImageBlock> m_upwellingWorkResult;
	size_t m_PhtonsEachProcess;
};


/* ==================================================================== */
/*                             Work processor                           */
/* ==================================================================== */

/**
* \brief Particle tracing worker -- looks for volume and surface interactions
* and tries to accumulate the resulting information at the image plane.
*/
class CapturePhotonWorker : public PhotonTracer {
public:
	inline CapturePhotonWorker(int maxDepth, int maxPathDepth,
		int rrDepth, bool bruteForce) : PhotonTracer(maxDepth, rrDepth, true),
		m_maxPathDepth(maxPathDepth), m_bruteForce(bruteForce) { }

	CapturePhotonWorker(Stream *stream, InstanceManager *manager);

	void serialize(Stream *stream, InstanceManager *manager) const;

	void prepare();
	ref<WorkProcessor> clone() const;
	ref<WorkResult> createWorkResult() const;
	void process(const WorkUnit *workUnit, WorkResult *workResult,
		const bool &stop);

	/**
	* \brief Handles particles emitted by a light source
	*
	* If a connection to the sensor is possible, compute the importance
	* and accumulate in the proper pixel of the accumulation buffer.
	*/
	void handleEmission(const PositionSamplingRecord &pRec,
		const Medium *medium, const Spectrum &weight);

	/**
	* \brief Handles particles interacting with a surface
	*
	* If a connection to the sensor is possible, compute the importance
	* and accumulate in the proper pixel of the accumulation buffer.
	*/
	void handleSurfaceInteraction(int depth, int nullInteractions, bool caustic,
		const Intersection &its, const Medium *medium,
		const Spectrum &weight);

	/**
	* \brief extended version of handleSurfaceInteraction
	*
	* If a connection to the sensor is possible, compute the importance
	* and accumulate in the proper pixel of the accumulation buffer.
	*/
	void handleSurfaceInteractionExt(int depth, int nullInteractions,
		bool delta, const Intersection &its, Point &previousPoint, const Medium *medium,
		const Spectrum &weight);

	/**
	* \brief Handles particles interacting with a medium
	*
	* If a connection to the sensor is possible, compute the importance
	* and accumulate in the proper pixel of the accumulation buffer.
	*/
	void handleMediumInteraction(int depth, int nullInteractions, bool caustic,
		const MediumSamplingRecord &mRec, const Medium *medium,
		const Vector &wi, const Spectrum &weight);

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~CapturePhotonWorker() { }
private:
	ref<const Sensor> m_sensor;
	ref<const ReconstructionFilter> m_rfilter;
	ref<CapturePhotonWorkResult> m_workResult;
	int m_maxPathDepth;
	bool m_bruteForce;
};


/* ==================================================================== */
/*                           Parallel process                           */
/* ==================================================================== */
/**
* Parallel particle tracing process - used to run this over
* a group of machines
*/
class CapturePhotonProcess : public PhotonProcess {
public:
	CapturePhotonProcess(const RenderJob *job, RenderQueue *queue,
		size_t sampleCount, size_t granularity, int maxDepth,
		int maxPathDepth, int rrDepth, bool bruteForce)
		: PhotonProcess(PhotonProcess::ETrace, sampleCount,
			granularity, "Simulating", job), m_job(job), m_queue(queue),
		m_maxDepth(maxDepth), m_maxPathDepth(maxPathDepth),
		m_rrDepth(rrDepth), m_bruteForce(bruteForce) {
	}

	void develop();

	/* ParallelProcess impl. */
	void processResult(const WorkResult *wr, bool cancelled);
	void bindResource(const std::string &name, int id);
	ref<WorkProcessor> createWorkProcessor() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~CapturePhotonProcess() { }
private:
	ref<const RenderJob> m_job;
	ref<RenderQueue> m_queue;
	ref<Film> m_film;
	//ref<ImageBlock> m_accum;
	int m_maxDepth;
	int m_maxPathDepth;
	int m_rrDepth;
	bool m_bruteForce;

	ref<Scene> m_scene;

	ref<ImageBlock> m_accum_downwell;

	ref<Film> m_film_upwell;
	ref<ImageBlock> m_accum_upwell;

	size_t m_totalPhotons;
};


MTS_NAMESPACE_END
#endif
