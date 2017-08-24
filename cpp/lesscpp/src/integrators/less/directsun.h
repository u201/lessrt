#if !defined(__LESS_DIRECT_SUN_H_)
#define __LESS_DIRECT_SUN_H_

#include <mitsuba/core/sched.h>
#include <mitsuba/render/renderjob.h>
#include <ctime>
#include <list>
MTS_NAMESPACE_BEGIN

/**
* \brief Work unit that specifies a rectangular region in an plane.
*
* Used for emiting paticles into scene for sun direct illumination
*
*/
class SunDirectWorkUnit : public WorkUnit {
public:
	inline SunDirectWorkUnit() {}

	/* WorkUnit implementation */
	void set(const WorkUnit *wu)
	{
		const SunDirectWorkUnit *sundirect = static_cast<const SunDirectWorkUnit *>(wu);
		m_offset = sundirect->m_offset;
		m_size = sundirect->m_size;
		m_particlePerM2 = sundirect->m_particlePerM2;
		m_sunPlaneHeight = sundirect->m_sunPlaneHeight;
	}
	void load(Stream *stream)
	{
		int data[6];
		stream->readIntArray(data, 6);
		m_offset.x = data[0];
		m_offset.y = data[1];
		m_size.x = data[2];
		m_size.y = data[3];
		m_extendedSceneSize.x = data[4];
		m_extendedSceneSize.y = data[5];
		m_particlePerM2 = stream->readInt();
		m_sunPlaneHeight = stream->readFloat();
	}
	void save(Stream *stream) const
	{
		int data[6];
		data[0] = m_offset.x;
		data[1] = m_offset.y;
		data[2] = m_size.x;
		data[3] = m_size.y;
		data[4] = m_extendedSceneSize.x;
		data[5] = m_extendedSceneSize.y;
		stream->writeIntArray(data, 6);
		stream->writeInt(m_particlePerM2);
		stream->writeFloat(m_sunPlaneHeight);
	}

	inline const Point2i &getOffset() const { return m_offset; }
	inline const Vector2i &getSize() const { return m_size; }
	inline const int &getParticleNum() const { return m_particlePerM2; }
	inline const float &getSunPlaneHeight() const { return m_sunPlaneHeight; }
	inline const Vector2i &getExtendedSceneSize() const { return m_extendedSceneSize; }

	inline void setOffset(const Point2i &offset) { m_offset = offset; }
	inline void setSize(const Vector2i &size) { m_size = size; }
	inline void setParticleNum(int num) { m_particlePerM2 = num; }
	inline void setSunPlaneHeight(float height) { m_sunPlaneHeight = height; }
	inline void setExtendedSceneSize(const Vector2i &extendedSize) { m_extendedSceneSize = extendedSize; }

	std::string toString() const
	{
		std::ostringstream oss;
		oss << "SunDirectWorkUnit[offset=" << m_offset.toString()
			<< ", size=" << m_size.toString() << "]";
		return oss.str();
	}
	MTS_DECLARE_CLASS()
private:
	Point2i m_offset;
	Vector2i m_size;
	Vector2i m_extendedSceneSize;//添加了边缘缓冲区的场景水平方向尺寸：宽度和高度
	int m_particlePerM2;
	float m_sunPlaneHeight; //发射太阳光线的平面的高度，一般为场景最高点加上一定缓冲
};





/**
* \brief Work result that specifies store the result of particle trace
*
* Used for emiting paticles into scene for sun direct illumination
*
*/
class SunDirectWorkResult :public WorkResult
{
public:

	void load(Stream *stream) {
		int data[5];
		stream->readIntArray(data, 5);
		m_offset.x = data[0];
		m_offset.y = data[1];
		m_size.x = data[2];
		m_size.y = data[3];
		m_pointNum = data[4];
		stream->readFloatArray(reinterpret_cast<float *>(point_cloud),
			m_pointNum * 3);
		m_particlePerM2 = stream->readInt();
	}

	void save(Stream *stream) const {
		int data[5];
		data[0] = m_offset.x;
		data[1] = m_offset.y;
		data[2] = m_size.x;
		data[3] = m_size.y;
		data[4] = m_pointNum;
		stream->writeIntArray(data, 5);
		stream->writeFloatArray(reinterpret_cast<const float *>(point_cloud),
			m_pointNum * 3);
		stream->writeInt(m_particlePerM2);
	}

	inline const Point2i &getOffset() const { return m_offset; }
	inline const Vector2i &getSize() const { return m_size; }
	inline float* getPointCloud() const { return point_cloud; }
	inline int getPointNum() const { return m_pointNum; }

	inline void setOffset(const Point2i &offset) { m_offset = offset; }
	inline void setSize(const Vector2i &size) { m_size = size; }
	inline void setParticleNum(const int num) { m_particlePerM2 = num; }
	inline void setPointCloud(float * pointcloud) { point_cloud = pointcloud; }
	inline void setPointNum(int num) { m_pointNum = num; }


	std::string toString() const {
		std::ostringstream oss;
		oss << "SunDirectWorkResult[" << endl
			<< "]";
		return oss.str();
	}
	MTS_DECLARE_CLASS()
private:
	float* point_cloud;
	int m_pointNum; //actual intersected number
	Point2i m_offset;
	Vector2i m_size;
	int m_particlePerM2;

};


#define N 999 //三位小数。
/**
* \brief Work processor for sun illumination
*
* Used for emiting paticles into scene for sun direct illumination
*
*/
class DirectSunWorkProcessor : public WorkProcessor {
public:
	/// Construct a new work processor
	DirectSunWorkProcessor() : WorkProcessor() { }
	/// Unserialize from a binary data stream (nothing to do in our case)
	DirectSunWorkProcessor(Stream *stream, InstanceManager *manager)
		: WorkProcessor(stream, manager) { }
	/// Serialize to a binary data stream (nothing to do in our case)
	void serialize(Stream *stream, InstanceManager *manager) const {
	}

	ref<WorkUnit> createWorkUnit() const {
		return new SunDirectWorkUnit();
	}
	ref<WorkResult> createWorkResult() const {
		return new SunDirectWorkResult();
	}
	ref<WorkProcessor> clone() const {
		return new DirectSunWorkProcessor(); // No state to clone in our case
	}
	/// No internal state, thus no preparation is necessary
	void prepare() {
		Scene *scene = static_cast<Scene *>(getResource("scene"));
		m_scene = new Scene(scene);
		m_sampler = static_cast<Sampler *>(getResource("sampler"));
		m_sensor = static_cast<Sensor *>(getResource("sensor"));
		m_scene->removeSensor(scene->getSensor());
		m_scene->addSensor(m_sensor);
		m_scene->setSensor(m_sensor);
		m_scene->initializeBidirectional();

		//m_scene->initialize();
	}
	/// Do the actual computation
	void process(const WorkUnit *workUnit, WorkResult *workResult,
		const bool &stop);
	MTS_DECLARE_CLASS()
private:
	ref<Scene> m_scene;
	ref<Sampler> m_sampler;
	ref<Sensor> m_sensor;
};


/**
* \brief Work parallel process
*
* Used for emiting paticles into scene for sun direct illumination
*
*/
class DirectSunProcess : public ParallelProcess {
public:
	DirectSunProcess(int particle_per_m2,
		int scene_w, int scene_h, int block_w, int block_h, float sunPlaneHeight,float Pe) :
		m_particlePerM2(particle_per_m2), m_currentOffset(0) {
		this->scene_w = scene_w;
		this->scene_h = scene_h;
		this->m_PePhoton = Pe;
		this->m_sunPlaneHeight = sunPlaneHeight;
		m_size = Vector2i(block_w, block_h);
		m_numBlocks = Vector2i(
			(int)std::ceil((float)scene_w / (float)block_w),
			(int)std::ceil((float)scene_h / (float)block_h));
		//m_all_point_cloud = std::vector<Point>();
		m_all_point_cloud = new float*[m_numBlocks.x*m_numBlocks.y];
		m_actual_num_of_each_block = new int[m_numBlocks.x*m_numBlocks.y];
	}
	//~DirectSunProcess()
	//{
	//	delete[] m_all_point_cloud;
	//}

	void develop()
	{
		//m_queue->signalRefresh(m_job);
	}

	void bindResource(const std::string &name, int id)
	{
		ParallelProcess::bindResource(name, id);
	}

	ref<WorkProcessor> createWorkProcessor() const {
		return new DirectSunWorkProcessor();
	}
	//std::vector<std::string> getRequiredPlugins() {
	//	std::vector<std::string> result;
	//	result.push_back("LessSim");
	//	return result;
	//}
	EStatus generateWork(WorkUnit *unit, int worker /* unused */) {
		if (m_currentOffset >= (int)m_numBlocks.x*m_numBlocks.y)
			return EFailure;
		SunDirectWorkUnit *wu = static_cast<SunDirectWorkUnit *>(unit);
		wu->setSize(m_size);
		int r = int(m_currentOffset / m_numBlocks.x);
		int c = m_currentOffset % m_numBlocks.x;
		Point2i offset = Point2i(r, c); // offset: 代表subregion的编号
		wu->setOffset(offset);
		wu->setParticleNum(m_particlePerM2);
		wu->setExtendedSceneSize(Vector2i(this->scene_w, this->scene_h));
		wu->setSunPlaneHeight(m_sunPlaneHeight);
		m_currentOffset++;
		return ESuccess;
	}
	void processResult(const WorkResult *result, bool cancelled) {
		if (cancelled) // indicates a work unit, which was
			return; // cancelled partly through its execution
		const SunDirectWorkResult *wr =
			static_cast<const SunDirectWorkResult *>(result);
			Point2i offset = wr->getOffset();
			m_actual_num_of_each_block[offset.x*m_numBlocks.x + offset.y] = wr->getPointNum();
			m_all_point_cloud[offset.x*m_numBlocks.x + offset.y] = wr->getPointCloud();

	}
	inline float** &getOutput() {
	return m_all_point_cloud;
	}
	inline int* &get_actual_point_num(int &size)
	{
		size = m_numBlocks.x*m_numBlocks.y;
		return m_actual_num_of_each_block;
	}

	MTS_DECLARE_CLASS()
public:
	int m_currentOffset;   //用于分配任务，到各个workprocessor
	Vector2i m_size;		//每个子区域的大小 单位为m
	int m_particlePerM2;    //每平米发射的点的光线条数
	float** m_all_point_cloud; //最后得到的交点数据
	int * m_actual_num_of_each_block;
	int scene_w, scene_h;//width and height of scene (m)
	Vector2i m_numBlocks; //讲场景平面划分成不同子区域，每个区域并行运算
	float m_sunPlaneHeight; //发射太阳光线的平面的高度，一般为场景最高点加上一定缓冲
	Point3f m_sceneCenter;// centor point of scene
	float m_PePhoton;//每个光子的能量，w
};



MTS_NAMESPACE_END

#endif 