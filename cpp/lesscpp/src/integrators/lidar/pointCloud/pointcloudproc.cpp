#include <boost/filesystem.hpp>

#include <mitsuba/core/sched.h>
#include <mitsuba/render/scene.h>
#include <vector>
#include <iostream>
#include <windows.h>
#include <cstdio>
#include <cstring>
//#include "lidarsampler.cpp"
#include "../lidarutils.h"
#include "../fftconv1d.cpp"

#include "../PulseGaussianFitting.h"
#include "../mpfit.h"
#include "../threadsafe_queue.h"



MTS_NAMESPACE_BEGIN

class PointCloudWorkUnit;
class PointCloudWorkResult;
class PointCloudWorkProcessor;
class PointCloudProcess;

struct DiscretePoint {
	Float x;
	Float y;
	Float z;
	Float a;
	int i;
	DiscretePoint() {}
	DiscretePoint(Float x, Float y, Float z, Float a, int i) : x(x), y(y), z(z), a(a), i(i) {}
};

struct WorkResultRecord {
	Float l;
	Spectrum w;
	WorkResultRecord() {}
	WorkResultRecord(Float l, Spectrum w) : l(l), w(w) {}
};

class PointCloudWorkUnit : public WorkUnit {
public:
	void set(const WorkUnit *workUnit) {
		const PointCloudWorkUnit *wu = static_cast<const PointCloudWorkUnit *>(workUnit);
		m_x = wu->m_x;
		m_y = wu->m_y;
		m_z = wu->m_z;
		m_u = wu->m_u;
		m_v = wu->m_v;
		m_w = wu->m_w;
		m_i = wu->m_i;
	}

	void load(Stream *stream) {
		m_x = stream->readDouble();
		m_y = stream->readDouble();
		m_z = stream->readDouble();
		m_u = stream->readDouble();
		m_v = stream->readDouble();
		m_w = stream->readDouble();
		m_i = stream->readInt();
	}

	void save(Stream *stream) const {
		stream->writeDouble(m_x);
		stream->writeDouble(m_y);
		stream->writeDouble(m_z);
		stream->writeDouble(m_u);
		stream->writeDouble(m_v);
		stream->writeDouble(m_w);
		stream->writeInt(m_i);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "PointCloudWorkUnit[" << endl
			<< "  x = " << m_x << endl
			<< "  y = " << m_y << endl
			<< "  z = " << m_z << endl
			<< "  u = " << m_u << endl
			<< "  v = " << m_v << endl
			<< "  w = " << m_w << endl
			<< "  i = " << m_i << endl
			<< "]";
		return oss.str();
	}

	inline Point getOrigin() const {
		return Point(m_x, m_y, m_z);
	}

	inline Vector getDirection() const {
		return Vector(m_u, m_v, m_w);
	}

	inline void set(Float x, Float y, Float z, Float u, Float v, Float w) {
		m_x = x;
		m_y = y;
		m_z = z;
		m_u = u;
		m_v = v;
		m_w = w;
	}

	MTS_DECLARE_CLASS()

public:

	double m_x;
	double m_y;
	double m_z;
	double m_u;
	double m_v;
	double m_w;
	int m_i;
};

/* ==================================================================== */
/*                           Work result impl.                          */
/* ==================================================================== */
class PointCloudWorkResult : public WorkResult {
public:

	void load(Stream *stream) {
		cout << "Work Result load" << endl;
		//stream->readDoubleArray(m_bin, m_numOfBins);
		//m_i = stream->readInt();
		//m_numOfBins = stream->readSize();
		//for (int i = 0; i < m_numOfBins; i++) {
		//	//cout << "PointCloudWorkResult::load i = " << i << endl;
		//	m_bin[i] = Spectrum(stream);
		//}

		//stream->readArray<Float>((Float *) m_lengths.data(), m_lengths.size());
		//stream->readArray<Float>((Float *) m_weights.data(), m_weights.size());
	}

	void save(Stream *stream) const {
		cout << "Work Result save" << endl;
		//stream->writeDoubleArray(m_bin, m_numOfBins);
		//stream->writeInt(m_i);
		//stream->writeSize(m_numOfBins);
		//for (int i = 0; i < m_numOfBins; i++) {
		//	//cout << "PointCloudWorkResult::save i = " << i << endl;
		//	m_bin[i].serialize(stream);
		//}

		//stream->writeFloatArray((Float *) m_lengths.data(), m_lengths.size());
		//stream->writeFloatArray((Float *) m_weights.data(), m_weights.size());
	}

	std::string toString() const {
		std::ostringstream oss;
		//oss << "PointCloudWorkResult[" << endl
		//	<< "  m_numOfBins = " << m_numOfBins << endl
		//	<< "  m_i = " << m_i << ", " << endl
		//	<< "]";
		return oss.str();
	}

	//void addToBin(int i, const Spectrum &w) {
	//	//cout << "i = " << i << ", num of bins = " << m_numOfBins << endl;
	//	if (i >= m_numOfBins) {
	//		return;
	//	}
	//	m_bin[i] += w;
	//	//cout << "addToBin m_bin[i] = " << m_bin[i].toString() << endl;
	//}

	void commit(Float length, Spectrum weight) {
		//m_lengths.push_back(length);
		//m_weights.push_back(weight);
		m_records.push_back(WorkResultRecord(length, weight));
	}

	//const std::vector<Float>& getLengths() const {
	//	return m_lengths;
	//}

	//const std::vector<Spectrum>& getWeights() const {
	//	return m_weights;
	//}

	MTS_DECLARE_CLASS()

public:
	
	//std::vector<Spectrum> m_bin;

	int m_i;

//private:
	//std::vector<Float> m_lengths;
	//std::vector<Spectrum> m_weights;
	std::vector<WorkResultRecord> m_records;
};

/* ==================================================================== */
/*                         Work processor impl.                         */
/* ==================================================================== */
class PointCloudWorkProcessor : public WorkProcessor {
public:
	PointCloudWorkProcessor() : WorkProcessor() {}

	PointCloudWorkProcessor(Stream *stream, InstanceManager *manager) : WorkProcessor(stream, manager) {}

	void serialize(Stream *stream, InstanceManager *manager) const {}

	ref<WorkUnit> createWorkUnit() const {
		return new PointCloudWorkUnit();
	}

	ref<WorkResult> createWorkResult() const {
		return new PointCloudWorkResult();
	}

	ref<WorkProcessor> clone() const {
		return new PointCloudWorkProcessor();
	}

	void prepare() {
		Scene *scene = static_cast<Scene *>(getResource("scene"));
		m_scene = new Scene(scene);
		m_random = new Random();
	}

	void process(const WorkUnit *workUnit, WorkResult *workResult, const bool &stop) {
		const PointCloudWorkUnit *wu = static_cast<const PointCloudWorkUnit *>(workUnit);
		PointCloudWorkResult *wr = static_cast<PointCloudWorkResult *>(workResult);

		//cout << wr << endl;

		wr->m_i = wu->m_i;
		//cout << "wr->m_i = " << wr->m_i << endl;
		wr->m_records.clear();
		//cout << "process wr->m_i = " << wr->m_i << ", clear" << endl;
		//wr->m_numOfBins = m_numOfBins;
		//wr->m_bin.clear();
		//wr->m_bin.resize(m_numOfBins);
		getWaveformResult(wu, wr);
		//cout << "process wr-> m_i = " << wr->m_i << ", size = " << wr->m_records.size() << endl;
	}

	// not override methods

	void getWaveformResult(const PointCloudWorkUnit *wu, PointCloudWorkResult *wr) {

		m_convergence = getConvergencePoint(wu);
		m_rotateRayTransform = getRotateRayTransform(wu);

		CircleBeamGridSampler sampler(m_axialDivision);

		while (sampler.hasNext()) {
			Vector2 s = sampler.next();
			Ray ray = generateRay(s, wu);
			Spectrum w = Spectrum(gaussian(s.length(), m_sigmaSquareOfBeam)) / m_weightTotal * m_pulseEnergy;
			trace(wu, wr, ray, w);
		}
	}

	void trace(const PointCloudWorkUnit *wu, PointCloudWorkResult *wr, Ray &ray, Spectrum w) {
		Intersection its;
		Float l = 0;
		Spectrum r;

		for (int depth = 0; depth < m_maxDepth && m_scene->rayIntersect(ray, its); depth++) {
			const BSDF *bsdf = its.getBSDF();
			r = bsdf->getDiffuseReflectance(its);

			record(wu, wr, its, l, w);

			// calculate new ray
			l += its.t;
			w = w * r;
			ray.setOrigin(its.p);
			ray.setDirection(generateRandomDirection(its));
		}
	}

	void record(const PointCloudWorkUnit *wu, PointCloudWorkResult *wr, const Intersection &its, const Float &l, const Spectrum &w) {
		Ray ray;
		ray.setOrigin(its.p);
		ray.setDirection(normalize(wu->getOrigin() - its.p));

		Intersection shadowIts;
		m_scene->rayIntersect(ray, shadowIts);
		// if it is not blocked and in fov
		Vector vectorFromConvergence = normalize(its.p - m_convergence);
		Vector wiWorld = its.shFrame.toWorld(its.wi);
		if (dot(wiWorld, its.shFrame.n) * dot(vectorFromConvergence, its.shFrame.n) < 0
			&& ((wu->getOrigin()) - its.p).length() < shadowIts.t
			&& dot(vectorFromConvergence, wu->getDirection()) > m_cosFov) {
			// calculate path length and weight
			Vector v = wu->getOrigin() - its.p;
			Float d = v.length();
			Float solidAngle = m_area / (d * d) * (dot(-wu->getDirection(), normalize(v)));

			Float length = l + its.t + d;
			Spectrum weight = solidAngle * its.getBSDF()->getDiffuseReflectance(its) * INV_PI * w * absDot(normalize(wiWorld), its.shFrame.n);

			wr->commit(length, weight);

			//int i = (int)((length - 2.0 * m_minRange) / (C * m_rate));
			//wr->addToBin(i, weight);

		}
	}

	Vector generateRandomDirection(const Intersection &its) {

		Float azimuth = 2.0 * M_PI * m_random->nextFloat();
		Float zenith = acos(m_random->nextFloat());

		Float x = sin(zenith) * cos(azimuth);
		Float y = sin(zenith) * sin(azimuth);
		Float z = cos(zenith);

		Vector n(its.shFrame.n);
		Vector wiWorld = its.shFrame.toWorld(its.wi);
		if (dot(wiWorld, n) < 0) {
			n = Vector(-n);
		}
		Vector direction = its.shFrame.s * x + its.shFrame.t * y + n * z;
		return direction;
	}

	Ray generateRay(const Vector2 &s, const PointCloudWorkUnit *wu) {
		const Float d = 10000;
		Float t = d * tan(m_fov);
		Point p = Point(s.x * t, -d, s.y * t);

		p = m_rotateRayTransform.transformAffine(p);

		Ray ray;
		ray.setOrigin(wu->getOrigin());
		ray.setDirection(normalize(Vector(p)));

		return ray;
	}

	Point getConvergencePoint(const PointCloudWorkUnit *wu) {
		Float r = sqrt(m_area * INV_PI);
		Point c = wu->getOrigin() - wu->getDirection() * (r / tan(m_fov));

		return c;
	}

	Transform getRotateRayTransform(const PointCloudWorkUnit *wu) {
		Transform t;
		Float degree;

		Vector z(0, -1, 0);
		Vector n = wu->getDirection();

		Vector axis = cross(z, n);
		if (abs(dot(axis, axis)) > 0) {
			degree = acos(dot(n, z));
			degree = radToDeg(degree);
			t = Transform::rotate(axis, degree);
		}

		return t;
	}

	MTS_DECLARE_CLASS()
public:
	ref<Scene> m_scene;

	Point m_convergence;
	Transform m_rotateRayTransform;

	ref<Random> m_random;

	// parameter

	Float m_sigmaSquareOfBeam;
	Spectrum m_weightTotal;
	Float m_cosFov;
	int m_numOfBins;

	// from XML

	int m_maxDepth;
	int m_axialDivision;

	Float m_fp;
	Float m_fov;
	Float m_pulseEnergy;
	Float m_rate;
	Float m_area;
	Float m_minRange;
	Float m_maxRange;
};

/* ==================================================================== */
/*                        Parallel process impl.                        */
/* ==================================================================== */
class PointCloudProcess : public ParallelProcess {
public:
	PointCloudProcess() : m_pos(0), m_numOfPulses(0) {}

	ref<WorkProcessor> createWorkProcessor() const {
		PointCloudWorkProcessor *processor = new PointCloudWorkProcessor();

		processor->m_sigmaSquareOfBeam = m_sigmaSquareOfBeam;
		processor->m_weightTotal = m_weightTotal;
		processor->m_cosFov = m_cosFov;
		processor->m_numOfBins = m_numOfBins;

		processor->m_maxDepth = m_maxDepth;
		processor->m_axialDivision = m_axialDivision;
		processor->m_fp = m_fp;
		processor->m_fov = m_fov;
		processor->m_pulseEnergy = m_pulseEnergy;
		processor->m_rate = m_rate;
		processor->m_area = m_area;
		processor->m_minRange = m_minRange;
		processor->m_maxRange = m_maxRange;

		return processor;
	}

	EStatus generateWork(WorkUnit *unit, int worker) {
		if (m_pos >= m_numOfPulses) {
			return EFailure;
		}
		PointCloudWorkUnit *wu = static_cast<PointCloudWorkUnit *>(unit);
		wu->set(m_x[m_pos], m_y[m_pos], m_z[m_pos], m_u[m_pos], m_v[m_pos], m_w[m_pos]);
		wu->m_i = m_pos;
		m_pos++;
		return ESuccess;
	}

	void processResult(const WorkResult *result, bool cancelled) {
		if (cancelled) {
			return;
		}

		const PointCloudWorkResult *wr = static_cast<const PointCloudWorkResult *>(result);

		// different from `waveformproc.cpp`
		//cout << "processResult start" << endl;
		//cout << "m_numOfBins = " << m_numOfBins << endl;
		std::vector<Spectrum> bin(m_numOfBins);
		int idx;
		//cout << "process Result wr->m_lengths.size = " << wr->m_lengths.size() << endl;
		//const std::vector<Float> &lengths = wr->getLengths();
		//const std::vector<Spectrum> &weights = wr->getWeights();
		//for (int i = 0; i < lengths.size(); ++i) {
		//	idx = (int)((lengths[i] - 2.0 * m_minRange) / (C * m_rate));
		//	addToBin(bin, idx, weights[i]);
		//}
		//for (int i = 0; i < wr->m_lengths.size(); ++i) {
		//	idx = (int)((wr->m_lengths[i] - 2.0 * m_minRange) / (C * m_rate));
		//	addToBin(bin, idx, wr->m_weights[i]);
		//}
		for (const WorkResultRecord& r : wr->m_records) {
			idx = (int)((r.l - 2.0 * m_minRange) / (C * m_rate));
			addToBin(bin, idx, r.w);
		}
		//cout << "add To bin finished" << endl;

		// BEGIN: get point cloud from accumulation

		// TODO: 1. convolve
		std::vector<double> accumulation(bin.size());
		//cout << "accum" << endl;
		for (int i = 0; i < bin.size(); ++i) {
			accumulation[i] = bin[i].eval(361);
			//cout << accumulation[i] << " ";
		}
		//cout << endl << "accumulation finished" << endl;
		std::vector<double> waveform = conv(accumulation, m_pulse);
		//cout << "wf" << endl;
		for (int i = 0; i < bin.size(); ++i) {
			//cout << waveform[i] << " ";
		}
		//cout << endl << "conv finished" << endl;

		// TODO: 2. gaussian decomposition
		std::vector<double> par;  // amp, center, sigma
		gaussianDecomposition(waveform, par);
		//cout << "par finished" << endl;

		// TODO: 3. output point cloud
		Float step = C * m_rate * 0.5;
		Float a;
		Float t;
		Point o = Point(m_x[wr->m_i], m_y[wr->m_i], m_z[wr->m_i]);
		Vector d = Vector(m_u[wr->m_i], m_v[wr->m_i], m_w[wr->m_i]);
		Point p;
		//std::ostringstream oss;
		//cout << "par.size = " << par.size() << endl;
		for (int i = 0; 3 * i < par.size(); i++) {
			a = par[i * 3] / (sqrt(2 * M_PI) * par[i * 3 + 2]);
			t = m_minRange + (par[i * 3 + 1] + 0.5) * step;
			//cout << "par center = " << par[i * 3 + 1] << endl;
			//cout << "step = " << step << endl;
			//cout << "t = " << t << endl;
			//cout << "o = " << o.toString() << endl;
			//cout << "d = " << d.toString() << endl;
			p = o + d * t;
			//cout << p.x << p.y << p.z << a << wr->m_i << endl;
			
			//oss << p.x << "\t" << p.y << "\t" << p.z << "\t" << a << "\t" << wr->m_i << endl;
			//m_pointCloudStringStream << oss.str();
			//oss.str("");

			//cout << p.x << "\t" << p.y << "\t" << p.z << "\t" << a << "\t" << wr->m_i << endl;
			DiscretePoint point(p.x, p.y, p.z, a, wr->m_i);
			m_pointCloud.push(point);
		}
		//cout << "output point cloud finished" << endl;

		// END: get point cloud from accumulation
	}

	void bindResource(const std::string &name, int id) {
		if (name == "scene") {
			m_scene = static_cast<Scene *>(Scheduler::getInstance()->getResource(id));
		}

		ParallelProcess::bindResource(name, id);
	}

	// not override methods

	void gaussianDecomposition(std::vector<double> &wf, std::vector<double> &par) {
		vector<int> peaks = peaksDetectFisrtOrderZeroCrosssing(wf, 0.005, 3); //0.005 is set to 1/200 of the maximum waveform amplitude to avoid too small peaks.
		//cout<<"peaks.size(): "<<peaks.size()<<endl;
		if (peaks.size() > 0) {
			vector<vector<int>> flexions = flexion_detect(wf, peaks);
			//cout << "after flexions_detect" << endl;

			//    vector<double> y_error;
				// Gaussian decomposition
			//vector<double> par; // amp, center, sigma
		//    if(intensityValueType != LidarProprietes::LIDAR_INTENSITY_NONE){ //0: only points(no need to create gaussian decomposition); 1: amplitude; 2: integral; 3: sigma; 4: solar signal; 5: all
			vector<double> y_error;
			vector<double> waveBinIndex;
			y_error.resize(wf.size());
			waveBinIndex.resize(wf.size());
			for (unsigned int i = 0; i < wf.size(); i++) {
				waveBinIndex[i] = i;
				//							*_wfdto.dStep;
				y_error[i] = 0.01;
			}

			for (unsigned int i = 0; i < peaks.size(); i++) {
				vector<double> sub_x(waveBinIndex.begin() + flexions[i][0], waveBinIndex.begin() + flexions[i][1]);
				vector<double> sub_y(wf.begin() + flexions[i][0], wf.begin() + flexions[i][1]);
				vector<double> comp_par = guess(sub_x, sub_y);
				//cout << "guess i = " << i << endl;
				par.push_back(comp_par[0]);
				par.push_back(comp_par[1]);
				par.push_back(comp_par[2]);
			}
			//cout << "before paramConstraints" << endl;

			mp_par *paramConstraints = new mp_par[peaks.size() * 3 * sizeof(mp_par)];
			memset(paramConstraints, 0, peaks.size() * 3 * sizeof(mp_par));
			for (unsigned int i = 0; i < peaks.size(); i++) {
				int idx = i * 3;

				//integral
				paramConstraints[idx].fixed = false;
				paramConstraints[idx].limited[0] = true;
				paramConstraints[idx].limits[0] = 0.0;

				//center
				paramConstraints[idx + 1].fixed = false;
				paramConstraints[idx + 1].limited[0] = true;
				paramConstraints[idx + 1].limited[1] = true;
				paramConstraints[idx + 1].limits[0] = waveBinIndex[flexions[i][0]];
				paramConstraints[idx + 1].limits[1] = waveBinIndex[flexions[i][1]];

				//sigma
				paramConstraints[idx + 2].fixed = false;
				paramConstraints[idx + 2].limited[0] = true;
				paramConstraints[idx + 2].limits[0] = 0.0;

				//						cout<<"range..."<<endl;
				//						cout<<paramConstraints[idx + 0].fixed<<paramConstraints[idx + 0].limited[0]<<paramConstraints[idx + 0].limited[1]<<endl;
				//						cout<<paramConstraints[idx + 1].fixed<<paramConstraints[idx + 1].limited[0]<<paramConstraints[idx + 1].limited[1]<<endl;
				//						cout<<paramConstraints[idx + 2].fixed<<paramConstraints[idx + 2].limited[0]<<paramConstraints[idx + 2].limited[1]<<endl;
			}
			XYData xydata;
			xydata.x = waveBinIndex.data();
			xydata.y = wf.data();
			xydata.y_error = y_error.data();
			//cout << "mpfit start" << endl;
			int status = mpfit(GaussianSum, waveBinIndex.size(), par.size(), par.data(), paramConstraints, 0, (void*)&xydata, 0);
			//cout << "mpfit over" << endl;
			if (status <= 0) {
				cout << "Failed to perform Gaussian decomposition." << endl;
			}
		}
	}

	void addToBin(std::vector<Spectrum> &bin, const int &i, const Spectrum &w) {
		if (i >= m_numOfBins) {
			//cout << "bin idx = " << i << endl;
			//cout << "numOfBins = " << m_numOfBins << endl;
			return;
		}
		//cout << w.toString() << endl;
		bin[i] += w;
	}

	void addGeometryConfiguration(Float x, Float y, Float z, Float u, Float v, Float w) {
		m_x.push_back(x);
		m_y.push_back(y);
		m_z.push_back(z);

		Vector d(u, v, w);
		d = normalize(d);
		m_u.push_back(d.x);
		m_v.push_back(d.y);
		m_w.push_back(d.z);

		//m_numOfPulses++;
	}

	//void outputToFile(const PointCloudWorkResult *wr) {

	//	std::string folderPath = m_scene->getDestinationFile().string();
	//	std::ostringstream oss;
	//	oss << folderPath << "\\" << wr->m_i << ".txt";
	//	std::string name = oss.str();
	//	std::ofstream fout(name);

	//	Float l = m_minRange * 2.0;
	//	Float d = C * m_rate;

	//	for (int i = 0; i < wr->m_numOfBins; i++) {
	//		//fout << l << "\t" << wr->m_bin[i].toString() << endl
	//		fout << l << "\t" << wr->m_bin[i].eval(361) << "\t" << wr->m_bin[i].eval(596) << "\t" << endl;  // TODO
	//		l += d;
	//	}
	//}

	//void addWaveform(const PointCloudWorkResult *wr) {
	//	//cout << "addWaveform" << endl;
	//	m_waveforms[wr->m_i] = std::vector<Spectrum>(wr->m_bin);
	//}

	//void outputWaveformToOneFile(std::string outName = "accumulation.txt") {
	//	//cout << "outputWaveformToOneFile" << endl;
	//	fs::path full_path(fs::initial_path());
	//	full_path = fs::system_complete(fs::path(m_scene->getDestinationFile().string(), fs::native));
	//	if (!fs::exists(full_path))
	//	{
	//		bool bRet = fs::create_directories(full_path);
	//		if (false == bRet)
	//		{
	//			cout << "no dir" << endl;
	//		}
	//	}

	//	//cout << m_scene->toString() << endl;
	//	std::string folderPath = m_scene->getDestinationFile().string();

	//	//cout << m_scene->getDestinationFile().string() << endl;

	//	//cout << folderPath + "\\waveform.txt" << endl;

	//	//std::ofstream fout(folderPath + "\\waveform.txt");
	//	std::ofstream fout(folderPath + "\\" + outName);

	//	//std::ostringstream oss;

	//	//cout << folderPath + "\\waveform.txt" << endl;

	//	Float l;
	//	Float d = C * m_rate;

	//	for (auto w : m_waveforms) {
	//		l = m_minRange * 2.0;
	//		for (auto s : w) {
	//			fout << l << "\t" << s.eval(361) << "\t" << s.eval(596) << "\t" << endl;  // TODO
	//			//oss << l << "\t" << s.eval(361) << "\t" << s.eval(596) << "\t" << endl;  // TODO
	//			l += d;
	//		}
	//	}

	//	fout.close();
	//	cout << "call writeFile" << endl;
	//	//writeFile(folderPath + "\\" + outName, oss);
	//}

	void outputPointCloudToOneFile(std::string outName = "cloud.txt") {
		//cout << "outputWaveformToOneFile" << endl;
		fs::path full_path(fs::initial_path());
		full_path = fs::system_complete(fs::path(m_scene->getDestinationFile().string(), fs::native));
		if (!fs::exists(full_path))
		{
			bool bRet = fs::create_directories(full_path);
			if (false == bRet)
			{
				cout << "output point cloud no dir" << endl;
			}
		}

		std::string folderPath = m_scene->getDestinationFile().string();
		std::ofstream fout(folderPath + "\\" + outName);

		//fout << m_pointCloudStringStream.str();

		//for (int i = 0; i < m_pointCloud.size(); i++) {
		//	fout << m_pointCloud[i].x << "\t"
		//		<< m_pointCloud[i].y << "\t"
		//		<< m_pointCloud[i].z << "\t"
		//		<< m_pointCloud[i].a << "\t"
		//		<< m_pointCloud[i].z << endl;
		//}

		DiscretePoint p;
		while (m_pointCloud.try_pop(p)) {
			fout << p.x << "\t" << p.y << "\t" << p.z << "\t" << p.a << "\t" << p.i << endl;
		}

		fout.close();
		cout << "call writeFile" << endl;
	}

	//void writeFile(const std::string &filename, const std::ostringstream &oss) {
	//	HANDLE hFile, hMapFile;
	//	LPVOID lpMapAddress;
	//	DWORD dFileSize = oss.str().size();
	//	cout << "string size = " << oss.str().size() << endl;

	//	hFile = CreateFile(filename.c_str(), /* file name */
	//		GENERIC_READ | GENERIC_WRITE, /* read/write access */
	//		FILE_SHARE_READ | FILE_SHARE_WRITE, /* no sharing of the file */
	//		NULL, /* default security */
	//		OPEN_ALWAYS, /* open new or existing file */
	//		FILE_ATTRIBUTE_NORMAL, /* routine file attributes */
	//		NULL); /* no file template */

	//	hMapFile = CreateFileMapping(hFile, /* file handle */
	//		NULL, /* default security */
	//		PAGE_READWRITE, /* read/write access to mapped pages */
	//		0, /* map entire file */
	//		dFileSize,
	//		TEXT("SharedObject")); /* named shared memory object */

	//	if (hFile == INVALID_HANDLE_VALUE)
	//	{
	//		puts("File handle invalid.");
	//	}

	//	lpMapAddress = MapViewOfFile(hMapFile, /* mapped object handle */
	//		FILE_MAP_ALL_ACCESS, /* read/write access */
	//		0, /* mapped view of entire file */
	//		0,
	//		dFileSize);

	//	if (lpMapAddress == NULL) {
	//		puts("Map failed");
	//	}

	//	/* write to shared memory */
	//	sprintf((char *)lpMapAddress, oss.str().c_str());
	//	UnmapViewOfFile(lpMapAddress);
	//	CloseHandle(hFile);
	//	CloseHandle(hMapFile);
	//}

	MTS_DECLARE_CLASS()
public:
	int m_pos;
	size_t m_numOfPulses;
	std::vector<Float> m_x;
	std::vector<Float> m_y;
	std::vector<Float> m_z;
	std::vector<Float> m_u;
	std::vector<Float> m_v;
	std::vector<Float> m_w;
	ref<Scene> m_scene;

	std::ostringstream m_pointCloudStringStream;
	gdface::mt::threadsafe_queue<DiscretePoint> m_pointCloud;
	//std::vector<std::vector<Spectrum> > m_waveforms;

	Float m_sigmaSquareOfBeam;
	Spectrum m_weightTotal;
	Float m_cosFov;
	size_t m_numOfBins;

	vector<double> m_pulse;

	// from XML

	int m_maxDepth;
	int m_axialDivision;

	Float m_fp;
	Float m_fov;
	Float m_pulseEnergy;
	Float m_rate;
	Float m_area;
	Float m_minRange;
	Float m_maxRange;

	std::string m_outputPath;

};

MTS_IMPLEMENT_CLASS(PointCloudWorkUnit, false, WorkUnit)
MTS_IMPLEMENT_CLASS(PointCloudWorkResult, false, WorkResult)
MTS_IMPLEMENT_CLASS_S(PointCloudWorkProcessor, false, WorkProcessor)
MTS_IMPLEMENT_CLASS(PointCloudProcess, false, ParallelProcess)

MTS_NAMESPACE_END
