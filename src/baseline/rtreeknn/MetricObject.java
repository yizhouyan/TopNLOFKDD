package baseline.rtreeknn;

import java.util.HashMap;

import util.SQConfig;

@SuppressWarnings("rawtypes")
public class MetricObject {
	private Object obj;
	private HashMap<Integer, Float> kNNINfo = new HashMap<Integer, Float>();

	private float kdist = 0;
	private float lrd = 0;
	private float lof = 0;

	public MetricObject(Object obj) {
		this.obj = obj;
	}

	public float computeKdist() {
		for (Float dist : kNNINfo.values()) {
			this.kdist = Math.max(this.kdist, dist);
		}
		return this.kdist;
	}

	public void addKNN(int neighborId, float distToNeighbor) {
		if (kNNINfo.size() >= SQConfig.K)
			return;
		kNNINfo.put(neighborId, distToNeighbor);
	}

	public float getLrd() {
		return lrd;
	}

	public void setLrd(float lrd) {
		this.lrd = lrd;
	}

	public float getLof() {
		return lof;
	}

	public void setLof(float lof) {
		this.lof = lof;
	}

	public Object getObj() {
		return obj;
	}

	public void setObj(Object obj) {
		this.obj = obj;
	}

	public float getKdist() {
		return kdist;
	}

	public void setKdist(float kdist) {
		this.kdist = kdist;
	}

	public HashMap<Integer, Float> getkNNINfo() {
		return kNNINfo;
	}

	public void setkNNINfo(HashMap<Integer, Float> kNNINfo) {
		this.kNNINfo = kNNINfo;
	}
}
