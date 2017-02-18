package baseline.pivotknn;
import java.util.HashMap;


@SuppressWarnings("rawtypes")
public class MetricObject implements Comparable {
	private Object obj;
	private HashMap<Long, Float> kNNINfo;
	private float distToPivot = 0.0f;

	private float kdist = 0;
	private float lrd = 0;
	private float lof = 0;
	
	public MetricObject(Object obj) {
		this.obj = obj;
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

	/**
	 * sort by the descending order
	 */
	@Override
	public int compareTo(Object o) {
		MetricObject other = (MetricObject) o;
		if (other.distToPivot > this.distToPivot)
			return 1;
		else if (other.distToPivot < this.distToPivot)
			return -1;
		else
			return 0;
	}

	public float getDistToPivot() {
		return distToPivot;
	}

	public void setDistToPivot(float distToPivot) {
		this.distToPivot = distToPivot;
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

	public HashMap<Long, Float> getkNNINfo() {
		return kNNINfo;
	}

	public void setkNNINfo(HashMap<Long, Float> kNNINfo) {
		this.kNNINfo = kNNINfo;
	}
}
