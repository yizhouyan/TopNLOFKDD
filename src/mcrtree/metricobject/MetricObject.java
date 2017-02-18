package mcrtree.metricobject;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import metricspace.*;
import metricspace.MetricFactory.L2Metric;
import util.PriorityQueue;
import util.SQConfig;

@SuppressWarnings("rawtypes")
public class MetricObject implements Comparable {
	private Object obj;
	private HashMap<Long, Float> kNNINfo;
	// public PriorityQueue pointPQ = new
	// PriorityQueue(PriorityQueue.SORT_ORDER_DESCENDING);
	private float distToPivot = 0.0f;

	private float kdist = 0;
	private float lrd = 0;
	private float lof = 0;
	private boolean prune = false;
	private boolean trueKNN = false;

	public MetricObject(Object obj) {
		this.obj = obj;
	}

	// public PriorityQueue getPointPQ() {
	// return pointPQ;
	// }
	//
	// public void setPointPQ(PriorityQueue pointPQ) {
	// this.pointPQ = pointPQ;
	// }

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

	public float maxDistToMicroCluster(MicroCluster mc, IMetricSpace metricSpace) throws IOException {
		float distToMC = 0.0f;
		distToMC = metricSpace.compDist(this.obj, mc.getClusterCoor());
		return distToMC + mc.getMaxRadius();
	}

	public float minDistToMicroClusterNoOverlap(MicroCluster mc, IMetricSpace metricSpace) throws IOException {
		float distToMC = 0.0f;
		distToMC = metricSpace.compDist(this.obj, mc.getClusterCoor());
		return distToMC - mc.getMaxRadius();
	}

	public float minDistToMicroClusterOverlap(float[] mc_1, float[] mc_2) {
		return getDistToHP(((Record) this.obj).getValue(), mc_1, mc_2);
	}

	/**
	 * Compute the distance from a point to the hyperplane between two pivots.
	 * 
	 * @param x
	 *            the query point.
	 * @param p_x
	 *            the pivot for x's partition.
	 * @param p_i
	 *            neighboring pivot.
	 * @return distance.
	 */
	public static float getDistToHP(float[] x, float[] p_x, float[] p_i) throws IllegalArgumentException {

		if (Arrays.equals(p_x, p_i)) {
			throw new IllegalArgumentException();
		}
		// find point on hyperplane
		float[] midpoint = new float[p_x.length];
		for (int i = 0; i < p_x.length; i++) {
			midpoint[i] = (p_x[i] + p_i[i]) / 2;
		}

		// find normal vector (a,b,c)
		float[] norm = new float[p_x.length];
		for (int i = 0; i < p_x.length; i++) {
			norm[i] = (p_i[i] - p_x[i]);
		}

		// find d term (plane in 3d = ax+by+cz+d =0)
		float d = 0f;
		for (int i = 0; i < p_x.length; i++) {
			d += norm[i] * midpoint[i];
		}
		d *= -1;

		// find distance
		// query point x = (x,y,z)
		// D = |ax + by + cz + d| / sqrt(a^2+b^2+c^2)
		float denom = 0f;

		for (int i = 0; i < p_x.length; i++) {
			denom += Math.pow(norm[i], 2);
		}
		denom = (float) Math.sqrt(denom);

		float dist = 0f;
		for (int i = 0; i < p_x.length; i++) {
			dist += norm[i] * x[i];
		}
		dist = Math.abs(dist + d);

		return Math.abs(dist / denom);
	}

	/**
	 * Compute the lower bound on the k-distance of p compute DistMin(p, MCs)
	 * and then sort, get the first
	 * 
	 * @return
	 * @throws IOException
	 */
	public float MinDistToMicroClusters(HashMap<Integer, MicroCluster> mcList, int currentMCId,
			HashMap<Integer, Boolean> overLapWithCurrentMC, IMetricSpace metricSpace) throws IOException {
		ArrayList<Float> minDistanceToMCs = new ArrayList<Float>();
		ArrayList<Integer> indexOfMCs = new ArrayList<Integer>();
		// ArrayList<Integer> minIndexOfMCs = new ArrayList<Integer>();
		// get min and max distance from the point to the corresponding Micro
		// Cluster
		for (Map.Entry<Integer, MicroCluster> entry : mcList.entrySet()) {
			if (entry.getKey() == currentMCId) {
				for (MetricObject mo : entry.getValue().getPointList()) {
					if (((Record) mo.getObj()).getRId() == ((Record) this.getObj()).getRId())
						continue;
					float tempMinDist = metricSpace.compDist(mo.getObj(), this.getObj());
					minDistanceToMCs.add(tempMinDist);
					indexOfMCs.add(entry.getKey());
				}
				// continue;
			} else {
				boolean overlap = overLapWithCurrentMC.get(entry.getKey());
				if (overlap) {
					float tempMinDist = this.minDistToMicroClusterOverlap(
							((Record) entry.getValue().getClusterCoor()).getValue(),
							((Record) mcList.get(currentMCId).getClusterCoor()).getValue());
					// System.out.println("OverLap Cluster " + tempMinDist);
					minDistanceToMCs.add(tempMinDist);
					indexOfMCs.add(entry.getKey());

				} else {
					float tempMinDist = this.minDistToMicroClusterNoOverlap(entry.getValue(), metricSpace);
					minDistanceToMCs.add(tempMinDist);
					indexOfMCs.add(entry.getKey());
				}
			}
		}
		float getKthMinDist = sortAndGetTopNMinDist(minDistanceToMCs, indexOfMCs, mcList, currentMCId, false);

		// System.out.println("K-th distance: " + getKthMinDist);
		return getKthMinDist;
	}

	public float MinDistToMicroClustersNoLocal(HashMap<Integer, MicroCluster> mcList, int currentMCId,
			HashMap<Integer, Boolean> overLapWithCurrentMC, IMetricSpace metricSpace,
			HashMap<Long, Float> localNeighbors) throws IOException {
		ArrayList<Float> minDistanceToMCs = new ArrayList<Float>();
		ArrayList<Integer> indexOfMCs = new ArrayList<Integer>();

		// load local info into
		for (Map.Entry<Long, Float> temp : localNeighbors.entrySet()) {
			minDistanceToMCs.add(temp.getValue());
			indexOfMCs.add(currentMCId);
		}

		for (Map.Entry<Integer, MicroCluster> entry : mcList.entrySet()) {
			if (entry.getKey() == currentMCId) {
				continue;
			} else {
				boolean overlap = overLapWithCurrentMC.get(entry.getKey());
				if (overlap) {
					float tempMinDist = this.minDistToMicroClusterOverlap(
							((Record) entry.getValue().getClusterCoor()).getValue(),
							((Record) mcList.get(currentMCId).getClusterCoor()).getValue());
					// System.out.println("OverLap Cluster " + tempMinDist);
					minDistanceToMCs.add(tempMinDist);
					indexOfMCs.add(entry.getKey());

				} else {
					float tempMinDist = this.minDistToMicroClusterNoOverlap(entry.getValue(), metricSpace);
					minDistanceToMCs.add(tempMinDist);
					indexOfMCs.add(entry.getKey());
				}
			}
		}

		float getKthMinDist = sortAndGetTopNMinDist(minDistanceToMCs, indexOfMCs, mcList, currentMCId, true);

		// System.out.println("K-th distance: " + getKthMinDist);
		return getKthMinDist;
	}

	public float sortAndGetTopNMinDist(ArrayList<Float> minDistanceToMCs, ArrayList<Integer> indexOfMCs,
			HashMap<Integer, MicroCluster> mcList, int currentMCId, boolean haveLocalKNN) {
		int countTotalNumPoints = 0;
		float totalMinDist = -1;
		boolean trueKNN = true;
		while (countTotalNumPoints < SQConfig.K) {
			int indexOfMin = -1;
			float minDist = Float.POSITIVE_INFINITY;
			// float minDist = -1;
			// find the minimal one from the list
			for (int i = 0; i < minDistanceToMCs.size(); i++) {
				// if(minDistanceToMCs.get(i)> minDist){
				if (minDistanceToMCs.get(i) < minDist) {
					indexOfMin = i;
					minDist = minDistanceToMCs.get(i);
				}
			}
			if (indexOfMCs.get(indexOfMin) == currentMCId)
				countTotalNumPoints += 1;
			else {
				countTotalNumPoints += mcList.get(indexOfMCs.get(indexOfMin)).getNumberPoints();
				mcList.get(currentMCId).addPossibleNeighbors(mcList.get(indexOfMCs.get(indexOfMin)));
				trueKNN = false;
			}
			minDistanceToMCs.remove(indexOfMin);
			indexOfMCs.remove(indexOfMin);
			totalMinDist = minDist;
		}
		if (haveLocalKNN)
			this.setTrueKNN(trueKNN);
		return totalMinDist;
	}

	public float sortAndGetTopNMaxDist(ArrayList<Float> maxDistanceToMCs, ArrayList<Integer> indexOfMCs,
			HashMap<Integer, MicroCluster> mcList, int currentMCId) {
		int countTotalNumPoints = 0;
		float totalMaxDist = -1;
		while (countTotalNumPoints < SQConfig.K) {
			int indexOfMax = -1;
			float minDist = Float.POSITIVE_INFINITY;
			// find the minimal one from the list
			for (int i = 0; i < maxDistanceToMCs.size(); i++) {
				if (maxDistanceToMCs.get(i) < minDist) {
					indexOfMax = i;
					minDist = maxDistanceToMCs.get(i);
				}
			}
			if (indexOfMCs.get(indexOfMax) == currentMCId)
				countTotalNumPoints += 1;
			else {
				countTotalNumPoints += mcList.get(indexOfMCs.get(indexOfMax)).getNumberPoints();
				mcList.get(currentMCId).addPossibleNeighbors(mcList.get(indexOfMCs.get(indexOfMax)));
			}

			maxDistanceToMCs.remove(indexOfMax);
			indexOfMCs.remove(indexOfMax);
			totalMaxDist = minDist;
		}
		return totalMaxDist;
	}

	public float MaxDistToMicroClustersNoLocal(HashMap<Integer, MicroCluster> mcList, int currentMCId,
			IMetricSpace metricSpace, HashMap<Long, Float> localNeighbors) throws IOException {
		ArrayList<Float> maxDistanceToMCs = new ArrayList<Float>();
		ArrayList<Integer> indexOfMCs = new ArrayList<Integer>();
		// load local info into
		for (Map.Entry<Long, Float> temp : localNeighbors.entrySet()) {
			maxDistanceToMCs.add(temp.getValue());
			indexOfMCs.add(currentMCId);
		}

		// get min and max distance from the point to the corresponding Micro
		// Cluster
		for (Map.Entry<Integer, MicroCluster> entry : mcList.entrySet()) {
			if (entry.getKey() == currentMCId) {
				continue;
			} else {
				float tempMaxDist = this.maxDistToMicroCluster(entry.getValue(), metricSpace);
				maxDistanceToMCs.add(tempMaxDist);
				indexOfMCs.add(entry.getKey());
			}
		}

		float getKthMaxDist = sortAndGetTopNMaxDist(maxDistanceToMCs, indexOfMCs, mcList, currentMCId);
		return getKthMaxDist;
	}

	public float MaxDistToMicroClusters(HashMap<Integer, MicroCluster> mcList, int currentMCId,
			IMetricSpace metricSpace) throws IOException {
		ArrayList<Float> maxDistanceToMCs = new ArrayList<Float>();
		ArrayList<Integer> indexOfMCs = new ArrayList<Integer>();
		// get min and max distance from the point to the corresponding Micro
		// Cluster
		for (Map.Entry<Integer, MicroCluster> entry : mcList.entrySet()) {
			if (entry.getKey() == currentMCId) {
				for (MetricObject mo : entry.getValue().getPointList()) {
					if (((Record) mo.getObj()).getRId() == ((Record) this.getObj()).getRId())
						continue;
					float tempMaxDist = metricSpace.compDist(mo.getObj(), this.getObj());
					maxDistanceToMCs.add(tempMaxDist);
					indexOfMCs.add(entry.getKey());
				}
				// continue;
			} else {
				float tempMaxDist = this.maxDistToMicroCluster(entry.getValue(), metricSpace);
				maxDistanceToMCs.add(tempMaxDist);
				indexOfMCs.add(entry.getKey());
			}
		}
		// System.out.println("Max Distance from Point to MCs: ");
		// for (int i = 0; i < maxDistanceToMCs.size(); i++) {
		// System.out.println("To MC " + maxIndexOfMCs.get(i) + " distance is "
		// + maxDistanceToMCs.get(i));
		// }
		float getKthMaxDist = sortAndGetTopNMaxDist(maxDistanceToMCs, indexOfMCs, mcList, currentMCId);
		return getKthMaxDist;
	}

	public boolean isPrune() {
		return prune;
	}

	public void setPrune(boolean prune) {
		this.prune = prune;
	}

	public HashMap<Long, Float> getkNNINfo() {
		return kNNINfo;
	}

	public void setkNNINfo(HashMap<Long, Float> kNNINfo) {
		this.kNNINfo = kNNINfo;
	}

	public boolean isTrueKNN() {
		return trueKNN;
	}

	public void setTrueKNN(boolean trueKNN) {
		this.trueKNN = trueKNN;
	}
}
