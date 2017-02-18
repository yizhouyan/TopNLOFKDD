package mcrtree.metricobject;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import metricspace.*;
import util.PriorityQueue;
import util.SQConfig;

public class MicroCluster {
	private Object clusterCoor;
	private ArrayList<MetricObject> pointList;
	private int numberPoints = 0;
	private float maxRadius = 0;
	private float kdistanceMin = Float.POSITIVE_INFINITY;
	private float kdistaneMax = -1;
	private float lofboundLower = Float.POSITIVE_INFINITY;
	private float lofboundUpper = -1;
	private float reachDistMin = Float.POSITIVE_INFINITY;
	private float reachDistMax = -1;
	private ArrayList<MicroCluster> possibleNeighbors = new ArrayList<MicroCluster>();
	private PriorityQueue pointPQ = new PriorityQueue(PriorityQueue.SORT_ORDER_DESCENDING);

	// private HashSet<Integer> neighborClusterSet;

	public float getPivotKDist() {
		if (pointPQ.size() == SQConfig.K)
			return pointPQ.getPriority();
		else
			return -1;
	}

	public Object getClusterCoor() {
		return clusterCoor;
	}

	public void addPossibleNeighbors(MicroCluster mc) {
		possibleNeighbors.add(mc);
	}

	public void setClusterCoor(Object clusterCoor) {
		this.clusterCoor = clusterCoor;
	}

	public ArrayList<MetricObject> getPointList() {
		return pointList;
	}

	public void setPointList(ArrayList<MetricObject> pointList) {
		this.pointList = pointList;
	}

	public int getNumberPoints() {
		return numberPoints;
	}

	public void setNumberPoints(int numberPoints) {
		this.numberPoints = numberPoints;
	}

	public float getMaxRadius() {
		return maxRadius;
	}

	public void setMaxRadius(float maxRadius) {
		this.maxRadius = maxRadius;
	}

	public MicroCluster(String objectInfo, IMetricSpace metricSpace) {
		clusterCoor = metricSpace.readObject(objectInfo, SQConfig.dims);
		this.pointList = new ArrayList<MetricObject>();
	}

	public void addPointToCurrentCluster(MetricObject mo, float distToPivot) {
		this.pointList.add(mo);
		this.maxRadius = Math.max(maxRadius, distToPivot);
		this.numberPoints++;
		if (pointPQ.size() < SQConfig.K) {
			pointPQ.insert(((Record) mo.getObj()).getRId(), distToPivot);
		} else {
			if (distToPivot < pointPQ.getPriority()) {
				pointPQ.pop();
				pointPQ.insert(((Record) mo.getObj()).getRId(), distToPivot);
			}
		}
		mo.setDistToPivot(distToPivot);
	}

	/**
	 * Compute Kdistance bound for this MC for each object p in this MC, do get
	 * kmin(p) and kmax(p)
	 * 
	 * @param mcList
	 * @param currentId
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void computeKdistanceBound(HashMap<Integer, MicroCluster> mcList, int currentId, IMetricSpace metricSpace,
			IMetric metric) throws IOException, InterruptedException {
		HashMap<Integer, Boolean> overLapWithCurrentMC = new HashMap<Integer, Boolean>();
		for (Map.Entry<Integer, MicroCluster> entry : mcList.entrySet()) {
			if (entry.getKey() == currentId)
				continue;
			// check if those two microcluster overlaps
			boolean overlap = overLapMCs(mcList.get(currentId), entry.getValue(), metricSpace);
			overLapWithCurrentMC.put(entry.getKey(), overlap);
		}
		// sort point list by distance to pivot, first traverse current to get
		// minDist and maxDist, then update with distance to other neighbor
		// clusters
		if (this.pointList.size() > SQConfig.K) {
			Collections.sort(pointList, new Comparator<MetricObject>() {
				public int compare(MetricObject map1, MetricObject map2) {
					if (map2.getDistToPivot() > map1.getDistToPivot())
						return 1;
					else if (map2.getDistToPivot() < map1.getDistToPivot())
						return -1;
					else
						return 0;
				}
			});

			for (int i = 0; i < pointList.size(); i++) {
				MetricObject o_S = pointList.get(i);
				HashMap<Long, Float> localNeighbors = findKNNForSingleObject(o_S, i, pointList, metric, metricSpace);
				float minDist = o_S.MinDistToMicroClustersNoLocal(mcList, currentId, overLapWithCurrentMC, metricSpace,
						localNeighbors);
				float maxDist = o_S.MaxDistToMicroClustersNoLocal(mcList, currentId, metricSpace, localNeighbors);
				if (minDist < this.kdistanceMin)
					this.kdistanceMin = minDist;
				if (maxDist > this.kdistaneMax)
					this.kdistaneMax = maxDist;
			}
		} else {
			for (MetricObject mo : this.pointList) {
				float minDist = mo.MinDistToMicroClusters(mcList, currentId, overLapWithCurrentMC, metricSpace);
				float maxDist = mo.MaxDistToMicroClusters(mcList, currentId, metricSpace);
				if (minDist < this.kdistanceMin)
					this.kdistanceMin = minDist;
				if (maxDist > this.kdistaneMax)
					this.kdistaneMax = maxDist;
			}
		}
	}

	private HashMap<Long, Float> findKNNForSingleObject(MetricObject o_R, int currentIndex,
			ArrayList<MetricObject> pointList, IMetric metric, IMetricSpace metricSpace)
			throws IOException, InterruptedException {
		float dist;
		PriorityQueue pq = new PriorityQueue(PriorityQueue.SORT_ORDER_DESCENDING);
		float theta = Float.POSITIVE_INFINITY;
		boolean kNNfound = false;
		int inc_current = currentIndex + 1;
		int dec_current = currentIndex - 1;
		float i = 0, j = 0; // i---increase j---decrease
		while ((!kNNfound) && ((inc_current < pointList.size()) || (dec_current >= 0))) {
			// System.out.println("increase: "+ inc_current+"; decrease:
			// "+dec_current);
			if ((inc_current > pointList.size() - 1) && (dec_current < 0))
				break;
			if (inc_current > pointList.size() - 1)
				i = Float.MAX_VALUE;
			if (dec_current < 0)
				j = Float.MAX_VALUE;
			if (i <= j) {
				MetricObject o_S = pointList.get(inc_current);
				dist = metric.dist(o_R.getObj(), o_S.getObj());
				if (pq.size() < SQConfig.K) {
					pq.insert(metricSpace.getID(o_S.getObj()), dist);
					theta = pq.getPriority();
				} else if (dist < theta) {
					pq.pop();
					pq.insert(metricSpace.getID(o_S.getObj()), dist);
					theta = pq.getPriority();
				}
				inc_current += 1;
				i = Math.abs(o_R.getDistToPivot() - o_S.getDistToPivot());
			} else {
				MetricObject o_S = pointList.get(dec_current);
				dist = metric.dist(o_R.getObj(), o_S.getObj());
				if (pq.size() < SQConfig.K) {
					pq.insert(metricSpace.getID(o_S.getObj()), dist);
					theta = pq.getPriority();
				} else if (dist < theta) {
					pq.pop();
					pq.insert(metricSpace.getID(o_S.getObj()), dist);
					theta = pq.getPriority();
				}
				dec_current -= 1;
				j = Math.abs(o_R.getDistToPivot() - o_S.getDistToPivot());
			}
			// System.out.println(pq.getPriority()+","+i+","+j);
			if (i > pq.getPriority() && j > pq.getPriority() && (pq.size() == SQConfig.K))
				kNNfound = true;
		}

		HashMap<Long, Float> kNNInfo = new HashMap<Long, Float>();
		o_R.setKdist(pq.getPriority());
		while (pq.size() > 0) {
			kNNInfo.put(pq.getValue(), pq.getPriority());
			pq.pop();
		}
		o_R.setkNNINfo(kNNInfo);

		return kNNInfo;
	}

	public void computeLofBoundForEachMC() {
		if (this.possibleNeighbors.size() == 0)
			return;
		float minReachLowerForNeighbors = Float.POSITIVE_INFINITY;
		float maxReachUpperForNeighbors = -1;
		for (MicroCluster mc : this.possibleNeighbors) {
			minReachLowerForNeighbors = Math.min(minReachLowerForNeighbors, mc.getReachDistMin());
			maxReachUpperForNeighbors = Math.max(maxReachUpperForNeighbors, mc.getReachDistMax());
		}
		// this.lofboundUpper = this.reachDistMax/minReachLowerForNeighbors;
		// this.lofboundLower =this.reachDistMin/maxReachUpperForNeighbors;
		this.lofboundUpper = Math.max(this.lofboundUpper, this.reachDistMax / minReachLowerForNeighbors);
		this.lofboundLower = Math.min(this.lofboundLower, this.reachDistMin / maxReachUpperForNeighbors);
	}

	public boolean overLapMCs(MicroCluster mc_1, MicroCluster mc_2, IMetricSpace metricSpace) {
		float dist = 0.0f;
		try {
			dist = metricSpace.compDist(mc_1.getClusterCoor(), mc_2.getClusterCoor());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if (dist < mc_1.getMaxRadius() + mc_2.getMaxRadius())
			return true;
		return false;
	}

	public void printMC() {
		System.out.println("Pivot Coor: " + ((Record) this.clusterCoor).toString());
		System.out.println("Number of points inside: " + this.numberPoints);
		System.out.println("Max Radius: " + this.maxRadius);
	}

	public void computeInternalMCBound() {
		if (this.numberPoints <= 1)
			return;
		float rmax = Math.max(2 * this.maxRadius, this.kdistaneMax);
		float rmin = kdistanceMin;
		this.reachDistMin = Math.min(rmin, this.reachDistMin);
		this.reachDistMax = Math.max(rmax, this.reachDistMax);
		this.lofboundUpper = (float) (rmax / rmin * 1.0);
		this.lofboundLower = (float) (rmin / rmax * 1.0);
	}

	public float maxDistBetweenTwoMCs(MicroCluster mc, IMetricSpace metricSpace) throws IOException {
		float dist = metricSpace.compDist(mc.clusterCoor, this.clusterCoor);
		return dist + this.maxRadius + mc.maxRadius;
	}

	public float minDistBetweenTwoMCs(MicroCluster mc, IMetricSpace metricSpace) throws IOException {
		float dist = metricSpace.compDist(mc.clusterCoor, this.clusterCoor);
		return (float) Math.max(0, dist - this.maxRadius - mc.maxRadius);
	}

	public void computeExternalMCBound(HashMap<Integer, MicroCluster> mcList, int currentId, IMetricSpace metricSpace)
			throws IOException {

		for (MicroCluster mc : this.possibleNeighbors) {
			float rmax = Math.max(this.maxDistBetweenTwoMCs(mc, metricSpace), mc.kdistaneMax);
			float rmin = Math.max(this.minDistBetweenTwoMCs(mc, metricSpace), mc.kdistanceMin);
			this.reachDistMin = Math.min(rmin, this.reachDistMin);
			this.reachDistMax = Math.max(rmax, this.reachDistMax);
			// if (this.lofboundUpper < rmax / rmin)
			// this.lofboundUpper = rmax / rmin;
			// if (this.lofboundLower > rmin / rmax)
			// this.lofboundLower = rmin / rmax;
		}
	}

	public float getKdistanceMin() {
		return kdistanceMin;
	}

	public void setKdistanceMin(float kdistanceMin) {
		this.kdistanceMin = kdistanceMin;
	}

	public float getKdistaneMax() {
		return kdistaneMax;
	}

	public void setKdistaneMax(float kdistaneMax) {
		this.kdistaneMax = kdistaneMax;
	}

	public float getLofboundLower() {
		return lofboundLower;
	}

	public void setLofboundLower(float lofboundLower) {
		this.lofboundLower = lofboundLower;
	}

	public float getLofboundUpper() {
		return lofboundUpper;
	}

	public void setLofboundUpper(float lofboundUpper) {
		this.lofboundUpper = lofboundUpper;
	}

	public float getReachDistMin() {
		return reachDistMin;
	}

	public void setReachDistMin(float reachDistMin) {
		this.reachDistMin = reachDistMin;
	}

	public float getReachDistMax() {
		return reachDistMax;
	}

	public void setReachDistMax(float reachDistMax) {
		this.reachDistMax = reachDistMax;
	}

	public ArrayList<MicroCluster> getPossibleNeighbors() {
		return possibleNeighbors;
	}

	public void setPossibleNeighbors(ArrayList<MicroCluster> possibleNeighbors) {
		this.possibleNeighbors = possibleNeighbors;
	}
}
