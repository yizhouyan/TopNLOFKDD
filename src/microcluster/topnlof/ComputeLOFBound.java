package microcluster.topnlof;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import microcluster.metricobject.MicroCluster;
import gnu.trove.procedure.TIntProcedure;
import metricspace.*;
import net.sf.jsi.Rectangle;
import net.sf.jsi.SpatialIndex;
import net.sf.jsi.rtree.RTree;
import util.PriorityQueue;
import util.SQConfig;

public class ComputeLOFBound {
	private HashMap<Integer, MicroCluster> mcList;
	SpatialIndex si;
	// private HashSet<Integer> prunedMCsIndex;

	public ComputeLOFBound(HashMap<Integer, MicroCluster> mcList) {
		this.mcList = mcList;
		this.buildRTreeForMCs();
	}

	private void buildRTreeForMCs() {
		this.si = new RTree();
		si.init(null);
		for (Map.Entry<Integer, MicroCluster> mc : mcList.entrySet()) {
			float[] coor = ((Record) mc.getValue().getClusterCoor()).getValue();
			float[] coorForR = { Math.max(0, coor[0] - mc.getValue().getMaxRadius()),
					Math.max(0, coor[1] - mc.getValue().getMaxRadius()),
					(float) Math.min(SQConfig.domainRange, coor[0] + mc.getValue().getMaxRadius()),
					(float) Math.min(SQConfig.domainRange, coor[1] + mc.getValue().getMaxRadius()) };
			Rectangle tempR = new Rectangle(coorForR[0], coorForR[1], coorForR[2], coorForR[3]);
			si.add(tempR, mc.getKey());
		}
	}

	public Set<Integer> ComputePossibleNeighborsForMC(MicroCluster mc) {
		if (mc.getNumberPoints() > SQConfig.K) {
			final HashSet<Integer> listOfNeighbors = new HashSet<Integer>();
			// find possible neighbors using 2*maxR+ kdistance
			float distRange = 2 * mc.getMaxRadius() + mc.getPivotKDist();
			float[] pivotCoor = ((Record) mc.getClusterCoor()).getValue();
			float[] recRange = { pivotCoor[0] - distRange, pivotCoor[1] - distRange, pivotCoor[0] + distRange,
					pivotCoor[1] + distRange };
			this.si.intersects(new Rectangle(recRange[0], recRange[1], recRange[2], recRange[3]), new TIntProcedure() {
				// be called with the results
				public boolean execute(int i) {
					listOfNeighbors.add(i);
					return true; // return true here to continue receiving
									// results
				}
			});

			return listOfNeighbors;
		} else
			return this.mcList.keySet();
	}

	public void ComputeKDistanceBoundForMCs(IMetricSpace metricSpace, IMetric metric) throws IOException, InterruptedException {
		for (Map.Entry<Integer, MicroCluster> mc : this.mcList.entrySet()) {
			// for each microcluster, traverse point inside
			mc.getValue().computeKdistanceBound(mcList, mc.getKey(), metricSpace, metric,
					this.ComputePossibleNeighborsForMC(mc.getValue()));
		}
	}

	public void ComputeFinalLOFBound(IMetricSpace metricSpace) throws IOException {
		// compute reachability distance for each micro cluster
		for (Map.Entry<Integer, MicroCluster> mc : this.mcList.entrySet()) {
			// for each microcluster, traverse point inside
			mc.getValue().computeInternalMCBound();
			mc.getValue().computeExternalMCBound(mcList, mc.getKey(), metricSpace);
			// System.out.println("MC: " + mc.getKey() + " LOF Upper Bound: " +
			// mc.getValue().getLofboundUpper()
			// + " LOF Lower Bound: " + mc.getValue().getLofboundLower());
		}
		// update lof bound for each micro cluster
		for (Map.Entry<Integer, MicroCluster> mc : this.mcList.entrySet()) {
			mc.getValue().computeLofBoundForEachMC();
		}
	}

	public void rankTopNOutliers(HashSet<Integer> prunedMCsIndex) {
		int numOfPrunedPoints = 0;
		PriorityQueue MicroClusterPQ = new PriorityQueue(PriorityQueue.SORT_ORDER_ASCENDING);
		for (Map.Entry<Integer, MicroCluster> mc : this.mcList.entrySet()) {
			if (MicroClusterPQ.size() < SQConfig.TOPN) {
				MicroClusterPQ.insert(mc.getKey(), mc.getValue().getLofboundLower());
			} else {
				if (mc.getValue().getLofboundUpper() < MicroClusterPQ.getPriority()){
					prunedMCsIndex.add(mc.getKey());
					numOfPrunedPoints += mc.getValue().getNumberPoints();
				}
				else if (mc.getValue().getLofboundLower() > MicroClusterPQ.getPriority()) {
					MicroClusterPQ.pop();
					MicroClusterPQ.insert(mc.getKey(), mc.getValue().getLofboundLower());
				}
			}
		}

		System.out.println("Number of Pruned Clusters: " + prunedMCsIndex.size());
		System.out.println("Number of Pruned Points: " + numOfPrunedPoints);
		// while(MicroClusterPQ.size()!=0){
		// System.out.println("MC: " + MicroClusterPQ.getValue() + "," +
		// MicroClusterPQ.getPriority());
		// MicroClusterPQ.pop();
		// }
	}
}
