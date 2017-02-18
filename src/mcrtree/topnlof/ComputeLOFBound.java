package mcrtree.topnlof;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import gnu.trove.procedure.TIntProcedure;

import mcrtree.metricobject.MicroCluster;
import metricspace.*;
import util.PriorityQueue;
import util.SQConfig;

public class ComputeLOFBound {
	private HashMap<Integer, MicroCluster> mcList;
	
	// private HashSet<Integer> prunedMCsIndex;

	public ComputeLOFBound(HashMap<Integer, MicroCluster> mcList) {
		this.mcList = mcList;
	}

	public void ComputeKDistanceBoundForMCs(IMetricSpace metricSpace, IMetric metric) throws IOException, InterruptedException {
		for (Map.Entry<Integer, MicroCluster> mc : this.mcList.entrySet()) {
			// for each microcluster, traverse point inside
			mc.getValue().computeKdistanceBound(mcList, mc.getKey(), metricSpace, metric);
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
