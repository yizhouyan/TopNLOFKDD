package cellpruning.lof.pruning.MultiThread;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import cellpruning.lof.pruning.partitionTreeNode;
import metricspace.IMetricSpace;
import metricspace.MetricObject;
import util.PriorityQueue;
import util.SQConfig;

public class InnerKNNSearch implements Runnable {

	private HashMap<Long, MetricObject> CanPrunePoints;
	private HashMap<Long, MetricObject> TrueKNNPoints;
	private HashMap<Long, MetricObject> lrdHM;
	private IMetricSpace metricSpace = null;
	private int indexStart;
	private int indexEnd;
	private partitionTreeNode ptn;
	private ArrayList<LargeCellStore> leaveNodes;
	private PriorityQueue topnLOF;

	public InnerKNNSearch(HashMap<Long, MetricObject> CanPrunePoints, HashMap<Long, MetricObject> TrueKNNPoints,
			HashMap<Long, MetricObject> lrdHM, int indexStart, int indexEnd, partitionTreeNode ptn,
			ArrayList<LargeCellStore> leaveNodes, PriorityQueue topnLOF, IMetricSpace metricSpace) {
		this.CanPrunePoints = CanPrunePoints;
		this.TrueKNNPoints = TrueKNNPoints;
		this.lrdHM = lrdHM;
		this.indexEnd = indexEnd;
		this.indexStart = indexStart;
		this.ptn = ptn;
		this.leaveNodes = leaveNodes;
		this.topnLOF = topnLOF;
		this.metricSpace = metricSpace;
	}

	@Override
	public void run() {
		
		for (int i = this.indexStart; i < this.indexEnd; i++) {
			try {
				InnerBucketKNNSearch(i);
			} catch (IOException | InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		
	}

	public void InnerBucketKNNSearch(int indexOfLeaveNodes) throws IOException, InterruptedException {
		HashMap<Long, MetricObject> TempCanPrunePoints = new HashMap<Long, MetricObject>();
		HashMap<Long, MetricObject> TempTrueKnnPoints = new HashMap<Long, MetricObject>();
		// start calculating LRD and LOF if possible
		// save those pruned points but need to recompute KNNs
		HashMap<Long, MetricObject> TempneedCalculatePruned = new HashMap<Long, MetricObject>();
		HashMap<Long, MetricObject> TemplrdHM = new HashMap<Long, MetricObject>();
		// save those cannot be pruned only by LRD value...
		HashMap<Long, MetricObject> TempneedCalLOF = new HashMap<Long, MetricObject>();
		// need more knn information, maybe knn is pruned...
		HashMap<Long, MetricObject> TempneedRecalLRD = new HashMap<Long, MetricObject>();
		leaveNodes.get(indexOfLeaveNodes).innerSearchWithEachLargeCell(TempCanPrunePoints, TempTrueKnnPoints, TemplrdHM,
				TempneedCalculatePruned, TempneedCalLOF, TempneedRecalLRD, indexOfLeaveNodes, SQConfig.thresholdLof,
				ptn, leaveNodes, topnLOF);
		
		synchronized (CanPrunePoints) {
			CanPrunePoints.putAll(TempCanPrunePoints);
		}
		synchronized (TrueKNNPoints) {
			TrueKNNPoints.putAll(TempTrueKnnPoints);
		}
		synchronized (lrdHM) {
			lrdHM.putAll(TemplrdHM);
		}

	}

}
