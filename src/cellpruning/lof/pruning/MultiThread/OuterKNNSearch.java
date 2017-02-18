package cellpruning.lof.pruning.MultiThread;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import cellpruning.lof.pruning.partitionTreeNode;
import metricspace.MetricObject;
import util.SQConfig;

public class OuterKNNSearch implements Runnable {
	private HashMap<Long, MetricObject> TrueKnnPoints;
	private ArrayList<LargeCellStore> large_cell_store;
	private partitionTreeNode ptn;
	private float[] partition_store;
	private int indexStart;
	private int indexEnd;

	public OuterKNNSearch(HashMap<Long, MetricObject> TrueKnnPoints, ArrayList<LargeCellStore> large_cell_store,
			int beginIndex, int endIndex, partitionTreeNode ptn, float[] partition_store) {
		this.TrueKnnPoints = TrueKnnPoints;
		this.large_cell_store = large_cell_store;
		this.ptn = ptn;
		this.partition_store = partition_store;
		this.indexEnd = endIndex;
		this.indexStart = beginIndex;
	}

	@Override
	public void run() {

		for (int i = this.indexStart; i < this.indexEnd; i++) {
			outerKNNSearch(i);
		}
	}

	public void outerKNNSearch(int indexOfLeaveNodes) {
		HashMap<Long, MetricObject> TempTrueKnnPoints = new HashMap<Long, MetricObject>();
		if (large_cell_store.get(indexOfLeaveNodes).isBreakIntoSmallCells()) {
			// find kNNs within the PR quad tree
			large_cell_store.get(indexOfLeaveNodes).findKnnsWithinPRTreeOutsideBucket(TempTrueKnnPoints,
					large_cell_store, indexOfLeaveNodes, ptn, partition_store, SQConfig.K, SQConfig.dims,
					SQConfig.indexOfIndependentDims);

		} else if (large_cell_store.get(indexOfLeaveNodes).getNumOfPoints() != 0) {
			// else find kNNs within the large cell
			large_cell_store.get(indexOfLeaveNodes).findKnnsForLargeCellOutsideBucket(TempTrueKnnPoints,
					large_cell_store, indexOfLeaveNodes, ptn, partition_store, SQConfig.K, SQConfig.dims,
					SQConfig.indexOfIndependentDims);
		}
		synchronized (TrueKnnPoints) {
			TrueKnnPoints.putAll(TempTrueKnnPoints);
		}
	}

}
