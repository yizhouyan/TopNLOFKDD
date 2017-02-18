package cellpruning.lof.pruning.MultiThread;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import metricspace.MetricObject;
import metricspace.Record;
import util.SQConfig;

public class ComputeLRD {
	/**
	 * Calculate LRD if possible
	 * 
	 * @param context
	 * @param o_S
	 * @param TrueKnnPoints
	 * @param CanPrunePoints
	 * @param needCalculatePruned
	 * @param lrdHM
	 * @param needCalLOF
	 * @param threshold
	 * @return
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static boolean CalLRDForSingleObject(MetricObject o_S, HashMap<Long, MetricObject> TrueKnnPoints,
			HashMap<Long, MetricObject> CanPrunePoints, HashMap<Long, MetricObject> needCalculatePruned,
			HashMap<Long, MetricObject> lrdHM, HashMap<Long, MetricObject> needCalLOF,
			HashMap<Long, MetricObject> needRecalLRD, float threshold, ArrayList<LargeCellStore> leaveNodes, int K)
			throws IOException, InterruptedException {
		float lrd_core = 0.0f;
		float reachdistMax = 0.0f;
		float minNNtoNN = Float.POSITIVE_INFINITY;
		// boolean canLRD = true;
		HashMap<Long, MetricObject> tempNeedCalculatePruned = new HashMap<Long, MetricObject>();

		long[] KNN_moObjectsID = o_S.getPointPQ().getValueSet();
		float[] moDistToKNN = o_S.getPointPQ().getPrioritySet();
		boolean canComLRD = true;
		for (int i = 0; i < moDistToKNN.length; i++) {
			long temp_kNNKey = KNN_moObjectsID[i];
			float temp_dist = moDistToKNN[i];
			if (!TrueKnnPoints.containsKey(temp_kNNKey)) {
				if (CanPrunePoints.containsKey(temp_kNNKey)) {
					tempNeedCalculatePruned.put(temp_kNNKey, CanPrunePoints.get(temp_kNNKey));
					reachdistMax = Math.max(temp_dist, reachdistMax);
					minNNtoNN = Math.min(minNNtoNN,
							leaveNodes.get(CanPrunePoints.get(temp_kNNKey).getIndexOfCPCellInList()).getCpDist());
				}
				canComLRD = false;
				continue;
			}
			float temp_reach_dist = Math.max(temp_dist, TrueKnnPoints.get(temp_kNNKey).getKdist());
			reachdistMax = Math.max(reachdistMax, temp_reach_dist);
			minNNtoNN = Math.min(minNNtoNN, TrueKnnPoints.get(temp_kNNKey).getNearestNeighborDist());
			lrd_core += temp_reach_dist;
		}
		// if (!canLRD) {
		// needRecalLRD.put(((Record) o_S.getObj()).getRId(), o_S);
		// return false;
		// }
		if (tempNeedCalculatePruned.isEmpty() && canComLRD) {
			lrd_core = 1.0f / (lrd_core / K * 1.0f);
			o_S.setLrdValue(lrd_core);
			o_S.setType('L');
			// calculate if this can prune? if can prune, then don't
			// calculate
			// lof
			float predictedLOF = reachdistMax / minNNtoNN;
			lrdHM.put(((Record) o_S.getObj()).getRId(), o_S);
			if (predictedLOF <= threshold) {
				o_S.setCanPrune(true);
				SQConfig.countPointBasedPruned++;
			} else
				needCalLOF.put(((Record) o_S.getObj()).getRId(), o_S);
			return true;
		} else {
			needCalculatePruned.putAll(tempNeedCalculatePruned);
			needRecalLRD.put(((Record) o_S.getObj()).getRId(), o_S);
			return false;
		}
	}

	/**
	 * Calculate LRD for some points that knns pruned
	 * 
	 * @param context
	 * @param o_S
	 * @param TrueKnnPoints
	 * @param needCalculatePruned
	 * @param lrdHM
	 * @param needCalLOF
	 * @param threshold
	 * @return
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static boolean ReCalLRDForSpecial(MetricObject o_S, HashMap<Long, MetricObject> TrueKnnPoints,
			HashMap<Long, MetricObject> needCalculatePruned, HashMap<Long, MetricObject> lrdHM,
			HashMap<Long, MetricObject> needCalLOF, HashMap<Long, MetricObject> needCalculateLRDPruned, float threshold,
			int K) throws IOException, InterruptedException {
		float lrd_core = 0.0f;
		float reachdistMax = 0.0f;
		float minNNtoNN = Float.POSITIVE_INFINITY;
		int countPruned = 0;

		long[] KNN_moObjectsID = o_S.getPointPQ().getValueSet();
		float[] moDistToKNN = o_S.getPointPQ().getPrioritySet();
		for (int i = 0; i < KNN_moObjectsID.length; i++) {
			long temp_kNNKey = KNN_moObjectsID[i];
			float temp_dist = moDistToKNN[i];
			float temp_reach_dist = 0.0f;
			if (!TrueKnnPoints.containsKey(temp_kNNKey)) {
				if (needCalculatePruned.containsKey(temp_kNNKey)) {
					temp_reach_dist = Math.max(temp_dist, needCalculatePruned.get(temp_kNNKey).getKdist());
					reachdistMax = Math.max(reachdistMax, temp_reach_dist);
					minNNtoNN = Math.min(minNNtoNN, needCalculatePruned.get(temp_kNNKey).getNearestNeighborDist());
				}
				else
					System.out.println("Error here");
			} else {
				temp_reach_dist = Math.max(temp_dist, TrueKnnPoints.get(temp_kNNKey).getKdist());
				reachdistMax = Math.max(reachdistMax, temp_reach_dist);
				minNNtoNN = Math.min(minNNtoNN, TrueKnnPoints.get(temp_kNNKey).getNearestNeighborDist());
			}
			lrd_core += temp_reach_dist;
		}
		lrd_core = 1.0f / (lrd_core / K * 1.0f);
		o_S.setLrdValue(lrd_core);
		o_S.setType('L');
		// calculate if this can prune? if can prune, then don't calculate
		// lof
		float predictedLOF = reachdistMax / minNNtoNN;
		lrdHM.put(((Record) o_S.getObj()).getRId(), o_S);
		if (predictedLOF <= threshold) {
			o_S.setCanPrune(true);
			SQConfig.countPointBasedPruned++;
		} else {
			needCalLOF.put(((Record) o_S.getObj()).getRId(), o_S);
			// long[] KNN_moObjectsID = o_S.getPointPQ().getValueSet();
			for (int i = 0; i < KNN_moObjectsID.length; i++) {
				long temp_kNNKey = KNN_moObjectsID[i];
				if (needCalculatePruned.containsKey(temp_kNNKey))
					needCalculateLRDPruned.put(((Record) needCalculatePruned.get(temp_kNNKey).getObj()).getRId(),
							needCalculatePruned.get(temp_kNNKey));
			}
		}
		// context.getCounter(Counters.LRDPrunedPoints).increment(countPruned);
		return true;
	}

}
