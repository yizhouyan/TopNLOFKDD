package microcluster.topnlof;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import metricspace.*;
import microcluster.metricobject.MetricObject;
import microcluster.metricobject.MicroCluster;
import util.PriorityQueue;
import util.SQConfig;

public class ComputeTopNLOF {
	private ArrayList<MetricObject> pointList;

	public ComputeTopNLOF(ArrayList<MetricObject> pointList) {
		this.pointList = pointList;
	}

	public HashMap<Long, MetricObject> ComputeKNNForUnPrunedPoints(IMetricSpace metricSpace, IMetric metric)
			throws IOException, InterruptedException {
		String centralPivotStr = "0";
		for (int j = 0; j < SQConfig.dims; j++) {
			centralPivotStr = centralPivotStr + "," + SQConfig.domainRange / 2.0;
		}

		Object cPivot = metricSpace.readObject(centralPivotStr, SQConfig.dims);
		for (int i = 0; i < pointList.size(); i++) {
			pointList.get(i).setDistToPivot(metric.dist(cPivot, pointList.get(i).getObj()));
		}
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
		HashMap<Long, Integer> MetricObjectIdToIndex = new HashMap<Long, Integer>();
		for (int i = 0; i < pointList.size(); i++) {
			pointList.get(i).setDistToPivot(metric.dist(cPivot, pointList.get(i).getObj()));
			MetricObjectIdToIndex.put(((Record) pointList.get(i).getObj()).getRId(), i);
		}

		HashMap<Long, MetricObject> kdistanceList = new HashMap<Long, MetricObject>();
		HashSet<Integer> morePointsKNN = new HashSet<Integer>();
		long begin = System.currentTimeMillis();
		System.out.println("...............Start Computing KNN for Unpruned Points............");
		// compute KNN for points that cannot be pruned
		for (int i = 0; i < pointList.size(); i++) {
			MetricObject o_S = pointList.get(i);
			if (o_S.isPrune()) {
				continue;
			} else {
				if (o_S.isTrueKNN()) {
					kdistanceList.put(metricSpace.getID(o_S.getObj()), o_S);
					continue;
				}
				o_S = findKNNForSingleObject(o_S, i, pointList, metric, metricSpace, kdistanceList, morePointsKNN,
						MetricObjectIdToIndex);
			}
		}

		System.out.println(
				"...............Start Computing KNN for Some Pruned Points............ " + morePointsKNN.size());
		// compute KNN for pruned point's KNN
		HashSet<Integer> morePointsKNNMore = new HashSet<Integer>();
		for (Integer i : morePointsKNN) {
			MetricObject o_S = pointList.get(i);
			o_S = findKNNForSingleObject(o_S, i, pointList, metric, metricSpace, kdistanceList, morePointsKNNMore,
					MetricObjectIdToIndex);
		}

		System.out.println(
				"...............Start Computing KNN for Some Pruned Points............ " + morePointsKNNMore.size());

		// compute KNN for pruned point's KNN's KNN
		morePointsKNN.clear();
		for (Integer i : morePointsKNNMore) {
			MetricObject o_S = pointList.get(i);
			if (kdistanceList.containsKey(metricSpace.getID(o_S.getObj())))
				continue;
			else
				o_S = findKNNForSingleObject(o_S, i, pointList, metric, metricSpace, kdistanceList, morePointsKNN,
						MetricObjectIdToIndex);
		}
		long end = System.currentTimeMillis();
		long second = (end - begin) / 1000;
		System.err.println("KNN computation time " + " takes " + second + " seconds");
		return kdistanceList;
	}

	public HashMap<Long, MetricObject> ComputeLRDForUnPrunedPoints(IMetricSpace metricSpace,
			HashMap<Long, MetricObject> kdistanceList) throws IOException, InterruptedException {
		HashMap<Long, MetricObject> lrdList = new HashMap<Long, MetricObject>();
		HashSet<MetricObject> prunedPoints = new HashSet<MetricObject>();
		System.out.println("...............Start Computing LRD for Unpruned Points............");
		for (MetricObject o_S : pointList) {
			if (o_S.isPrune())
				continue;
			CalLRDForSingleObject(o_S, kdistanceList, lrdList, prunedPoints, metricSpace);
		}
		System.out.println(
				"...............Start Computing LRD for Some Pruned Points............ " + prunedPoints.size());
		// compute lrd for some pruned points
		HashSet<MetricObject> moreprunedPoints = new HashSet<MetricObject>();
		for (MetricObject o_S : prunedPoints) {
			CalLRDForSingleObject(o_S, kdistanceList, lrdList, moreprunedPoints, metricSpace);
		}
		return lrdList;
	}

	private void CalLRDForSingleObject(MetricObject o_S, HashMap<Long, MetricObject> kdistanceList,
			HashMap<Long, MetricObject> lrdList, HashSet<MetricObject> prunedPoints, IMetricSpace metricSpace)
			throws IOException, InterruptedException {
		float lrd_core = 0.0f;
		for (Map.Entry<Long, Float> tempKnn : o_S.getkNNINfo().entrySet()) {
			float temp_dist = tempKnn.getValue();
			float temp_reach_dist = Math.max(temp_dist, kdistanceList.get(tempKnn.getKey()).getKdist());
			lrd_core += temp_reach_dist;
			if (kdistanceList.get(tempKnn.getKey()).isPrune()) {
				prunedPoints.add(kdistanceList.get(tempKnn.getKey()));
			}
		}
		lrd_core = 1.0f / (lrd_core / SQConfig.K * 1.0f);
		// System.out.println(lrd_core);
		o_S.setLrd(lrd_core);
		lrdList.put(metricSpace.getID(o_S.getObj()), o_S);
	}

	public void ComputeLOFAndTopN(HashMap<Long, MetricObject> lrdList, IMetricSpace metricSpace)
			throws IOException, InterruptedException {
		PriorityQueue TopNPQ = new PriorityQueue(PriorityQueue.SORT_ORDER_ASCENDING);
		System.out.println("...............Start Computing LOF for Unpruned Points............");
		for (MetricObject o_S : pointList) {
			if (o_S.isPrune())
				continue;
			CalLOFForSingleObject(o_S, lrdList);
			if (TopNPQ.size() < SQConfig.TOPN) {
				TopNPQ.insert(metricSpace.getID(o_S.getObj()), o_S.getLof());
			} else if (o_S.getLof() > TopNPQ.getPriority()) {
				// System.out.println("Update: " + TopNPQ.getPriority() + " to "
				// + o_S.getLof());
				TopNPQ.pop();
				TopNPQ.insert(metricSpace.getID(o_S.getObj()), o_S.getLof());
			}
		}

		// write top-n lof to file
		BufferedWriter out;
		try {
			out = new BufferedWriter(new FileWriter(new File(SQConfig.outputFile)));
			while (TopNPQ.size() > 0) {
				out.write(TopNPQ.getValue() + "," + TopNPQ.getPriority());
				out.newLine();
				TopNPQ.pop();
			}
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void CalLOFForSingleObject(MetricObject o_S, HashMap<Long, MetricObject> lrdList)
			throws IOException, InterruptedException {
		float lof_core = 0.0f;
		if (o_S.getLrd() == 0)
			lof_core = 0;
		else {
			for (Map.Entry<Long, Float> tempKnn : o_S.getkNNINfo().entrySet()) {
				float temp_lrd = lrdList.get(tempKnn.getKey()).getLrd();
				if (temp_lrd == 0 || o_S.getLrd() == 0) {
					lof_core = lof_core;
				} else
					lof_core += temp_lrd / o_S.getLrd() * 1.0f;
			}
			lof_core = lof_core / SQConfig.K * 1.0f;
		}
		if (Float.isNaN(lof_core) || Float.isInfinite(lof_core))
			lof_core = 0;
		o_S.setLof(lof_core);
	}

	/**
	 * find kNN using pivot based index
	 * 
	 * @return MetricObject with kdistance and knns
	 * @throws InterruptedException
	 */
	private MetricObject findKNNForSingleObject(MetricObject o_R, int currentIndex, ArrayList<MetricObject> pointList,
			IMetric metric, IMetricSpace metricSpace, HashMap<Long, MetricObject> kdistanceList,
			HashSet<Integer> morePointsKNN, HashMap<Long, Integer> MetricObjectIdToIndex)
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
					if (o_S.isPrune())
						morePointsKNN.add(inc_current);
				} else if (dist < theta) {
					pq.pop();
					pq.insert(metricSpace.getID(o_S.getObj()), dist);
					theta = pq.getPriority();
					if (o_S.isPrune())
						morePointsKNN.add(inc_current);
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
		o_R.setKdist(pq.getPriority());
		kdistanceList.put(metricSpace.getID(o_R.getObj()), o_R);
		HashMap<Long, Float> kNNInfo = new HashMap<Long, Float>();
		while (pq.size() > 0) {
			kNNInfo.put(pq.getValue(), pq.getPriority());
			if (pointList.get(MetricObjectIdToIndex.get(pq.getValue())).isPrune())
				morePointsKNN.add(MetricObjectIdToIndex.get(pq.getValue()));
			pq.pop();
		}
		o_R.setkNNINfo(kNNInfo);
		return o_R;
	}

}
