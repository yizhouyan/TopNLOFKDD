package cellpruning.lof.pruning;

import java.io.IOException;
import java.util.HashMap;

import metricspace.MetricObject;

public class ComputeLOF {
	public static void CalLOFForSingleObject(MetricObject o_S, HashMap<Long, MetricObject> lrdHm, int K,
			float thresholdLof) throws IOException, InterruptedException {
		float lof_core = 0.0f;
		if (o_S.getLrdValue() == 0)
			lof_core = 0;
		else {
			long[] KNN_moObjectsID = o_S.getPointPQ().getValueSet();
			for (int i = 0; i < KNN_moObjectsID.length; i++) {
				long temp_kNNKey = KNN_moObjectsID[i];
				if (!lrdHm.containsKey(temp_kNNKey)) {
//					System.out.println(temp_kNNKey);
					return;
				}
				float temp_lrd = lrdHm.get(temp_kNNKey).getLrdValue();
				if (temp_lrd == 0 || o_S.getLrdValue() == 0)
					continue;
				else
					lof_core += temp_lrd / o_S.getLrdValue() * 1.0f;
			}
			lof_core = lof_core / K * 1.0f;
		}
		if (Float.isNaN(lof_core) || Float.isInfinite(lof_core))
			lof_core = 0;
		o_S.setLofValue(lof_core);
		o_S.setType('O'); // calculated LOF
		if (lof_core <= thresholdLof) {
			o_S.setCanPrune(true);
		}
	}
	
	public static void CalLOFFinal(MetricObject o_S, HashMap<Long, MetricObject> lrdHm, int K,
			float thresholdLof) throws IOException, InterruptedException {
		float lof_core = 0.0f;
		if (o_S.getLrdValue() == 0)
			lof_core = 0;
		else {
			long[] KNN_moObjectsID = o_S.getPointPQ().getValueSet();
			for (int i = 0; i < KNN_moObjectsID.length; i++) {
				long temp_kNNKey = KNN_moObjectsID[i];
				if (!lrdHm.containsKey(temp_kNNKey)) {
					System.out.println(temp_kNNKey);
					return;
				}
				float temp_lrd = lrdHm.get(temp_kNNKey).getLrdValue();
				if (temp_lrd == 0 || o_S.getLrdValue() == 0)
					continue;
				else
					lof_core += temp_lrd / o_S.getLrdValue() * 1.0f;
			}
			lof_core = lof_core / K * 1.0f;
		}
		if (Float.isNaN(lof_core) || Float.isInfinite(lof_core))
			lof_core = 0;
		o_S.setLofValue(lof_core);
		o_S.setType('O'); // calculated LOF
		if (lof_core <= thresholdLof) {
			o_S.setCanPrune(true);
		}
	}
}
