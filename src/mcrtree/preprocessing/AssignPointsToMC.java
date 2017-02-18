package mcrtree.preprocessing;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import gnu.trove.procedure.TIntProcedure;

import mcrtree.metricobject.MetricObject;
import mcrtree.metricobject.MicroCluster;
import metricspace.*;
import util.SQConfig;

public class AssignPointsToMC {
	private HashMap<Integer, MicroCluster> mcList;
	private IMetricSpace metricSpace = null;
	private IMetric metric = null;

	public AssignPointsToMC(HashMap<Integer, MicroCluster> mcList) {
		this.mcList = mcList;
		try {
			this.readMetricAndMetricSpace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * get MetricSpace and metric from configuration
	 * 
	 * @param conf
	 * @throws IOException
	 */
	public void readMetricAndMetricSpace() throws IOException {
		try {
			metricSpace = MetricSpaceUtility.getMetricSpace(SQConfig.strMetricSpace);
			metric = MetricSpaceUtility.getMetric(SQConfig.strMetric);
			metricSpace.setMetric(metric);
		} catch (InstantiationException e) {
			throw new IOException("InstantiationException");
		} catch (IllegalAccessException e) {
			e.printStackTrace();
			throw new IOException("IllegalAccessException");
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
			throw new IOException("ClassNotFoundException");
		}
	}

	public void findClosestMCPivotAndAdd(Object obj) {
		int closestMCId = -1;
		float minDist = Float.MAX_VALUE;
		for (Map.Entry<Integer, MicroCluster> entry : this.mcList.entrySet()) {
			float currentDist = 0.0f;
			try {
				currentDist = metric.dist(entry.getValue().getClusterCoor(), obj);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			if (currentDist < minDist) {
				closestMCId = entry.getKey();
				minDist = currentDist;
			}
		}
		if (closestMCId == -1) {
			System.out.println("Cannot find closest micro cluster...");
		} else {
			this.addPointToClosestMC(closestMCId, minDist, obj);
		}
	}
	
	public void addPointToClosestMC(int closestMCId, float minDist, Object currentPoint) {
		MetricObject mo = new MetricObject(currentPoint);
		this.mcList.get(closestMCId).addPointToCurrentCluster(mo, minDist);
	}
	
	public void assignPointsToMC() {
		// Read one instace at a time from the dataset
		// Dataset format: each line contain a set of value id, v1, v2, v3...
		// separated by comma
		BufferedReader in;
		try {
			in = new BufferedReader(new FileReader(SQConfig.dataset));
			String line = null;
			while ((line = in.readLine()) != null) {
				Object currentPoint = metricSpace.readObject(line, SQConfig.dims);
				this.findClosestMCPivotAndAdd(currentPoint);
			}
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
}
