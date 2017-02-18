package microcluster.preprocessing;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import microcluster.metricobject.MetricObject;

import microcluster.metricobject.MicroCluster;
import gnu.trove.procedure.TIntProcedure;
import metricspace.*;
import net.sf.jsi.Point;
import net.sf.jsi.Rectangle;
import net.sf.jsi.SpatialIndex;
import net.sf.jsi.rtree.RTree;
import util.SQConfig;

public class AssignPointsToMC {
	private HashMap<Integer, MicroCluster> mcList;
	private IMetricSpace metricSpace = null;
	private IMetric metric = null;

	public SpatialIndex buildRTreeWithMCList() {
		System.out.println("...................Build a R-Tree for Micro Clusters...........");
		// Create and initialize an rtree
		SpatialIndex si = new RTree();
		si.init(null);
		for (Map.Entry<Integer, MicroCluster> mc : mcList.entrySet()) {
			float[] coor = ((Record) mc.getValue().getClusterCoor()).getValue();
			Rectangle tempR = new Rectangle(coor[0], coor[1], coor[0], coor[1]);
			si.add(tempR, mc.getKey());
		}
		return si;
	}

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

	public void findClosestMCPivotAndAdd(final Object obj, SpatialIndex si) {
//		int closestMCId = -1;
//		float minDist = Float.MAX_VALUE;
		float [] coord = ((Record)obj).getValue();
		si.nearest(new Point(coord[0], coord[1]), new TIntProcedure() { 
			public boolean execute(int i) {
				float currentDist = 0.0f;
				try {
					currentDist = metric.dist(mcList.get(i).getClusterCoor(), obj);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				this.addPointToClosestMC(i, currentDist, obj);
				return false; 
			}

			private void addPointToClosestMC(int i, float currentDist, Object obj) {
				// TODO Auto-generated method stub
				MetricObject mo = new MetricObject(obj);
				mcList.get(i).addPointToCurrentCluster(mo, currentDist);
			}
		}, Float.POSITIVE_INFINITY);
	}

	public void assignPointsToMC() {
		SpatialIndex si = buildRTreeWithMCList();
		// Read one instace at a time from the dataset
		// Dataset format: each line contain a set of value id, v1, v2, v3...
		// separated by comma
		BufferedReader in;
		try {
			in = new BufferedReader(new FileReader(SQConfig.dataset));
			String line = null;
			while ((line = in.readLine()) != null) {
				Object currentPoint = metricSpace.readObject(line, SQConfig.dims);
				this.findClosestMCPivotAndAdd(currentPoint, si);
			}
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
}
