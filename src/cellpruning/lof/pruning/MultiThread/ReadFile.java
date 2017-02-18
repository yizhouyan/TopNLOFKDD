package cellpruning.lof.pruning.MultiThread;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import metricspace.IMetricSpace;
import metricspace.MetricObject;
import util.SQConfig;

public class ReadFile implements Runnable {

	private static BufferedReader in = null;
	private ArrayList<MetricObject> pointList;
	private IMetricSpace metricSpace = null;

	public ReadFile(IMetricSpace metricSpace, ArrayList<MetricObject> pointList) {
		this.metricSpace = metricSpace;
		this.pointList = pointList;
	}

	static {
		try {
			in = new BufferedReader(new FileReader(SQConfig.dataset), 1000);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
//		System.out.println("Start Reading Points");
		String line = null;
		int count = 0;
		while (true) {
			ArrayList<String> strList = new ArrayList<String>();
			synchronized (in) {
				try {
					while ((line = in.readLine()) != null) {
						if (count < 30) {
							strList.add(line);
							count++;
						} else {
							strList.add(line);
							count = 0;
							break;
						}
					}
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			ArrayList<MetricObject> objList = new ArrayList<MetricObject>();
			for (String str : strList) {
				Object currentPoint = metricSpace.readObject(str, SQConfig.dims);
				MetricObject mo = new MetricObject(currentPoint);
				objList.add(mo);
			}
			synchronized (pointList) {
				this.pointList.addAll(objList);
			}
			if (line == null)
				break;
		}
//		System.out.println(pointList.size() + " Points Added!");
	}

	public static void main(String[] args) throws IOException {
//		long begin = System.currentTimeMillis();
//		IMetricSpace metricSpace = null;
//		IMetric metric = null;
//		ArrayList<MetricObject> pointList = new ArrayList<MetricObject>();
//		try {
//			metricSpace = MetricSpaceUtility.getMetricSpace(SQConfig.strMetricSpace);
//			metric = MetricSpaceUtility.getMetric(SQConfig.strMetric);
//			metricSpace.setMetric(metric);
//		} catch (InstantiationException e) {
//			throw new IOException("InstantiationException");
//		} catch (IllegalAccessException e) {
//			e.printStackTrace();
//			throw new IOException("IllegalAccessException");
//		} catch (ClassNotFoundException e) {
//			e.printStackTrace();
//			throw new IOException("ClassNotFoundException");
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		List<Thread> list = new ArrayList<Thread>();
//
//		for (int i = 0; i < 20; i++) {
//			Thread thread = new Thread(new ReadFile(metricSpace,pointList), i + "");
//			thread.start();
//			list.add(thread);
//		}
//		try {
//			for (Thread thread : list)
//				thread.join();
//		} catch (InterruptedException e) {
//			e.printStackTrace();
//		}
//		long end = System.currentTimeMillis();
//		long second = (end - begin) / 1000;
//		System.err.println("Computing Top-N LOF with Pruning takes " + second + " seconds");
//		System.out.println(pointList.size());
	}

}
