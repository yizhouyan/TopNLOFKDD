package baseline.pivotknn;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import baseline.pivotknn.MetricObject;
import util.PriorityQueue;
import metricspace.IMetric;
import metricspace.IMetricSpace;
import metricspace.MetricSpaceUtility;
import metricspace.Record;
import util.SQConfig;

public class ComputeTopNLOF {
	private ArrayList<MetricObject> pointList;
	private IMetricSpace metricSpace = null;
	private IMetric metric = null;

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

	public ComputeTopNLOF() throws IOException {
		readMetricAndMetricSpace();
		ReadInputFile();
	}

	public void ReadInputFile() {
		pointList = new ArrayList<MetricObject>();
		BufferedReader in;
		int count = 0;
		try {
			in = new BufferedReader(new FileReader(SQConfig.dataset));
			String line = null;
			while ((line = in.readLine()) != null) {
				Object currentPoint = metricSpace.readObject(line, SQConfig.dims);
				MetricObject mo = new MetricObject(currentPoint);
				this.pointList.add(mo);
				count++;
				if (count % 10000000 == 0)
					System.out.println(pointList.size() + " Points Added!");
			}
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println(pointList.size() + " Points Added!");
	}

	public HashMap<Long, Float> ComputeKNN() throws IOException, InterruptedException {
		System.out.println("Start Computing KNN....");
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
		System.out.println("End sorting....");
		HashMap<Long, Float> kdistanceList = new HashMap<Long, Float>();
		long begin = System.currentTimeMillis();
		System.out.println("...............Start Computing KNN ...............");
		// compute KNN for points that cannot be pruned
		for (int i = 0; i < pointList.size(); i++) {
			MetricObject o_S = pointList.get(i);
			o_S = findKNNForSingleObject(o_S, i, pointList, kdistanceList);
			if (i % 10000000 == 0)
				System.out.println(i + " points finished!");
		}
		long end = System.currentTimeMillis();
		long second = (end - begin) / 1000;
		System.err.println("KNN computation time " + " takes " + second + " seconds");
		return kdistanceList;
	}

	public HashMap<Long, Float> ComputeLRD(HashMap<Long, Float> kdistanceList)
			throws IOException, InterruptedException {
		HashMap<Long, Float> lrdList = new HashMap<Long, Float>();
		System.out.println("...............Start Computing LRD ............");
		for (MetricObject o_S : pointList) {
			CalLRDForSingleObject(o_S, kdistanceList, lrdList);
		}
		return lrdList;
	}

	private void CalLRDForSingleObject(MetricObject o_S, HashMap<Long, Float> kdistanceList,
			HashMap<Long, Float> lrdList) throws IOException, InterruptedException {
		float lrd_core = 0.0f;

		for (Map.Entry<Long, Float> tempKnn : o_S.getkNNINfo().entrySet()) {
			float temp_dist = tempKnn.getValue();
			float temp_reach_dist = Math.max(temp_dist, kdistanceList.get(tempKnn.getKey()));
			lrd_core += temp_reach_dist;
		}
		lrd_core = 1.0f / (lrd_core / SQConfig.K * 1.0f);
		// System.out.println(lrd_core);
		o_S.setLrd(lrd_core);
		lrdList.put(metricSpace.getID(o_S.getObj()), lrd_core);
	}

	public void ComputeLOFAndTopN(HashMap<Long, Float> lrdList) throws IOException, InterruptedException {
		PriorityQueue TopNPQ = new PriorityQueue(PriorityQueue.SORT_ORDER_ASCENDING);
		System.out.println("...............Start Computing LOF ............");
		HashMap<Long, MetricObject> idToPoint = new HashMap<Long, MetricObject>();

		for (MetricObject o_S : pointList) {
			CalLOFForSingleObject(o_S, lrdList);
			idToPoint.put(((Record) o_S.getObj()).getRId(), o_S);
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
				// System.out.println(TopNPQ.getValue() + "," +
				// TopNPQ.getPriority());
				out.newLine();
				TopNPQ.pop();
			}
			// for(MetricObject mo: pointList){
			// out.write("ID: " + ((Record)mo.getObj()).getRId());
			// out.newLine();
			// for (Map.Entry<Long, Float> e : mo.getkNNINfo().entrySet()) {
			// out.write("KNN: " + e.getKey() + "," + e.getValue());
			// out.newLine();
			// }
			// out.write("LRD: " + mo.getLrd() + "," + mo.getLof());
			// out.newLine();
			// }
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void CalLOFForSingleObject(MetricObject o_S, HashMap<Long, Float> lrdList)
			throws IOException, InterruptedException {
		float lof_core = 0.0f;
		if (o_S.getLrd() == 0)
			lof_core = 0;
		else {
			for (Map.Entry<Long, Float> tempKnn : o_S.getkNNINfo().entrySet()) {
				float temp_lrd = lrdList.get(tempKnn.getKey());
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
			HashMap<Long, Float> kdistanceList) throws IOException, InterruptedException {
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
				} else if (dist < theta) {
					pq.pop();
					pq.insert(metricSpace.getID(o_S.getObj()), dist);
					theta = pq.getPriority();
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
		kdistanceList.put(metricSpace.getID(o_R.getObj()), pq.getPriority());
		HashMap<Long, Float> kNNInfo = new HashMap<Long, Float>();
		while (pq.size() > 0) {
			kNNInfo.put(pq.getValue(), pq.getPriority());
			pq.pop();
		}
		o_R.setkNNINfo(kNNInfo);
		return o_R;
	}

	public static void main(String[] args) {
		System.out.println("Set Arguments.......");
		Options options = new Options();

		Option paraK = new Option("k", true, "K for KNN Search");
		paraK.setRequired(true);
		options.addOption(paraK);

		Option paraTopN = new Option("n", true, "Top-N Number");
		paraTopN.setRequired(true);
		options.addOption(paraTopN);

		Option inputFilePath = new Option("i", true, "Input File Path");
		inputFilePath.setRequired(true);
		options.addOption(inputFilePath);

		Option paraDim = new Option("d", true, "Dimensions");
		paraDim.setRequired(true);
		options.addOption(paraDim);

		Option paraDomain = new Option("r", true, "Domain Range");
		paraDomain.setRequired(true);
		options.addOption(paraDomain);

		Option outputFilePath = new Option("o", true, "Output File Path");
		outputFilePath.setRequired(true);
		options.addOption(outputFilePath);
		
		CommandLineParser parser = new DefaultParser();
		HelpFormatter formatter = new HelpFormatter();
		CommandLine cmd;

		try {
			cmd = parser.parse(options, args);
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			formatter.printHelp("utility-name", options);

			System.exit(1);
			return;
		}
		SQConfig.K = Integer.parseInt(cmd.getOptionValue("k"));
		SQConfig.TOPN = Integer.parseInt(cmd.getOptionValue("n"));
		SQConfig.dataset = cmd.getOptionValue("i");
		SQConfig.dims = Integer.parseInt(cmd.getOptionValue("d"));
		float tempDomainRange = Float.parseFloat(cmd.getOptionValue("r"));
		SQConfig.domainRange = tempDomainRange;
		SQConfig.outputFile = cmd.getOptionValue("o");
		
		System.out.println("K = " + SQConfig.K);
		System.out.println("Top-N = " + SQConfig.TOPN);
		System.out.println("Input File Path =  " + SQConfig.dataset);
		System.out.println("Output File Path = " + SQConfig.outputFile);
		System.out.println("Dim =  " + SQConfig.dims);
		System.out.println("Domain Range =  " + tempDomainRange);

		long begin = System.currentTimeMillis();
		try {
			ComputeTopNLOF topn = new ComputeTopNLOF();
			System.out.println("Reading Dataset takes " + ((System.currentTimeMillis() - begin) / 1000) + " seconds");
			HashMap<Long, Float> kdistList = topn.ComputeKNN();
			HashMap<Long, Float> lrdList = topn.ComputeLRD(kdistList);
			topn.ComputeLOFAndTopN(lrdList);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		long end = System.currentTimeMillis();
		long second = (end - begin) / 1000;
		System.err.println("Total computation time " + " takes " + second + " seconds");

	}
}
