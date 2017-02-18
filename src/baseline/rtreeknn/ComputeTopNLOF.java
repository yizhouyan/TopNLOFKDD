package baseline.rtreeknn;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import gnu.trove.procedure.TIntProcedure;
import util.PriorityQueue;
import metricspace.IMetric;
import metricspace.IMetricSpace;
import baseline.rtreeknn.MetricObject;
import metricspace.MetricSpaceUtility;

import metricspace.Record;
import net.sf.jsi.Point;
import net.sf.jsi.Rectangle;
import net.sf.jsi.SpatialIndex;
import net.sf.jsi.rtree.RTree;
import util.SQConfig;

public class ComputeTopNLOF {
	private HashMap<Integer, MetricObject> pointList;
	private IMetricSpace metricSpace = null;
	private IMetric metric = null;
	private SpatialIndex si;

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
		pointList = new HashMap<Integer, MetricObject>();
		BufferedReader in;
		int count = 0;
		try {
			in = new BufferedReader(new FileReader(SQConfig.dataset));
			String line = null;
			while ((line = in.readLine()) != null) {
				Object currentPoint = metricSpace.readObject(line, SQConfig.dims);
				MetricObject mo = new MetricObject(currentPoint);
				this.pointList.put((int) ((Record) mo.getObj()).getRId(), mo);
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

	public void BuildRTreeForPoints() {
		si = new RTree();
		si.init(null);
		for (Map.Entry<Integer, MetricObject> entry : pointList.entrySet()) {
			float[] coor = ((Record) entry.getValue().getObj()).getValue();
			Rectangle tempR = new Rectangle(coor[0], coor[1], coor[0], coor[1]);
			si.add(tempR, entry.getKey());
		}
	}

	public HashMap<Integer, Float> ComputeKNN() throws IOException, InterruptedException {
		System.out.println("Start Computing KNN....");
		int count = 0;
		for (Map.Entry<Integer, MetricObject> entry : pointList.entrySet()) {
			final MetricObject tempMO = entry.getValue();
			float[] coord = ((Record) tempMO.getObj()).getValue();
			si.nearestN(new Point(coord[0], coord[1]), new TIntProcedure() {
				public boolean execute(int i) {
					float currentDist = 0.0f;
					try {
						currentDist = metric.dist(pointList.get(i).getObj(), tempMO.getObj());
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					tempMO.addKNN(i, currentDist);
					return true;
				}
			}, SQConfig.K, Float.POSITIVE_INFINITY);
			count++;
			if (count % 10000 == 0) {
				System.out.println(count + " Points KNN Search Finished !!!!!!");
			}
		}

		System.out.println("End KNN Search....");

		HashMap<Integer, Float> kdistanceList = new HashMap<Integer, Float>();
		System.out.println("...............Start Computing Kdistance ...............");
		// compute KNN for points that cannot be pruned
		for (Map.Entry<Integer, MetricObject> entry : pointList.entrySet()) {
			if (entry.getValue().getkNNINfo().size() != SQConfig.K)
				System.out.println("No KNN Found!!!!");
			else {
				kdistanceList.put(entry.getKey(), entry.getValue().computeKdist());
			}
		}
		System.out.println("...............End Computing Kdistance ...............");
		return kdistanceList;
	}

	public HashMap<Integer, Float> ComputeLRD(HashMap<Integer, Float> kdistanceList)
			throws IOException, InterruptedException {
		HashMap<Integer, Float> lrdList = new HashMap<Integer, Float>();
		System.out.println("...............Start Computing LRD ............");
		for (MetricObject o_S : pointList.values()) {
			CalLRDForSingleObject(o_S, kdistanceList, lrdList);
		}
		return lrdList;
	}

	private void CalLRDForSingleObject(MetricObject o_S, HashMap<Integer, Float> kdistanceList,
			HashMap<Integer, Float> lrdList) throws IOException, InterruptedException {
		float lrd_core = 0.0f;

		for (Map.Entry<Integer, Float> tempKnn : o_S.getkNNINfo().entrySet()) {
			float temp_dist = tempKnn.getValue();
			float temp_reach_dist = Math.max(temp_dist, kdistanceList.get(tempKnn.getKey()));
			lrd_core += temp_reach_dist;
		}
		lrd_core = 1.0f / (lrd_core / SQConfig.K * 1.0f);
		// System.out.println(lrd_core);
		o_S.setLrd(lrd_core);
		lrdList.put((int) metricSpace.getID(o_S.getObj()), lrd_core);
	}

	public void ComputeLOFAndTopN(HashMap<Integer, Float> lrdList) throws IOException, InterruptedException {
		PriorityQueue TopNPQ = new PriorityQueue(PriorityQueue.SORT_ORDER_ASCENDING);
		System.out.println("...............Start Computing LOF ............");

		for (MetricObject o_S : pointList.values()) {
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
				// System.out.println(TopNPQ.getValue() + "," +
				// TopNPQ.getPriority());
				out.newLine();
				TopNPQ.pop();
			}
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void CalLOFForSingleObject(MetricObject o_S, HashMap<Integer, Float> lrdList)
			throws IOException, InterruptedException {
		float lof_core = 0.0f;
		if (o_S.getLrd() == 0)
			lof_core = 0;
		else {
			for (Map.Entry<Integer, Float> tempKnn : o_S.getkNNINfo().entrySet()) {
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
			topn.BuildRTreeForPoints();
			HashMap<Integer, Float> kdistList = topn.ComputeKNN();
			HashMap<Integer, Float> lrdList = topn.ComputeLRD(kdistList);
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
