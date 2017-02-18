package microcluster.topnlof;

import java.io.IOException;
import java.util.ArrayList;
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

import metricspace.*;
import microcluster.metricobject.MetricObject;
import microcluster.metricobject.MicroCluster;
import microcluster.preprocessing.AssignPointsToMC;
import microcluster.preprocessing.BuildCFTree;
import util.SQConfig;

/**
 * The main function of the project
 * 
 * @author yizhouyan
 *
 */
public class TopNLOFDetection {
	private HashMap<Integer, MicroCluster> mcList;
	private HashSet<Integer> prunedMCsIndex;
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

	public TopNLOFDetection() {
		this.mcList = new HashMap<Integer, MicroCluster>();
	}

	public void preprocessing() {
		// build a CF-Tree
		BuildCFTree cfTree = new BuildCFTree();
		System.out.println(".............Start training CF-Tree..............");
		cfTree.startTrainingCFTree();
		// cfTree.printGeneratedTree();
		System.out.println(".............Start Generating Micro Clusters..............");
		this.mcList = cfTree.generateMicroClusters(metricSpace);
		System.out.println("Size of Micro Clusters: " + mcList.size());
		
		AssignPointsToMC assignPoints = new AssignPointsToMC(this.mcList);
		System.out.println(".............Start Assigning Points to Micro Clusters..............");
		assignPoints.assignPointsToMC();
//		for(Map.Entry<Integer, MicroCluster> mc: mcList.entrySet()){
//			System.out.println("Cluster " + mc.getKey() + "," + mc.getValue().getNumberPoints());
//		}
	}

	public void computeLOFBound() throws IOException, InterruptedException {
		ComputeLOFBound computeBound = new ComputeLOFBound(this.mcList);
		computeBound.ComputeKDistanceBoundForMCs(metricSpace, metric);
		System.out.println("..............End Computing Kdistance Bounds..............");
		computeBound.ComputeFinalLOFBound(metricSpace);
		System.out.println("..............End Computing LOF Bounds..............");
		this.prunedMCsIndex = new HashSet<Integer>();
		computeBound.rankTopNOutliers(this.prunedMCsIndex);
		System.out.println("..............End Removing Micro Clusters..............");
	}

	public void printMCList() {
		System.out.println("Current MC List Size: " + mcList.size());
		for (MicroCluster mc : mcList.values()) {
			mc.printMC();
		}
	}

	public void ComputeTopNLOFValues() throws IOException, InterruptedException {
		ArrayList<MetricObject> pointList = new ArrayList<MetricObject>();
		for (Map.Entry<Integer, MicroCluster> entry : mcList.entrySet()) {
			if (this.prunedMCsIndex.contains(entry.getKey())) {
				// points can be pruned
				for (MetricObject mo : entry.getValue().getPointList()) {
					mo.setPrune(true);
					pointList.add(mo);
				}
			} else {
				for (MetricObject mo : entry.getValue().getPointList()) {
					pointList.add(mo);
				}
			}
		}
		ComputeTopNLOF topN = new ComputeTopNLOF(pointList);
		HashMap<Long, MetricObject> kdistanceList = topN.ComputeKNNForUnPrunedPoints(metricSpace, metric);
		HashMap<Long, MetricObject> lrdList = topN.ComputeLRDForUnPrunedPoints(metricSpace, kdistanceList);
		topN.ComputeLOFAndTopN(lrdList, metricSpace);
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

		Option paraRadius = new Option("c", true, "Cluster Radius");
		paraRadius.setRequired(true);
		options.addOption(paraRadius);

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
		SQConfig.domainRange = Float.parseFloat(cmd.getOptionValue("r"));
		SQConfig.clusterRadiusRate = Double.parseDouble(cmd.getOptionValue("c"));
		SQConfig.distThreshold = SQConfig.domainRange * SQConfig.clusterRadiusRate;
		SQConfig.outputFile = cmd.getOptionValue("o");
		
		System.out.println("K = " + SQConfig.K);
		System.out.println("Top-N = " + SQConfig.TOPN);
		System.out.println("Input File Path =  " + SQConfig.dataset);
		System.out.println("Output File Path = " + SQConfig.outputFile);
		System.out.println("Dim =  " + SQConfig.dims);
		System.out.println("Domain Range =  " + SQConfig.domainRange);
		System.out.println("Cluster Radius Rate =  " + SQConfig.clusterRadiusRate);

		long begin = System.currentTimeMillis();
		TopNLOFDetection topn = new TopNLOFDetection();
		try {
			topn.readMetricAndMetricSpace();
			topn.preprocessing();
			System.out.println("..............End Preprocessing..............");
			System.out.println("..............Preprocessing takes " + (System.currentTimeMillis()-begin)/1000 + " seconds.............");
			// topn.printMCList();
			topn.computeLOFBound();
			topn.ComputeTopNLOFValues();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		long end = System.currentTimeMillis();
		long second = (end - begin)/1000;
		System.err.println("Total computation time " + " takes " + second + " seconds");
	}
}
