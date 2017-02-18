package cellpruning.lof.pruning.MultiThread;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import cellpruning.lof.pruning.partitionTreeNode;
import cellpruning.lof.pruning.firstknn.prQuadTree.prQuadLeaf;
import metricspace.*;
import util.PriorityQueue;
import util.SQConfig;

public class ComputeTopNLOFWithPruning {
	private ArrayList<MetricObject> pointList = new ArrayList<MetricObject>();
	private IMetricSpace metricSpace = null;
	private IMetric metric = null;
	private float thresholdLof = 0.0f;
	private PriorityQueue topnLOF = new PriorityQueue(PriorityQueue.SORT_ORDER_ASCENDING);

	public ComputeTopNLOFWithPruning() throws IOException {
		this.readMetricAndMetricSpace();
		this.thresholdLof = SQConfig.thresholdLof;
		this.pointList = new ArrayList<MetricObject>();
		this.ReadInputFile();
	}

	public void ReadInputFile() {
		System.out.println("Start Reading Points");
		List<Thread> list = new ArrayList<Thread>();

		for (int i = 0; i < SQConfig.numThreads; i++) {
			Thread thread = new Thread(new ReadFile(metricSpace, pointList), i + "");
			thread.start();
			list.add(thread);
		}
		try {
			for (Thread thread : list)
				thread.join();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		// pointList = new ArrayList<MetricObject>();
		// BufferedReader in;
		//// long count = 0;
		// try {
		// in = new BufferedReader(new FileReader(SQConfig.dataset));
		// String line = null;
		// while ((line = in.readLine()) != null) {
		// Object currentPoint = metricSpace.readObject(line, SQConfig.dims);
		// MetricObject mo = new MetricObject(currentPoint);
		// this.pointList.add(mo);
		//// count++;
		//// if (count % 10000000 == 0)
		//// System.out.println(pointList.size() + " Points Added!");
		// }
		// in.close();
		// } catch (IOException e) {
		// // TODO Auto-generated catch block
		// e.printStackTrace();
		// }
		System.out.println(pointList.size() + " Points Added!");
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

	public static float[] maxOfTwoFloatArray(float[] x, float[] y) {
		float[] newArray = new float[x.length];
		for (int i = 0; i < x.length; i++) {
			newArray[i] = Math.max(x[i], y[i]);
		}
		return newArray;
	}

	public void ComputeTopNLOFAndPrune() throws IOException, InterruptedException {
		System.out.println("Start Computation!");
		long beginPreprocessing = System.currentTimeMillis();
		ArrayList<LargeCellStore> leaveNodes = new ArrayList<LargeCellStore>();
		ClosestPair cpObj = new ClosestPair(SQConfig.domainSpace);
		float[] independentCoordinates = new float[SQConfig.num_independentDim * 2];
		for (int i = 0; i < SQConfig.num_independentDim; i++) {
			independentCoordinates[i * 2] = SQConfig.domainSpace[SQConfig.indexOfIndependentDims[i] * 2];
			independentCoordinates[i * 2 + 1] = SQConfig.domainSpace[SQConfig.indexOfIndependentDims[i] * 2 + 1];
		}

		partitionTreeNode ptn = cpObj.divideAndConquer(pointList, independentCoordinates,
				SQConfig.indexOfIndependentDims, leaveNodes, SQConfig.K, metric, metricSpace);

		// save points that can prune
		HashMap<Long, MetricObject> CanPrunePoints = new HashMap<Long, MetricObject>();

		// save points that can find exact knns
		HashMap<Long, MetricObject> TrueKnnPoints = new HashMap<Long, MetricObject>();
		HashMap<Long, MetricObject> lrdHM = new HashMap<Long, MetricObject>();

		System.out.println("..........End Preprocessing............ ");
		System.out.println("..........Preprocessing takes " + (System.currentTimeMillis() - beginPreprocessing) / 1000
				+ " seconds............ ");
		System.out.println("Leave Node Size: " + leaveNodes.size());

		int leaveNodesCountForEachThread = (int) Math.ceil(leaveNodes.size() / SQConfig.thresholdLof);
		List<Thread> list = new ArrayList<Thread>();

		for (int i = 0; i < SQConfig.numThreads; i++) {
			Thread thread = new Thread(
					new InnerKNNSearch(CanPrunePoints, TrueKnnPoints, lrdHM, i * leaveNodesCountForEachThread,
							Math.min((i + 1) * leaveNodesCountForEachThread, leaveNodes.size()), ptn, leaveNodes,
							topnLOF, metricSpace),
					i + "");
			thread.start();
			list.add(thread);
		}
		try {
			for (Thread thread : list)
				thread.join();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		System.out.println("................ End Inner KNN Search ....................");
		System.out.println("Can Prune " + CanPrunePoints.size() + " Points!");

		list.clear();
		for (int i = 0; i < SQConfig.numThreads; i++) {
			Thread thread = new Thread(new OuterKNNSearch(TrueKnnPoints, leaveNodes, i * leaveNodesCountForEachThread,
					Math.min((i + 1) * leaveNodesCountForEachThread, leaveNodes.size()), ptn, independentCoordinates),
					i + "");
			thread.start();
			list.add(thread);
		}
		try {
			for (Thread thread : list)
				thread.join();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		
		System.out.println("................ End Outer KNN Search ....................");

		// save those pruned points but need to recompute KNNs
		HashMap<Long, MetricObject> needCalculatePruned = new HashMap<Long, MetricObject>();
		// save those cannot be pruned only by LRD value...
		HashMap<Long, MetricObject> needCalLOF = new HashMap<Long, MetricObject>();
		// need more knn information, maybe knn is pruned...
		HashMap<Long, MetricObject> needRecalLRD = new HashMap<Long, MetricObject>();

		// calculate LRD for points that can not be pruned
		for (MetricObject mo : pointList) {
			if (mo.isCanPrune())
				continue;
			if (mo.getType() == 'T')
				ComputeLRD.CalLRDForSingleObject(mo, TrueKnnPoints, CanPrunePoints, needCalculatePruned, lrdHM,
						needCalLOF, needRecalLRD, thresholdLof, leaveNodes, SQConfig.K);
		}

		// for those pruned by cell-based pruning, find kNNs for these
		// points
		for (MetricObject mo : needCalculatePruned.values()) {
			int tempIndex = mo.getIndexOfCPCellInList();
			prQuadLeaf curLeaf = leaveNodes.get(tempIndex).findLeafWithSmallCellIndex(
					leaveNodes.get(tempIndex).getRootForPRTree(), mo.getIndexForSmallCell(),
					SQConfig.indexOfIndependentDims);
			if (!mo.isInsideKNNfind()) {
				leaveNodes.get(tempIndex).findKnnsForOnePointInsideBucket(TrueKnnPoints, mo, curLeaf,
						leaveNodes.get(tempIndex), SQConfig.K, SQConfig.indexOfIndependentDims);
			}
			leaveNodes.get(tempIndex).findKnnsForOnePointOutsideBucket(TrueKnnPoints, mo, leaveNodes,
					leaveNodes.get(tempIndex), ptn, independentCoordinates, SQConfig.K, SQConfig.dims,
					SQConfig.indexOfIndependentDims);

		}

		// knn's knn is pruned...
		HashMap<Long, MetricObject> needCalculateLRDPruned = new HashMap<Long, MetricObject>();
		// calculate LRD for some points again....
		for (MetricObject mo : needRecalLRD.values()) {
			ComputeLRD.ReCalLRDForSpecial(mo, TrueKnnPoints, needCalculatePruned, lrdHM, needCalLOF,
					needCalculateLRDPruned, thresholdLof, SQConfig.K);
		}

		// (needs implementation) calculate LRD for points that needs
		// Calculate LRD (deal with needCalculateLRDPruned)
		for (MetricObject mo : needCalculateLRDPruned.values()) {
			float lrd_core = 0.0f;
			boolean canCalLRD = true;
			long[] KNN_moObjectsID = mo.getPointPQ().getValueSet();
			float[] moDistToKNN = mo.getPointPQ().getPrioritySet();

			for (int i = 0; i < KNN_moObjectsID.length; i++) {
				long knn_mo = KNN_moObjectsID[i];
				// first point out which large cell it is in
				// (tempIndexX, tempIndexY)
				float kdistknn = 0.0f;
				if (TrueKnnPoints.containsKey(knn_mo))
					kdistknn = TrueKnnPoints.get(knn_mo).getKdist();
				else if (CanPrunePoints.containsKey(knn_mo) && (!needCalculatePruned.containsKey(knn_mo))) {
					MetricObject newKnnFind = CanPrunePoints.get(knn_mo);
					int tempIndex = newKnnFind.getIndexOfCPCellInList();
					prQuadLeaf curLeaf = leaveNodes.get(tempIndex).findLeafWithSmallCellIndex(
							leaveNodes.get(tempIndex).getRootForPRTree(), newKnnFind.getIndexForSmallCell(),
							SQConfig.indexOfIndependentDims);
					if (!newKnnFind.isInsideKNNfind()) {
						leaveNodes.get(tempIndex).findKnnsForOnePointInsideBucket(TrueKnnPoints, newKnnFind, curLeaf,
								leaveNodes.get(tempIndex), SQConfig.K, SQConfig.indexOfIndependentDims);
					}
					leaveNodes.get(tempIndex).findKnnsForOnePointOutsideBucket(TrueKnnPoints, newKnnFind, leaveNodes,
							leaveNodes.get(tempIndex), ptn, independentCoordinates, SQConfig.K, SQConfig.dims,
							SQConfig.indexOfIndependentDims);
					if (TrueKnnPoints.containsKey(knn_mo))
						kdistknn = TrueKnnPoints.get(knn_mo).getKdist();
					else {
						canCalLRD = false;
						break;
					}
				} else {
					canCalLRD = false;
					break;
				}
				float temp_reach_dist = Math.max(moDistToKNN[i], kdistknn);
				lrd_core += temp_reach_dist;
				// System.out.println("Found KNNs for pruning point: " +
				// mo.getKdist());
			}
			if (canCalLRD) {
				lrd_core = 1.0f / (lrd_core / SQConfig.K * 1.0f);
				mo.setLrdValue(lrd_core);
				mo.setType('L');
				lrdHM.put(((Record) mo.getObj()).getRId(), mo);
			}
			// else{
			// System.out.println("Cannot compute LRD");
			// }
		} // end calculate LRD for pruned points
			////////////////////////////////////////////////////////////////////////////////////////
			// calculate LOF for points that can calculate
		for (MetricObject mo : pointList) {
			if (!mo.isCanPrune() && mo.getType() != 'O') {
				ComputeLOF.CalLOFFinal(mo, lrdHM, SQConfig.K, thresholdLof);
				if (mo.getType() == 'O' && mo.getLofValue() > thresholdLof) {
					float tempLofValue = mo.getLofValue();
					if (topnLOF.size() < SQConfig.TOPN) {
						topnLOF.insert(metricSpace.getID(mo.getObj()), tempLofValue);

					} else if (tempLofValue > topnLOF.getPriority()) {
						topnLOF.pop();
						topnLOF.insert(metricSpace.getID(mo.getObj()), tempLofValue);
						if (thresholdLof < topnLOF.getPriority())
							thresholdLof = topnLOF.getPriority();
					}
				}
			}
		}
		if (topnLOF.size() == SQConfig.TOPN && thresholdLof < topnLOF.getPriority())
			thresholdLof = topnLOF.getPriority();

		System.out.println("Point-based Pruning: " + SQConfig.countPointBasedPruned);
		System.err.println("computation finished");
		// write top-n lof to file
		BufferedWriter out;
		try {
			out = new BufferedWriter(new FileWriter(new File(SQConfig.outputFile)));
			while (topnLOF.size() > 0) {
				// System.out.println(topnLOF.getValue() + "," +
				// topnLOF.getPriority());
				out.write(topnLOF.getValue() + "," + topnLOF.getPriority());
				out.newLine();
				topnLOF.pop();
			}

			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		for (MetricObject mo : this.pointList) {
			if (mo.isCanPrune() || mo.getType() == 'O')
				continue;
			else {
				if (needRecalLRD.containsKey(((Record) mo.getObj()).getRId())) {
					System.out.println("Contains in need re cal lof");
				} else
					System.out.println(mo.getType());
			}
		}
	} // end
		// function

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

		Option outputFilePath = new Option("o", true, "Output File Path");
		outputFilePath.setRequired(true);
		options.addOption(outputFilePath);
		
		Option paraDim = new Option("d", true, "Dimensions");
		paraDim.setRequired(true);
		options.addOption(paraDim);

		Option paraDomain = new Option("r", true, "Domain Range");
		paraDomain.setRequired(true);
		options.addOption(paraDomain);
		
		Option paraThreadsNum = new Option("t", true, "Thread Number");
		paraThreadsNum.setRequired(true);
		options.addOption(paraThreadsNum);

		Option paraThreshold = new Option("h", true, "Start Threshold");
		paraThreshold.setRequired(true);
		options.addOption(paraThreshold);
		
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
		SQConfig.domains = new float[2];
		SQConfig.domains[0] = 0.0f;
		SQConfig.domains[1] = tempDomainRange;
		SQConfig.domainSpace = new float[2 * SQConfig.dims];
		for (int i = 0; i < SQConfig.dims; i++) {
			SQConfig.domainSpace[i * 2] = 0.0f;
			SQConfig.domainSpace[i * 2 + 1] = tempDomainRange;
		}
		SQConfig.numThreads = Integer.parseInt(cmd.getOptionValue("t"));
		SQConfig.thresholdLof = Float.parseFloat(cmd.getOptionValue("h"));
		SQConfig.outputFile = cmd.getOptionValue("o");
		
		System.out.println("K = " + SQConfig.K);
		System.out.println("Top-N = " + SQConfig.TOPN);
		System.out.println("Input File Path =  " + SQConfig.dataset);
		System.out.println("Output File Path = " + SQConfig.outputFile);
		System.out.println("Dim =  " + SQConfig.dims);
		System.out.println("Domain Range =  " + tempDomainRange);
		System.out.println("Threads Num =  " + SQConfig.numThreads);
		System.out.println("Start threshold = " + SQConfig.thresholdLof);
		
		long begin = System.currentTimeMillis();
		try {
			ComputeTopNLOFWithPruning computeTopNLOFWithPruning = new ComputeTopNLOFWithPruning();
			System.out.println("Reading Dataset takes " + ((System.currentTimeMillis() - begin) / 1000) + " seconds");
			computeTopNLOFWithPruning.ComputeTopNLOFAndPrune();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		long end = System.currentTimeMillis();
		long second = (end - begin) / 1000;
		System.err.println("Computing Top-N LOF with Pruning takes " + second + " seconds");
	}
}
