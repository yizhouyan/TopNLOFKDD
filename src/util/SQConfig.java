package util;

import net.cftree.mc.CFTree;

public class SQConfig {
	/** ============================ basic threshold ================ */
	/** number of K */
	public static int K = 6;
	/** N for Top-N */
	public static int TOPN = 10;
	/** input file path */
	public static String dataset = "small_mass.csv";
	// public static final String dataset = "InputFile";
	/** number of dimensions */
	public static int dims = 2;
	/** domain range */
	public static float[] domains = { 0.0f, 10000.0f };
	/** domain Space */
	public static float[] domainSpace = { 0.0f, 10000.0f, 0.0f, 10000.0f };
	/** metric space */
	public static final String strMetricSpace = "vector";
	/** metric */
	public static final String strMetric = "L2Metric";

	/** current threshold of Top-N LOF */
	public static float thresholdLof = 0.0f;
	/** output file path */
	public static String outputFile = "result.csv";
	public static int countPointBasedPruned = 0;
	/** Number of threads */
	public static int numThreads = 50;

	/**
	 * ============================ parameters for Multi-dimensional Pruning
	 * ==================
	 */

	/**
	 * index of independent dims, have to change when running multidimensional
	 * data
	 */
	public static final int[] indexOfIndependentDims = { 0, 1 };
	/**
	 * number of independent dims, have to change when running multidimensional
	 * data
	 */
	public static final int num_independentDim = 2;
	/**
	 * number of correlated dims, have to change when running multidimensional
	 * data
	 */
	public static final int num_correlatedDim = 0;
	/**
	 * index of correlated dims, have to change when running multidimensional
	 * data
	 */
	public static final int[] indexOfCorrelatedDims = {};

	/**
	 * ============================ parameters for CF-Tree ==================
	 */
	/**
	 * number of max entries in each node, set to be large since we don't limit
	 * the max entry size
	 */
	public static final int maxNodeEntries = 400000000;
	/** initial distance threshold (= sqrt(radius)) */
	public static double clusterRadiusRate = 0.1;

	/** domain range of the given dataset [0,domainRange] */
	public static double domainRange = 10000;

	public static double distThreshold = domainRange * clusterRadiusRate;
	/** distance function of CF-Tree */
	public static final int distFunction = CFTree.D0_DIST;

	/** ============================= seperators ================ */
	/** seperator for items of every record in the index */
	public static final String sepStrForRecord = ",";
	public static final String sepStrForKeyValue = "\t";
	public static final String sepStrForIDDist = "|";
	public static final String sepSplitForIDDist = "\\|";

}
