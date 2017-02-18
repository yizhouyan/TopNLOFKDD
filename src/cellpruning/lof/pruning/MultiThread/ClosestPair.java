package cellpruning.lof.pruning.MultiThread;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;

import cellpruning.lof.pruning.partitionTreeInternal;
import cellpruning.lof.pruning.partitionTreeNode;
import metricspace.*;
import util.SQConfig;

public class ClosestPair {
	private static float[] partitionArea;

	public ClosestPair(float[] partitionCoor) {
		this.partitionArea = partitionCoor.clone();
	}

	public static class Pair {
		public IMetric metric = null;
		public double distance = 0.0;
		public boolean canMerge = true;

		public partitionTreeNode ptn = null;

		public Pair() {
		}

		public Pair(partitionTreeNode ptn, double distance, boolean canMerge) {
			this.ptn = ptn;
			this.distance = distance;
			this.canMerge = canMerge;
		}

		public Pair(MetricObject point1, MetricObject point2, IMetric metric) {
			this.metric = metric;
			calcDistance(point1, point2);
		}

		public void update(double distance) {
			this.distance = distance;
		}

		public void calcDistance(MetricObject point1, MetricObject point2) {
			this.distance = distance(point1, point2, metric);
		}
	}

	public static double distance(MetricObject p1, MetricObject p2, IMetric metric) {
		try {
			return metric.dist(p1.getObj(), p2.getObj());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return 0.0;
		}
	}

	public static Pair bruteForce(ArrayList<MetricObject> points, float[] coordinates, IMetric metric) {
		int numPoints = points.size();
		if (numPoints < 2)
			return null;
		Pair pair = new Pair(points.get(0), points.get(1), metric);
		if (numPoints > 2) {
			for (int i = 0; i < numPoints - 1; i++) {
				MetricObject point1 = points.get(i);
				for (int j = i + 1; j < numPoints; j++) {
					MetricObject point2 = points.get(j);
					double distance = distance(point1, point2, metric);
					if (distance < pair.distance)
						pair.update(distance);
				}
			}
		}
		return pair;
	}

	/**
	 * Sort all points based on given dimension
	 **/
	public static void sortByDim(ArrayList<MetricObject> points, final int dim) {
		Collections.sort(points, new Comparator<MetricObject>() {
			public int compare(MetricObject point1, MetricObject point2) {
				if (((Record) point1.getObj()).getValue()[dim] < ((Record) point2.getObj()).getValue()[dim])
					return -1;
				if (((Record) point1.getObj()).getValue()[dim] > ((Record) point2.getObj()).getValue()[dim])
					return 1;
				return 0;
			}
		});
	}

	public static int getLongerDimToDeal(float[] coordinates) {
		int dim = -1;
		float maxCoor = -1;
		for (int i = 0; i < coordinates.length / 2; i++) {
			float diff = coordinates[i * 2 + 1] - coordinates[i * 2];
			if (diff > maxCoor) {
				dim = i;
				maxCoor = diff;
			}
		}
		return dim;
	}

	public static boolean checkOverlapDims(float[] partitionArea, float[] checkedArea, int[] independentDims) {
		for (int i = 0; i < independentDims.length; i++) {
			if (partitionArea[2 * independentDims[i]] == checkedArea[2 * i])
				return true;
			if (partitionArea[2 * independentDims[i] + 1] == checkedArea[2 * i + 1])
				return true;
		}
		return false;
	}

	public static partitionTreeNode divideAndConquer(ArrayList<MetricObject> points, float[] coordinates,
			int[] independentDims, ArrayList<LargeCellStore> leaveNodes, int K, IMetric metric,
			IMetricSpace metricspace) {
		HashMap<Integer, ArrayList<MetricObject>> pointsSortedByDims = new HashMap<Integer, ArrayList<MetricObject>>();
		// sort by each independentDims
		for (int i = 0; i < independentDims.length; i++) {
			ArrayList<MetricObject> pointsSorted = new ArrayList<MetricObject>(points);
			sortByDim(pointsSorted, independentDims[i]);
			pointsSortedByDims.put(i, pointsSorted);
		}

		int dimToDeal = getLongerDimToDeal(coordinates);
		Pair tempRes = divideAndConquerByDim(pointsSortedByDims, independentDims, dimToDeal, coordinates, leaveNodes, K,
				metric, metricspace);

		partitionTreeNode ptn;
		if (tempRes.canMerge) {
			// System.out.println("Closest Pair: " + tempRes.distance);
			// create a new Large Cell Store
			ptn = new LargeCellStore(coordinates, points, (float) (tempRes.distance), metric, metricspace);
			((LargeCellStore) ptn).computePriorityForLargeCell(true);
			leaveNodes.add((LargeCellStore) ptn);
		} else {
			ptn = tempRes.ptn;
		}
		return ptn;
	}

	private static Pair divideAndConquerByDim(HashMap<Integer, ArrayList<MetricObject>> pointsSortedByDims,
			int[] independentDims, int dimToDeal, float[] coordinates, ArrayList<LargeCellStore> leaveNodes, int K,
			IMetric metric, IMetricSpace metricspace) {
		// System.out.println("Deal With coordinate:" + dimToDeal);
		int numPoints = pointsSortedByDims.get(dimToDeal).size();
		if (numPoints < 20 * K) {
			Pair pairRes = SecondaryDivideAndConquerByDim(pointsSortedByDims, independentDims, dimToDeal, coordinates,
					K, metric, metricspace);
			// pairRes.closestPairList.add(pairRes.distance);
			return pairRes;
			// return bruteForce(pointsSortedByX, metric);
		}

		// divide the dataset into left and right based on dimension "dimToDeal"
		int dividingIndex = numPoints >>> 1;
		ArrayList<MetricObject> leftOfCenter = new ArrayList<MetricObject>(
				pointsSortedByDims.get(dimToDeal).subList(0, dividingIndex));
		ArrayList<MetricObject> rightOfCenter = new ArrayList<MetricObject>(
				pointsSortedByDims.get(dimToDeal).subList(dividingIndex, numPoints));
		float centerDim = ((Record) rightOfCenter.get(0).getObj()).getValue()[independentDims[dimToDeal]];

		float[] leftCoordinates = coordinates.clone();
		leftCoordinates[dimToDeal * 2 + 1] = centerDim;

		float[] rightCoordinates = coordinates.clone();
		rightCoordinates[dimToDeal * 2] = centerDim;

		// System.out.println("Coordinates:" + coordinates[0] + "," +
		// coordinates[1]
		// + "," + coordinates[2] + "," + coordinates[3] + "," + centerX);

		// deal with left first (Call divide and Conquer)
		HashMap<Integer, ArrayList<MetricObject>> newPointsSortedByDims = new HashMap<Integer, ArrayList<MetricObject>>();
		// sort by each independentDims
		for (int i = 0; i < independentDims.length; i++) {
			if (i == dimToDeal) {
				newPointsSortedByDims.put(i, leftOfCenter);
			} else {
				ArrayList<MetricObject> tempList = new ArrayList<MetricObject>(leftOfCenter);
				sortByDim(tempList, independentDims[i]);
				newPointsSortedByDims.put(i, tempList);
			}
		}

		int newDimToDeal = getLongerDimToDeal(leftCoordinates);
		Pair closestPairLeft = divideAndConquerByDim(newPointsSortedByDims, independentDims, newDimToDeal,
				leftCoordinates, leaveNodes, K, metric, metricspace);

		// deal with right (Call divide and Conquer)
		newPointsSortedByDims.clear();
		// sort by each independentDims
		for (int i = 0; i < independentDims.length; i++) {
			if (i == dimToDeal) {
				newPointsSortedByDims.put(i, rightOfCenter);
			} else {
				ArrayList<MetricObject> tempList = new ArrayList<MetricObject>(rightOfCenter);
				sortByDim(tempList, independentDims[i]);
				newPointsSortedByDims.put(i, tempList);
			}
		}
		newDimToDeal = getLongerDimToDeal(rightCoordinates);
		Pair closestPairRight = divideAndConquerByDim(newPointsSortedByDims, independentDims, newDimToDeal,
				rightCoordinates, leaveNodes, K, metric, metricspace);

		int secondDim = AnotherDimNotDeal(independentDims, dimToDeal);
		return dealTwoCPs(closestPairLeft, closestPairRight, pointsSortedByDims.get(secondDim), dimToDeal, secondDim,
				independentDims, centerDim, coordinates, leftCoordinates, rightCoordinates, leftOfCenter, rightOfCenter,
				leaveNodes, metric, metricspace);
	}

	/**
	 * find another dim other than current "dimToDeal" suppose independentDims
	 * has more than 2 dims to deal
	 * 
	 * @param independentDims
	 * @param dealToDeal
	 * @return
	 */
	private static int AnotherDimNotDeal(int[] independentDims, int dimToDeal) {
		for (int i = 0; i < independentDims.length; i++) {
			if (i != dimToDeal)
				return i;
		}
		return -1;
	}

	private static Pair dealTwoCPs(Pair closestPairLeft, Pair closestPairRight,
			ArrayList<MetricObject> pointsSortedByDims, int dimToDeal, int secondDim, int[] independentDims,
			double centerDim, float[] coordinates, float[] leftCoordinates, float[] rightCoordinates,
			ArrayList<MetricObject> leftOfCenter, ArrayList<MetricObject> rightOfCenter,
			ArrayList<LargeCellStore> leaveNodes, IMetric metric, IMetricSpace metricspace) {
		// check if these can be combined
		if (closestPairLeft.canMerge && closestPairRight.canMerge) {
			// if the two closest pair is not that different
			double minCP = Math.min(closestPairLeft.distance, closestPairRight.distance);
			double maxCP = Math.max(closestPairLeft.distance, closestPairRight.distance);
			// System.out.println("Both can merge");
			if ((closestPairLeft.distance <= 1e-10 && closestPairRight.distance <= 1e-10)
					|| (minCP >= 1e-10 && maxCP / minCP <= 10)) {
				closestPairLeft.distance = minCP;
				ArrayList<MetricObject> tempList = new ArrayList<MetricObject>();
				double shortestDistance = closestPairLeft.distance;

				for (MetricObject point : pointsSortedByDims)
					if (Math.abs(centerDim
							- ((Record) point.getObj()).getValue()[independentDims[dimToDeal]]) < shortestDistance)
						tempList.add(point);
				for (int i = 0; i < tempList.size() - 1; i++) {
					MetricObject point1 = tempList.get(i);
					for (int j = i + 1; j < tempList.size(); j++) {
						MetricObject point2 = tempList.get(j);
						if ((((Record) point2.getObj()).getValue()[independentDims[secondDim]]
								- ((Record) point1.getObj())
										.getValue()[independentDims[secondDim]]) >= shortestDistance)
							break;
						double distance = distance(point1, point2, metric);
						if (distance < closestPairLeft.distance) {
							closestPairLeft.update(distance);
							shortestDistance = distance;
						}
					}
				}
				// System.out.println("Update CP: " + closestPairLeft.distance);
				return closestPairLeft;
			} // end if compare two cp
			else { // both can merge but in fact these two cannot merge because
					// of the cps
					// create two leave nodes
					// System.out.println("CP different, cannot merge");
				LargeCellStore leftLargeCell = new LargeCellStore(leftCoordinates, leftOfCenter,
						(float) (closestPairLeft.distance), metric, metricspace);
				// System.out.println("New bucket: " + leftCoordinates[0] + ","
				// + leftCoordinates[1] + "," +
				// leftCoordinates[2] + "," +
				// leftCoordinates[3] + "," + leftCoordinates[4] + "," +
				// leftCoordinates[5] + "," +
				// leftCoordinates[6] + "," +
				// leftCoordinates[7] + "," + closestPairLeft.distance);
				leftLargeCell
						.computePriorityForLargeCell(checkOverlapDims(partitionArea, leftCoordinates, independentDims));
				leaveNodes.add(leftLargeCell);
				LargeCellStore rightLargeCell = new LargeCellStore(rightCoordinates, rightOfCenter,
						(float) (closestPairRight.distance), metric, metricspace);
				// System.out.println("New bucket: " + rightCoordinates[0] + ","
				// + rightCoordinates[1] + "," +
				// rightCoordinates[2] + "," +
				// rightCoordinates[3] + "," + rightCoordinates[4] + "," +
				// rightCoordinates[5] + "," +
				// rightCoordinates[6] + "," +
				// rightCoordinates[7] + "," + closestPairRight.distance);
				rightLargeCell.computePriorityForLargeCell(
						checkOverlapDims(partitionArea, rightCoordinates, independentDims));
				leaveNodes.add(rightLargeCell);
				// then create one internal node and return this internal node
				partitionTreeInternal pti = new partitionTreeInternal(coordinates);
				leftLargeCell.setParentNode(pti);
				rightLargeCell.setParentNode(pti);
				pti.addNewChild(leftLargeCell);
				pti.addNewChild(rightLargeCell);
				return new Pair(pti, 0, false);
			}
		} // end if both can merge
		else if ((!closestPairLeft.canMerge) && closestPairRight.canMerge) {
			// System.out.println("Left can not merge, Right can merge ");
			// change the can merge one to a leave node
			LargeCellStore rightLargeCell = new LargeCellStore(rightCoordinates, rightOfCenter,
					(float) (closestPairRight.distance), metric, metricspace);
			// System.out.println("New bucket: " + rightCoordinates[0] + "," +
			// rightCoordinates[1] + "," +
			// rightCoordinates[2] + "," +
			// rightCoordinates[3] + "," + rightCoordinates[4] + "," +
			// rightCoordinates[5] + "," +
			// rightCoordinates[6] + "," +
			// rightCoordinates[7] + "," + closestPairRight.distance);
			rightLargeCell
					.computePriorityForLargeCell(checkOverlapDims(partitionArea, rightCoordinates, independentDims));
			leaveNodes.add(rightLargeCell);
			// then create one internal node and return this internal node
			partitionTreeInternal pti = new partitionTreeInternal(coordinates);
			closestPairLeft.ptn.setParentNode(pti);
			rightLargeCell.setParentNode(pti);
			pti.addNewChild(closestPairLeft.ptn);
			pti.addNewChild(rightLargeCell);
			return new Pair(pti, 0, false);
		} else if (closestPairLeft.canMerge && (!closestPairRight.canMerge)) {
			// System.out.println("Left can merge, Right can not merge ");
			LargeCellStore leftLargeCell = new LargeCellStore(leftCoordinates, leftOfCenter,
					(float) (closestPairLeft.distance), metric, metricspace);
			// System.out.println("New bucket: " + leftCoordinates[0] + "," +
			// leftCoordinates[1] + "," +
			// leftCoordinates[2] + "," +
			// leftCoordinates[3] + "," + leftCoordinates[4] + "," +
			// leftCoordinates[5] + "," +
			// leftCoordinates[6] + "," +
			// leftCoordinates[7] + "," +closestPairLeft.distance);
			leftLargeCell
					.computePriorityForLargeCell(checkOverlapDims(partitionArea, leftCoordinates, independentDims));
			leaveNodes.add(leftLargeCell);
			// then create one internal node and return this internal node
			partitionTreeInternal pti = new partitionTreeInternal(coordinates);
			leftLargeCell.setParentNode(pti);
			closestPairRight.ptn.setParentNode(pti);
			pti.addNewChild(leftLargeCell);
			pti.addNewChild(closestPairRight.ptn);
			return new Pair(pti, 0, false);
		} else { // if both cannot merge
			// System.out.println("Both can not merge");
			partitionTreeInternal pti = new partitionTreeInternal(coordinates);
			closestPairLeft.ptn.setParentNode(pti);
			closestPairRight.ptn.setParentNode(pti);
			pti.addNewChild(closestPairLeft.ptn);
			pti.addNewChild(closestPairRight.ptn);
			return new Pair(pti, 0, false);
		}
	}

	public static Pair SecondaryDivideAndConquerByDim(HashMap<Integer, ArrayList<MetricObject>> pointsSortedByDims,
			int[] independentDims, int dimToDeal, float[] coordinates, int K, IMetric metric,
			IMetricSpace metricspace) {
		// System.out.println("Deal with dimension :" + dimToDeal);
		int numPoints = pointsSortedByDims.get(dimToDeal).size();
		if (numPoints <= 3)
			return bruteForce(pointsSortedByDims.get(dimToDeal), coordinates, metric);

		// divide the dataset into left and right
		int dividingIndex = numPoints >>> 1;
		ArrayList<MetricObject> leftOfCenter = new ArrayList<MetricObject>(
				pointsSortedByDims.get(dimToDeal).subList(0, dividingIndex));
		ArrayList<MetricObject> rightOfCenter = new ArrayList<MetricObject>(
				pointsSortedByDims.get(dimToDeal).subList(dividingIndex, numPoints));

		float centerDim = ((Record) rightOfCenter.get(0).getObj()).getValue()[independentDims[dimToDeal]];
		float[] leftCoordinates = coordinates.clone();
		leftCoordinates[dimToDeal * 2 + 1] = centerDim;

		float[] rightCoordinates = coordinates.clone();
		rightCoordinates[dimToDeal * 2] = centerDim;

		// deal with left (Call divide and conquer)
		HashMap<Integer, ArrayList<MetricObject>> newPointsSortedByDims = new HashMap<Integer, ArrayList<MetricObject>>();
		// sort by each independentDims
		for (int i = 0; i < independentDims.length; i++) {
			if (i == dimToDeal) {
				newPointsSortedByDims.put(i, leftOfCenter);
			} else {
				ArrayList<MetricObject> tempList = new ArrayList<MetricObject>(leftOfCenter);
				sortByDim(tempList, independentDims[i]);
				newPointsSortedByDims.put(i, tempList);
			}
		}
		int newDimToDeal = getLongerDimToDeal(leftCoordinates);
		Pair closestPairLeft = SecondaryDivideAndConquerByDim(newPointsSortedByDims, independentDims, newDimToDeal,
				leftCoordinates, K, metric, metricspace);

		// deal with right (Call divide and Conquer)
		newPointsSortedByDims.clear();
		// sort by each independentDims
		for (int i = 0; i < independentDims.length; i++) {
			if (i == dimToDeal) {
				newPointsSortedByDims.put(i, rightOfCenter);
			} else {
				ArrayList<MetricObject> tempList = new ArrayList<MetricObject>(rightOfCenter);
				sortByDim(tempList, independentDims[i]);
				newPointsSortedByDims.put(i, tempList);
			}
		}

		newDimToDeal = getLongerDimToDeal(rightCoordinates);
		Pair closestPairRight = SecondaryDivideAndConquerByDim(newPointsSortedByDims, independentDims, newDimToDeal,
				rightCoordinates, K, metric, metricspace);

		double minCP = Math.min(closestPairLeft.distance, closestPairRight.distance);
		closestPairLeft.distance = minCP;
		ArrayList<MetricObject> tempList = new ArrayList<MetricObject>();
		double shortestDistance = closestPairLeft.distance;
		int secondDim = AnotherDimNotDeal(independentDims, dimToDeal);

		for (MetricObject point : pointsSortedByDims.get(secondDim))
			if (Math.abs(
					centerDim - ((Record) point.getObj()).getValue()[independentDims[dimToDeal]]) < shortestDistance)
				tempList.add(point);
		for (int i = 0; i < tempList.size() - 1; i++) {
			MetricObject point1 = tempList.get(i);
			for (int j = i + 1; j < tempList.size(); j++) {
				MetricObject point2 = tempList.get(j);
				if ((((Record) point2.getObj()).getValue()[independentDims[secondDim]]
						- ((Record) point1.getObj()).getValue()[independentDims[secondDim]]) >= shortestDistance)
					break;
				double distance = distance(point1, point2, metric);
				if (distance < closestPairLeft.distance) {
					closestPairLeft.update(distance);
					shortestDistance = distance;
				}
			}
		}
		return closestPairLeft;
	}

	public static int[] CalculateIndexPerDim(int number, int numIndependentDim) {
		int[] results = new int[numIndependentDim];
		for (int i = numIndependentDim - 1; i >= 0; i--) {
			int temp = (int) Math.pow(2, i);
			results[i] = (int) Math.floor(number / temp);
			number = number - temp * results[i];
		}
		return results;
	}

	public static void main(String[] args) throws IOException, InterruptedException {
		// IMetricSpace metricSpace = null;
		// IMetric metric = null;
		// try {
		// metricSpace = MetricSpaceUtility.getMetricSpace("vector");
		// metric = MetricSpaceUtility.getMetric("L2Metric");
		// metricSpace.setMetric(metric);
		// } catch (Exception e) {
		// System.out.println("Exception caught");
		// }
		// // read in objects
		// ArrayList<MetricObject> moList = new ArrayList<MetricObject>();
		// BufferedReader currentReader = null;
		// try {
		// currentReader = new BufferedReader(new FileReader("./InputFile"));
		// String line;
		// while ((line = currentReader.readLine()) != null) {
		// /** parse line */
		// MetricObject mo = parseObject(1, line, metricSpace, 2);
		// moList.add(mo);
		// }
		// } catch (FileNotFoundException e) {
		// // TODO Auto-generated catch block
		// e.printStackTrace();
		// } catch (IOException e) {
		// // TODO Auto-generated catch block
		// e.printStackTrace();
		// }
		// // System.out.println(bruteForce(moList, metric).distance);
		// float[] partitionArea = { 0, 10000, 0, 10000 };
		// float[] coordinates = { 0, 10000, 0, 10000 };
		// ArrayList<LargeCellStore> leaveNodes = new
		// ArrayList<LargeCellStore>();
		// int[] independentDims = { 0, 1 };
		// ClosestPair cpObj = new ClosestPair(partitionArea);
		// partitionTreeNode ptn = cpObj.divideAndConquer(moList, coordinates,
		// independentDims, leaveNodes, 3, metric,
		// metricSpace);
		//
		// System.out.println("Total number of large buckets: " +
		// leaveNodes.size());
		// // for (LargeCellStore lcs : leaveNodes) {
		// // System.out.println("Coordinates: " + lcs.getCoordinates()[0] + ","
		// +
		// // lcs.getCoordinates()[1] + ","
		// // + lcs.getCoordinates()[2] + "," + lcs.getCoordinates()[3] + "," +
		// // "CP:" + lcs.getCpDist());
		// // }
		// HashMap<Long, MetricObject> CanPrunePoints = new HashMap<Long,
		// MetricObject>();
		// int K = 3;
		//
		// // seperate each large bucket into a QuadTree
		// for (int i = 0; i < leaveNodes.size(); i++) {
		// if (leaveNodes.get(i).getNumOfPoints() == 0)
		// continue;
		// else if (leaveNodes.get(i).getNumOfPoints() > K * 5) {
		// leaveNodes.get(i).seperateToSmallCells(CanPrunePoints, i, 10, K,
		// independentDims, null, 2);
		// if (leaveNodes.get(i).isBreakIntoSmallCells())
		// // System.out.println("Build PRQuad Tree successful!");
		// if (!leaveNodes.get(i).isBreakIntoSmallCells() &&
		// leaveNodes.get(i).getNumOfPoints() > K * 20) {
		// // System.out.println("Build another tree........");
		// leaveNodes.get(i).seperateLargeNoPrune(K, i, independentDims, null);
		// }
		// }
		// System.out.println("Large Bucket coordinates: " +
		// leaveNodes.get(i).getCoordinates()[0] + ","
		// + leaveNodes.get(i).getCoordinates()[1] + "," +
		// leaveNodes.get(i).getCoordinates()[2] + ","
		// + leaveNodes.get(i).getCoordinates()[3]);
		// System.out.println("Closest Pair: " + leaveNodes.get(i).getCpDist() +
		// ", small cell size:"
		// + leaveNodes.get(i).getSmallCellSize());
		// System.out.println("Break into small cells? " +
		// leaveNodes.get(i).isBreakIntoSmallCells());
		//
		// if (leaveNodes.get(i).isBreakIntoSmallCells()) {
		// System.out.println("Number of small cells per dim :" +
		// leaveNodes.get(i).getNumSmallCells()[0] + ","
		// + leaveNodes.get(i).getNumSmallCells()[1]);
		// System.out.println("Leave Nodes Count: " +
		// leaveNodes.get(i).getPrLeaves().size());
		// }
		// }

	}

}