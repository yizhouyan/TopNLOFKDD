package cellpruning.lof.pruning.firstknn.prQuadTree;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Stack;

import metricspace.*;

public abstract class prQuadNode {
	private float[] coordinates;
	private int[] indexInSmallCell;
	private prQuadNode parentNode = null;
	private int[] numSmallCells;
	private float smallCellSize = 0.0f;

	public float[] getCoordinates() {
		return coordinates;
	}

	public void setCoordinates(float[] coordinates) {
		this.coordinates = coordinates;
	}

	public int[] getNumSmallCells() {
		return numSmallCells;
	}

	public void setNumSmallCells(int[] numSmallCells) {
		this.numSmallCells = numSmallCells;
	}

	public prQuadNode(float[] coordinates, int[] indexInSmallCell, prQuadNode parentNode, int[] numSmallCells,
			float smallCellSize) {
		this.coordinates = coordinates;
		this.indexInSmallCell = indexInSmallCell;
		this.parentNode = parentNode;
		this.numSmallCells = numSmallCells;
		this.smallCellSize = smallCellSize;
	}

	/**
	 * check if these two ranges overlap?
	 * 
	 * @param expectedRange
	 * @param checkedRange
	 * @return
	 */
	boolean checkRange(double[] expectedRange, double[] checkedRange) {
		// check (x1,y2) and (x2,y1)
		if (expectedRange[0] > checkedRange[1] || checkedRange[0] > expectedRange[1])
			return false;
		if (expectedRange[3] < checkedRange[2] || checkedRange[3] < expectedRange[2])
			return false;
		return true;
	}

	public void generateChilden(HashMap<Long, MetricObject> CanPrunePoints, prQuadInternal curPRNode,
			Stack<prQuadInternal> prQuadTree,
			HashMap<prQuadInternal, ArrayList<MetricObject>> mapQuadInternalWithPoints, ArrayList<prQuadLeaf> prLeaves,
			float[] largeCellCoor, int[] numSmallCells, int indexOfLeaveNodes, int[] independentDims,
			int[] correlatedDims, int K, boolean withPrune) {
		int count = 0;
		int[][] numCellsPerDim = new int[independentDims.length][2];
		for (int i = 0; i < independentDims.length; i++) {
			numCellsPerDim[i][0] = (int) Math.ceil(curPRNode.getNumSmallCells()[i] / 2);
			numCellsPerDim[i][1] = curPRNode.getNumSmallCells()[i] - numCellsPerDim[i][0];
		}

		// init small cell index (the beginning index of the second one)
		// if we have index(3,9), then the first one contains (3,6) the second
		// one contains(7,9)
		int[] midSmallIndex = new int[independentDims.length];
		for (int i = 0; i < independentDims.length; i++) {
			midSmallIndex[i] = curPRNode.getIndexInSmallCell()[2 * i] + numCellsPerDim[i][0];
		}

		int[][] newIndexInSmall = new int[independentDims.length][4];
		for (int i = 0; i < independentDims.length; i++) {
			newIndexInSmall[i][0] = curPRNode.getIndexInSmallCell()[2 * i];
			newIndexInSmall[i][1] = midSmallIndex[i] - 1;
			newIndexInSmall[i][2] = midSmallIndex[i];
			newIndexInSmall[i][3] = curPRNode.getIndexInSmallCell()[2 * i + 1];
		}

		// new Coordinate Generated
		float[] midCoor = new float[independentDims.length];
		for (int i = 0; i < independentDims.length; i++) {
			midCoor[i] = largeCellCoor[2 * i] + curPRNode.getSmallCellSize() * midSmallIndex[i];
		}
		float[][] newCoor = new float[independentDims.length][3];
		for (int i = 0; i < independentDims.length; i++) {
			newCoor[i][0] = curPRNode.getCoordinates()[2 * i];
			newCoor[i][1] = midCoor[i];
			newCoor[i][2] = curPRNode.getCoordinates()[2 * i + 1];
		}

		// Init each small bucket (dims * 2 in total)
		int numNewBucket = (int) Math.pow(2, independentDims.length);
		ArrayList[] pointsForS = new ArrayList[numNewBucket];
		for (int i = 0; i < numNewBucket; i++)
			pointsForS[i] = new ArrayList<MetricObject>();

		ArrayList<MetricObject> listOfPoints = mapQuadInternalWithPoints.get(curPRNode);
		for (MetricObject mo : listOfPoints) {
			int[] tempIndex = mo.getIndexForSmallCell();
			int result = 0;
			for (int i = 0; i < independentDims.length; i++) {
				int temp = (tempIndex[i] < midSmallIndex[i]) ? 0 : 1;
				result = result + (int) Math.pow(2, i) * temp;
			}
			pointsForS[result].add(mo);
		}

		// seperate each small bucket (previous i = x j = y)
		for (int i = 0; i < numNewBucket; i++) {
			int[] tempIndexEachDim = CalculateIndexPerDim(i, independentDims.length); // 0
																						// or
																						// 1
																						// for
																						// each
																						// dim
			int[] numCellsPerDimPerBucket = new int[independentDims.length];
			for (int j = 0; j < independentDims.length; j++) {
				numCellsPerDimPerBucket[j] = numCellsPerDim[j][tempIndexEachDim[j]];
			}
			for (int j = 0; j < independentDims.length; j++) {
				if (numCellsPerDimPerBucket[j] <= 0)
					continue;
			}
			// check size and num of points inside, if size too small or num
			// of points small, don't divide
			boolean reachMinSizeEachDim = true;
			for (int j = 0; j < independentDims.length; j++) {
				if (numCellsPerDimPerBucket[j] != 1)
					reachMinSizeEachDim = false;
			}
			float[] newConstructed = new float[independentDims.length * 2];
			for (int j = 0; j < independentDims.length; j++) {
				newConstructed[2 * j] = newCoor[j][tempIndexEachDim[j]];
				newConstructed[2 * j + 1] = newCoor[j][tempIndexEachDim[j] + 1];
			}
			boolean canPrune = false;
			int flag = 0;
			int[] newIndexInsideSmall = new int[independentDims.length * 2];
			for (int j = 0; j < independentDims.length; j++) {
				newIndexInsideSmall[2 * j] = newIndexInSmall[j][2 * tempIndexEachDim[j]];
				newIndexInsideSmall[2 * j + 1] = newIndexInSmall[j][2 * tempIndexEachDim[j] + 1];
				if (newIndexInsideSmall[2 * j] < 4 || newIndexInsideSmall[2 * j] >= numSmallCells[j] - 5) {
					flag = 1;
				}
			}
			if (reachMinSizeEachDim) { // will generate a leave node
				if (withPrune && pointsForS[i].size() > K && flag == 0) {
					// check correlatedDims, save max and min, compare to
					// smallCellSize, check if can prune
					boolean canPruneWithCorrelatedDims = true;
					for (int dd = 0; dd < correlatedDims.length; dd++) {
						float minValue = Float.MAX_VALUE;
						float maxValue = Float.MIN_VALUE;
						for (int mm = 0; mm < pointsForS[i].size(); mm++) {
							float tempValue = ((Record) ((MetricObject) pointsForS[i].get(mm)).getObj())
									.getValue()[correlatedDims[dd]];
							minValue = Math.min(minValue, tempValue);
							maxValue = Math.max(maxValue, tempValue);
						}
						if (maxValue - minValue > curPRNode.getSmallCellSize()) {
							canPruneWithCorrelatedDims = false;
							break;
						}
					} // end traverse all points all dims
						// set all points to be can prune
					if (canPruneWithCorrelatedDims) {
						canPrune = true;
						for (int mm = 0; mm < pointsForS[i].size(); mm++) {
							((MetricObject) pointsForS[i].get(mm)).setCanPrune(true);
							((MetricObject) pointsForS[i].get(mm)).setIndexOfCPCellInList(indexOfLeaveNodes);
							CanPrunePoints.put(((Record) ((MetricObject) pointsForS[i].get(mm)).getObj()).getRId(),
									((MetricObject) pointsForS[i].get(mm)));
//							System.out.println("Prune :"+((MetricObject) pointsForS[i].get(mm)).getIndexForSmallCell()[0] + "," + ((MetricObject) pointsForS[i].get(mm)).getIndexForSmallCell()[1]);
						}
					} // end set all points to be can prune
				} // end if (withPrune && pointsForS[i].size() > K) {
				
				prQuadLeaf newLeaf = new prQuadLeaf(newConstructed, newIndexInsideSmall, curPRNode,
						numCellsPerDimPerBucket, curPRNode.getSmallCellSize(), pointsForS[i], pointsForS[i].size(),
						canPrune);
				if (pointsForS[i].size() > 0) {
					curPRNode.addNewChild(newLeaf);
					// context.getCounter(Counters.SmallCellsNum).increment(1);
				}
				if ((!canPrune) && (pointsForS[i].size() > 0))
					prLeaves.add(newLeaf);
				count = count + pointsForS[i].size();
//				System.out.println("Get node: " + newLeaf.printQuadLeaf());
			} // end if(reachMinSizeEachDim cell)
			else if (pointsForS[i].size() < K + 1) {
				// create a leaf no matter the size of the bucket
				prQuadLeaf newLeaf = new prQuadLeaf(newConstructed, newIndexInsideSmall, curPRNode,
						numCellsPerDimPerBucket, curPRNode.getSmallCellSize(), pointsForS[i], pointsForS[i].size(),
						false);
				if (pointsForS[i].size() > 0) {
					curPRNode.addNewChild(newLeaf);
					// context.getCounter(Counters.SmallCellsNum).increment(1);
				}
				if (pointsForS[i].size() > 0)
					prLeaves.add(newLeaf);
				count = count + pointsForS[i].size();
//				System.out.println("Get node: " + newLeaf.printQuadLeaf());
			} else {
				// create a internal node and add the arraylist to the
				// hashmap structure
				prQuadInternal NewLeaf = new prQuadInternal(newConstructed, newIndexInsideSmall, curPRNode,
						numCellsPerDimPerBucket, curPRNode.getSmallCellSize());
				curPRNode.addNewChild(NewLeaf);
				mapQuadInternalWithPoints.put(NewLeaf, pointsForS[i]);
				count = count+ pointsForS[i].size();
				prQuadTree.push(NewLeaf);
//				System.out.println("Get node: " + NewLeaf.printQuadInternal());
			} // end else
		} // end for
			// mapQuadInternalWithPoints.remove(curPRNode);
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

	public boolean areaInsideSafeArea(float[] safeArea, float[] extendedArea) {
		if (extendedArea[0] >= safeArea[0] && extendedArea[1] <= safeArea[1] && extendedArea[2] >= safeArea[2]
				&& extendedArea[3] <= safeArea[3])
			return true;
		else
			return false;
	}

	public int[] getIndexInSmallCell() {
		return indexInSmallCell;
	}

	public void setIndexInSmallCell(int[] indexInSmallCell) {
		this.indexInSmallCell = indexInSmallCell;
	}

	public boolean checkIfParentNull() {
		if (parentNode == null)
			return true;
		else
			return false;
	}

	public prQuadNode getParentNode() {
		return parentNode;
	}

	public void setParentNode(prQuadNode parentNode) {
		this.parentNode = parentNode;
	}

	public float getSmallCellSize() {
		return smallCellSize;
	}

	public void setSmallCellSize(float smallCellSize) {
		this.smallCellSize = smallCellSize;
	}

	public String printPRQuadNode() {
		String str = "";
		str = str + "Point Coordinates: ";
		for (int i = 0; i < coordinates.length; i++) {
			str += coordinates[i] + ",";
		}
		if (parentNode == null)
			str = str + "NULL ";
		else
			str = str + "Not empty";
		for (int i = 0; i < numSmallCells.length; i++) {
			str += "numSmallCell " + i + " :" + numSmallCells[i] + ",   ";
		}
		str = str + "smallCellSize: " + smallCellSize;
		if (this.getClass().getName().endsWith("prQuadLeaf"))
			str = str + ",   #of points inside" + ((prQuadLeaf) this).getNumPoints();
		return str;
	}
}
