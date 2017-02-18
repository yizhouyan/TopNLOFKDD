package cellpruning.lof.pruning;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Stack;

import cellpruning.lof.pruning.firstknn.prQuadTree.*;
import metricspace.*;
import util.PriorityQueue;

public class LargeCellStore extends partitionTreeNode {

	/**
	 * Coordinates of the large cell, save independent dims, for example,
	 * dims[0,1,2,3], independentDims [2,3] then save coordinates for dims [2,3]
	 */
	private float[] coordinates;

	/** list of points in the Large grid cell */
	private ArrayList<MetricObject> listOfPoints;

	/** Closest pair distance inside the large cell */
	private float cpDist;

	/** Number of points inside the Large cell */
	private int numOfPoints;

	/** Break up into small cells, only break up those independent dims */
	private float smallCellSize;

	/**
	 * Keep track of number of small cells per dim, only for independent dims
	 */
	private int[] numSmallCells;

	private IMetric metric = null;

	private IMetricSpace metricSpace = null;

	private boolean breakIntoSmallCells = false;

	private prQuadInternal rootForPRTree;
	// save leaves that cannot be pruned
	ArrayList<prQuadLeaf> prLeaves;

	/** priority of the bucket, used for sorting */
	private double bucketPriority = 1;

	public LargeCellStore(float[] coordinates, ArrayList<MetricObject> listOfPoints, float cpDist, IMetric metric,
			IMetricSpace metricspace) {
		this.coordinates = coordinates.clone();
		this.listOfPoints = listOfPoints;
		this.numOfPoints = this.listOfPoints.size();
		this.cpDist = cpDist;
		// this.setIndexForLeaveNodesList(indexForLeaveNodesList);
		this.metric = metric;
		this.metricSpace = metricspace;
	}

	public int[] getNumSmallCells() {
		return numSmallCells;
	}

	public void setNumSmallCells(int[] numSmallCells) {
		this.numSmallCells = numSmallCells.clone();
	}

	public void addPoints(MetricObject newpoint) {
		listOfPoints.add(newpoint);
		numOfPoints++;
	}

	public ArrayList<MetricObject> getListOfPoints() {
		return listOfPoints;
	}

	public int getTotalNumOfDim() {
		return coordinates.length / 2;
	}

	/**
	 * compute priority for each large cell
	 * 
	 * @param isOnBoundary
	 */
	public void computePriorityForLargeCell(boolean isOnBoundary) {
		if (isOnBoundary) {
			this.bucketPriority = 0;
			return;
		} else
			this.bucketPriority = this.bucketPriority * Math.log(this.numOfPoints);
	}

	public float innerSearchWithEachLargeCell(HashMap<Long, MetricObject> CanPrunePoints,
			HashMap<Long, MetricObject> TrueKNNPoints, HashMap<Long, MetricObject> lrdHM, int K, int indexOfLeaveNodes,
			float thresholdLof, partitionTreeNode ptn, ArrayList<LargeCellStore> leaveNodes, PriorityQueue topnLOF,
			int topNNumber, int[] independentDims, int[] dependentDims, int totalDim)
			throws IOException, InterruptedException {

		HashMap<Long, MetricObject> TempCanPrunePoints = new HashMap<Long, MetricObject>();
		HashMap<Long, MetricObject> TempTrueKnnPoints = new HashMap<Long, MetricObject>();
		// start calculating LRD and LOF if possible
		// save those pruned points but need to recompute KNNs
		HashMap<Long, MetricObject> TempneedCalculatePruned = new HashMap<Long, MetricObject>();
		HashMap<Long, MetricObject> TemplrdHM = new HashMap<Long, MetricObject>();
		// save those cannot be pruned only by LRD value...
		HashMap<Long, MetricObject> TempneedCalLOF = new HashMap<Long, MetricObject>();
		// need more knn information, maybe knn is pruned...
		HashMap<Long, MetricObject> TempneedRecalLRD = new HashMap<Long, MetricObject>();

		// build index for the LargeCellStore
		if (this.numOfPoints == 0) {
			return thresholdLof;
		} else if (this.numOfPoints > 5 * K) {
			this.seperateToSmallCells(TempCanPrunePoints, indexOfLeaveNodes, thresholdLof, K, independentDims,
					dependentDims, totalDim);
			if (!this.breakIntoSmallCells && this.numOfPoints > K * 20) {
				this.seperateLargeNoPrune(K, indexOfLeaveNodes, independentDims, dependentDims);
			}
		}

		CanPrunePoints.putAll(TempCanPrunePoints);

		// Inner bucket KNN search
		if (this.breakIntoSmallCells) {
			this.findKnnsWithinPRTreeInsideBucket(TempTrueKnnPoints, this, K, independentDims);
		} else if (this.numOfPoints != 0) {
			// else find kNNs within the large cell
			this.findKnnsForLargeCellInsideBucket(TempTrueKnnPoints, this, K, independentDims);
		}

		TrueKNNPoints.putAll(TempTrueKnnPoints);
		// Inner bucket LRD computation
		for (MetricObject mo : TempTrueKnnPoints.values()) {
			ComputeLRD.CalLRDForSingleObject(mo, TempTrueKnnPoints, TempCanPrunePoints, TempneedCalculatePruned,
					TemplrdHM, TempneedCalLOF, TempneedRecalLRD, thresholdLof, leaveNodes, K);
		}

		// Inner bucket LOF computation
		for (MetricObject mo : TempneedCalLOF.values()) {
			ComputeLOF.CalLOFForSingleObject(mo, TemplrdHM, K, thresholdLof);
			if (mo.getType() == 'O' && mo.getLofValue() > thresholdLof) {
				float tempLofValue = mo.getLofValue();
				if (topnLOF.size() < topNNumber) {
					topnLOF.insert(metricSpace.getID(mo.getObj()), tempLofValue);

				} else if (tempLofValue > topnLOF.getPriority()) {
					topnLOF.pop();
					topnLOF.insert(metricSpace.getID(mo.getObj()), tempLofValue);
					if (thresholdLof < topnLOF.getPriority())
						thresholdLof = topnLOF.getPriority();
					// System.out.println("Threshold updated: " +
					// thresholdLof);
				}
			}
		}

		if (topnLOF.size() == topNNumber && thresholdLof < topnLOF.getPriority())
			thresholdLof = topnLOF.getPriority();
		lrdHM.putAll(TemplrdHM);
		return thresholdLof;
	}

	public float multiplyArray(float[] array) {
		float result = 1;
		for (float element : array) {
			result *= element;
		}
		return result;
	}

	// partitionTreeNode ptn, float[] partition_store,
	public void seperateToSmallCells(HashMap<Long, MetricObject> CanPrunePoints, int indexOfLeaveNodes, float threshold,
			int K, int[] independentDims, int[] correlatedDims, int totalDim) {
		// save independent dim range
		float[] LargeCellRange = new float[independentDims.length];
		for (int i = 0; i < independentDims.length; i++) {
			LargeCellRange[i] = coordinates[i * 2 + 1] - coordinates[i * 2];
		}

		// smallCellSize using closest pair and threshold
		smallCellSize = (float) (threshold * this.cpDist / (2 * Math.sqrt(totalDim)));

		// cell calculated as average cell size
		float smallCellSize_predict = (float) Math.pow(multiplyArray(LargeCellRange) * K / numOfPoints,
				1.0 / LargeCellRange.length);

		// if the small cell size too small, then don't use this size to build
		// PRQuadTree
		if (smallCellSize < smallCellSize_predict / 10) {
			return;
		}
		if (smallCellSize > smallCellSize_predict / 5) {
			smallCellSize = smallCellSize_predict / 5;
		}

		for (float LargeCellRangeTemp : LargeCellRange) {
			if (smallCellSize >= LargeCellRangeTemp)
				return;
		}
		numSmallCells = new int[independentDims.length];
		// calculate how many small cells for each partition per dimension
		for (int i = 0; i < LargeCellRange.length; i++) {
			numSmallCells[i] = (int) Math.ceil(LargeCellRange[i] / smallCellSize);
		}

		for (int i = 0; i < LargeCellRange.length; i++) {
			if (numSmallCells[i] < 10)
				return;
		}

		breakIntoSmallCells = true;

		for (MetricObject mo : listOfPoints) {
			Record record = (Record) mo.getObj();
			int[] indexInSmall = new int[independentDims.length];
			for (int i = 0; i < independentDims.length; i++) {
				float tempValue = record.getValue()[independentDims[i]];
				indexInSmall[i] = (int) (Math.floor(tempValue - coordinates[2 * i])
						/ (smallCellSize + Float.MIN_VALUE));
			}
			mo.setIndexForSmallCell(indexInSmall);
		}
		// build up PR quadtree
		float[] largeCellCoor = coordinates.clone();
		buildPRQuadTree(CanPrunePoints, numSmallCells, smallCellSize, listOfPoints, numOfPoints, largeCellCoor,
				indexOfLeaveNodes, independentDims, correlatedDims, K, true);
	}

	public boolean QuerySurroundingBucketsForCP(partitionTreeNode ptn, float[] expectedRange) {
		Stack<partitionTreeInternal> partitionTree = new Stack<partitionTreeInternal>();
		partitionTree.push((partitionTreeInternal) ptn);
		while (!partitionTree.isEmpty()) {
			partitionTreeInternal tempInternal = partitionTree.pop();
			ArrayList<partitionTreeNode> tempChildNodes = tempInternal.getChildNodes();
			// check children
			for (int i = 0; i < tempChildNodes.size(); i++) {
				if (tempChildNodes.get(i).getClass().getName().endsWith("partitionTreeInternal")) {
					float[] tempCoordinates = ((partitionTreeInternal) tempChildNodes.get(i)).getCoordinates();
					if (checkRange(expectedRange, tempCoordinates)) {
						partitionTree.push((partitionTreeInternal) tempChildNodes.get(i));
					}
				} else if (tempChildNodes.get(i).getClass().getName().endsWith("LargeCellStore")) {
					float[] tempCoordinates = ((LargeCellStore) tempChildNodes.get(i)).getCoordinates();
					if (checkRange(expectedRange, tempCoordinates)
							&& ((LargeCellStore) tempChildNodes.get(i)).getCpDist() < this.cpDist) {
						return false;
					}
				} else {
					System.out.println("Unknown Bucket Node Type!");
				}
			}
		}
		return true;
	}

	public void seperateLargeNoPrune(int K, int indexOfLeaveNodes, int[] independentDims, int[] correlatedDims) {
		float[] LargeCellRange = new float[independentDims.length];
		for (int i = 0; i < independentDims.length; i++) {
			LargeCellRange[i] = coordinates[i * 2 + 1] - coordinates[i * 2];
		}

		smallCellSize = (float) Math.pow(multiplyArray(LargeCellRange) * 4 * K / numOfPoints,
				1.0 / LargeCellRange.length);
		for (float LargeCellRangeTemp : LargeCellRange) {
			if (smallCellSize >= LargeCellRangeTemp)
				return;
		}
		numSmallCells = new int[independentDims.length];
		// calculate how many small cells for each partition per dimension
		for (int i = 0; i < LargeCellRange.length; i++) {
			numSmallCells[i] = (int) Math.ceil(LargeCellRange[i] / smallCellSize);
		}
		breakIntoSmallCells = true;
		for (MetricObject mo : listOfPoints) {
			Record record = (Record) mo.getObj();
			int[] indexInSmall = new int[independentDims.length];
			for (int i = 0; i < independentDims.length; i++) {
				float tempValue = record.getValue()[independentDims[i]];
				indexInSmall[i] = (int) (Math.floor(tempValue - coordinates[2 * i])
						/ (smallCellSize + Float.MIN_VALUE));
			}
			mo.setIndexForSmallCell(indexInSmall);
		}
		// build up PR quadtree
		float[] largeCellCoor = coordinates.clone();
		buildPRQuadTree(null, numSmallCells, smallCellSize, listOfPoints, numOfPoints, largeCellCoor, indexOfLeaveNodes,
				independentDims, correlatedDims, K, false);
	}

	/**
	 * build up PR quad tree with information from the large cell
	 * 
	 * @param numSmallCells
	 * @param smallCellSize
	 * @param listOfPoints
	 * @param numOfPoints
	 * @param largeCellCoor
	 * @return root of PR QuadTree
	 */
	public void buildPRQuadTree(HashMap<Long, MetricObject> CanPrunePoints, int[] numSmallCells, float smallCellSize,
			ArrayList<MetricObject> listOfPoints, int numOfPoints, float[] largeCellCoor, int indexOfLeaveNodes,
			int[] independentDims, int[] correlatedDims, int K, boolean withPrune) {
		// init root
		int[] indexRangeInSmallCell = new int[numSmallCells.length * 2];
		for (int i = 0; i < numSmallCells.length; i++) {
			indexRangeInSmallCell[2 * i] = 0;
			indexRangeInSmallCell[2 * i + 1] = numSmallCells[i] - 1;
		}

		rootForPRTree = new prQuadInternal(largeCellCoor, indexRangeInSmallCell, null, numSmallCells, smallCellSize);
		Stack<prQuadInternal> prQuadTree = new Stack<prQuadInternal>();
		HashMap<prQuadInternal, ArrayList<MetricObject>> mapQuadInternalWithPoints = new HashMap<prQuadInternal, ArrayList<MetricObject>>();
		mapQuadInternalWithPoints.put(rootForPRTree, listOfPoints);
		// save leaves in the pr tree
		prLeaves = new ArrayList<prQuadLeaf>();
		prQuadTree.push(rootForPRTree);
		int count = 0;
		while (!prQuadTree.empty()) {
			/**
			 * pop up the quad node and divide into 4 parts check if each part
			 * contains enough points if contains K+1 points, create a
			 * prQuadInternal and push to stack if contains less than K points,
			 * create a prQuadLeaf and save the pointer if can not divide (reach
			 * minimum size), create prQuadLeaf and save the pointer
			 */
			prQuadInternal curPRNode = prQuadTree.pop();
			curPRNode.generateChilden(CanPrunePoints, curPRNode, prQuadTree, mapQuadInternalWithPoints, prLeaves,
					largeCellCoor, numSmallCells, indexOfLeaveNodes, independentDims, correlatedDims, K, withPrune);
		}
	}

	public void traverseLargeCell(MetricObject curPoint, LargeCellStore large_cell_store, int K) {
		// traverse points
		float dist = 0.0f;
		float theta;
		if (curPoint.pointPQ.size() > 0)
			theta = curPoint.pointPQ.getPriority();
		else
			theta = Float.POSITIVE_INFINITY;
		for (int i = 0; i < large_cell_store.getNumOfPoints(); i++) {
			MetricObject o_S = large_cell_store.getListOfPoints().get(i);
			if (((Record) o_S.getObj()).getRId() == ((Record) curPoint.getObj()).getRId()) {
				continue;
			} else if (o_S.getType() == 'C')
				continue;
			else {
				try {
					dist = metric.dist(curPoint.getObj(), o_S.getObj());
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if (curPoint.pointPQ.size() < K) {
					curPoint.pointPQ.insert(metricSpace.getID(o_S.getObj()), dist);
					theta = curPoint.pointPQ.getPriority();
				} else if (dist < theta) {
					curPoint.pointPQ.pop();
					curPoint.pointPQ.insert(metricSpace.getID(o_S.getObj()), dist);
					theta = curPoint.pointPQ.getPriority();
				}
			}
		}
	}

	public void findKnns(MetricObject curPoint, prQuadLeaf curCheckLeaf, int K) {
		// traverse points
		float dist = 0.0f;
		float theta;
		if (curPoint.pointPQ.size() > 0)
			theta = curPoint.pointPQ.getPriority();
		else
			theta = Float.POSITIVE_INFINITY;
		for (int i = 0; i < curCheckLeaf.getNumPoints(); i++) {
			MetricObject o_S = curCheckLeaf.getListOfPoints().get(i);
			if (((Record) o_S.getObj()).getRId() == ((Record) curPoint.getObj()).getRId()) {
				continue;
			} else if (o_S.getType() == 'C')
				continue;
			else {
				try {
					dist = metric.dist(curPoint.getObj(), o_S.getObj());
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if (curPoint.pointPQ.size() < K) {
					curPoint.pointPQ.insert(metricSpace.getID(o_S.getObj()), dist);
					theta = curPoint.pointPQ.getPriority();
				} else if (dist < theta) {
					curPoint.pointPQ.pop();
					curPoint.pointPQ.insert(metricSpace.getID(o_S.getObj()), dist);
					theta = curPoint.pointPQ.getPriority();
				}
			}
		}
	}

	static boolean checkRange(float[] expectedRange, float[] checkedRange) {
		for (int i = 0; i < expectedRange.length / 2; i++) {
			if (expectedRange[2 * i] > checkedRange[2 * i + 1] || checkedRange[2 * i] > expectedRange[2 * i + 1])
				return false;
		}
		return true;
	}

	/**
	 * check if the extended area is in checked area
	 * 
	 * @param checkedArea
	 * @param extendedArea
	 * @return
	 */
	boolean insideCheckedArea(float[] checkedArea, float[] extendedArea) {
		for (int i = 0; i < checkedArea.length / 2; i++) {
			if (extendedArea[2 * i] < checkedArea[2 * i] || extendedArea[2 * i + 1] > checkedArea[2 * i + 1])
				return false;
		}
		return true;
	}

	public void findKnnsWithinOneCell(MetricObject curPoint, prQuadLeaf curLeaf, LargeCellStore large_cell_store, int K,
			int[] independentDims) {
		float[] curPointCoor = ((Record) curPoint.getObj()).getValue();
		float kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
		Stack<prQuadNode> checkWithinOneTree = new Stack<prQuadNode>();
		prQuadNode tempCheckBreakNode = curLeaf;
		float[] extendArea = new float[independentDims.length * 2];
		float[] largeCellCoor = large_cell_store.coordinates;
		for (int i = 0; i < independentDims.length; i++) {
			extendArea[2 * i] = Math.max(largeCellCoor[2 * i], curPointCoor[independentDims[i]] - kdist);
			extendArea[2 * i + 1] = Math.min(largeCellCoor[2 * i + 1], curPointCoor[independentDims[i]] + kdist);
		}
		float[] checkedCoordinates = tempCheckBreakNode.getCoordinates();
		if (insideCheckedArea(checkedCoordinates, extendArea)) {
			return;
		}
		while (tempCheckBreakNode.getParentNode() != null) {
			// first add the parent node
			checkWithinOneTree.push(tempCheckBreakNode.getParentNode());
			// then add brothers
			for (prQuadNode brother : ((prQuadInternal) tempCheckBreakNode.getParentNode()).getChildNodes()) {
				if (brother.equals(tempCheckBreakNode)) /////////
					continue;
				else if (checkRange(extendArea, brother.getCoordinates())) {
					checkWithinOneTree.push(brother);
				} else
					continue;
			}
			// traverse the stack until has only one element(the parent) left,
			// all brother traversed
			while (checkWithinOneTree.size() > 1) {
				prQuadNode tempNode = checkWithinOneTree.pop();
				if (tempNode.getClass().getName().endsWith("prQuadLeaf")) {
					if (tempNode.equals(curLeaf)) {
						continue;
					}
					findKnns(curPoint, (prQuadLeaf) tempNode, K);
					kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
					// update extend area
					float[] newExtendArea = new float[independentDims.length * 2];
					for (int i = 0; i < independentDims.length; i++) {
						newExtendArea[2 * i] = Math.max(largeCellCoor[2 * i], curPointCoor[independentDims[i]] - kdist);
						newExtendArea[2 * i + 1] = Math.min(largeCellCoor[2 * i + 1],
								curPointCoor[independentDims[i]] + kdist);
					}
					extendArea = newExtendArea;
					if (insideCheckedArea(checkedCoordinates, extendArea)) {
						return;
					}
				} else {
					// add children
					for (prQuadNode children : ((prQuadInternal) tempNode).getChildNodes()) {
						if (children.equals(curLeaf))
							continue;
						else if (checkRange(extendArea, children.getCoordinates())) {
							checkWithinOneTree.push(children);
						} else
							continue;
					}
				} // end else
			} // end while
				// stack size == 1 only one parent left, all brothers traversed
			tempCheckBreakNode = checkWithinOneTree.pop(); // let parent be the
															// new node
			checkedCoordinates = tempCheckBreakNode.getCoordinates();
			if (insideCheckedArea(checkedCoordinates, extendArea)) {
				return;
			}
		} // end while
	}

	public void savePriorityQueueToKNN(MetricObject curPoint, int K) {
		if (curPoint.pointPQ.size() != K) {
			System.out.println("Less than K points in Priority Queue");
			// curPoint.setCanPrune(true);
			return;
		}
		curPoint.setKdist(curPoint.pointPQ.getPriority());
		float NNDist = Float.MAX_VALUE;
		float[] allPriority = curPoint.pointPQ.getPrioritySet();
		for (int i = 0; i < allPriority.length; i++) {
			NNDist = Math.min(NNDist, allPriority[i]);
		}
		curPoint.setNearestNeighborDist(NNDist);
		// curPoint.setExpandDist(expandDist);
		curPoint.setType('T');
	}

	/**
	 * search the partitionTreeNode and find supporting large cells
	 * 
	 * @param ExtendArea
	 *            the coordinate of the point's extended area
	 * @param ptn
	 *            the root node of the partition tree (binary tree)
	 * @param currentCell
	 *            current LargeCellStore, supportingLargeCells will not include
	 *            current cell
	 * @return the supporting cells for the current point
	 */
	public ArrayList<LargeCellStore> searchSupportingLargeCells(float[] ExtendArea, partitionTreeNode ptn,
			LargeCellStore currentCell) {
		ArrayList<LargeCellStore> supportingLargeCells = new ArrayList<LargeCellStore>();
		Stack<partitionTreeInternal> stackOfInternals = new Stack<partitionTreeInternal>();
		if (ptn.getClass().getName().endsWith("LargeCellStore")) {
			return supportingLargeCells;
		} else if (ptn.getClass().getName().endsWith("partitionTreeInternal")) {
			stackOfInternals.push((partitionTreeInternal) ptn);
		}
		while (!stackOfInternals.isEmpty()) {
			// check the coordinates of each child node
			ArrayList<partitionTreeNode> tempChildNodes = stackOfInternals.pop().getChildNodes();
			for (int i = 0; i < tempChildNodes.size(); i++) {
				partitionTreeNode tempPTN = tempChildNodes.get(i);
				if (tempPTN.getClass().getName().endsWith("LargeCellStore")
						&& checkRange(ExtendArea, ((LargeCellStore) tempPTN).getCoordinates())
						&& !(tempPTN == currentCell)) {
					// if(tempPTN == currentCell){
					// System.out.println("This is exactly the same cell~ not a
					// support");
					// }
					supportingLargeCells.add((LargeCellStore) tempPTN);
				} else if (tempPTN.getClass().getName().endsWith("partitionTreeInternal")
						&& checkRange(ExtendArea, ((partitionTreeInternal) tempPTN).getCoordinates())) {
					stackOfInternals.push((partitionTreeInternal) tempPTN);
				}
			}
		}
		return supportingLargeCells;
	}

	public void findKnnsForOnePointInsideBucket(HashMap<Long, MetricObject> TrueKnnPoints, MetricObject curPoint,
			prQuadLeaf curLeaf, LargeCellStore currentLeafNode, int K, int[] independentDims) {
		curPoint.setInsideKNNfind(true);
		float kdist = Float.POSITIVE_INFINITY;

		// first find kNNs within the large cell and bound a partition area for
		// largeCell
		// first find kNNs within the leaf
		findKnns(curPoint, curLeaf, K);
		kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;

		// then find KNNs within the large cell
		findKnnsWithinOneCell(curPoint, curLeaf, currentLeafNode, K, independentDims);
		kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
		// System.out.println("old kdistance: " + kdist);
		// check if kNNs exceeds the large cell
		curPoint.setLargeCellExpand(calExtendDistance(currentLeafNode, curPoint, kdist, independentDims));
		// if not exceed the large cell, don't need to traverse other large
		// cells
		if (curPoint.getLargeCellExpand() <= 1e-9) {
			savePriorityQueueToKNN(curPoint, K);
			TrueKnnPoints.put(((Record) curPoint.getObj()).getRId(), curPoint);
		}
	}

	// prQuadLeaf curLeaf
	public void findKnnsForOnePointOutsideBucket(HashMap<Long, MetricObject> TrueKnnPoints, MetricObject curPoint,
			ArrayList<LargeCellStore> large_cell_store, LargeCellStore currentLeafNode, partitionTreeNode ptn,
			float[] partition_store, int K, int num_dims, int[] independentDims) {
		float kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;

		// if exceed the large cell, traverse nearby large cells
		// include more supporting cells
		float[] curPointCoor = ((Record) curPoint.getObj()).getValue();
		float[] ExtendArea = new float[independentDims.length * 2];
		for (int i = 0; i < independentDims.length; i++) {
			ExtendArea[2 * i] = (float) Math.max(partition_store[2 * i], curPointCoor[independentDims[i]] - kdist);
			ExtendArea[2 * i + 1] = (float) Math.min(partition_store[2 * i + 1],
					curPointCoor[independentDims[i]] + kdist);
		}

		ArrayList<LargeCellStore> supportingLargeCells = searchSupportingLargeCells(ExtendArea, ptn, currentLeafNode);
		// System.out.println("Size of Support:" + supportingLargeCells.size());
		// for each supporting cell, traverse until not exceed the checked area
		for (LargeCellStore supportingCell : supportingLargeCells) {
			if (checkRange(ExtendArea, supportingCell.getCoordinates())) {
				// find a leaf to start
				if (supportingCell.breakIntoSmallCells) {
					findOverLapTreeNodesAndSearch(supportingCell.getRootForPRTree(), ExtendArea, curPoint, K);
					// prQuadLeaf tempLeaf =
					// RangeQuery(supportingCell.getRootForPRTree(),
					// ExtendArea);
					// if (tempLeaf != null) {
					// // then find KNNs within the large cell
					// findKnnsWithinOneCell(curPoint, tempLeaf, supportingCell,
					// K, independentDims);
					// kdist = curPoint.pointPQ.size() == K ?
					// curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
					// float[] newExtendArea = new float[independentDims.length
					// * 2];
					// for (int i = 0; i < independentDims.length; i++) {
					// newExtendArea[2 * i] = (float) Math.max(partition_store[2
					// * i],
					// curPointCoor[independentDims[i]] - kdist);
					// newExtendArea[2 * i + 1] = (float)
					// Math.min(partition_store[2 * i + 1],
					// curPointCoor[independentDims[i]] + kdist);
					// }
					// ExtendArea = newExtendArea;
					//
					// }
//					traverseLargeCell(curPoint, supportingCell, K);
					kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
					float[] newExtendArea = new float[independentDims.length * 2];
					for (int i = 0; i < independentDims.length; i++) {
						newExtendArea[2 * i] = (float) Math.max(partition_store[2 * i],
								curPointCoor[independentDims[i]] - kdist);
						newExtendArea[2 * i + 1] = (float) Math.min(partition_store[2 * i + 1],
								curPointCoor[independentDims[i]] + kdist);
					}
					ExtendArea = newExtendArea;
				} // end if
				else { // traverse Large cell
					traverseLargeCell(curPoint, supportingCell, K);
					kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
					float[] newExtendArea = new float[independentDims.length * 2];
					for (int i = 0; i < independentDims.length; i++) {
						newExtendArea[2 * i] = (float) Math.max(partition_store[2 * i],
								curPointCoor[independentDims[i]] - kdist);
						newExtendArea[2 * i + 1] = (float) Math.min(partition_store[2 * i + 1],
								curPointCoor[independentDims[i]] + kdist);
					}
					ExtendArea = newExtendArea;
				}
			} // end if(checkRange(ExtendArea, supportingCell.getCoordinates()))
		} // end for
			// bound supporting area for the partition
		savePriorityQueueToKNN(curPoint, K);
		TrueKnnPoints.put(((Record) curPoint.getObj()).getRId(), curPoint);
	}

	public void findOverLapTreeNodesAndSearch(prQuadInternal prRoot, float[] expectedRange, MetricObject curPoint,
			int K) {
		ArrayList<prQuadLeaf> prQuadLeafList = new ArrayList<prQuadLeaf>();
		Stack<prQuadInternal> prQuadTreeInternal = new Stack<prQuadInternal>();
		prQuadTreeInternal.push(prRoot);
		while (!prQuadTreeInternal.empty()) {
			prQuadInternal curPRNode = prQuadTreeInternal.pop();
			// traverse 4 childs and save to the stack if inside the expecting
			// range
			for (prQuadNode tempNode : curPRNode.getChildNodes()) {
				// check range
				if (!checkRange(expectedRange, tempNode.getCoordinates()))
					continue;
				if (tempNode.getClass().getName().endsWith("prQuadInternal")) {
					prQuadTreeInternal.push((prQuadInternal) tempNode);
				} else { // leaf
					prQuadLeafList.add((prQuadLeaf) tempNode);
				}
			}
		}
		// traverse each leaf nodes to find kNN
		float dist = 0.0f;
		float theta;
		if (curPoint.pointPQ.size() == K)
			theta = curPoint.pointPQ.getPriority();
		else
			theta = Float.POSITIVE_INFINITY;
		for (prQuadLeaf currentLeaf : prQuadLeafList) {
			for (MetricObject o_S : currentLeaf.getListOfPoints()) {
				if (((Record) o_S.getObj()).getRId() == ((Record) curPoint.getObj()).getRId()) {
					continue;
				} else if (o_S.getType() == 'C')
					continue;
				else {
					try {
						dist = metric.dist(curPoint.getObj(), o_S.getObj());
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					if (curPoint.pointPQ.size() < K) {
						curPoint.pointPQ.insert(metricSpace.getID(o_S.getObj()), dist);
						theta = curPoint.pointPQ.getPriority();
					} else if (dist < theta) {
						curPoint.pointPQ.pop();
						curPoint.pointPQ.insert(metricSpace.getID(o_S.getObj()), dist);
						theta = curPoint.pointPQ.getPriority();
					}
				}
			}
		}
	}

	public prQuadLeaf findLeafWithSmallCellIndex(prQuadInternal prRoot, int[] indexForSmallCells,
			int[] independentDims) {
		Stack<prQuadInternal> prQuadTreeInternal = new Stack<prQuadInternal>();
		prQuadTreeInternal.push(prRoot);
		// System.out.println(xx + "," + yy);
		while (!prQuadTreeInternal.empty()) {
			prQuadInternal curPRNode = prQuadTreeInternal.pop();
			// traverse childs and save to the stack if inside the expecting
			// range
			for (prQuadNode tempNode : curPRNode.getChildNodes()) {
				boolean tempFlag = true;
				int[] indexRangeOfSmallCell = tempNode.getIndexInSmallCell();
				for (int i = 0; i < independentDims.length; i++) {
					if (indexRangeOfSmallCell[2 * i] > indexForSmallCells[i]
							|| indexForSmallCells[i] > indexRangeOfSmallCell[2 * i + 1]) {
						tempFlag = false;
						break;
					}
				}
				if (tempFlag) {
					if (tempNode.getClass().getName().endsWith("prQuadInternal")) {
						prQuadTreeInternal.push((prQuadInternal) tempNode);
					} else { // leaf
						return (prQuadLeaf) tempNode;
					}
				} else
					continue;
			}
		}
		return null;
	}

	public prQuadLeaf RangeQuery(prQuadInternal prRoot, float[] expectedRange) {
		Stack<prQuadInternal> prQuadTreeInternal = new Stack<prQuadInternal>();
		prQuadTreeInternal.push(prRoot);
		while (!prQuadTreeInternal.empty()) {
			prQuadInternal curPRNode = prQuadTreeInternal.pop();
			// traverse 4 childs and save to the stack if inside the expecting
			// range
			for (prQuadNode tempNode : curPRNode.getChildNodes()) {
				// check range
				if (!checkRange(expectedRange, tempNode.getCoordinates()))
					continue;
				if (tempNode.getClass().getName().endsWith("prQuadInternal")) {
					prQuadTreeInternal.push((prQuadInternal) tempNode);
				} else { // leaf
					return (prQuadLeaf) tempNode;
				}
			}
		}
		return null;
	}

	public float calExtendDistance(LargeCellStore large_cell_store, MetricObject curPoint, float kdist,
			int[] independentDims) {
		float extendDist = 0.0f;
		float[] curPointCoor = ((Record) curPoint.getObj()).getValue();
		float[] largeCellCoor = large_cell_store.coordinates;
		for (int i = 0; i < independentDims.length; i++) {
			extendDist = (float) Math.max(extendDist,
					largeCellCoor[2 * i] - (curPointCoor[independentDims[i]] - kdist));
			extendDist = (float) Math.max(extendDist,
					(curPointCoor[independentDims[i]] + kdist) - largeCellCoor[2 * i + 1]);
		}
		return extendDist;
	}

	public void findKnnsWithinPRTreeInsideBucket(HashMap<Long, MetricObject> TempTrueKnnPoints,
			LargeCellStore large_cell_store, int K, int[] independentDims) {
		// find kNNs for each point in the prLeaves
		for (prQuadLeaf curLeaf : this.prLeaves) {
			for (MetricObject curPoint : curLeaf.getListOfPoints()) {
				findKnnsForOnePointInsideBucket(TempTrueKnnPoints, curPoint, curLeaf, large_cell_store, K,
						independentDims);
			}
		}
	}

	public void findKnnsWithinPRTreeOutsideBucket(HashMap<Long, MetricObject> TrueKnnPoints,
			ArrayList<LargeCellStore> large_cell_store, int xx, partitionTreeNode ptn, float[] partition_store, int K,
			int num_dims, int[] independentDims) {
		// find kNNs for each point in the prLeaves
		for (prQuadLeaf curLeaf : this.prLeaves) {
			for (MetricObject curPoint : curLeaf.getListOfPoints()) {
				if (curPoint.getType() == 'F' && !curPoint.isCanPrune())
					findKnnsForOnePointOutsideBucket(TrueKnnPoints, curPoint, large_cell_store,
							large_cell_store.get(xx), ptn, partition_store, K, num_dims, independentDims);
			}
		}
	}

	public void findKnnsForLargeCellInsideBucket(HashMap<Long, MetricObject> TrueKnnPoints,
			LargeCellStore large_cell_store, int K, int[] independentDims) {
		// find knns for each point in the large cell
		for (MetricObject curPoint : large_cell_store.getListOfPoints()) {
			findKnnsForOnePointInLargeCellInsideBucket(TrueKnnPoints, curPoint, large_cell_store, K, independentDims);
		}
	}

	public void findKnnsForLargeCellOutsideBucket(HashMap<Long, MetricObject> TrueKnnPoints,
			ArrayList<LargeCellStore> large_cell_store, int xx, partitionTreeNode ptn, float[] partition_store, int K,
			int num_dims, int[] independentDims) {
		// find knns for each point in the large cell
		for (MetricObject curPoint : large_cell_store.get(xx).getListOfPoints()) {
			if (curPoint.getType() == 'F')
				findKnnsForOnePointOutsideBucket(TrueKnnPoints, curPoint, large_cell_store, large_cell_store.get(xx),
						ptn, partition_store, K, num_dims, independentDims);
		}
	}

	public void findKnnsForOnePointInLargeCellInsideBucket(HashMap<Long, MetricObject> TrueKnnPoints,
			MetricObject curPoint, LargeCellStore large_cell_store, int K, int[] independentDims) {
		curPoint.setInsideKNNfind(true);
		float kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
		// first traverse the large cell
		traverseLargeCell(curPoint, large_cell_store, K);
		kdist = curPoint.pointPQ.size() == K ? curPoint.pointPQ.getPriority() : Float.POSITIVE_INFINITY;
		// check if kNNs exceeds the large cell
		curPoint.setLargeCellExpand(calExtendDistance(large_cell_store, curPoint, kdist, independentDims));
		if (curPoint.getLargeCellExpand() <= 1e-9) {
			savePriorityQueueToKNN(curPoint, K);
			TrueKnnPoints.put(((Record) curPoint.getObj()).getRId(), curPoint);
		}
	}

	public boolean isBreakIntoSmallCells() {
		return breakIntoSmallCells;
	}

	public void setBreakIntoSmallCells(boolean breakIntoSmallCells) {
		this.breakIntoSmallCells = breakIntoSmallCells;
	}

	/** print Large cell store */
	public String printCellStoreWithSupport() {
		String str = "";
		for (float x : coordinates) {
			str += x + ",";
		}
		str += "\n";
		str += "number of points:" + numOfPoints + "\n" + "closest pair distance: " + cpDist + "\n"
				+ "Points in detail: ";
		for (Iterator<MetricObject> itr = listOfPoints.iterator(); itr.hasNext();) {
			str = str + itr.next().getObj().toString() + "\n";

		}
		return str.substring(0, str.length());
	}

	public String printCellStoreDetailedInfo() {
		String str = "";
		for (float x : coordinates) {
			str += x + ",";
		}
		str += "\n";
		str += "number of points:" + numOfPoints + "\n" + "closest pair distance: " + cpDist + "\n";
		str += "isbreakup?" + breakIntoSmallCells + "\n";
		str += "small cell size = " + smallCellSize;
		return str.substring(0, str.length());
	}

	public float getCpDist() {
		return cpDist;
	}

	public void setCpDist(float cpDist) {
		this.cpDist = cpDist;
	}

	public float getSmallCellSize() {
		return smallCellSize;
	}

	public void setSmallCellSize(float smallCellSize) {
		this.smallCellSize = smallCellSize;
	}

	public IMetric getMetric() {
		return metric;
	}

	public void setMetric(IMetric metric) {
		this.metric = metric;
	}

	public int getNumOfPoints() {
		return numOfPoints;
	}

	public void setNumOfPoints(int numOfPoints) {
		this.numOfPoints = numOfPoints;
	}

	public float[] getCoordinates() {
		return this.coordinates;
	}

	public prQuadInternal getRootForPRTree() {
		return rootForPRTree;
	}

	public void setRootForPRTree(prQuadInternal rootForPRTree) {
		this.rootForPRTree = rootForPRTree;
	}

	public ArrayList<prQuadLeaf> getPrLeaves() {
		return prLeaves;
	}

	public void setPrLeaves(ArrayList<prQuadLeaf> prLeaves) {
		this.prLeaves = prLeaves;
	}

	public double getBucketPriority() {
		return bucketPriority;
	}

	public void setBucketPriority(double bucketPriority) {
		this.bucketPriority = bucketPriority;
	}

}
