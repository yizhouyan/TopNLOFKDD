package metricspace;

import java.util.HashMap;
import java.util.Map;

import util.PriorityQueue;

@SuppressWarnings("rawtypes")
public class MetricObject {

	private Object obj;
	// private Map<Long,Float> knnInDetail = new HashMap<Long,Float>();
	public PriorityQueue pointPQ = new PriorityQueue(PriorityQueue.SORT_ORDER_DESCENDING);
	// private Map<Long, coreInfoKNNs> knnMoreDetail = new HashMap<Long,
	// coreInfoKNNs>();
	private float kdist = -1;
	private float lrdValue = -1;
	private float lofValue = -1;
	private char type = 'F';

	private float nearestNeighborDist = Float.MAX_VALUE;
	private boolean canPrune = false;
	private int[] indexForSmallCell;
	private int indexOfCPCellInList = -1; // index of which cell it is in
	private float largeCellExpand = 0.0f;
	private boolean insideKNNfind = false;

	public MetricObject(Object obj) {
		this.obj = obj;
	}

	public float getNearestNeighborDist() {
		return nearestNeighborDist;
	}

	public void setNearestNeighborDist(float nearestNeighborDist) {
		this.nearestNeighborDist = nearestNeighborDist;
	}

	public char getType() {
		return type;
	}

	public void setType(char type) {
		this.type = type;
	}

	public float getLrdValue() {
		return lrdValue;
	}

	public void setLrdValue(float lrdValue) {
		this.lrdValue = lrdValue;
	}

	public float getLofValue() {
		return lofValue;
	}

	public void setLofValue(float lofValue) {
		this.lofValue = lofValue;
	}

	public int[] getIndexForSmallCell() {
		return indexForSmallCell;
	}

	public void setIndexForSmallCell(int[] indexForSmallCell) {
		this.indexForSmallCell = indexForSmallCell;
	}

	public boolean isCanPrune() {
		return canPrune;
	}

	public void setCanPrune(boolean canPrune) {
		this.canPrune = canPrune;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		Record r = (Record) obj;
		sb.append(", Knn in detail: ");
		return sb.toString();
	}

	public Object getObj() {
		return obj;
	}

	public void setObj(Object obj) {
		this.obj = obj;
	}

	public float getKdist() {
		return kdist;
	}

	public void setKdist(float kdist) {
		this.kdist = kdist;
	}

	@SuppressWarnings("unchecked")
	public static void main(String[] args) {

	}

	public int getIndexOfCPCellInList() {
		return indexOfCPCellInList;
	}

	public void setIndexOfCPCellInList(int indexOfCPCellInList) {
		this.indexOfCPCellInList = indexOfCPCellInList;
	}

	public float getLargeCellExpand() {
		return largeCellExpand;
	}

	public void setLargeCellExpand(float largeCellExpand) {
		this.largeCellExpand = largeCellExpand;
	}

	public PriorityQueue getPointPQ() {
		return pointPQ;
	}

	public void setPointPQ(PriorityQueue pointPQ) {
		this.pointPQ = pointPQ;
	}

	public boolean isInsideKNNfind() {
		return insideKNNfind;
	}

	public void setInsideKNNfind(boolean insideKNNfind) {
		this.insideKNNfind = insideKNNfind;
	}
}