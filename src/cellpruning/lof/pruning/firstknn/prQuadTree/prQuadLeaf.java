package cellpruning.lof.pruning.firstknn.prQuadTree;

import java.util.ArrayList;

import metricspace.MetricObject;

public class prQuadLeaf extends prQuadNode {
	private int numPoints;
	private ArrayList<MetricObject> listOfPoints;
	private boolean canPrune = false;

	public int getNumPoints() {
		return numPoints;
	}

	public void setNumPoints(int numPoints) {
		this.numPoints = numPoints;
	}

	public ArrayList<MetricObject> getListOfPoints() {
		return listOfPoints;
	}

	public void setListOfPoints(ArrayList<MetricObject> listOfPoints) {
		this.listOfPoints = listOfPoints;
	}

	public boolean isCanPrune() {
		return canPrune;
	}

	public void setCanPrune(boolean canPrune) {
		this.canPrune = canPrune;
	}

	public prQuadLeaf(float[] coordinates, int[] indexInSmallCell, prQuadNode parentNode, int[] numSmallCells,
			float smallCellSize, ArrayList<MetricObject> listOfPoints, int numPoints, boolean canPrune) {
		super(coordinates, indexInSmallCell, parentNode, numSmallCells, smallCellSize);
		this.listOfPoints = listOfPoints;
		this.numPoints = numPoints;
		this.canPrune = canPrune;
	}

	public String printQuadLeaf() {
		String str = super.printPRQuadNode();
		str = str + "Num of Points in the leaf: " + numPoints + "   , if can prune: " + canPrune;
		return str;
	}
}
