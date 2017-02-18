package cellpruning.lof.pruning;

import java.util.ArrayList;

public class partitionTreeInternal extends partitionTreeNode {
	private ArrayList<partitionTreeNode> childNodes;
	private float[] coordinates;

	public partitionTreeInternal(float[] coordinates) {
		this.coordinates = new float[coordinates.length];
		this.coordinates = coordinates;
		childNodes = new ArrayList<partitionTreeNode>();
	}

	public void addNewChild(partitionTreeNode newChild) {
		childNodes.add(newChild);
	}

	public ArrayList<partitionTreeNode> getChildNodes() {
		return childNodes;
	}

	public void setChildNodes(ArrayList<partitionTreeNode> childNodes) {
		this.childNodes = childNodes;
	}

	public float[] getCoordinates() {
		return coordinates;
	}

	public String printQuadInternal() {
		String str = "";
		str = str + "Point Coordinates: ";
		str = str + coordinates[0] + "," + coordinates[1] + "," + coordinates[2] + "," + coordinates[3] + ",";
		if (this.parentNode == null)
			str = str + "NULL ";
		else
			str = str + "Not empty";
		str = str + "Num of Childs: " + childNodes.size();
		return str;
	}

}
