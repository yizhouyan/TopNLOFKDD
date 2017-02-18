package mcrtree.preprocessing;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import mcrtree.metricobject.MicroCluster;
import metricspace.*;
import net.cftree.mcrtree.CFTree;
/**
 * Call Build CF-Tree function
 * @author yizhouyan
 * 
 */
import util.SQConfig;

public class BuildCFTree {
	private CFTree birchTree;
	
	public BuildCFTree() {
		this.birchTree = new CFTree(SQConfig.maxNodeEntries, SQConfig.distThreshold, SQConfig.distFunction, false);
	}

	public void startTrainingCFTree() {
		// Read one instace at a time from the dataset
		// Dataset format: each line contain a set of value id, v1, v2, v3...
		// separated by comma
		BufferedReader in;
		try {
			in = new BufferedReader(new FileReader(SQConfig.dataset));
			String line = null;
			while ((line = in.readLine()) != null) {
				String[] tmp = line.split(SQConfig.sepStrForRecord);

				float[] x = new float[tmp.length - 1];
				for (int i = 1; i <= x.length; i++) {
					x[i - 1] = Float.parseFloat(tmp[i]);
				}

				// training birch, one instance at a time...
				boolean inserted = birchTree.insertEntry(x);
				if (!inserted) {
					System.err.println("ERROR: NOT INSERTED!");
					System.exit(1);
				}
			}
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		birchTree.finishedInsertingData();
	}

	/**
	 * Generate micro clusters based on the built tree, just generate center
	 * point coordinates
	 */
	public HashMap<Integer, MicroCluster> generateMicroClusters(IMetricSpace metricSpace) {
		return this.birchTree.generateMCFromCFTree(metricSpace);
	}

	public void printGeneratedTree() {
		// System.out.println("****************** CF-Tree *******************");
		// birchTree.printCFTree();
		System.out.println("****************** LEAVES *******************");
		birchTree.printLeafEntries();
	}

}
