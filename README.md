Experimental Code for KDD 2017 
For paper: Scalable Top-N Local Outlier Detection 
If you find any problems, please contact yyan2@wpi.edu.

Includes three parts: 

1. Baseline Methods 

	Two different indexing methods for KNN search are implemented. 

	--- Pivot-based KNN Search : 

		Main Class: baseline.pivotknn.ComputeTopNLOF

		From Paper: Bhaduri, Kanishka, Bryan L. Matthews, and Chris R. Giannella. "Algorithms for speeding up distance-based outlier detection." Proceedings of the 17th ACM SIGKDD international conference on Knowledge Discovery and Data Mining. ACM, 2011.

	--- R-Tree based KNN Search: 

		Main Class: baseline.rtreeknn.ComputeTopNLOF

		From Paper: Guttman, Antonin. R-trees: a dynamic index structure for spatial searching. Vol. 14. No. 2. ACM, 1984.

2. MicroCluster --- State-of-the-art method: 

		Main Class: microcluster.topnlof.TopNLOFDetection

		From Paper: Jin, Wen, Anthony KH Tung, and Jiawei Han. "Mining top-n local outliers in large databases." Proceedings of the seventh ACM SIGKDD international conference on Knowledge discovery and data mining. ACM, 2001.
3. TOLF --- Proposed method
		
		Two versions of TOLF: Single-thread and Multi-thread version. 

		Main Class for Single-Thread Version: cellpruning.lof.pruning.ComputeTopNLOFWithPruning

		Main Class for Multi-Thread Version: cellpruning.lof.pruning.MultiThread.ComputeTopNLOFWithPruning

		From Paper: Scalable Top-N Local Outlier Detection (In Submission)

