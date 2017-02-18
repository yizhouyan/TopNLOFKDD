package net.sf.jsi;

import gnu.trove.procedure.TIntProcedure;
import net.sf.jsi.rtree.RTree;

public class TestRTree {
	public static void main(String[] args) {
		// Create and initialize an rtree
		SpatialIndex si = new RTree();
		si.init(null);

		final Rectangle[] rects = new Rectangle[100];
		rects[0] = new Rectangle(0, 10, 0, 10);
		rects[1] = new Rectangle(0, 11, 1, 20);
		rects[2] = new Rectangle(30, 40, 30, 50);
		rects[3] = new Rectangle(60, 70, 80, 90);
		si.add(rects[0], 0);
		si.add(rects[1], 1);
		si.add(rects[2], 2);
		si.add(rects[3], 3);
//		si.nearestN(new Point(36.3f, 84.3f), // the point for which we want to
//												// find nearby rectangles
//				new TIntProcedure() { // a procedure whose execute() method will
//										// be called with the results
//					public boolean execute(int i) {
//						System.out.println("Rectangle " + i + " " + rects[i]);
//						// log.info("Rectangle " + i + " " + rects[i] + ",
//						// distance=" + rects[i].distance(p));
//						return true; // return true here to continue receiving
//										// results
//					}
//				}, 1, // the number of nearby rectangles to find
//				Float.MAX_VALUE // Don't bother searching further than this.
//								// MAX_VALUE means search everything
//		);
		si.intersects(new Rectangle(20, 20, 100, 100), new TIntProcedure() {
			// be called with the results
			public boolean execute(int i) {
				System.out.println("Rectangle " + i + " " + rects[i]);
				// log.info("Rectangle " + i + " " + rects[i] + ",
				// distance=" + rects[i].distance(p));
				return true; // return true here to continue receiving
								// results
			}
		});

	}
}
