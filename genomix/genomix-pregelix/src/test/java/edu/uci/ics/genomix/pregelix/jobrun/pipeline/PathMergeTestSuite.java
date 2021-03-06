package edu.uci.ics.genomix.pregelix.jobrun.pipeline;

import junit.framework.Test;
import edu.uci.ics.genomix.pregelix.jobrun.BasicGraphCleanTestSuite;

public class PathMergeTestSuite extends BasicGraphCleanTestSuite {
    
    public static Test suite() throws Exception {
        String pattern = "MERGE";
        String testSet[] = {"SimpleBubble"}; 
//                "2", "3", "4", "5", "6", "7", "8", "9", "Head_5", "Head_6", "Head_7", "Head_8", "Head_9",
//                "Head_10", "LeftAdj", "RightAdj", "SimplePath", "ThreeDuplicate", "SimpleBridgePath", "SimpleTreePath",
//                "RingPath", "CyclePath", "SelfTandemRepeat", "TandemRepeatWithMergeEdge",
//                "TandemRepeatWithUnmergeEdge", "ComplexTandemRepeat", "MultiTandemRepeat", "MultiTandemRepeat2",
//                "MultiTandemRepeat3", "TandemRepeatWithThreeNodes", "TandemRepeatAndCycle", "ThreeNodesCycle",
//                "SameGeneCodeWithEdge", "SameGeneCodeWithoutEdge", "synthetic" };
        init(pattern, testSet);
        BasicGraphCleanTestSuite testSuite = new BasicGraphCleanTestSuite();
        return makeTestSuite(testSuite);
    }
}
