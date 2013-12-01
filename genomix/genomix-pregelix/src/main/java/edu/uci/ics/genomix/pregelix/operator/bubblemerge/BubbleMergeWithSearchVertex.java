package edu.uci.ics.genomix.pregelix.operator.bubblemerge;

import java.util.Iterator;
import java.util.logging.Logger;

import edu.uci.ics.genomix.config.GenomixJobConf;
import edu.uci.ics.genomix.pregelix.client.Client;
import edu.uci.ics.genomix.pregelix.io.BubbleMergeWithSearchVertexValueWritable;
import edu.uci.ics.genomix.pregelix.io.VertexValueWritable.State;
import edu.uci.ics.genomix.pregelix.io.common.ArrayListWritable;
import edu.uci.ics.genomix.pregelix.io.message.BubbleMergeWithSearchMessage;
import edu.uci.ics.genomix.pregelix.operator.DeBruijnGraphCleanVertex;
import edu.uci.ics.genomix.pregelix.operator.aggregator.StatisticsAggregator;
import edu.uci.ics.genomix.type.EDGETYPE;
import edu.uci.ics.genomix.type.VKmer;
import edu.uci.ics.genomix.type.VKmerList;

/**
 * Graph clean pattern: Simple Bubble Merge
 */
public class BubbleMergeWithSearchVertex extends
        DeBruijnGraphCleanVertex<BubbleMergeWithSearchVertexValueWritable, BubbleMergeWithSearchMessage> {

    private static final Logger LOG = Logger.getLogger(BubbleMergeWithSearchVertex.class.getName());

    private float MAX_BFS_LENGTH = -1;

    /**
     * initiate kmerSize, maxIteration
     */
    @Override
    public void initVertex() {
        super.initVertex();
        if (MAX_BFS_LENGTH < 0)
            MAX_BFS_LENGTH = Float.parseFloat(getContext().getConfiguration().get(
                    GenomixJobConf.BUBBLE_MERGE_MAX_DISSIMILARITY));
        if (outgoingMsg == null)
            outgoingMsg = new BubbleMergeWithSearchMessage();
        StatisticsAggregator.preGlobalCounters.clear();
        //        else
        //            StatisticsAggregator.preGlobalCounters = BasicGraphCleanVertex.readStatisticsCounterResult(getContext().getConfiguration());
        counters.clear();
        getVertexValue().getCounters().clear();
    }

    @Override
    public void compute(Iterator<BubbleMergeWithSearchMessage> msgIterator) {
        initVertex();
        if (getSuperstep() == 1) {
            /** begin BFS in source vertices **/
            BubbleMergeWithSearchVertexValueWritable vertex = new BubbleMergeWithSearchVertexValueWritable();
            int internalKmerLength = vertex.getInternalKmer().getKmerLetterLength();
        } else if (getSuperstep() == 2) {
            /** continue BFS **/
            BubbleMergeWithSearchVertexValueWritable vertex = new BubbleMergeWithSearchVertexValueWritable();
            while (msgIterator.hasNext()) {
                BubbleMergeWithSearchMessage incomingMsg = msgIterator.next();

                // get msg flag
                byte flag = (byte) (vertex.getState() & State.BUBBLE_WITH_SEARCH_FLAG_MASK);

                if (flag == State.UPDATE_PATH_IN_NEXT) {
                    int internalKmerLength = vertex.getInternalKmer().getKmerLetterLength();
                    VKmer source = incomingMsg.getPathList().getPosition(0);
                    if (internalKmerLength + incomingMsg.getPreKmerLength() > MAX_BFS_LENGTH) {
                        // send back to source vertex (pathList and internalKmer)
                        outgoingMsg.reset();
                        outgoingMsg.setPathList(incomingMsg.getPathList());
                        outgoingMsg.setInternalKmer(incomingMsg.getInternalKmer());
                        outgoingMsg.setEdgeTypeList(incomingMsg.getEdgeTypeList());
                        sendMsg(source, outgoingMsg);
                    } else {
                        // if numBranches > 1, update numBranches
                        if (vertex.outDegree() > 1) {
                            outgoingMsg.reset();
                            outgoingMsg.setFlag(State.UPDATE_BRANCH_IN_SRC);
                            outgoingMsg.setNumBranches(vertex.outDegree());
                            sendMsg(source, outgoingMsg);
                        }

                        // send to next (pathList and internalKmer)  
                        for (EDGETYPE et : EDGETYPE.OUTGOING) {
                            for (VKmer dest : vertex.getEdgeMap(et).keySet()) {
                                outgoingMsg.reset();

                                // set flag and source vertex
                                outFlag &= EDGETYPE.CLEAR;
                                outFlag |= et.mirror().get();
                                outgoingMsg.setFlag(outFlag);
                                outgoingMsg.setSourceVertexId(getVertexId());

                                // update pathList
                                VKmerList pathList = incomingMsg.getPathList();
                                pathList.append(getVertexId());
                                outgoingMsg.setPathList(pathList);

                                // update internalKmer
                                VKmer internalKmer = incomingMsg.getInternalKmer();
                                internalKmer.mergeWithKmerInDir(et, Integer.parseInt(GenomixJobConf.KMER_LENGTH),
                                        vertex.getInternalKmer());
                                outgoingMsg.setInternalKmer(internalKmer);

                                // update edgeTypes
                                ArrayListWritable<EDGETYPE> edgeTypes = incomingMsg.getEdgeTypeList();
                                edgeTypes.add(et);
                                outgoingMsg.setEdgeTypeList(edgeTypes);

                                sendMsg(dest, outgoingMsg);
                            }
                        }
                    }
                } else if (flag == State.UPDATE_BRANCH_IN_SRC) {
                    // update numBranches in src
                    vertex.setNumBranches(vertex.getNumBranches() + incomingMsg.getNumBranches() - 1);
                } else if (flag == State.END_NOTICE_IN_SRC) {
                    // update pathList
                    vertex.getArrayOfPathList().add(incomingMsg.getPathList());

                    // update internalKmer
                    vertex.getArrayOfInternalKmer().add(incomingMsg.getInternalKmer());

                    // update numBranches
                    int numBranches = vertex.getNumBranches();
                    numBranches--;
                    if (numBranches == 0) {
                        /* process in src */
                        // step1: figure out which path to keep
                        int k = 2;
                        VKmerList pathList = vertex.getArrayOfPathList().get(k);
                        VKmer mergeKmer = vertex.getArrayOfInternalKmer().get(k);
                        ArrayListWritable<EDGETYPE> edgeTypes = vertex.getArrayOfEdgeTypes().get(k);

                        // step2: replace internalKmer with mergeKmer and clear edges towards path nodes
                        vertex.setInternalKmer(mergeKmer);
                        vertex.getEdgeMap(vertex.getArrayOfEdgeTypes().get(k).get(0)).remove(pathList.getPosition(1));

                        // step3: send kill message to path nodes
                        for (int i = 0; i < pathList.size(); i++) {
                            VKmer dest = pathList.getPosition(i);
                            if (dest.equals(getVertexId()))
                                continue;
                            outgoingMsg.reset();
                            outgoingMsg.setFlag(State.KILL_MESSAGE_FROM_SOURCE);
                            outgoingMsg.setSourceVertexId(pathList.getPosition(i - 1));
                            if (i + 1 < pathList.size())
                                outgoingMsg.setInternalKmer(pathList.getPosition(i + 1)); // use internalKmer field to store next node 

                            // store edgeType in msg.ArrayListWritable<EDGETYPE>(0)
                            outgoingMsg.getEdgeTypeList().add(edgeTypes.get(i - 1));
                            outgoingMsg.getEdgeTypeList().add(edgeTypes.get(i));
                        }
                    }
                } else if (flag == State.KILL_MESSAGE_FROM_SOURCE) {
                    // send msg to delete edges except the path nodes
                    for (EDGETYPE et : EDGETYPE.values()) {
                        for (VKmer dest : vertex.getEdgeMap(et).keySet()) {
                            if ((et == incomingMsg.getEdgeTypeList().get(0) && dest.equals(incomingMsg
                                    .getSourceVertexId()))
                                    || (et == incomingMsg.getEdgeTypeList().get(1) && dest.equals(incomingMsg
                                            .getInternalKmer())))
                                continue;
                            outgoingMsg.reset();
                            outFlag |= State.PRUNE_DEAD_EDGE;
                            outFlag |= et.mirror().get();
                            outgoingMsg.setFlag(outFlag);
                            outgoingMsg.setSourceVertexId(getVertexId());
                        }
                    }
                } else if (flag == State.PRUNE_DEAD_EDGE) {
                    EDGETYPE meToNeighborEdgetype = EDGETYPE.fromByte(incomingMsg.getFlag());
                    vertex.getEdgeMap(meToNeighborEdgetype).remove(incomingMsg.getSourceVertexId());
                }
            }
            voteToHalt();
        } else if (getSuperstep() == 3) {
        } else if (getSuperstep() == 4) {
        }
        voteToHalt();
    }

    public static void main(String[] args) throws Exception {
        Client.run(args, getConfiguredJob(null, BubbleMergeWithSearchVertex.class));
    }

}
