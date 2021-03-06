package edu.uci.ics.genomix.pregelix.operator.bridgeremove;

import java.io.IOException;
import java.util.Iterator;

import edu.uci.ics.genomix.data.config.GenomixJobConf;
import edu.uci.ics.genomix.data.types.DIR;
import edu.uci.ics.genomix.data.types.EDGETYPE;
import edu.uci.ics.genomix.data.types.VKmer;
import edu.uci.ics.genomix.pregelix.base.DeBruijnGraphCleanVertex;
import edu.uci.ics.genomix.pregelix.base.MessageWritable;
import edu.uci.ics.genomix.pregelix.base.VertexValueWritable;
import edu.uci.ics.genomix.pregelix.types.GraphMutations;

/**
 * Graph clean pattern: Remove Bridge
 * 
 */
public class BridgeRemoveVertex extends DeBruijnGraphCleanVertex<VertexValueWritable, MessageWritable> {

    private int MIN_LENGTH_TO_KEEP = -1;

    /**
     * initiate kmerSize, maxIteration
     */
    @Override
    public void initVertex() {
        super.initVertex();
        if (MIN_LENGTH_TO_KEEP == -1)
            MIN_LENGTH_TO_KEEP = Integer.parseInt(getContext().getConfiguration().get(
                    GenomixJobConf.BRIDGE_REMOVE_MAX_LENGTH));
        if (outgoingMsg == null)
            outgoingMsg = new MessageWritable();
    }

    /**
     * step 1: detect neighbor of bridge vertex
     */
    public void detectBridgeNeighbor() {
        //detect neighbor of bridge vertex
        VertexValueWritable vertex = getVertexValue();
        if (vertex.getDegree() == 3) {
            for (DIR d : DIR.values()) {
                //only 1 incoming and 2 outgoing || 2 incoming and 1 outgoing are valid 
                if (vertex.degree(d) == 2) {
                    for (EDGETYPE et : d.edgeTypes()) {
                        for (VKmer dest : vertex.getEdges(et)) {
                            sendMsg(dest, outgoingMsg);
                        }
                    }
                }
            }
        }
        voteToHalt();
    }

    /**
     * step2: remove bridge pattern
     */
    public void removeBridge(Iterator<MessageWritable> msgIterator) {
        VertexValueWritable vertex = getVertexValue();
        //only the vertex which has and only has 2 degree can be bridge vertex
        if (vertex.getDegree() == 2) {
            //count #receivedMsg
            int count = 0;
            while (msgIterator.hasNext()) {
                msgIterator.next();
                if (count == 3)
                    break;
                count++;
            }
            //remove bridge
            if (count == 2) { //I'm bridge vertex
                if (vertex.getKmerLength() < MIN_LENGTH_TO_KEEP) {
                    broadcastKillself();
                    deleteVertex(getVertexId());
                    getCounters().findCounter(GraphMutations.Num_RemovedBridges).increment(1);
                }
            }
        }
    }

    @Override
    public void compute(Iterator<MessageWritable> msgIterator) throws IOException {
        initVertex();
        if (getSuperstep() == 1) {
            detectBridgeNeighbor();
        } else if (getSuperstep() == 2) {
            removeBridge(msgIterator);
        } else if (getSuperstep() == 3) {
            pruneDeadEdges(msgIterator);
        }
        voteToHalt();
    }

}
