package edu.uci.ics.genomix.pregelix.checker;

import java.util.Iterator;

import edu.uci.ics.genomix.pregelix.io.VertexValueWritable;
import edu.uci.ics.genomix.pregelix.io.VertexValueWritable.State;
import edu.uci.ics.genomix.pregelix.io.message.MessageWritable;
import edu.uci.ics.genomix.pregelix.operator.BasicGraphCleanVertex;
import edu.uci.ics.genomix.pregelix.operator.aggregator.StatisticsAggregator;
import edu.uci.ics.genomix.pregelix.type.MessageFlag;
import edu.uci.ics.genomix.type.VKmerBytesWritable;

public class SymmetryCheckerVertex extends
    BasicGraphCleanVertex<VertexValueWritable, MessageWritable> {

    @Override
    public void initVertex() {
        super.initVertex();
        if(incomingMsg == null)
            incomingMsg = new MessageWritable();
        if(outgoingMsg == null)
            outgoingMsg = new MessageWritable();
        else
            outgoingMsg.reset();
        if(destVertexId == null)
            destVertexId = new VKmerBytesWritable();
        if(getSuperstep() == 1)
            StatisticsAggregator.preGlobalCounters.clear();
//        else
//            StatisticsAggregator.preGlobalCounters = BasicGraphCleanVertex.readStatisticsCounterResult(getContext().getConfiguration());
        counters.clear();
        getVertexValue().getCounters().clear();
        outFlag = 0;
    }
    
    @Override
    public void compute(Iterator<MessageWritable> msgIterator) throws Exception {
        initVertex();
        if(getSuperstep() == 1){
            sendSettledMsgToAllNeighborNodes(getVertexValue());
            voteToHalt();
        } else if(getSuperstep() == 2){
            //check if the corresponding edge exists
            while(msgIterator.hasNext()){
                incomingMsg = msgIterator.next();
                byte meToNeighborDir = (byte) (incomingMsg.getFlag() & MessageFlag.DIR_MASK);
                byte neighborToMeDir = mirrorDirection(meToNeighborDir);
                boolean exist = getVertexValue().getEdgeList(neighborToMeDir).contains(incomingMsg.getSourceVertexId());
                if(!exist){
                    getVertexValue().setState(State.IS_ERROR);
                }
            }
            voteToHalt();
        }
    }

}