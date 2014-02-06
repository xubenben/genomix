package edu.uci.ics.genomix.pregelix.operator.prune;

import java.util.Iterator;

import edu.uci.ics.genomix.data.types.DIR;
import edu.uci.ics.genomix.data.types.EDGETYPE;
import edu.uci.ics.genomix.data.types.VKmer;
import edu.uci.ics.genomix.pregelix.base.DeBruijnGraphCleanVertex;

public class PruneVertex extends DeBruijnGraphCleanVertex<PruneValue, PruneMessage> {
	
	
	public void sendPruneGraphMsg(){		
		PruneValue vertex = getVertexValue();
		PruneMessage msg = new PruneMessage();
		DIR dir = DIR.FORWARD;
		if(vertex.getVisitedList().size() > 0){
			for (EDGETYPE et : dir.edgeTypes()) {
                for (VKmer next : vertex.getEdges(et)) {
                	if (vertex.isPrunable(et, vertex.neighborsLastChar(et, next, kmerSize))){
                		 msg.setToPruneEdge(et.mirror());
                		 msg.setToPruneId(getVertexId());
                		 sendMsg(next, msg);
                		 vertex.getEdges(et).remove(next);
                	}
                   
                }
            }
		}
	}
	
	public void recievePruneGraphMsg(PruneMessage msg){
		PruneValue vertex = new PruneValue();
		if(vertex.getVisitedList().size() > 0){
             vertex.getEdges(msg.getToPruneEdge()).remove(msg.getToPruneId());
 
        }
	}
	
	
	@Override
	public void compute(Iterator<PruneMessage> msgIterator) throws Exception {
		// TODO Auto-generated method stub
		
	}

}
