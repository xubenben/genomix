package edu.uci.ics.genomix.pregelix.operator.prune;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

import edu.uci.ics.genomix.data.types.EDGETYPE;
import edu.uci.ics.genomix.data.types.VKmer;
import edu.uci.ics.genomix.pregelix.base.MessageWritable;
public class PruneMessage extends MessageWritable {
	
    private VKmer toPruneId = null; 
    private EDGETYPE toPruneEdge = null; 
    
    protected class FIELDS extends MESSAGE_FIELDS {
        public static final byte EDGE_TO_PRUNE = 1 << 1;
    }
    
    public void readFields(DataInput in) throws IOException {
        super.readFields(in);
        if ((messageFields & FIELDS.EDGE_TO_PRUNE) != 0) {
        	getToPruneId().readFields(in);
        	toPruneEdge = EDGETYPE.fromByte(in.readByte());
        }

    }
    
    @Override
    public void write(DataOutput out) throws IOException {
        super.write(out);
        if (toPruneId != null && toPruneId.getKmerLetterLength() > 0) {
        	toPruneId.write(out);
        	out.writeByte(toPruneEdge.get());
        }
    }

    
    @Override
    protected byte getActiveMessageFields() {
        byte fields = super.getActiveMessageFields();
        if (toPruneId != null && toPruneId.getKmerLetterLength() > 0) {
            fields |= FIELDS.EDGE_TO_PRUNE;
        }
        return fields;
    }
    
    
    public VKmer getToPruneId() {
        if (toPruneId == null)
        	toPruneId = new VKmer();
        return toPruneId;
    }

    public void setToPruneId(VKmer toPruneId) {
        this.toPruneId = toPruneId;
    }

    public EDGETYPE getToPruneEdge() {
        return this.toPruneEdge;
    }

    public void setToPruneEdge(EDGETYPE toPruneEdge) {
        this.toPruneEdge = toPruneEdge;
    }
    
    
    @Override
    public void reset() {
        super.reset();
        toPruneId = null;
        toPruneEdge = null;

    }
}
