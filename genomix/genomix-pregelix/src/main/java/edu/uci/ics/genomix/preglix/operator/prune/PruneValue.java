package edu.uci.ics.genomix.preglix.operator.prune;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import org.apache.hadoop.io.WritableUtils;

import edu.uci.ics.genomix.data.types.EDGETYPE;
import edu.uci.ics.genomix.data.types.VKmer;
import edu.uci.ics.genomix.pregelix.base.VertexValueWritable;

public class PruneValue extends VertexValueWritable{
	private static final long serialVersionUID = 1L;
	
    List<VKmer> visitedList;
    /**
     * This BitSet is supposed to keep track of those neighbors, each node has chosen at least on one of the passing walk.
     */
    BitSet visitedNeighbors ;
    
    protected static class FIELDS {
        public static final byte VISITED_LIST = 0b01;
        public static final byte VISITED_NEIGHBORS = 0b10;
    }
    
    @Override
    public void readFields(DataInput in) throws IOException {
        super.readFields(in);
        
        if ((state & FIELDS.VISITED_LIST) != 0) {
            getVisitedList().clear();
            int count = in.readInt();
            for (int i = 0; i < count; i++) {
                VKmer m = new VKmer();
                m.readFields(in);
                visitedList.add(m);
            }
        }
        
        if ((state & FIELDS.VISITED_NEIGHBORS) != 0) {
        	readBitSet(in, 8, visitedNeighbors);
        }
        
    }
    
    @Override
    public void write(DataOutput out) throws IOException {
        state = 0;
        super.write(out);
        if (visitedList != null && visitedList.size() > 0) {
            out.writeInt(visitedList.size());
            for (VKmer m : visitedList) {
                m.write(out);
            }
        }
        if (visitedNeighbors != null){
        	writeBitSet(out, 8 , visitedNeighbors);
        }
    }
    
    @Override
    public void reset() {
        super.reset();
        visitedList = null;
        visitedNeighbors = null;
    }
    
    public List<VKmer> getVisitedList() {
        if (visitedList == null) {
            visitedList = new ArrayList<>();
        }
        return visitedList;
    }

    public void setVisitedList(List<VKmer> visitedList) {
        this.visitedList = visitedList;
    }
    
    private static final void readBitSet(DataInput stream, int nbits, 
    		BitSet bitSet) throws IOException {
    	bitSet.clear();
    	long initialBits = WritableUtils.readVLong(stream);
    	long last = 0L;
    	while (0L != initialBits) {
    		last = Long.lowestOneBit(initialBits);
    		initialBits ^= last;
    		bitSet.set(Long.numberOfTrailingZeros(last));
    	}
    		    
    	for (int offset=Long.SIZE; offset < nbits; offset+=Byte.SIZE) {
    		byte bits = stream.readByte();
    		while (0 != bits) {
    			last = Long.lowestOneBit(bits);
    			bits ^= last;
    			bitSet.set(Long.numberOfTrailingZeros(last) + offset);
    		}
    	}
    }
    
    private static final void writeBitSet(DataOutput stream, int nbits,
    		BitSet bitSet) throws IOException {
    	long bits = 0L;
    
    	int bitSetIndex = bitSet.nextSetBit(0);
    	for (;bitSetIndex >= 0 && bitSetIndex < Long.SIZE;
    			bitSetIndex=bitSet.nextSetBit(bitSetIndex+1)) {
    		bits |= 1L << bitSetIndex;
    	}
    	WritableUtils.writeVLong(stream,bits);

    	if (nbits > Long.SIZE) {
    		bits = 0L;
    		for (int lastWordWritten = 0; bitSetIndex >= 0 && bitSetIndex < nbits; 
    				bitSetIndex = bitSet.nextSetBit(bitSetIndex+1)) {
    			int bitsIndex = bitSetIndex % Byte.SIZE;
    			int word = (bitSetIndex-Long.SIZE) / Byte.SIZE;
    			if (word > lastWordWritten) {
    				stream.writeByte((byte)bits);
    				bits = 0L;
    				for (lastWordWritten++;lastWordWritten<word;lastWordWritten++) {
    					stream.writeByte((byte)bits);
    				}
    			}
    			bits |= 1L << bitsIndex;
    		}
    		stream.writeByte((byte)bits);
    	}
    }
    
    public BitSet getVisitedNeighbors() {
        if (visitedNeighbors == null) {
            visitedNeighbors = new BitSet();
        }
        return visitedNeighbors;
    }

    public void setVisitedNeighbors(BitSet visitedNeighbors) {
        this.visitedNeighbors = visitedNeighbors;
    }
    
	public char neighborsLastChar(EDGETYPE edge, VKmer neighbor, int kmerSize){
		if(edge.equals(EDGETYPE.FF) || edge.equals(EDGETYPE.RF)){
			return neighbor.toString().charAt(kmerSize -1);
		}
		
		return neighbor.reverse().toString().charAt(kmerSize -1);
		
	}
	
	public boolean isPrunable(EDGETYPE edge , char last){
		int index = -1;
		if(edge.equals(EDGETYPE.FF) || edge.equals(EDGETYPE.FR)){
			switch(last){
			case 'A':  index = 0;
				break;
			case 'T':  index = 1;
				break;
			case 'C':  index = 2;
        		break;
			case 'G':  index = 3;
				break;
		}
			if (visitedNeighbors.get(index)){
				return false;
			}else {
				return true;
			}
		}else {
			switch(last){
			case 'A':  index = 4;
				break;
			case 'T':  index = 5;
				break;
			case 'C':  index = 6;
        		break;
			case 'G':  index = 7;
				break;
		}
			if (visitedNeighbors.get(index)){
				return false;
			}else {
				return true;
			}
		}
		
	}
}
