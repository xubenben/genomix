package edu.uci.ics.genomix.pregelix.operator.scaffolding2;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import org.apache.hadoop.io.WritableUtils;

import edu.uci.ics.genomix.data.types.EDGETYPE;
import edu.uci.ics.genomix.data.types.Kmer;
import edu.uci.ics.genomix.data.types.ReadHeadInfo;
import edu.uci.ics.genomix.data.types.VKmer;
import edu.uci.ics.genomix.pregelix.base.VertexValueWritable;
public class RayValue extends VertexValueWritable {
    private static final long serialVersionUID = 1L;
    
    Boolean flippedFromInitialDirection = null;
    //boolean visited = false;
    List<VKmer> visitedList;
    boolean intersection = false;
    boolean stopSearch = false;
    /**
     * This BitSet is supposed to keep track of those neighbors, each node has chosen at least on one of the passing walk.
     */
    BitSet visitedNeighbors ;

    Integer pendingCandidateBranches = null;
    ArrayList<RayMessage> candidateMsgs = null;

    protected static class FIELDS {
        public static final byte DIR_FLIPPED_VS_INITIAL = 0b01;
        public static final byte DIR_SAME_VS_INITIAL = 0b10;
        public static final byte VISITED_LIST = 0b1 << 2;
        public static final byte INTERSECTION = 0b1 << 3;
        public static final byte STOP_SEARCH = 0b1 << 4;
        public static final byte PENDING_CANDIDATE_BRANCHES = 1 << 5;
        public static final byte CANDIDATE_MSGS = 1 << 6;
        public static final short VISITED_NEIGHBORS = 1 << 7;
    }

    @Override
    public void readFields(DataInput in) throws IOException {
        super.readFields(in);
        if (((state) & FIELDS.DIR_FLIPPED_VS_INITIAL) != 0) {
            flippedFromInitialDirection = true;
        } else if (((state) & FIELDS.DIR_SAME_VS_INITIAL) != 0) {
            flippedFromInitialDirection = false;
        } else {
            flippedFromInitialDirection = null;
        }
        //visited = ((state & FIELDS.VISITED) != 0);        
        intersection = ((state & FIELDS.INTERSECTION) != 0);
        stopSearch = ((state & FIELDS.STOP_SEARCH) != 0);
        
        if ((state & FIELDS.VISITED_LIST) != 0) {
            getVisitedList().clear();
            int count = in.readInt();
            for (int i = 0; i < count; i++) {
                VKmer m = new VKmer();
                m.readFields(in);
                visitedList.add(m);
            }
        }
        
        if ((state & FIELDS.PENDING_CANDIDATE_BRANCHES) != 0) {
            pendingCandidateBranches = in.readInt();
        }
        if ((state & FIELDS.CANDIDATE_MSGS) != 0) {
            getCandidateMsgs().clear();
            int count = in.readInt();
            for (int i = 0; i < count; i++) {
                RayMessage m = new RayMessage();
                m.readFields(in);
                candidateMsgs.add(m);
            }
        }
        if ((state & FIELDS.VISITED_NEIGHBORS) != 0) {
        	readBitSet(in, 8, visitedNeighbors);
        }
        
    }

    @Override
    public void write(DataOutput out) throws IOException {
        state = 0;
        if (flippedFromInitialDirection != null) {
            state |= flippedFromInitialDirection ? FIELDS.DIR_FLIPPED_VS_INITIAL : FIELDS.DIR_SAME_VS_INITIAL;
        }
        if (intersection) {
            state |= FIELDS.INTERSECTION;
        }
        if (stopSearch) {
            state |= FIELDS.STOP_SEARCH;
        }
        
        if (visitedList != null && visitedList.size() > 0) {
            state |= FIELDS.VISITED_LIST;
        } 
        
        if (pendingCandidateBranches != null) {
            state |= FIELDS.PENDING_CANDIDATE_BRANCHES;
        }
        if (candidateMsgs != null && candidateMsgs.size() > 0) {
            state |= FIELDS.CANDIDATE_MSGS;
        }
        super.write(out);
        
        if (visitedList != null && visitedList.size() > 0) {
            out.writeInt(visitedList.size());
            for (VKmer m : visitedList) {
                m.write(out);
            }
        }
        
        if (pendingCandidateBranches != null) {
            out.writeInt(pendingCandidateBranches);
        }
        if (candidateMsgs != null && candidateMsgs.size() > 0) {
            out.writeInt(candidateMsgs.size());
            for (RayMessage m : candidateMsgs) {
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
        flippedFromInitialDirection = null;
        intersection = false;
        stopSearch = false;
        visitedList = null;
        pendingCandidateBranches = null;
        candidateMsgs = null;
        visitedNeighbors = null;
    }

    /**
     * @return whether or not I have any readids that **could** contribute to the current walk
     */
    public boolean isOutOfRange(int myOffset, int walkLength, int maxDist) {
        int numBasesToSkip = Math.max(0, walkLength - maxDist - myOffset);
        int myLength = getKmerLength() - Kmer.getKmerLength() + 1;
        if (!flippedFromInitialDirection) {
            // TODO fix max offset to be distance
            // cut off the beginning
            if (numBasesToSkip > myLength) {
                // start > offsets I contain-- no valid readids
                return true;
            }
            return getUnflippedReadIds().getOffSetRange(numBasesToSkip, ReadHeadInfo.MAX_OFFSET_VALUE).isEmpty();
        } else {
            // cut off the end 
            if (myLength - numBasesToSkip < 0) {
                // my max is negative-- no valid readids
                return true;
            }
            //FIXME
            return getFlippedReadIds().getOffSetRange(0, Math.max(0, getKmerLength() - numBasesToSkip)).isEmpty();
            //return getFlippedReadIds().getOffSetRange(0, myLength - numBasesToSkip).isEmpty();
        }
    }

    public ArrayList<RayMessage> getCandidateMsgs() {
        if (candidateMsgs == null) {
            candidateMsgs = new ArrayList<>();
        }
        return candidateMsgs;
    }

    public void setCandidateMsgs(ArrayList<RayMessage> candidateMsgs) {
        this.candidateMsgs = candidateMsgs;
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
    
	public void setNeighbors(EDGETYPE edge, VKmer neighborId, int kmerSize){
		int index = -1;
		String nId = neighborId.toString();
		String reversenId = neighborId.reverse().toString();
		char nextNuc = nId.charAt(kmerSize -1) ;
		char reverseNextNuc = reversenId.charAt(kmerSize - 1);
		if(edge.equals(EDGETYPE.FF)){
			switch(index){
				case 0:  nextNuc = 'A';
					break;
				case 1:  nextNuc = 'T';
					break;
				case 2:  nextNuc = 'C';
            		break;
				case 3:  nextNuc = 'G';
					break;
			}
			visitedNeighbors.set(index);
		}
		
		if(edge.equals(EDGETYPE.FR)){
			switch(index){
				case 0:  reverseNextNuc = 'A';
					break;
				case 1:  reverseNextNuc = 'T';
					break;
				case 2:  reverseNextNuc = 'C';
            		break;
				case 3:  reverseNextNuc = 'G';
					break;
			}
			visitedNeighbors.set(index);
		}
		
		if(edge.equals(EDGETYPE.RF)){
			switch(index){
				case 4:  nextNuc = 'A';
					break;
				case 5:  nextNuc = 'T';
					break;
				case 6:  nextNuc = 'C';
            		break;
				case 7:  nextNuc = 'G';
					break;
			}
			visitedNeighbors.set(index);
		}
		
		if(edge.equals(EDGETYPE.RR)){
			switch(index){
				case 4:  reverseNextNuc = 'A';
					break;
				case 5:  reverseNextNuc = 'T';
					break;
				case 6:  reverseNextNuc = 'C';
            		break;
				case 7:  reverseNextNuc = 'G';
					break;
			}
			visitedNeighbors.set(index);
		}
		
	}
}
