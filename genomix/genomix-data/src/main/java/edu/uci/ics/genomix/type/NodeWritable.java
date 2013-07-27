package edu.uci.ics.genomix.type;

import java.io.ByteArrayOutputStream;
import java.io.DataInput;
import java.io.DataOutput;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.Serializable;

import org.apache.hadoop.io.DataOutputBuffer;
import org.apache.hadoop.io.WritableComparable;


public class NodeWritable implements WritableComparable<NodeWritable>, Serializable{

    private static final long serialVersionUID = 1L;
    public static final NodeWritable EMPTY_NODE = new NodeWritable();
    
    private PositionListWritable nodeIdList;
    private KmerListWritable forwardForwardList;
    private KmerListWritable forwardReverseList;
    private KmerListWritable reverseForwardList;
    private KmerListWritable reverseReverseList;
    private VKmerBytesWritable kmer;
    
    // merge/update directions
    public static class DirectionFlag {
        public static final byte DIR_FF = 0b00 << 0;
        public static final byte DIR_FR = 0b01 << 0;
        public static final byte DIR_RF = 0b10 << 0;
        public static final byte DIR_RR = 0b11 << 0;
        public static final byte DIR_MASK = 0b11 << 0;
    }
    
    public NodeWritable() {
        nodeIdList = new PositionListWritable();
        forwardForwardList = new KmerListWritable();
        forwardReverseList = new KmerListWritable();
        reverseForwardList = new KmerListWritable();
        reverseReverseList = new KmerListWritable();
        kmer = new VKmerBytesWritable();  // in graph construction - not set kmerlength Optimization: VKmer
    }
    
    public NodeWritable(PositionListWritable nodeIdList, KmerListWritable FFList, KmerListWritable FRList,
            KmerListWritable RFList, KmerListWritable RRList, VKmerBytesWritable kmer) {
        this();
        set(nodeIdList, FFList, FRList, RFList, RRList, kmer);
    }
    
    public void set(NodeWritable node){
        set(node.nodeIdList, node.forwardForwardList, node.forwardReverseList, node.reverseForwardList, 
                node.reverseReverseList, node.kmer);
    }
    
    public void set(PositionListWritable nodeIdList, KmerListWritable FFList, KmerListWritable FRList,
            KmerListWritable RFList, KmerListWritable RRList, VKmerBytesWritable kmer2) {
        this.nodeIdList.set(nodeIdList);
        this.forwardForwardList.setCopy(FFList);
        this.forwardReverseList.setCopy(FRList);
        this.reverseForwardList.setCopy(RFList);
        this.reverseReverseList.setCopy(RRList);
        this.kmer.setAsCopy(kmer2);
    }

    public void reset() {
        this.nodeIdList.reset();
        this.forwardForwardList.reset();
        this.forwardReverseList.reset();
        this.reverseForwardList.reset();
        this.reverseReverseList.reset();
        this.kmer.reset(0);
    }
    
    
    public PositionListWritable getNodeIdList() {
        return nodeIdList;
    }

    public void setNodeIdList(PositionListWritable nodeIdList) {
        this.nodeIdList.set(nodeIdList);
    }

    public VKmerBytesWritable getKmer() {
        return kmer;
    }

    public void setKmer(VKmerBytesWritable kmer) {
        this.kmer.setAsCopy(kmer);
    }
    
    public int getKmerLength() {
        return kmer.getKmerLetterLength();
    }
    
    public KmerListWritable getFFList() {
        return forwardForwardList;
    }

    public KmerListWritable getFRList() {
        return forwardReverseList;
    }

    public KmerListWritable getRFList() {
        return reverseForwardList;
    }

    public KmerListWritable getRRList() {
        return reverseReverseList;
    }
    
	public void setFFList(KmerListWritable forwardForwardList) {
		this.forwardForwardList.setCopy(forwardForwardList);
	}

	public void setFRList(KmerListWritable forwardReverseList) {
		this.forwardReverseList.setCopy(forwardReverseList);
	}

	public void setRFList(KmerListWritable reverseForwardList) {
		this.reverseForwardList.setCopy(reverseForwardList);
	}

	public void setRRList(KmerListWritable reverseReverseList) {
		this.reverseReverseList.setCopy(reverseReverseList);
	}

	public KmerListWritable getListFromDir(byte dir) {
        switch (dir & DirectionFlag.DIR_MASK) {
            case DirectionFlag.DIR_FF:
                return getFFList();
            case DirectionFlag.DIR_FR:
                return getFRList();
            case DirectionFlag.DIR_RF:
                return getRFList();
            case DirectionFlag.DIR_RR:
                return getRRList();
            default:
                throw new RuntimeException("Unrecognized direction in getListFromDir: " + dir);
        }
    }
	
	/**
	 * Returns the length of the byte-array version of this node
	 */
	public int getSerializedLength() {
	    return nodeIdList.getLength() + forwardForwardList.getLength() + forwardReverseList.getLength() + 
	            reverseForwardList.getLength() + reverseReverseList.getLength() + kmer.getLength();
	}
	
	/**
     * Return this Node's representation as a new byte array 
     */
    public byte[] marshalToByteArray() throws IOException {
        ByteArrayOutputStream baos = new ByteArrayOutputStream(getSerializedLength());
        DataOutputStream out = new DataOutputStream(baos);
        write(out);
        return baos.toByteArray();
    }
    
    public void setAsCopy(byte[] data, int offset) {
        int curOffset = offset;
        nodeIdList.set(data, curOffset);
        
        curOffset += nodeIdList.getLength();
        forwardForwardList.setCopy(data, curOffset);
        curOffset += forwardForwardList.getLength();
        forwardReverseList.setCopy(data, curOffset);
        curOffset += forwardReverseList.getLength();
        reverseForwardList.setCopy(data, curOffset);
        curOffset += reverseForwardList.getLength();
        reverseReverseList.setCopy(data, curOffset);
        
        curOffset += reverseReverseList.getLength();
        kmer.setAsCopy(data, curOffset);
    }
	
    @Override
    public void write(DataOutput out) throws IOException {
        this.nodeIdList.write(out);
        this.forwardForwardList.write(out);
        this.forwardReverseList.write(out);
        this.reverseForwardList.write(out);
        this.reverseReverseList.write(out);
        this.kmer.write(out);
    }

    @Override
    public void readFields(DataInput in) throws IOException {
        reset();
        this.nodeIdList.readFields(in);
        this.forwardForwardList.readFields(in);
        this.forwardReverseList.readFields(in);
        this.reverseForwardList.readFields(in);
        this.reverseReverseList.readFields(in);
        this.kmer.readFields(in);
    }

    @Override
    public int compareTo(NodeWritable other) {
        return this.kmer.compareTo(other.kmer);
    }
    
    @Override
    public int hashCode() {
        return this.kmer.hashCode();
    }
    
    @Override
    public boolean equals(Object o) {
        if (o instanceof NodeWritable) {
            NodeWritable nw = (NodeWritable) o;
            return (this.nodeIdList.equals(nw.nodeIdList)
                    && this.forwardForwardList.equals(nw.forwardForwardList)
                    && this.forwardReverseList.equals(nw.forwardReverseList)
                    && this.reverseForwardList.equals(nw.reverseForwardList)
                    && this.reverseReverseList.equals(nw.reverseReverseList) && this.kmer.equals(nw.kmer));
        }
        return false;
    }
    
    @Override
    public String toString() {
        StringBuilder sbuilder = new StringBuilder();
        sbuilder.append('{');
        sbuilder.append(nodeIdList.toString()).append('\t');
        sbuilder.append(forwardForwardList.toString()).append('\t');
        sbuilder.append(forwardReverseList.toString()).append('\t');
        sbuilder.append(reverseForwardList.toString()).append('\t');
        sbuilder.append(reverseReverseList.toString()).append('\t');
        sbuilder.append(kmer.toString()).append('}');
        return sbuilder.toString();
    }

    public void mergeForwardNext(final NodeWritable nextNode, int initialKmerSize) {
        this.forwardForwardList.setCopy(nextNode.forwardForwardList);
        this.forwardReverseList.setCopy(nextNode.forwardReverseList);
        kmer.mergeWithFFKmer(initialKmerSize, nextNode.getKmer());
    }

    public void mergeForwardPre(final NodeWritable preNode, int initialKmerSize) {
        this.reverseForwardList.setCopy(preNode.reverseForwardList);
        this.reverseReverseList.setCopy(preNode.reverseReverseList);
        kmer.mergeWithRRKmer(initialKmerSize, preNode.getKmer());
    }
    
    public int inDegree() {
        return reverseReverseList.getCountOfPosition() + reverseForwardList.getCountOfPosition();
    }

    public int outDegree() {
        return forwardForwardList.getCountOfPosition() + forwardReverseList.getCountOfPosition();
    }

    /*
     * Return if this node is a "path" compressible node, that is, it has an in-degree and out-degree of 1 
     */
    public boolean isPathNode() {
        return inDegree() == 1 && outDegree() == 1;
    }

    public boolean isSimpleOrTerminalPath() {
        return isPathNode() || (inDegree() == 0 && outDegree() == 1) || (inDegree() == 1 && outDegree() == 0);
    }
}
