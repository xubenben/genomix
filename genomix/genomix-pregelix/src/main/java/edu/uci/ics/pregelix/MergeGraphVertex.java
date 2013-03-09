package edu.uci.ics.pregelix;

import java.io.IOException;
import java.util.Iterator;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.ByteWritable;
import org.apache.hadoop.io.BytesWritable;
import org.apache.hadoop.io.NullWritable;

import edu.uci.ics.pregelix.api.graph.Vertex;
import edu.uci.ics.pregelix.api.job.PregelixJob;
import edu.uci.ics.pregelix.api.util.BspUtils;
import edu.uci.ics.pregelix.dataflow.util.IterationUtils;
import edu.uci.ics.pregelix.example.client.Client;
import edu.uci.ics.pregelix.example.io.MessageWritable;
import edu.uci.ics.pregelix.hdfs.HDFSOperation;

/*
 * vertexId: BytesWritable
 * vertexValue: ByteWritable
 * edgeValue: NullWritable
 * message: MessageWritable
 * 
 * DNA:
 * A: 00
 * C: 01
 * G: 10
 * T: 11
 * 
 * succeed node
 *  A 00000001 1
 *  G 00000010 2
 *  C 00000100 4
 *  T 00001000 8
 * precursor node
 *  A 00010000 16
 *  G 00100000 32
 *  C 01000000 64
 *  T 10000000 128
 *  
 * For example, ONE LINE in input file: 00,01,10	0001,0010,
 * That means that vertexId is ACG, its succeed node is A and its precursor node is C.
 * The succeed node and precursor node will be stored in vertexValue and we don't use edgeValue.
 * The details about message are in edu.uci.ics.pregelix.example.io.MessageWritable. 
 */
public class MergeGraphVertex extends Vertex<BytesWritable, ByteWritable, NullWritable, MessageWritable>{
	
	private byte[] tmpSourceVertextId;
	private byte[] tmpDestVertexId;
	private byte[] tmpChainVertexId;
	private byte[] tmpNeighberBytes = new byte[1];
	private byte tmpVertexValue;
	private MessageWritable tmpMsg = new MessageWritable();
	public static final int k = 3; //kmer, k = 3
	/**
	 * Naive Algorithm for merge graph
	 */
	@Override
	public void compute(Iterator<MessageWritable> msgIterator) {
		if (getSuperstep() == 1) {
			if(GraphVertexOperation.isHead(getVertexValue())){
				tmpSourceVertextId = getVertexId().getBytes(); 
				tmpDestVertexId = GraphVertexOperation.getDestVertexId(tmpSourceVertextId,
						getVertexValue().get());
				tmpMsg.setSourceVertexIdOrNeighberInfo(tmpSourceVertextId);
				tmpChainVertexId = new byte[0];
				tmpMsg.setChainVertexId(tmpChainVertexId);
				sendMsg(new BytesWritable(tmpDestVertexId),tmpMsg);
				}
		}
		//path node sends message back to head node
		else if(getSuperstep()%2 == 0){
			 if(msgIterator.hasNext()){
				tmpMsg = msgIterator.next();
				if(!tmpMsg.isRear()){
					if(GraphVertexOperation.isPathVertex(getVertexValue())){
						tmpSourceVertextId = tmpMsg.getSourceVertexIdOrNeighberInfo();
					    //GraphVertexOperation.getDestVertexId(getVertexId().getBytes(), getVertexValue().get());
						tmpNeighberBytes[0] = getVertexValue().get();
						tmpMsg.setSourceVertexIdOrNeighberInfo(tmpNeighberBytes); //set neighber
						tmpChainVertexId = tmpMsg.getChainVertexId();
						if(tmpChainVertexId.length == 0){
							tmpMsg.setChainVertexId(getVertexId().getBytes());
							tmpMsg.setLengthOfChain(k);
						}
						else{
							tmpMsg.setChainVertexId(GraphVertexOperation.updateChainVertexId(tmpChainVertexId,
									tmpMsg.getLengthOfChain(),getVertexId().getBytes()));
							tmpMsg.incrementLength();
							//deleteVertex(getVertexId());
						}
						sendMsg(new BytesWritable(tmpSourceVertextId),tmpMsg);
						
					}
					else if(GraphVertexOperation.isRear(getVertexValue())){
						tmpSourceVertextId = tmpMsg.getSourceVertexIdOrNeighberInfo();
						tmpMsg.setSourceVertexIdOrNeighberInfo(getVertexId().getBytes());
						tmpMsg.setRear(true);
						sendMsg(new BytesWritable(tmpSourceVertextId),tmpMsg);
					}
					voteToHalt();
				}
				else{
					tmpVertexValue = GraphVertexOperation.updateRightNeighberByVertexId(getVertexValue().get(),
							tmpMsg.getSourceVertexIdOrNeighberInfo());
					setVertexValue(new ByteWritable(tmpVertexValue));
					setVertexId(new BytesWritable(tmpMsg.getChainVertexId()));
					signalTerminate();
				}
			}
			 else voteToHalt();
		}
		//head node sends message to path node
		else if(getSuperstep()%2 == 1){
			if (msgIterator.hasNext()){
				tmpMsg = msgIterator.next();
				tmpNeighberBytes = tmpMsg.getSourceVertexIdOrNeighberInfo();
				tmpChainVertexId = tmpMsg.getChainVertexId();
				if(!tmpMsg.isRear()){
					byte[] lastKmer = GraphVertexOperation.getLastKmer(tmpChainVertexId, 
							tmpMsg.getLengthOfChain());
					tmpDestVertexId = GraphVertexOperation.getDestVertexId(lastKmer, tmpNeighberBytes[0]);
					tmpSourceVertextId = getVertexId().getBytes();
					tmpMsg.setSourceVertexIdOrNeighberInfo(tmpSourceVertextId);
					sendMsg(new BytesWritable(tmpDestVertexId),tmpMsg);
				}
				else{	
					tmpDestVertexId = GraphVertexOperation.getDestVertexId(getVertexId().getBytes(),
							getVertexValue().get());
					sendMsg(new BytesWritable(tmpDestVertexId),tmpMsg);
					voteToHalt();
				}
			}
			voteToHalt();
		}
	}
	
   private void signalTerminate() {
        Configuration conf = getContext().getConfiguration();
        try {
            IterationUtils.writeForceTerminationState(conf, BspUtils.getJobId(conf));
            writeMergeGraphResult(conf, true);
        } catch (Exception e) {
            throw new IllegalStateException(e);
        }
    }
   
   private void writeMergeGraphResult(Configuration conf, boolean terminate) {
       try {
           FileSystem dfs = FileSystem.get(conf);
           String pathStr = IterationUtils.TMP_DIR + BspUtils.getJobId(conf) + "MergeGraph";
           Path path = new Path(pathStr);
           if (!dfs.exists(path)) {
               FSDataOutputStream output = dfs.create(path, true);
               output.writeBoolean(terminate);
               output.flush();
               output.close();
           }
       } catch (IOException e) {
           throw new IllegalStateException(e);
       }
   }

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
        PregelixJob job = new PregelixJob(MergeGraphVertex.class.getSimpleName());
        job.setVertexClass(MergeGraphVertex.class);
        /**
         * BinaryInput and BinaryOutput
         */
        job.setVertexInputFormatClass(BinaryLoadGraphInputFormat.class); 
        job.setVertexOutputFormatClass(BinaryLoadGraphOutputFormat.class); 
        job.setOutputKeyClass(BytesWritable.class);
        job.setOutputValueClass(ByteWritable.class);
        Client.run(args, job);
	}
}