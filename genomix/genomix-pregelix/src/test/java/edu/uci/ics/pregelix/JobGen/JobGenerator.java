package edu.uci.ics.pregelix.JobGen;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.ByteWritable;
import org.apache.hadoop.io.BytesWritable;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;

import edu.uci.ics.pregelix.BinaryLoadGraphInputFormat;
import edu.uci.ics.pregelix.BinaryLoadGraphOutputFormat;
import edu.uci.ics.pregelix.LoadGraphVertex;
import edu.uci.ics.pregelix.MergeGraphVertex;
import edu.uci.ics.pregelix.LoadGraphVertex.SimpleLoadGraphVertexOutputFormat;
import edu.uci.ics.pregelix.TestLoadGraphVertex;
import edu.uci.ics.pregelix.TextLoadGraphInputFormat;
import edu.uci.ics.pregelix.api.job.PregelixJob;


public class JobGenerator {

    private static String outputBase = "src/test/resources/jobs/";
    private static String HDFS_INPUTPATH = "/webmap";
    private static String HDFS_OUTPUTPAH = "/result";
    
    private static void generateLoadGraphJob(String jobName, String outputPath) throws IOException {
    	PregelixJob job = new PregelixJob(jobName);
    	job.setVertexClass(LoadGraphVertex.class);
    	job.setVertexInputFormatClass(TextLoadGraphInputFormat.class);
        job.setVertexOutputFormatClass(SimpleLoadGraphVertexOutputFormat.class);
        FileInputFormat.setInputPaths(job, HDFS_INPUTPATH);
        FileOutputFormat.setOutputPath(job, new Path(HDFS_OUTPUTPAH));
        job.getConfiguration().setLong(PregelixJob.NUM_VERTICE, 20);
        job.getConfiguration().writeXml(new FileOutputStream(new File(outputPath)));
    }
    
    private static void genLoadGraph() throws IOException {
    	generateLoadGraphJob("LoadGraph", outputBase + "LoadGraph.xml");
    }
    
    private static void generateBinaryLoadGraphJob(String jobName, String outputPath) throws IOException {
    	PregelixJob job = new PregelixJob(jobName);
    	job.setVertexClass(MergeGraphVertex.class);
    	job.setVertexInputFormatClass(BinaryLoadGraphInputFormat.class);
        job.setVertexOutputFormatClass(BinaryLoadGraphOutputFormat.class);
        job.setOutputKeyClass(BytesWritable.class);
        job.setOutputValueClass(ByteWritable.class);
        FileInputFormat.setInputPaths(job, HDFS_INPUTPATH);
        FileOutputFormat.setOutputPath(job, new Path(HDFS_OUTPUTPAH));
        job.getConfiguration().writeXml(new FileOutputStream(new File(outputPath)));
    }
    
    private static void genBinaryLoadGraph() throws IOException {
    	generateBinaryLoadGraphJob("BinaryLoadGraph", outputBase + "BinaryLoadGraph.xml");
    }
    
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		genBinaryLoadGraph();
		//genBasicBinaryLoadGraph();
	}

}