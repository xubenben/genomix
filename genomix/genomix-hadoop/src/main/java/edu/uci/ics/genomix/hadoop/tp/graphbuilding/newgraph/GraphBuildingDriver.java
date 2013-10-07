package edu.uci.ics.genomix.hadoop.tp.graphbuilding.newgraph;

import java.io.IOException;

import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapred.FileInputFormat;
import org.apache.hadoop.mapred.FileOutputFormat;
import org.apache.hadoop.mapred.JobClient;
import org.apache.hadoop.mapred.JobConf;
import org.apache.hadoop.mapred.SequenceFileOutputFormat;
import org.apache.hadoop.mapred.TextInputFormat;
import org.apache.hadoop.mapred.TextOutputFormat;
import org.apache.hadoop.mapred.lib.NLineInputFormat;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.uci.ics.genomix.type.KmerBytesWritable;
import edu.uci.ics.genomix.type.NodeWritable;
import edu.uci.ics.genomix.type.VKmerBytesWritable;

@SuppressWarnings("deprecation")
public class GraphBuildingDriver {
    
    private static class Options {
        @Option(name = "-inputpath", usage = "the input path", required = true)
        public String inputPath;

        @Option(name = "-outputpath", usage = "the output path", required = true)
        public String outputPath;

        @Option(name = "-num-reducers", usage = "the number of reducers", required = true)
        public int numReducers;

        @Option(name = "-kmer-size", usage = "the size of kmer", required = true)
        public int sizeKmer;
        
        @Option(name = "-num-lines-per-map", usage = "the number of lines per map", required = false)
        public int linesPerMap;
    }
    
    public void run(String inputPath, String outputPath, int numReducers, int sizeKmer,
            int linesPerMap, boolean seqOutput, String defaultConfPath) throws IOException{
        JobConf conf = new JobConf(GraphBuildingDriver.class);
        conf.setInt("sizeKmer", sizeKmer);
        if (defaultConfPath != null) {
            conf.addResource(new Path(defaultConfPath));
        }

        conf.setJobName("Genomix Graph Building");
        conf.setMapperClass(GraphBuildingMapper.class);
        conf.setReducerClass(GraphBuildingReducer.class);

        conf.setMapOutputKeyClass(KmerBytesWritable.class);
        conf.setMapOutputValueClass(NodeWritable.class);
        
        //InputFormat and OutputFormat for Reducer
//        conf.setInputFormat(NLineInputFormat.class);
//        conf.setInt("mapred.line.input.format.linespermap", linesPerMap);
//        conf.setInt("io.sort.mb", 150);
        conf.setInputFormat(TextInputFormat.class);
        if (seqOutput == true)
            conf.setOutputFormat(SequenceFileOutputFormat.class);
        else
            conf.setOutputFormat(TextOutputFormat.class);
        
        //Output Key/Value Class
        conf.setOutputKeyClass(VKmerBytesWritable.class);
        conf.setOutputValueClass(NodeWritable.class);
        
        FileInputFormat.setInputPaths(conf, new Path(inputPath));
        FileOutputFormat.setOutputPath(conf, new Path(outputPath));
        conf.setNumReduceTasks(numReducers);

        FileSystem dfs = FileSystem.get(conf);
        dfs.delete(new Path(outputPath), true);
        JobClient.runJob(conf);
    }
    
    public static void main(String[] args) throws Exception {
        Options options = new Options();
        CmdLineParser parser = new CmdLineParser(options);
        parser.parseArgument(args);
        GraphBuildingDriver driver = new GraphBuildingDriver();
        driver.run(options.inputPath, options.outputPath, options.numReducers, 
                options.sizeKmer, options.linesPerMap, true, null);
    }
}