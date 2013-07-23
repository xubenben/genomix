package edu.uci.ics.genomix.hadoop.contrailgraphbuilding;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hdfs.MiniDFSCluster;
import org.apache.hadoop.mapred.JobConf;
import org.apache.hadoop.mapred.MiniMRCluster;
import org.junit.Test;

@SuppressWarnings("deprecation")
public class GraphBuildingTest {

    private JobConf conf = new JobConf();
    private static final String ACTUAL_RESULT_DIR = "actual";
    private static final String HADOOP_CONF_PATH = ACTUAL_RESULT_DIR + File.separator + "conf.xml";
    private static final String DATA_PATH = "data/webmap/test.txt";
    private static final String HDFS_PATH = "/webmap";
    private static final String RESULT_PATH = "/result";
    
//    private static final int COUNT_REDUCER = 2;
    private static final int SIZE_KMER = 5;
    private static final int READ_LENGTH = 8;
    
    private MiniDFSCluster dfsCluster;
    private MiniMRCluster mrCluster;
    private FileSystem dfs;
    
    @Test
    public void test() throws Exception {
        FileUtils.forceMkdir(new File(ACTUAL_RESULT_DIR));
        FileUtils.cleanDirectory(new File(ACTUAL_RESULT_DIR));
        startHadoop();
        TestMapKmerToNode();
        cleanupHadoop();
    }
    
    public void TestMapKmerToNode() throws Exception {
        GenomixDriver driver = new GenomixDriver();
        driver.run(HDFS_PATH, RESULT_PATH, 0, SIZE_KMER, READ_LENGTH, false, HADOOP_CONF_PATH);
        dumpResult();
    }
    
    private void startHadoop() throws IOException {
        FileSystem lfs = FileSystem.getLocal(new Configuration());
        lfs.delete(new Path("build"), true);
        System.setProperty("hadoop.log.dir", "logs");
        dfsCluster = new MiniDFSCluster(conf, 2, true, null);
        dfs = dfsCluster.getFileSystem();
        mrCluster = new MiniMRCluster(4, dfs.getUri().toString(), 2);

        Path src = new Path(DATA_PATH);
        Path dest = new Path(HDFS_PATH + "/");
        dfs.mkdirs(dest);
        dfs.copyFromLocalFile(src, dest);

        DataOutputStream confOutput = new DataOutputStream(new FileOutputStream(new File(HADOOP_CONF_PATH)));
        conf.writeXml(confOutput);
        confOutput.flush();
        confOutput.close();
    }

    private void cleanupHadoop() throws IOException {
        mrCluster.shutdown();
        dfsCluster.shutdown();
    }
    
    private void dumpResult() throws IOException {
        Path src = new Path(RESULT_PATH);
        Path dest = new Path(ACTUAL_RESULT_DIR);
        dfs.copyToLocalFile(src, dest);
    }
}
