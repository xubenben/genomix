/*
 * Copyright 2009-2012 by The Regents of the University of California
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * you may obtain a copy of the License from
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package edu.uci.ics.genomix.hadoop.gage;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hdfs.MiniDFSCluster;
import org.apache.hadoop.io.SequenceFile;
import org.apache.hadoop.io.SequenceFile.Writer;
import org.apache.hadoop.mapred.Counters;
import org.apache.hadoop.mapred.JobConf;
import org.apache.hadoop.mapred.MiniMRCluster;
import org.junit.Test;

import edu.uci.ics.genomix.config.GenomixJobConf;
import edu.uci.ics.genomix.hadoop.graph.GraphStatistics;
import edu.uci.ics.genomix.hadoop.graphbuilding.checkingtool.ResultsCheckingDriver;
import edu.uci.ics.genomix.hadoop.pmcommon.HadoopMiniClusterTest;
import edu.uci.ics.genomix.type.Node;
import edu.uci.ics.genomix.type.VKmer;

public class GetFastaStatsTest {
    private JobConf conf = new JobConf();
    private static final String ACTUAL_RESULT_DIR = "actual";
    private static final String BINSOURCE_PATH = "GAGETest.fasta";
    private static final String TXTSOURCE_PATH = "gageTxt.fasta";
    private static final String HDFS_PATH = "/gage";
    private static final String HADOOP_CONF_PATH = ACTUAL_RESULT_DIR + File.separator + "conf.xml";
    private static final String RESULT_PATH = "/gageresult";

    private static final int COUNT_REDUCER = 4;
    private MiniDFSCluster dfsCluster;
    private MiniMRCluster mrCluster;
    private FileSystem dfs;

    public void writeTestCaseToFile() throws IOException {
        dfs = FileSystem.getLocal(conf);
        //        String targetFileName = "gageTest";
        Path targetPath = new Path(BINSOURCE_PATH);

        SequenceFile.Writer writer = new SequenceFile.Writer(dfs, conf, targetPath, VKmer.class, Node.class);
        String[] keyList = { "AGCTG", "CTGAC", "GTACT", "TCAGA", "CATCT", "GACTC" };
        String[] valueList = { "AGTCTCTCTCTCTAGACTCTCTCTTTCTAGAGCTCTCAG", "AGTCTCTCTCTCTAGACTCTCTCTTTCTAGAGCTCTCAGCGT",
                "AGTCTCTCTCTCTAGACTCTCTCTTTCT", "AGTCTCTCTCTCTAGACTCTCTCTTTCTATCGA", "ATCGC",
                "CTTTCTAGAGCTCTCTCTCTCTAGACTCTCTCTCTC" };
        VKmer outputKey = new VKmer();
        Node outputValue = new Node();
        VKmer tempInternalKmer = new VKmer();
        for (int i = 0; i < keyList.length; i++) {
            outputKey.setAsCopy(keyList[i]);
            tempInternalKmer.setAsCopy(valueList[i]);
            outputValue.setInternalKmer(tempInternalKmer);
            writer.append(outputKey, outputValue);
        }
        writer.close();

        String textTargetFileName = TXTSOURCE_PATH;
        BufferedWriter txtWriter = new BufferedWriter(new FileWriter(textTargetFileName));
        for (int i = 0; i < keyList.length; i++) {
            txtWriter.write(">node_" + keyList[i]);
            txtWriter.newLine();
            txtWriter.write(valueList[i]);
            txtWriter.newLine();
        }
        txtWriter.close();
    }

    @SuppressWarnings("deprecation")
    @Test
    public void Test() throws Exception {
        writeTestCaseToFile();
        FileUtils.forceMkdir(new File(ACTUAL_RESULT_DIR));
        FileUtils.cleanDirectory(new File(ACTUAL_RESULT_DIR));
        startHadoop();
        GetFastaStatsJob driver = new GetFastaStatsJob();
        JobConf conf = new JobConf(HADOOP_CONF_PATH);
        conf.setInt(GenomixJobConf.MIN_CONTIG_LENGTH, 25);
        conf.setInt(GenomixJobConf.EXPECTED_GENOME_SIZE, 150);
        conf.setBoolean(GenomixJobConf.USE_BAYLOR_FORMAT, false);
        conf.setBoolean(GenomixJobConf.OLD_STYLE, true);
        Counters counters = GraphStatistics.run(HDFS_PATH, RESULT_PATH, conf);
        GraphStatistics.getFastaStatsForGage(RESULT_PATH, counters, conf);
//        driver.run(HDFS_PATH, RESULT_PATH, COUNT_REDUCER, conf);
        dumpResult();
        cleanupHadoop();
        if (!compareWithGageSourceCodeResults("actual/metrics.txt"))
            throw new Exception("the results are not same!");
    }

    public boolean compareWithGageSourceCodeResults(String mapreducePath) throws Exception {
        String[] args = { "-o", "-min", "25", "-genomeSize", "150", TXTSOURCE_PATH };
        String src = gageDriver(args);
        BufferedReader br = new BufferedReader(new FileReader(new File(mapreducePath)));
        StringBuffer target = new StringBuffer();
        String line;
        while ((line = br.readLine()) != null) {
            target.append(line);
        }
        String pureTarget = target.toString().replaceAll("\\n", "");
        String pureSrc = src.replaceAll("\\n", "");
        if (pureTarget.equals(pureSrc))
            return true;
        else
            return false;
    }

    public static String gageDriver(String[] args) throws Exception {
        if (args.length < 1) {
            //            printUsage();
            System.exit(1);
        }

        boolean useBaylorFormat = false;
        boolean oldStyle = false;
        long genomeSize = 0;
        int initialVal = 0;
        while (args[initialVal].startsWith("-")) {
            if (args[initialVal].trim().equalsIgnoreCase("-b")) {
                useBaylorFormat = true;
            } else if (args[initialVal].trim().equalsIgnoreCase("-min")) {
                GageSourceGetFastaStats.MIN_LENGTH = Integer.parseInt(args[++initialVal]);
            } else if (args[initialVal].trim().equalsIgnoreCase("-o")) {
                oldStyle = true;
            } else if (args[initialVal].trim().equalsIgnoreCase("-genomeSize")) {
                initialVal++;
                genomeSize = Long.parseLong(args[initialVal]);
                System.err.println("Found genome size at position " + initialVal + " with value " + args[initialVal]
                        + " aka " + genomeSize);
            } else {
                System.err.println("Unknown parameter " + args[initialVal]
                        + " specified, please specify -b for Baylor-style output.");
                System.exit(1);
            }
            initialVal++;
        }

        //        for (int i = initialVal; i < args.length; i++) {
        String assemblyTitle = args[initialVal].trim().split("/")[0];
        GageSourceGetFastaStats f = new GageSourceGetFastaStats(useBaylorFormat, oldStyle, genomeSize);
        String[] splitLine = args[initialVal].trim().split(",");

        for (int j = 0; j < splitLine.length; j++) {
            f.processFile(splitLine[j]);
        }
        return (f.toString(true, assemblyTitle));
        //        }

    }

    public void startHadoop() throws IOException {
        FileSystem lfs = FileSystem.getLocal(new Configuration());
        lfs.delete(new Path("build"), true);
        System.setProperty("hadoop.log.dir", "logs");
        dfsCluster = new MiniDFSCluster(conf, 1, true, null);
        dfs = dfsCluster.getFileSystem();
        mrCluster = new MiniMRCluster(1, dfs.getUri().toString(), 1);

        Path src = new Path(BINSOURCE_PATH);
        Path dest = new Path(HDFS_PATH + "/");
        dfs.mkdirs(dest);
        dfs.copyFromLocalFile(src, dest);

        DataOutputStream confOutput = new DataOutputStream(new FileOutputStream(new File(HADOOP_CONF_PATH)));
        conf.writeXml(confOutput);
        confOutput.flush();
        confOutput.close();
    }

    public void cleanupHadoop() throws IOException {
        mrCluster.shutdown();
        dfsCluster.shutdown();
    }

    private void dumpResult() throws IOException {
        Path src = new Path(RESULT_PATH);
        Path dest = new Path(ACTUAL_RESULT_DIR);
        dfs.copyToLocalFile(src, dest);
        HadoopMiniClusterTest.copyResultsToLocal(RESULT_PATH, "actual/metrics.txt", true, conf, true, dfs);
    }
}
