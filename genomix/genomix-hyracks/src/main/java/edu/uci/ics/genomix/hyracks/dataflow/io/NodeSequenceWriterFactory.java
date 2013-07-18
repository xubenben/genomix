/*
 * Copyright 2009-2013 by The Regents of the University of California
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
package edu.uci.ics.genomix.hyracks.dataflow.io;

import java.io.DataOutput;
import java.io.IOException;

import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.SequenceFile;
import org.apache.hadoop.io.SequenceFile.CompressionType;
import org.apache.hadoop.io.SequenceFile.Writer;
import org.apache.hadoop.mapred.JobConf;

import edu.uci.ics.genomix.data.Marshal;
import edu.uci.ics.genomix.hyracks.dataflow.MapReadToNodeOperator;
import edu.uci.ics.genomix.hyracks.job.GenomixJobConf;
import edu.uci.ics.genomix.oldtype.NodeWritable;
import edu.uci.ics.genomix.oldtype.PositionWritable;
import edu.uci.ics.hyracks.api.context.IHyracksTaskContext;
import edu.uci.ics.hyracks.api.exceptions.HyracksDataException;
import edu.uci.ics.hyracks.dataflow.common.data.accessors.ITupleReference;
import edu.uci.ics.hyracks.hdfs.api.ITupleWriter;
import edu.uci.ics.hyracks.hdfs.api.ITupleWriterFactory;
import edu.uci.ics.hyracks.hdfs.dataflow.ConfFactory;

@SuppressWarnings("deprecation")
public class NodeSequenceWriterFactory implements ITupleWriterFactory {

    /**
     * 
     */
    private static final long serialVersionUID = 1L;

    public static final int InputNodeIDField = MapReadToNodeOperator.OutputNodeIDField;
    public static final int InputCountOfKmerField = MapReadToNodeOperator.OutputCountOfKmerField;
    public static final int InputFFField = MapReadToNodeOperator.OutputForwardForwardField;
    public static final int InputFRField = MapReadToNodeOperator.OutputForwardReverseField;
    public static final int InputRFField = MapReadToNodeOperator.OutputReverseForwardField;
    public static final int InputRRField = MapReadToNodeOperator.OutputReverseReverseField;

    public static final int InputKmerBytesField = MapReadToNodeOperator.OutputKmerBytesField;

    private ConfFactory confFactory;
    private final int kmerlength;

    public NodeSequenceWriterFactory(JobConf conf) throws HyracksDataException {
        this.confFactory = new ConfFactory(conf);
        this.kmerlength = conf.getInt(GenomixJobConf.KMER_LENGTH, GenomixJobConf.DEFAULT_KMERLEN);
    }

    public class TupleWriter implements ITupleWriter {

        public TupleWriter(ConfFactory confFactory) {
            this.cf = confFactory;
        }

        ConfFactory cf;
        Writer writer = null;
        NodeWritable node = new NodeWritable(kmerlength);

        @Override
        public void open(DataOutput output) throws HyracksDataException {
            try {
                writer = SequenceFile.createWriter(cf.getConf(), (FSDataOutputStream) output, NodeWritable.class,
                        NullWritable.class, CompressionType.NONE, null);
            } catch (IOException e) {
                throw new HyracksDataException(e);
            }
        }

        @Override
        public void write(DataOutput output, ITupleReference tuple) throws HyracksDataException {
            node.getNodeID().setNewReference(tuple.getFieldData(InputNodeIDField),
                    tuple.getFieldStart(InputNodeIDField));
            node.getFFList().setNewReference(tuple.getFieldLength(InputFFField) / PositionWritable.LENGTH,
                    tuple.getFieldData(InputFFField), tuple.getFieldStart(InputFFField));
            node.getFRList().setNewReference(tuple.getFieldLength(InputFRField) / PositionWritable.LENGTH,
                    tuple.getFieldData(InputFRField), tuple.getFieldStart(InputFRField));
            node.getRFList().setNewReference(tuple.getFieldLength(InputRFField) / PositionWritable.LENGTH,
                    tuple.getFieldData(InputRFField), tuple.getFieldStart(InputRFField));
            node.getRRList().setNewReference(tuple.getFieldLength(InputRRField) / PositionWritable.LENGTH,
                    tuple.getFieldData(InputRRField), tuple.getFieldStart(InputRRField));

            node.getKmer().setNewReference(
                    Marshal.getInt(tuple.getFieldData(NodeSequenceWriterFactory.InputCountOfKmerField),
                            tuple.getFieldStart(NodeSequenceWriterFactory.InputCountOfKmerField)),
                    tuple.getFieldData(InputKmerBytesField), tuple.getFieldStart(InputKmerBytesField));

            try {
                writer.append(node, NullWritable.get());
            } catch (IOException e) {
                throw new HyracksDataException(e);
            }
        }

        @Override
        public void close(DataOutput output) throws HyracksDataException {
        }

    }

    /**
     * Input schema:
     * (Position, LengthCount, InComingPosList, OutgoingPosList, Kmer)
     */
    @Override
    public ITupleWriter getTupleWriter(IHyracksTaskContext ctx) throws HyracksDataException {
        return new TupleWriter(confFactory);
    }

}
