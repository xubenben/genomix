package edu.uci.ics.genomix.hyracks.dataflow.aggregators;

import java.io.DataOutput;
import java.io.IOException;
import java.nio.ByteBuffer;

import edu.uci.ics.hyracks.api.comm.IFrameTupleAccessor;
import edu.uci.ics.hyracks.api.context.IHyracksTaskContext;
import edu.uci.ics.hyracks.api.dataflow.value.RecordDescriptor;
import edu.uci.ics.hyracks.api.exceptions.HyracksDataException;
import edu.uci.ics.hyracks.data.std.util.ArrayBackedValueStorage;
import edu.uci.ics.hyracks.dataflow.common.comm.io.ArrayTupleBuilder;
import edu.uci.ics.hyracks.dataflow.std.group.AggregateState;
import edu.uci.ics.hyracks.dataflow.std.group.IAggregatorDescriptor;
import edu.uci.ics.hyracks.dataflow.std.group.IAggregatorDescriptorFactory;

public class MergeReadIDAggregateFactory implements IAggregatorDescriptorFactory {

    /**
     * 
     */
    private static final long serialVersionUID = 1L;

    private final int ValidPosCount;

    public MergeReadIDAggregateFactory(int readLength, int kmerLength) {
        ValidPosCount = getPositionCount(readLength, kmerLength);
    }

    public static int getPositionCount(int readLength, int kmerLength) {
        return readLength - kmerLength + 1;
    }

    public static final int InputReadIDField = AggregateReadIDAggregateFactory.OutputReadIDField;
    public static final int InputPositionListField = AggregateReadIDAggregateFactory.OutputPositionListField;

    public static final int BYTE_SIZE = 1;
    public static final int INTEGER_SIZE = 4;

    /**
     * (ReadID, {(PosInRead,{OtherPositoin..},Kmer) ...} to Aggregate as
     * (ReadID, Storage[posInRead]={PositionList,Kmer})
     */
    @Override
    public IAggregatorDescriptor createAggregator(IHyracksTaskContext ctx, RecordDescriptor inRecordDescriptor,
            RecordDescriptor outRecordDescriptor, int[] keyFields, int[] keyFieldsInPartialResults)
            throws HyracksDataException {
        return new IAggregatorDescriptor() {

            class PositionArray {
                public ArrayBackedValueStorage[] forwardStorages;
                public ArrayBackedValueStorage[] reverseStorages;
                public int count;

                public PositionArray() {
                    forwardStorages = new ArrayBackedValueStorage[ValidPosCount];
                    reverseStorages = new ArrayBackedValueStorage[ValidPosCount];
                    for (int i = 0; i < ValidPosCount; i++) {
                        forwardStorages[i] = new ArrayBackedValueStorage();
                        reverseStorages[i] = new ArrayBackedValueStorage();
                    }
                    count = 0;
                }

                public void reset() {
                    for (int i = 0; i < ValidPosCount; i++) {
                        forwardStorages[i].reset();
                        reverseStorages[i].reset();
                    }
                    count = 0;
                }
            }

            @Override
            public AggregateState createAggregateStates() {

                return new AggregateState(new PositionArray());
            }

            @Override
            public void init(ArrayTupleBuilder tupleBuilder, IFrameTupleAccessor accessor, int tIndex,
                    AggregateState state) throws HyracksDataException {
                PositionArray positionArray = (PositionArray) state.state;
                positionArray.reset();

                pushIntoStorage(accessor, tIndex, positionArray);

                // make fake fields
                for (int i = 0; i < ValidPosCount * 2; i++) {
                    tupleBuilder.addFieldEndOffset();
                }
            }

            private void pushIntoStorage(IFrameTupleAccessor accessor, int tIndex, PositionArray positionArray)
                    throws HyracksDataException {
                int leadbyte = accessor.getTupleStartOffset(tIndex) + accessor.getFieldSlotsLength();
                int fieldOffset = leadbyte + accessor.getFieldStartOffset(tIndex, InputPositionListField);
                ByteBuffer fieldBuffer = accessor.getBuffer();

                while (fieldOffset < leadbyte + accessor.getFieldEndOffset(tIndex, InputPositionListField)) {
                    byte posInRead = fieldBuffer.get(fieldOffset);

                    ArrayBackedValueStorage[] storage = positionArray.forwardStorages;
                    boolean hasKmer = true;
                    if (posInRead < 0) {
                        storage = positionArray.reverseStorages;
                        posInRead = (byte) -posInRead;
                        hasKmer = false;
                    }
                    if (storage[posInRead - 1].getLength() > 0) {
                        throw new IllegalArgumentException("Reentering into an exist storage");
                    }
                    fieldOffset += BYTE_SIZE;

                    // read poslist
                    fieldOffset += writeBytesToStorage(storage[posInRead - 1], fieldBuffer, fieldOffset);
                    // read Kmer
                    if (hasKmer) {
                        fieldOffset += writeBytesToStorage(storage[posInRead - 1], fieldBuffer, fieldOffset);
                    }

                    positionArray.count += 1;
                }

            }

            private int writeBytesToStorage(ArrayBackedValueStorage storage, ByteBuffer fieldBuffer, int fieldOffset)
                    throws HyracksDataException {
                int lengthPosList = fieldBuffer.getInt(fieldOffset);
                try {
                    storage.getDataOutput().writeInt(lengthPosList);
                    fieldOffset += INTEGER_SIZE;
                    if (lengthPosList > 0) {
                        storage.getDataOutput().write(fieldBuffer.array(), fieldOffset, lengthPosList);
                    }
                } catch (IOException e) {
                    throw new HyracksDataException("Failed to write into temporary storage");
                }
                return lengthPosList + INTEGER_SIZE;
            }

            @Override
            public void reset() {
                // TODO Auto-generated method stub

            }

            @Override
            public void aggregate(IFrameTupleAccessor accessor, int tIndex, IFrameTupleAccessor stateAccessor,
                    int stateTupleIndex, AggregateState state) throws HyracksDataException {
                PositionArray positionArray = (PositionArray) state.state;
                pushIntoStorage(accessor, tIndex, positionArray);
            }

            @Override
            public void outputPartialResult(ArrayTupleBuilder tupleBuilder, IFrameTupleAccessor accessor, int tIndex,
                    AggregateState state) throws HyracksDataException {
                throw new IllegalStateException("partial result method should not be called");
            }

            @Override
            public void outputFinalResult(ArrayTupleBuilder tupleBuilder, IFrameTupleAccessor accessor, int tIndex,
                    AggregateState state) throws HyracksDataException {
                PositionArray positionArray = (PositionArray) state.state;

                if (positionArray.count != ValidPosCount * 2) {
                    throw new IllegalStateException("Final aggregate position number is invalid");
                }
                DataOutput fieldOutput = tupleBuilder.getDataOutput();
                try {
                    for (int i = 0; i < ValidPosCount; i++) {
                        fieldOutput.write(positionArray.forwardStorages[i].getByteArray(),
                                positionArray.forwardStorages[i].getStartOffset(),
                                positionArray.forwardStorages[i].getLength());
                        tupleBuilder.addFieldEndOffset();

                        fieldOutput.write(positionArray.reverseStorages[i].getByteArray(),
                                positionArray.reverseStorages[i].getStartOffset(),
                                positionArray.reverseStorages[i].getLength());
                        tupleBuilder.addFieldEndOffset();
                    }
                } catch (IOException e) {
                    throw new HyracksDataException("I/O exception when writing aggregation to the output buffer.");
                }
            }

            @Override
            public void close() {
                // TODO Auto-generated method stub

            }

        };
    }
}