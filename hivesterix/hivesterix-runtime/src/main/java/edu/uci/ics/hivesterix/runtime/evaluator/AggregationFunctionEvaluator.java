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
package edu.uci.ics.hivesterix.runtime.evaluator;

import java.io.DataOutput;
import java.io.IOException;
import java.util.List;

import org.apache.hadoop.hive.ql.exec.ExprNodeEvaluator;
import org.apache.hadoop.hive.ql.metadata.HiveException;
import org.apache.hadoop.hive.ql.plan.ExprNodeDesc;
import org.apache.hadoop.hive.ql.udf.generic.GenericUDAFCount;
import org.apache.hadoop.hive.ql.udf.generic.GenericUDAFEvaluator;
import org.apache.hadoop.hive.ql.udf.generic.GenericUDAFEvaluator.AggregationBuffer;
import org.apache.hadoop.hive.serde2.SerDe;
import org.apache.hadoop.hive.serde2.SerDeException;
import org.apache.hadoop.hive.serde2.objectinspector.ObjectInspector;
import org.apache.hadoop.hive.serde2.objectinspector.PrimitiveObjectInspector;
import org.apache.hadoop.hive.serde2.typeinfo.TypeInfo;
import org.apache.hadoop.io.BytesWritable;

import edu.uci.ics.hivesterix.serde.lazy.LazyObject;
import edu.uci.ics.hyracks.algebricks.common.exceptions.AlgebricksException;
import edu.uci.ics.hyracks.algebricks.runtime.base.ICopyAggregateFunction;
import edu.uci.ics.hyracks.dataflow.common.data.accessors.IFrameTupleReference;

@SuppressWarnings("deprecation")
public class AggregationFunctionEvaluator implements ICopyAggregateFunction {

    /**
     * the mode of aggregation function
     */
    private GenericUDAFEvaluator.Mode mode;

    /**
     * an array of evaluators
     */
    private ExprNodeEvaluator[] evaluators;

    /**
     * udaf evaluator partial
     */
    private GenericUDAFEvaluator udafPartial;

    /**
     * udaf evaluator complete
     */
    private GenericUDAFEvaluator udafComplete;

    /**
     * cached parameter objects
     */
    private Object[] cachedParameters;

    /**
     * cached row objects
     */
    private LazyObject<? extends ObjectInspector> cachedRowObject;

    /**
     * the output channel
     */
    private DataOutput out;

    /**
     * aggregation buffer
     */
    private AggregationBuffer aggBuffer;

    /**
     * we only use lazy serde to do serialization
     */
    private SerDe lazySer;

    /**
     * the output object inspector for this aggregation function
     */
    private ObjectInspector outputInspector;

    /**
     * the output object inspector for this aggregation function
     */
    private ObjectInspector outputInspectorPartial;

    /**
     * parameter inspectors
     */
    private ObjectInspector[] parameterInspectors;

    /**
     * output make sure the aggregation functio has least object creation
     * 
     * @param desc
     * @param oi
     * @param output
     */
    public AggregationFunctionEvaluator(List<ExprNodeDesc> inputs, List<TypeInfo> inputTypes, String genericUDAFName,
            GenericUDAFEvaluator.Mode aggMode, boolean distinct, ObjectInspector oi, DataOutput output,
            ExprNodeEvaluator[] evals, ObjectInspector[] pInspectors, Object[] parameterCache, SerDe serde,
            LazyObject<? extends ObjectInspector> row, GenericUDAFEvaluator udafunctionPartial,
            GenericUDAFEvaluator udafunctionComplete, ObjectInspector outputOi, ObjectInspector outputOiPartial) {
        // shared object across threads
        this.out = output;
        this.mode = aggMode;
        this.parameterInspectors = pInspectors;

        // thread local objects
        this.evaluators = evals;
        this.cachedParameters = parameterCache;
        this.cachedRowObject = row;
        this.lazySer = serde;
        this.udafPartial = udafunctionPartial;
        this.udafComplete = udafunctionComplete;
        this.outputInspector = outputOi;
        this.outputInspectorPartial = outputOiPartial;
    }

    @Override
    public void init() throws AlgebricksException {
        try {
            aggBuffer = udafPartial.getNewAggregationBuffer();
        } catch (HiveException e) {
            throw new AlgebricksException(e);
        }
    }

    @Override
    public void step(IFrameTupleReference tuple) throws AlgebricksException {
        readIntoCache(tuple);
        processRow();
    }

    private void processRow() throws AlgebricksException {
        try {
            // get values by evaluating them
            for (int i = 0; i < cachedParameters.length; i++) {
                cachedParameters[i] = evaluators[i].evaluate(cachedRowObject);
            }
            processAggregate();
        } catch (HiveException e) {
            throw new AlgebricksException(e);
        }
    }

    private void processAggregate() throws HiveException {
        /**
         * accumulate the aggregation function
         */
        switch (mode) {
            case PARTIAL1:
            case COMPLETE:
                udafPartial.iterate(aggBuffer, cachedParameters);
                break;
            case PARTIAL2:
            case FINAL:
                if (udafPartial instanceof GenericUDAFCount.GenericUDAFCountEvaluator) {
                    Object parameter = ((PrimitiveObjectInspector) parameterInspectors[0])
                            .getPrimitiveWritableObject(cachedParameters[0]);
                    udafPartial.merge(aggBuffer, parameter);
                } else
                    udafPartial.merge(aggBuffer, cachedParameters[0]);
                break;
            default:
                break;
        }
    }

    /**
     * serialize the result
     * 
     * @param result
     *            the evaluation result
     * @throws IOException
     * @throws AlgebricksException
     */
    private void serializeResult(Object result, ObjectInspector oi) throws IOException, AlgebricksException {
        try {
            BytesWritable outputWritable = (BytesWritable) lazySer.serialize(result, oi);
            out.write(outputWritable.getBytes(), 0, outputWritable.getLength());
        } catch (SerDeException e) {
            throw new AlgebricksException(e);
        }
    }

    /**
     * bind the tuple reference to the cached row object
     * 
     * @param r
     */
    private void readIntoCache(IFrameTupleReference r) {
        cachedRowObject.init(r);
    }

    @Override
    public void finish() throws AlgebricksException {
        // aggregator
        try {
            Object result = null;
            result = udafPartial.terminatePartial(aggBuffer);
            if (mode == GenericUDAFEvaluator.Mode.COMPLETE || mode == GenericUDAFEvaluator.Mode.FINAL) {
                result = udafComplete.terminate(aggBuffer);
                serializeResult(result, outputInspector);
            } else {
                serializeResult(result, outputInspectorPartial);
            }
        } catch (HiveException e) {
            throw new AlgebricksException(e);
        } catch (IOException e) {
            throw new AlgebricksException(e);
        }
    }

    @Override
    public void finishPartial() throws AlgebricksException {
        // aggregator.
        try {
            Object result = null;
            // get aggregations
            result = udafPartial.terminatePartial(aggBuffer);
            serializeResult(result, outputInspectorPartial);
        } catch (HiveException e) {
            throw new AlgebricksException(e);
        } catch (IOException e) {
            throw new AlgebricksException(e);
        }
    }
}
