/* 
 * MIT License
 * 
 * Copyright (c) 2024 Francesco Da Dalt
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package FlinkBFESTest;

import org.apache.flink.api.common.eventtime.WatermarkStrategy;
import org.apache.flink.api.common.functions.FlatMapFunction;
import org.apache.flink.api.common.serialization.SimpleStringSchema;
import org.apache.flink.streaming.api.datastream.DataStream;
import org.apache.flink.streaming.api.environment.StreamExecutionEnvironment;
import org.apache.flink.streaming.api.functions.KeyedProcessFunction;
import org.apache.flink.util.Collector;
import org.apache.flink.api.java.tuple.Tuple2;
import org.apache.flink.util.OutputTag;
import java.lang.Integer;

import java.util.Properties;
import java.util.ArrayList;
import java.util.List;

public class FlinkBFESTest {

    public static void main(String[] args) throws Exception {
        // Set up the streaming execution environment
        StreamExecutionEnvironment env = StreamExecutionEnvironment.getExecutionEnvironment();

        List<String> elements = new ArrayList<>();
        for (int i = 0; i < 1000000; i++) {
            elements.add("key" + (i % 1000) + ":" + (i % 42));
        }

        // Create a source (this example uses a flatMap for generating dummy data;
        // replace with a real Kafka source or other)
        @SuppressWarnings("deprecation")
        DataStream<String> source = env.fromCollection(elements);

        // Transform the stream into tuples of key-value pairs
        DataStream<Tuple2<String, Integer>> keyValueStream = source.flatMap(
                new FlatMapFunction<String, Tuple2<String, Integer>>() {
                    @Override
                    public void flatMap(String value, Collector<Tuple2<String, Integer>> out) {
                        String[] parts = value.split(":");
                        out.collect(new Tuple2<String, Integer>(parts[0], Integer.parseInt(parts[1])));
                    }
                });

        // Apply the BFESFunction
        @SuppressWarnings("deprecation")
        DataStream<BFES_Snapshot> result = keyValueStream.keyBy(t -> t.f0)
                .process(new BFESFunction<String, Integer>(100, 2));
        
        // Output
        result.print();

        // Execute the Flink job
        env.execute("BFES Function Example");
    }
}
