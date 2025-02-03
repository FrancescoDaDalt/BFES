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
        DataStream<BFES_Snapshot> result = keyValueStream.keyBy(t -> t.f0)
                .process(new BFESFunction<String, Integer>(100, 2));

        // Output
        result.print();

        // Execute the Flink job
        env.execute("BFES Function Example");
    }
}
