import org.apache.flink.api.common.functions.RuntimeContext;
import org.apache.flink.api.java.tuple.Tuple2;
import org.apache.flink.streaming.api.datastream.DataStream;
import org.apache.flink.streaming.api.environment.StreamExecutionEnvironment;
import org.apache.flink.streaming.api.functions.KeyedProcessFunction;
import org.apache.flink.streaming.api.windowing.time.Time;
import org.apache.flink.streaming.api.collector.selector.OutputTag;
import org.apache.flink.streaming.api.datastream.KeyedStream;

public class FlinkBFESTest {

    public static void main(String[] args) throws Exception {
        // Set up the execution environment
        StreamExecutionEnvironment env = StreamExecutionEnvironment.getExecutionEnvironment();

        // Create a sample toy data stream of Key-Value pairs (e.g., String as key and Integer as value)
        DataStream<Tuple2<String, Integer>> sourceStream = env.fromElements(
            Tuple2.of("key1", 1),
            Tuple2.of("key1", 2),
            Tuple2.of("key2", 3),
            Tuple2.of("key1", 4),
            Tuple2.of("key2", 5),
            Tuple2.of("key3", 6)
        );

        // Key the stream by the key (first element of Tuple2)
        KeyedStream<Tuple2<String, Integer>, String> keyedStream = sourceStream.keyBy(value -> value.f0);

        // Apply the BFESFunction
        DataStream<BFES_Snapshot> resultStream = keyedStream.process(new BFESFunction<String, Integer>());

        // Print the result stream (you can replace this with a sink to save the output)
        resultStream.print();

        // Execute the Flink job
        env.execute("BFES Function Test");
    }
}
