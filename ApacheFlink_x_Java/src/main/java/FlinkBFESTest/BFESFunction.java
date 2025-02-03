package FlinkBFESTest;

import org.apache.flink.api.common.functions.RichMapFunction;
import org.apache.flink.streaming.api.functions.ProcessFunction;
import org.apache.flink.util.Collector;
import org.apache.flink.api.java.tuple.Tuple2;
import org.apache.flink.configuration.Configuration;

// It takes in Key-value pairs and returns at the end the sketch datastructure (bytearray containing the all relevant sketch data and a string describing metadata) which must be processed in C++
public class BFESFunction<K, V extends Number> extends ProcessFunction<Tuple2<K, V>, BFES_Snapshot> {

    private BFES<K, V> bfes;

    private int num_counters;
    private int depth;

    private boolean timerSet = false;

    public BFESFunction(int num_counters, int depth) {
        this.num_counters = num_counters;
        this.depth = depth;
    }

    @Override
    public void open(Configuration parameters) {
        this.bfes = new BFES(num_counters, depth);
    }

    @Override
    public void processElement(Tuple2<K, V> key_value, Context ctx, Collector<BFES_Snapshot> out) throws Exception {
        K key = key_value.f0;
        V value = key_value.f1;
        bfes.insert(key, value);
        if (!this.timerSet) {
            long tm = ctx.timerService().currentProcessingTime() + 10;
            ctx.timerService().registerProcessingTimeTimer(tm);
            this.timerSet = true;
        }
    }

    @Override
    public void onTimer(long timestamp, OnTimerContext ctx, Collector<BFES_Snapshot> out) throws Exception {
        BFES_Snapshot snapshot = bfes.getSnapshot();
        snapshot.name = String.valueOf(timestamp) + "_" + snapshot.name;
        out.collect(snapshot);
        this.timerSet = false;
        this.bfes.clear();
    }
}
