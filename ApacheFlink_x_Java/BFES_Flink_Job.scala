import org.apache.flink.api.scala.StreamExecutionEnvironment;
import org.apache.flink.streaming.api.scala.DataStream;
import org.apache.flink.streaming.api.scala._

object BFESJob {

    def main(args: Array[String]): Unit = {
        // Set up the Flink execution environment
        val env = StreamExecutionEnvironment.getExecutionEnvironment

        // Example input data stream
        val inputStream: DataStream[String] = env.fromElements("apple", "banana", "apple", "apple", "orange")

        // Apply the CountMinSketch function
        val resultStream = inputStream
            .keyBy(value => value) // Key by the value (for per-key frequency tracking)
            .process(new CountMinSketchFunction(100, 5)) // Create and apply the Count-Min Sketch function

        // Print the result to stdout
        resultStream.print()

        // Execute the Flink job
        env.execute("CountMinSketch Job")
    }
}
