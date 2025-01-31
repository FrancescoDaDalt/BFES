import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import org.apache.datasketches.hll.HllSketch;

import MurmurHash3;

public class BFES_Snapshot {
    public byte[] bytearray;
    public String name;

    public BFES_Snapshot (byte[] bytearray, String name) {
        this.bytearray = bytearray;
        this.name = name;
    }
}

public class BFES <K, V> {
    private final int[] axes;
    private final int num_counters;
    private final int num_nodes;
    private final int depth; 
    private final V[] sketch;
    private V total;
    private HllSketch hll;

    final int lgK = 10;
    
    public BFES(int num_counters, int depth) {
        this.axes = new int[depth];
        this.sketch = new float[num_counters];
        this.depth = depth;
        this.num_counters = num_counters;
        this.total = 0;
        this.axes[0] = num_counters - 2 * (depth - 1);
        for (int i = 1; i < depth; i++) {
            this.axes[i] = 2;
        }
        for (int i = 0; i < num_counters; i++) {
            this.sketch[i] = 0;
        }
        this.num_nodes = compute_num_nodes();
        this.hll = new HllSketch(lgK);
    }

    private int compute_num_nodes() {
		int length = 1;
		for (int a: this.axes) {
			length *= a;
		}
		return length;
	}

    public void insert(K raw_key, V value) {
        this.total += value;

        int int_key = raw_key.hashCode();

        byte[] keyAuxBytes = ByteBuffer.allocate(Long.BYTES)
                                        .order(ByteOrder.LITTLE_ENDIAN) 
                                        .putLong(int_key)
                                        .array();

        long[] extradest = MurmurHash3.murmurHash3_x64_128(keyAuxBytes, 42);
        long hashed_key = extradest[0];
        this.hll.insert(hashed_key);

        final int lower = 1;
		final int upper = this.depth;
		final int loopsize = upper - lower;
		int running_key = hashed_key;
		int offset = 0;
		final int axis = this.axes[0];
		final int hash = running_key % axis;
		this.sketch[hash + offset] += value;
		offset += axis;
		running_key /= axis;
		for (int idx = 0; idx < loopsize; idx++) {
			final int i = idx + lower;
			final int axis = this.axes[i];
			final int hash = running_key % axis;
			running_key /= axis;
			if (hash > 0) {
				this.sketch[hash - 1 + offset] += value;
			}
			offset += axis - 1;
		}
    }

    public String getDescriptor () {
        String out = "int" + "_" + total.getTypeName() + "_" + num_counters.toString() + "_" + depth.toString();
        return out;
    }

    public byte[] toByteArray () throws IOException {
        try (ByteArrayOutputStream byteArrayOutputStream = new ByteArrayOutputStream();
             ObjectOutputStream objectOutputStream = new ObjectOutputStream(byteArrayOutputStream)) {

            int num_distinct_est = (int) hll.getEstimate();
            objectOutputStream.writeObject(num_distinct_est);
            objectOutputStream.writeObject(total);
            objectOutputStream.writeObject(sketch);

            return byteArrayOutputStream.toByteArray();
        }
    }

    public BFES_Snapshot getSnapshot () {
        return BFES_Snapshot(toByteArray(), getDescriptor());
    } 
}

