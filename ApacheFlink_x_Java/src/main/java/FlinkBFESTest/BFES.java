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

import java.util.ArrayList;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import org.apache.datasketches.hll.HllSketch;

public class BFES<K, V extends Number> {
    private final int[] axes;
    private final int num_counters;
    private final int num_nodes;
    private final int depth;
    private final ArrayList<V> sketch;
    private V total;
    private HllSketch hll;

    final int lgK = 10;

    public BFES(int num_counters, int depth) {
        this.axes = new int[depth];
        this.sketch = new ArrayList<V>(num_counters);
        this.depth = depth;
        this.num_counters = num_counters;
        this.total = (V) (Number) 0;
        this.axes[0] = num_counters - 2 * (depth - 1);
        for (int i = 1; i < depth; i++) {
            this.axes[i] = 2;
        }
        for (int i = 0; i < num_counters; i++) {
            this.sketch.add((V) (Number) 0.0);
        }
        this.num_nodes = compute_num_nodes();
        this.hll = new HllSketch(lgK);
    }

    private int compute_num_nodes() {
        int length = 1;
        for (int a : this.axes) {
            length *= a;
        }
        return length;
    }

    public void clear () {
        for (int i = 0; i < num_counters; i++) {
            this.sketch.set(i, (V) (Number) 0.0);
        }
        this.total = (V) (Number) 0.0;
        this.hll = new HllSketch(lgK);
    }

    public void insert(K raw_key, V value) {
        this.total = (V) (Number) (this.total.doubleValue() + value.doubleValue());

        int int_key = raw_key.hashCode();

        byte[] keyAuxBytes = ByteBuffer.allocate(Long.BYTES)
                .order(ByteOrder.LITTLE_ENDIAN)
                .putLong(int_key)
                .array();

        long[] extradest = MurmurHash3.murmurHash3_x64_128(keyAuxBytes, 42);
        long hashed_key = extradest[0];
        this.hll.update((int) hashed_key);

        final int lower = 1;
        final int upper = this.depth;
        final int loopsize = upper - lower;
        int running_key = (int) hashed_key;
        int offset = 0;
        int axis = this.axes[0];
        int hash = (((running_key % axis) + axis) % axis);
        this.sketch.set(hash + offset,
                (V) (Number) (this.sketch.get(hash + offset).doubleValue() + value.doubleValue()));
        offset += axis;
        running_key /= axis;
        for (int idx = 0; idx < loopsize; idx++) {
            final int i = idx + lower;
            axis = this.axes[i];
            hash = running_key % axis;
            running_key /= axis;
            if (hash > 0) {
                this.sketch.set(hash - 1 + offset,
                        (V) (Number) (this.sketch.get(hash - 1 + offset).doubleValue() + value.doubleValue()));
            }
            offset += axis - 1;
        }
    }

    public String getDescriptor() {
        String out = "int" + "_" + total.getClass().getTypeName() + "_" + String.valueOf(num_counters) + "_"
                + String.valueOf(depth);
        return out;
    }

    public byte[] toByteArray() throws IOException {
        try (ByteArrayOutputStream byteArrayOutputStream = new ByteArrayOutputStream();
                ObjectOutputStream objectOutputStream = new ObjectOutputStream(byteArrayOutputStream)) {

            int num_distinct_est = (int) hll.getEstimate();
            objectOutputStream.writeObject(num_distinct_est);
            objectOutputStream.writeObject(total);
            objectOutputStream.writeObject(sketch);

            return byteArrayOutputStream.toByteArray();
        }
    }

    public BFES_Snapshot getSnapshot() throws IOException {
        return new BFES_Snapshot(toByteArray(), getDescriptor());
    }
}
