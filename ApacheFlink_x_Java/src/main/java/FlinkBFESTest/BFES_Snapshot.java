package FlinkBFESTest;

import java.nio.charset.StandardCharsets;

public class BFES_Snapshot {
    public byte[] bytearray;
    public String name;

    public BFES_Snapshot(byte[] bytearray, String name) {
        this.bytearray = bytearray;
        this.name = name;
    }

    @Override
    public String toString() {
        return this.name + "_" + new String(this.bytearray, StandardCharsets.UTF_8);
    }

}