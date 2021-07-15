package de.tu_darmstadt.rs.synbio.mapping.compatibility;

import java.util.HashMap;
import java.util.Map;

public class CompatibilityMatrix<T> {

    Map<T, Map<T, Boolean>> compatibility;

    public CompatibilityMatrix() {
        compatibility = new HashMap<>();
    }

    public void addEntry(T source, T destination, Boolean isCompatible) {

        if (!compatibility.containsKey(source)) {
            compatibility.put(source, new HashMap<>());
        }
        compatibility.get(source).put(destination, isCompatible);
    }

    public boolean isCompatible(T source, T destination) {
        return compatibility.get(source).get(destination);
    };

}
