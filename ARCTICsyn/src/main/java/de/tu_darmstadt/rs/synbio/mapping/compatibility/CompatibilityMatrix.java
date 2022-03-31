package de.tu_darmstadt.rs.synbio.mapping.compatibility;

import java.util.HashMap;
import java.util.Map;

public class CompatibilityMatrix<T> {

    /* key order: source, destination, 2nd source */
    Map<T, Map<T, Map<T, Boolean>>> compatibility;

    public CompatibilityMatrix() {
        compatibility = new HashMap<>();
    }

    public void addEntry(T source, T destination, T secondSource, Boolean isCompatible) {

        if (!compatibility.containsKey(source)) {
            compatibility.put(source, new HashMap<>());
        }

        if (!compatibility.get(source).containsKey(destination)) {
            compatibility.get(source).put(destination, new HashMap<>());
        }

        compatibility.get(source).get(destination).put(secondSource, isCompatible);
    }

    public boolean isCompatible(T source, T destination, T secondSource) {
        return compatibility.get(source).get(destination).get(secondSource);
    };

}
