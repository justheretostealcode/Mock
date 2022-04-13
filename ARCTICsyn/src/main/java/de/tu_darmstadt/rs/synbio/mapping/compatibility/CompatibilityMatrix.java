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
    }

    public String getCompatibilityInfo() {

        long totalPairs = 0;
        long totalTriples = 0;

        long compatiblePairs = 0;
        long compatibleTriples = 0;

        for (T source : compatibility.keySet()) {
            for (T dest : compatibility.get(source).keySet()) {

                if (dest.equals("input_1") || dest.equals("input_2") ||dest.equals("input_3"))
                    continue;

                if (source.equals(dest))
                    continue;

                totalPairs ++;

                if (compatibility.get(source).get(dest).get(null))
                    compatiblePairs ++;

                for (T second : compatibility.get(source).get(dest).keySet()) {

                    if (source.equals(second) || dest.equals(second))
                        continue;

                    totalTriples ++;

                    if (compatibility.get(source).get(dest).get(second))
                        compatibleTriples ++;

                }
            }
        }

        return "pairs: " + (double) compatiblePairs/totalPairs + " triple: " + (double) compatibleTriples/totalTriples;
    }
}
