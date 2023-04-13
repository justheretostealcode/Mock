package de.tu_darmstadt.rs.synbio.synthesis.enumeration;

import de.tu_darmstadt.rs.synbio.common.LogicType;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

public class BuildWorker implements Callable<List<TreeCircuit>> {

    private final EnumeratorFast enumerator;
    private final TreeCircuit startCircuit;
    private final int level;
    private final int maxDepth;

    List<TreeCircuit> results;

    public BuildWorker(EnumeratorFast enumerator, TreeCircuit startCircuit, int level, int maxDepth) {
        this.enumerator = enumerator;
        this.startCircuit = startCircuit;
        this.level = level;
        this.maxDepth = maxDepth;
    }

    @Override
    public List<TreeCircuit> call() throws Exception {

        results = new ArrayList<>();

        buildCircuits(startCircuit, level);

        return results;
    }

    private void buildCircuits(TreeCircuit circuit, int level) {

        // if max depth reached --> abort
        if (level >= maxDepth)
            return;

        int numberOfInputs = circuit.getNumOpenInputs(false);

        // iterate over rows with corresponding number of gates/entries
        for (int i = 0; i < enumerator.combinations.get(numberOfInputs - 1).size(); i++) {

            if (enumerator.combinations.get(numberOfInputs - 1).get(i).contains(LogicType.OUTPUT_OR2)) // limit or gate to output row
                continue;

            TreeCircuit newCircuit = new TreeCircuit(circuit);
            boolean noRedundantInverters = newCircuit.addLevel(enumerator.combinations.get(numberOfInputs - 1).get(i));

            if (!noRedundantInverters)
                continue;

            if (enumerator.circuitAllowedByPreFilter(newCircuit))
                results.add(newCircuit);

            buildCircuits(newCircuit, level + 1);
        }
    }
}
