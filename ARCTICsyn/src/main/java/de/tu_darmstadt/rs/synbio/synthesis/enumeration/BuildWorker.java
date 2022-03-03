package de.tu_darmstadt.rs.synbio.synthesis.enumeration;

import de.tu_darmstadt.rs.synbio.common.LogicType;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

public class BuildWorker implements Callable<List<PrimitiveCircuit>> {

    private final EnumeratorFast enumerator;
    private final PrimitiveCircuit startCircuit;
    private final int level;
    private final int maxDepth;

    List<PrimitiveCircuit> results;

    public BuildWorker(EnumeratorFast enumerator, PrimitiveCircuit startCircuit, int level, int maxDepth) {
        this.enumerator = enumerator;
        this.startCircuit = startCircuit;
        this.level = level;
        this.maxDepth = maxDepth;
    }

    @Override
    public List<PrimitiveCircuit> call() throws Exception {

        results = new ArrayList<>();

        buildCircuits(startCircuit, level);

        return results;
    }

    private void buildCircuits(PrimitiveCircuit circuit, int level) {

        // if max depth reached --> abort
        if (level >= maxDepth)
            return;

        // if circuit is empty --> build start circuits and recurse
        if (level == 0) {

            for (int i = 0; i < enumerator.combinations.get(0).size(); i++) {
                PrimitiveCircuit newCircuit = new PrimitiveCircuit(circuit);
                newCircuit.insertEntry(0, i);
                buildCircuits(newCircuit, 1);
            }

            // if circuit is not empty --> extend by next row
        } else {

            // get number of inputs of lower level
            PrimitiveCircuit.Entry entry = circuit.getEntry(level - 1);
            int numberOfInputs = enumerator.getNumberOfInputs(enumerator.mapEntryToRow(entry));

            // iterate over rows with corresponding number of gates/entries
            for (int i = 0; i < enumerator.combinations.get(numberOfInputs - 1).size(); i++) {

                if (enumerator.combinations.get(numberOfInputs - 1).get(i).contains(LogicType.OUTPUT_OR2)) // limit or gate to output row
                    continue;

                PrimitiveCircuit newCircuit = new PrimitiveCircuit(circuit);
                newCircuit.addEntry(level, numberOfInputs - 1, i);

                if (enumerator.circuitAllowedByPreFilter(newCircuit))
                    results.add(newCircuit);

                buildCircuits(newCircuit, level + 1);
            }
        }
    }
}
