package de.tu_darmstadt.rs.synbio.synthesis.enumeration;

import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import org.apache.commons.lang3.StringUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

public class BuildWorker implements Callable<List<TreeCircuit>> {

    private final EnumeratorFast enumerator;
    private final GateLibrary library;
    private final TreeCircuit startCircuit;
    private final int level;
    private final int maxDepth;
    private final int maxWeight;

    private final List<LogicType> gateTypes;

    List<TreeCircuit> results;

    public BuildWorker(EnumeratorFast enumerator, GateLibrary library, TreeCircuit startCircuit, List<LogicType> gateTypes, int level, int maxDepth, int maxWeight) {
        this.enumerator = enumerator;
        this.library = library;
        this.startCircuit = startCircuit;
        this.gateTypes = gateTypes;
        this.level = level;
        this.maxDepth = maxDepth;
        this.maxWeight = maxWeight;
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
        for (int i = 0; i < (int) Math.pow(gateTypes.size(), numberOfInputs); i++) {

            List<LogicType> row = generateRow(numberOfInputs, i);

            if (row == null)
                continue;

            TreeCircuit newCircuit = new TreeCircuit(circuit);
            boolean noRedundantInverters = newCircuit.addLevel(row);

            if (!noRedundantInverters)
                continue;

            if (enumerator.circuitAllowedByPreFilter(newCircuit))
                results.add(newCircuit);

            buildCircuits(newCircuit, level + 1);
        }
    }

    private List<LogicType> generateRow(int length, int index) {

        List<LogicType> row = new ArrayList<>();

        String combination = Integer.toString(index, gateTypes.size());
        combination = StringUtils.leftPad(combination, length, '0');

        for (int j = 0; j < combination.length(); j++) {
            row.add(j, gateTypes.get(Character.getNumericValue(combination.charAt(j))));
        }

        if (fulfillsMaxWeight(row) && isNotEmpty(row) && isCoveredByLibrary(row))
            if (!(length != 1 && (row.contains(LogicType.OUTPUT_BUFFER) || row.contains(LogicType.OUTPUT_OR2)))) // limit output gates to output row
                return row;

        return null;
    }

    private boolean fulfillsMaxWeight(List<LogicType> row) {

        int weight = 0;

        for (LogicType type : row) {
            weight += type.getWeight();
        }

        return weight <= maxWeight;
    }

    private boolean isNotEmpty(List<LogicType> row) {

        for (LogicType type : row) {
            if (type != LogicType.EMPTY)
                return true;
        }

        return false;
    }

    boolean isCoveredByLibrary(List<LogicType> gates) {

        // check if enough gate groups are available
        int numGates = 0;

        for (LogicType element : gates) {
            if (element != LogicType.EMPTY)
                numGates ++;
        }

        if (numGates > library.getGroups().size())
            return false;

        // check availability of gates of each type
        for (LogicType type : gateTypes) {

            if (type == LogicType.EMPTY)
                continue;

            int occurrences = 0;

            for (LogicType gateType : gates) {
                if (gateType == type)
                    occurrences++;
            }

            if (occurrences > library.getNumAvailableGates(type))
                return false;
        }

        return true;
    }
}
