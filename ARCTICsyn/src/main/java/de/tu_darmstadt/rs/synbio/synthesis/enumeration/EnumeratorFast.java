package de.tu_darmstadt.rs.synbio.synthesis.enumeration;

import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.common.TruthTable;
import de.tu_darmstadt.rs.synbio.common.circuit.*;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.synthesis.SynthesisConfiguration;
import de.tu_darmstadt.rs.synbio.synthesis.util.ExpressionParser;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.jgrapht.Graphs;
import org.logicng.formulas.Variable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.Collectors;

public class EnumeratorFast {

    private static final Logger logger = LoggerFactory.getLogger(EnumeratorFast.class);

    /* Configuration */
    private final SynthesisConfiguration synConfig;
    private TruthTable targetTruthTable;
    private final GateLibrary gateLib;
    private final int feasibility;

    private int currentMinSize;

    private final int numGroupsInLib;
    private final Circuit templateCircuit;
    private final List<LogicType> gateTypes;

    private HashMap<TruthTable, Circuit> resultCircuits;

    /* enumeration buffers */
    private final List<String> gateInputNames;
    private final List<String> inputVars;
    private final List<String> intermediateVariables;
    private final List<InputGate> inputGates;

    /* constructor for guided enumeration */

    public EnumeratorFast(GateLibrary gateLib, TruthTable targetTruthTable, SynthesisConfiguration synConfig) {

        this(gateLib, targetTruthTable.getSupportSize(), synConfig);
        currentMinSize = Integer.MAX_VALUE - synConfig.getWeightRelaxation();

        this.targetTruthTable = targetTruthTable;
    }

    /* constructor for free enumeration */

    public EnumeratorFast(GateLibrary gateLib, int feasibility, SynthesisConfiguration synConfig) {

        this.gateLib = gateLib;
        this.feasibility = feasibility;
        this.synConfig = synConfig;

        this.gateTypes = new ArrayList<>();
        this.gateTypes.add(LogicType.EMPTY);
        this.gateTypes.addAll(gateLib.getGateTypes());

        this.numGroupsInLib = gateLib.getGroups().size();

        /* init final fields */
        templateCircuit = new Circuit();
        Gate output = new OutputGate("O");
        templateCircuit.addVertex(output);

        gateInputNames = new ArrayList<>();
        gateInputNames.add("x");
        gateInputNames.add("y");

        inputVars = new ArrayList<>();
        inputGates = new ArrayList<>();
        char varName = 'a';
        for (int i = 0; i < feasibility; i ++) {
            inputVars.add(""+varName);
            inputGates.add(new InputGate(ExpressionParser.parse(Character.toString(varName)), Character.toString(varName)));
            varName ++;
        }

        intermediateVariables = new ArrayList<>();
        for (int i = 0; i < 100; i++) {
            intermediateVariables.add("m" + String.format("%03d", i)); //TODO: make dynamically grow
        }
    }

    /* main algorithm */

    List<List<List<LogicType>>> combinations;

    public void enumerate() {

        logger.info("setting up enumeration...");

        int availableProcessors = Runtime.getRuntime().availableProcessors() - 1;

        List<PrimitiveCircuit> rawCircuits = new ArrayList<>();
        resultCircuits = new HashMap<>();

        OptionalInt maxGateFeasibility = gateTypes.stream().mapToInt(LogicType::getNumInputs).max();
        if (maxGateFeasibility.isEmpty())
            return;

        /* generate valid combinations of gates */

        // max number of gates in uppermost level
        int maxRowLength = (int) Math.pow(maxGateFeasibility.getAsInt(), synConfig.getMaxDepth() - 1);

        // lists containing combinations
        combinations = new ArrayList<>(maxRowLength);

        for (int i = 0; i < maxRowLength; i ++) {
            combinations.add(i, new ArrayList<>());
        }

        // generate rows (combinations of gates and empty slots)
        for (int rowLength = 1; rowLength <= maxRowLength; rowLength ++) {

            List<List<LogicType>> lengthCombinations = combinations.get(rowLength - 1);

            for (int i = 0; i < (int) Math.pow(gateTypes.size(), rowLength); i ++) {

                String combination = Integer.toString(i, gateTypes.size());
                combination = StringUtils.leftPad(combination, rowLength, '0');

                ArrayList<LogicType> row = new ArrayList<>(rowLength);

                for (int j = 0; j < combination.length(); j++) {
                    row.add(j, gateTypes.get(Character.getNumericValue(combination.charAt(j))));
                }

                if (isCoveredByLibrary(row) && isNotEmpty(row))
                    if (!(rowLength != 1 && row.contains(LogicType.OR2))) // limit OR gate to output row
                        lengthCombinations.add(row);
            }
        }

        /* call recursive circuit build function */

        List<BuildWorker> buildWorkers = new ArrayList<>();

        PrimitiveCircuit emptyCircuit = new PrimitiveCircuit(synConfig.getMaxDepth());

        /* initialize build workers with differing start gates */
        for (int i = 0; i < combinations.get(0).size(); i++) {
            PrimitiveCircuit newCircuit = new PrimitiveCircuit(emptyCircuit);
            newCircuit.insertEntry(0, i);
            buildWorkers.add(new BuildWorker(this, newCircuit, 1, synConfig.getMaxDepth()));
        }

        int numThreads = synConfig.getSynLimitThreadsNum() == 0 ? availableProcessors : Math.min(synConfig.getSynLimitThreadsNum(), availableProcessors);
        logger.info("possible parallelism during building: " + buildWorkers.size());
        logger.info("number of threads resulting from hardware and synthesis settings: " + numThreads);

        ExecutorService executor = Executors.newFixedThreadPool(numThreads);
        List<Future<List<PrimitiveCircuit>>> buildResults = Collections.emptyList();

        try {
            buildResults = executor.invokeAll(buildWorkers);
        } catch (InterruptedException e) {
            logger.error(e.getMessage());
        }

        for (Future<List<PrimitiveCircuit>> resultList : buildResults) {

            try {
                rawCircuits.addAll(resultList.get());
            }  catch (Exception e) {
                logger.error(e.getMessage());
            }
        }

        buildResults.clear();
        executor.shutdownNow();

        logger.info("found " + rawCircuits.size() + " pre-filtered circuits.");

        /* filter structurally equivalent circuits */

        rawCircuits = filterRedundantCircuits(rawCircuits);

        logger.info("found " + rawCircuits.size() + " structurally different circuits.");

        /* evaluate circuits in multiple threads */

        List<EnumerationWorker> workers = new ArrayList<>();

        // generate input mapping (primary inputs --> unbounded circuit inputs)
        for (PrimitiveCircuit circuit : rawCircuits){
            workers.add(new EnumerationWorker(this, circuit, feasibility, targetTruthTable));
        }

        executor = Executors.newFixedThreadPool(numThreads);

        List<Future<List<EnumerationResult>>> results = Collections.emptyList();

        try {
            results = executor.invokeAll(workers);
        } catch (InterruptedException e) {
            logger.error(e.getMessage());
        }

        /* build graphs of result circuits and perform final checks */

        for (Future<List<EnumerationResult>> resultList : results) {

            try {
                for (EnumerationResult result : resultList.get()) {
                    buildCircuitGraph(result.getCircuit(), result.getInputMapping());
                }
            }  catch (Exception e) {
                logger.error(e.getMessage());
            }
        }

        executor.shutdownNow();
    }

    /*private void buildCircuits(PrimitiveCircuit circuit, int level) {

        // if max depth reached --> abort
        if (level >= synConfig.getMaxDepth())
            return;

        // if circuit is empty --> build start circuits and recurse
        if (level == 0) {

            for (int i = 0; i < combinations.get(0).size(); i++) {
                PrimitiveCircuit newCircuit = new PrimitiveCircuit(circuit);
                newCircuit.insertEntry(0, i);
                buildCircuits(newCircuit, 1);
            }

            // if circuit is not empty --> extend by next row
        } else {

            // get number of inputs of lower level
            PrimitiveCircuit.Entry entry = circuit.getEntry(level - 1);
            int numberOfInputs = getNumberOfInputs(mapEntryToRow(entry));

            // iterate over rows with corresponding number of gates/entries
            for (int i = 0; i < combinations.get(numberOfInputs - 1).size(); i++) {

                if (combinations.get(numberOfInputs - 1).get(i).contains(LogicType.OR2)) // limit OR gate to output row
                    continue;

                PrimitiveCircuit newCircuit = new PrimitiveCircuit(circuit);
                newCircuit.addEntry(level, numberOfInputs - 1, i);

                if (circuitAllowedByPreFilter(newCircuit))
                    intermediateCircuits.add(newCircuit);

                buildCircuits(newCircuit, level + 1);
            }
        }
    }*/

    boolean circuitAllowedByPreFilter(PrimitiveCircuit circuit) {

        // check if circuit is implementable with library
        if (!isCoveredByLibrary(circuit.getList().stream().map(this::mapEntryToRow).flatMap(Collection::stream).collect(Collectors.toCollection(ArrayList::new))))
            return false;

        // check is feasibility is met
        if (getNumberOfUnboundInputs(circuit) < feasibility)
            return false;

        // check circuit weight. if it exceeds maximum --> continue
        if (getCircuitWeight(circuit) > synConfig.getMaxWeight())
            return false;

        // check if circuit contains redundant inverters
        if (hasRedundantInverters(circuit))
            return false;

        return true;
    }

    private void buildCircuitGraph(PrimitiveCircuit candidate, String inputMapping) {

        /* build circuit */
        Circuit circuit = new Circuit();
        Graphs.addGraph(circuit, templateCircuit);

        // add first logic gate after output
        LogicType type = mapEntryToRow(candidate.getEntry(0)).get(0);
        int gateCounter = 0;
        Gate newGate = new LogicGate(type.name() + "_" + gateCounter, type);
        circuit.addVertex(newGate);
        circuit.addEdge(newGate, circuit.getOutputBuffer(), new Wire(circuit.getOutputBuffer().getExpression().variables().first()));

        // add rows of gates and wire them
        List<Gate> prevRow = new ArrayList<>();
        List<Gate> currentRow = new ArrayList<>();
        prevRow.add(newGate);

        int inputCount = 0;

        for (int level = 1; level <= candidate.getDepth(); level ++) {

            int column = 0;

            for (Gate anchorGate : prevRow) {

                for (Variable anchorVariable : anchorGate.getExpression().variables()) {

                    if (level < candidate.getDepth())
                        type = mapEntryToRow(candidate.getEntry(level)).get(column);
                    else
                        type = LogicType.EMPTY;

                    if (type != LogicType.EMPTY) {

                        gateCounter++;
                        newGate = new LogicGate(type.name() + "_" + gateCounter, type);
                        circuit.addVertex(newGate);
                        circuit.addEdge(newGate, anchorGate, new Wire(anchorVariable));

                        currentRow.add(newGate);
                    } else {

                        InputGate inputGate = inputGates.get(Integer.parseInt(String.valueOf(inputMapping.charAt(inputCount))));

                        if (!circuit.containsVertex(inputGate))
                            circuit.addVertex(inputGate);

                        if (circuit.containsEdge(inputGate, anchorGate))
                            return;

                        circuit.addEdge(inputGate, anchorGate, new Wire(anchorVariable));

                        inputCount ++;
                    }

                    column++;
                }
            }

            prevRow = new ArrayList<>(currentRow);
            currentRow.clear();
        }

        // remove redundant gates
        if (!circuit.removeRedundantGates())
            return;

        // compare to current min size
        int weight = circuit.getWeight();
        if (weight > currentMinSize + synConfig.getWeightRelaxation())
            return;

        TruthTable truthTable = new TruthTable(circuit.getExpression());

        if (!truthTable.equalsLogically(targetTruthTable)) {
            circuit.print(new File("test.dot"));
            evaluatePrimitiveCircuit(candidate, inputMapping);
            logger.error("Synthesis error: " + truthTable.toString());
        }

        // filter out structurally equivalent circuits
        for (TruthTable libTruthTable : resultCircuits.keySet()) {

            if (targetTruthTable != null || libTruthTable.equalsLogically(truthTable)) {

                if (circuit.isEquivalent(resultCircuits.get(libTruthTable))) {
                    return;
                }
            }
        }

        if (weight < currentMinSize) {
            currentMinSize = weight;

            // update result circuits
            resultCircuits.entrySet().removeIf(e -> e.getValue().getWeight() > (currentMinSize + synConfig.getWeightRelaxation()));
        }

        resultCircuits.put(truthTable, circuit);
    }

    /* helper functions for handling of primitive circuits */

    List<LogicType> mapEntryToRow(PrimitiveCircuit.Entry entry) {
        return combinations.get(entry.combId).get(entry.index);
    }

    private boolean isCoveredByLibrary(List<LogicType> row) {

        // check if enough gate groups are available
        int numGates = 0;

        for (LogicType element : row) {
            if (element != LogicType.EMPTY)
                numGates ++;
        }

        if (numGates > numGroupsInLib)
            return false;

        // check availability of gates of each type
        for (LogicType type : gateTypes) {

            if (type == LogicType.EMPTY)
                continue;

            int occurrences = 0;

            for (LogicType gateType : row) {
                if (gateType == type)
                    occurrences++;
            }

            if (occurrences > gateLib.getNumAvailableGates(type))
                return false;
        }

        return true;
    }

    int getNumberOfInputs(List<LogicType> row) {

        int numInputs = 0;

        for (LogicType gateType : row) {
            if (gateType != LogicType.EMPTY)
                numInputs += gateType.getNumInputs();
        }

        return numInputs;
    }

    int getNumberOfUnboundInputs(PrimitiveCircuit circuit) {

        // get total number of gate inputs
        int totalInputs = circuit.getList().stream().map(this::mapEntryToRow).mapToInt(this::getNumberOfInputs).sum();

        // get number of gates above level 0 (bound inputs)
        int boundInputs = 0;

        for (int level = 1; level < circuit.getDepth(); level ++) {

            PrimitiveCircuit.Entry entry = circuit.getEntry(level);

            for(LogicType gateType : mapEntryToRow(entry)) {
                if (gateType != LogicType.EMPTY)
                    boundInputs ++;
            }
        }

        return totalInputs - boundInputs;
    }

    private int getCircuitWeight(PrimitiveCircuit circuit) {

        int weight = 0;

        for (PrimitiveCircuit.Entry entry : circuit.getList()) {

            for (LogicType gateType : mapEntryToRow(entry)) {
                weight += gateType.getWeight();
            }
        }

        return weight;
    }

    private boolean hasRedundantInverters(PrimitiveCircuit circuit) {

        for (int level = 0; level < (circuit.getDepth() - 1); level ++) {

            int upperColumn = 0;

            PrimitiveCircuit.Entry entry = circuit.getEntry(level);

            for (LogicType gateType : mapEntryToRow(entry)) {

                if (gateType == LogicType.NOT) {

                    PrimitiveCircuit.Entry upperEntry = circuit.getEntry(level + 1);
                    List<LogicType> upperRow = mapEntryToRow(upperEntry);

                    if (upperRow.get(upperColumn) == LogicType.NOT)
                        return true;
                }

                upperColumn += gateType.getNumInputs();
            }
        }
        return false;
    }

    private boolean isNotEmpty(List<LogicType> row) {

        for (LogicType type : row) {
            if (type != LogicType.EMPTY)
                return true;
        }

        return false;
    }

    String evaluatePrimitiveCircuit(PrimitiveCircuit circuit, String inputMapping) {

        int varCount = 0;       /* counter for instantiating intermediate variables */
        Integer gateCount = 0;  /* counter for identifying gates */

        /* list of currently unconnected variables */
        List<String> currentVars = new ArrayList<>();

        /* map of variables to gate IDs */
        Map<String, Integer> varsToGates = new HashMap<>();

        /* list of variables to be connected to primary inputs */
        List<String> finalVars = new ArrayList<>();

        /* temporary list of replaced variables */
        List<String> replacedVars = new ArrayList<>();

        String varName;
        String expression = "";

        for (int row = 0; row < circuit.getDepth(); row ++) {

            List<LogicType> currentRow = mapEntryToRow(circuit.getEntry(row));

            /* update current vars */
            currentVars.removeAll(replacedVars);
            currentVars.removeAll(finalVars);
            Collections.sort(currentVars);

            replacedVars.clear();

            for (int column = 0; column < currentRow.size(); column ++) {

                /* empty type means gate input will be connected to circuit input */
                if (row > 0 && currentRow.get(column) == LogicType.EMPTY) {
                    finalVars.add(currentVars.get(column));
                    continue;
                }

                /* insert expression of gate */
                String gateExpression = currentRow.get(column).getExpressionString();

                if (expression.isBlank()) {   /* output gate: add expression */
                    expression = gateExpression;
                } else {                    /* every other gate: replace variable by expression */
                    String varToReplace = currentVars.get(column);
                    expression = expression.replaceAll(varToReplace, gateExpression);
                    replacedVars.add(varToReplace);
                }

                /* replace gate inputs by intermediate variables */
                for (String inputName : gateInputNames) {
                    if (expression.contains(inputName)) {
                        varName = intermediateVariables.get(varCount);
                        expression = expression.replaceAll(inputName, varName);
                        currentVars.add(varName);
                        varsToGates.put(varName, gateCount);
                        varCount++;
                    }
                }

                gateCount ++;
            }
        }

        /* update current variables. all left current and final variables are connected to inputs */
        currentVars.removeAll(finalVars);
        currentVars.removeAll(replacedVars);

        List<String> allVars = new ArrayList<>();
        allVars.addAll(currentVars);
        allVars.addAll(finalVars);
        Collections.sort(allVars);

        /* map of gate IDs to substituted inputs */
        Map<Integer, String> gatesToSubstitutions = new HashMap<>(); // TODO: for 3 or more input gates: replace values by sets of formulas

        int substCount = 0;
        for (String var : allVars) {

            /* get input to substitute according to input mapping */
            String substitution = inputVars.get(Character.getNumericValue(inputMapping.charAt(substCount)));

            /* check if input is already connected to other variable of gate and abort */
            if (gatesToSubstitutions.containsKey(varsToGates.get(var)) && gatesToSubstitutions.get(varsToGates.get(var)).equals(substitution)) {
                return null;
            }

            if (varsToGates.containsKey(var))
                gatesToSubstitutions.put(varsToGates.get(var), substitution);

            expression = expression.replaceAll(var, substitution);
            substCount ++;
        }

        return expression;
    }

    /* structural equivalence check */

    private List<HashMap<Coordinates, List<List<LogicType>>>> pathDB;

    private List<PrimitiveCircuit> filterRedundantCircuits(List<PrimitiveCircuit> inputCircuits) {

        // calculate all paths
        pathDB = new ArrayList<>();

        for (PrimitiveCircuit circuit : inputCircuits) {
            pathDB.add(getAllPaths(circuit));
        }

        // filter structurally equivalent circuits
        List<Integer> resultCircuitIndices = new ArrayList<>();
        resultCircuitIndices.add(0);

        for (int i = 0; i < inputCircuits.size(); i ++) {

            boolean isUnique = true;

            for (Integer cmpIndex : resultCircuitIndices) {

                if (structurallyEquivalent(inputCircuits.get(i), i, inputCircuits.get(cmpIndex), cmpIndex)) {
                    isUnique = false;
                    break;
                }
            }

            if (isUnique) {
                resultCircuitIndices.add(i);
                //logger.info("circuit " + i + " added");
            }
        }

        // build list of filtered circuits
        List<PrimitiveCircuit> filteredCircuits = new ArrayList<>();

        for (Integer index : resultCircuitIndices) {
            filteredCircuits.add(inputCircuits.get(index));
        }

        return filteredCircuits;
    }

    private boolean structurallyEquivalent(PrimitiveCircuit circuit1, int circuit1Index, PrimitiveCircuit circuit2, int circuit2Index) {

        if (circuit1.getDepth() != circuit2.getDepth())
            return false;

        if (getCircuitWeight(circuit1) != getCircuitWeight(circuit2))
            return false;

        // get paths
        HashMap<Coordinates, List<List<LogicType>>> paths1 = pathDB.get(circuit1Index);
        HashMap<Coordinates, List<List<LogicType>>> paths2 = pathDB.get(circuit2Index);

        // create mark list
        List<List<Boolean>> circuit2Marks = new ArrayList<>();

        for (PrimitiveCircuit.Entry entry : circuit2.getList()) {

            ArrayList<Boolean> markRow = new ArrayList<>();
            circuit2Marks.add(markRow);

            for (LogicType ignored : mapEntryToRow(entry)) {
                markRow.add(false);
            }
        }

        // iterate over elements of circuit 1 and search equivalents in circuit 2
        for (int level = 0; level < circuit1.getDepth(); level ++) {

            if (mapEntryToRow(circuit1.getEntry(level)).size() != mapEntryToRow(circuit2.getEntry(level)).size())
                return false;

            int numColumns = mapEntryToRow(circuit1.getEntry(level)).size();

            for (int column = 0; column < numColumns; column ++) {

                LogicType element = mapEntryToRow(circuit1.getEntry(level)).get(column);
                Coordinates coords1 = new Coordinates(level, column);
                List<List<LogicType>> element1Paths = paths1.get(coords1);

                boolean equivalentFound = false;

                // search equivalent gate in circuit 2
                for (int column2 = 0; column2 < numColumns; column2 ++) {

                    // if element is same type and is not yet marked
                    if (mapEntryToRow(circuit2.getEntry(level)).get(column2) == element && !circuit2Marks.get(level).get(column2)) {

                        Coordinates coords2 = new Coordinates(level, column2);
                        List<List<LogicType>> element2Paths = paths2.get(coords2);

                        if (element2Paths.containsAll(element1Paths) && element1Paths.containsAll(element2Paths)) {
                            circuit2Marks.get(level).set(column2, true);
                            equivalentFound = true;
                            break;
                        }
                    }
                }

                if (!equivalentFound)
                    return false;
            }
        }

        // check if every element of circuit 2 is marked
        for (List<Boolean> level : circuit2Marks) {
            for (Boolean marked : level) {
                if (!marked)
                    return false;
            }
        }

        return true;
    }

    private HashMap<Coordinates, List<List<LogicType>>> getAllPaths(PrimitiveCircuit circuit) {

        HashMap<Coordinates, List<List<LogicType>>> paths = new HashMap<>();

        for (int level = circuit.getDepth() - 1; level >= 0; level --) {

            int rowLength = mapEntryToRow(circuit.getEntry(level)).size();

            for (int column = 0; column < rowLength; column ++) {

                Coordinates currentNode = new Coordinates(level, column);
                List<Coordinates> connectedNodes = getInputNodes(circuit, currentNode);

                paths.putIfAbsent(currentNode, new ArrayList<>());

                if (connectedNodes.isEmpty()) {

                    ArrayList<LogicType> path = new ArrayList<>();
                    path.add(mapEntryToRow(circuit.getEntry(level)).get(column));
                    paths.get(currentNode).add(path);
                }

                for (Coordinates node : connectedNodes) {

                    for (List<LogicType> path : paths.get(node)) {

                        ArrayList<LogicType> extendedPath = new ArrayList<>(path);
                        extendedPath.add(mapEntryToRow(circuit.getEntry(level)).get(column));
                        paths.get(currentNode).add(extendedPath);
                    }

                }
            }
        }

        return paths;
    }

    private List<Coordinates> getInputNodes(PrimitiveCircuit circuit, Coordinates coords) {

        List<Coordinates> inputNodes = new ArrayList<>();

        LogicType element = mapEntryToRow(circuit.getEntry(coords.level)).get(coords.column);
        int numInputs = element.getNumInputs();

        // if element is on uppermost level --> return empty list
        if (coords.level == (circuit.getDepth() - 1)) {
            return inputNodes;
        }

        int upperInput = 0;

        for (int currentColumn = 0; currentColumn < coords.column; currentColumn ++) {

            LogicType currentElement = mapEntryToRow(circuit.getEntry(coords.level)).get(currentColumn);

            if (currentElement != LogicType.EMPTY)
                upperInput += currentElement.getNumInputs();
        }

        for (int input = upperInput; input < (upperInput + numInputs); input ++) {

            inputNodes.add(new Coordinates(coords.level + 1, input));
        }

        return inputNodes;

    }

    static private class Coordinates {

        public final Integer level;
        public final Integer column;

        public Coordinates(Integer level, Integer column) {
            this.level = level;
            this.column = column;
        }

        @Override
        public boolean equals(Object o) {
            if (!(o instanceof Coordinates)) {
                return false;
            }
            Coordinates cmp  = (Coordinates) o;
            EqualsBuilder builder = new EqualsBuilder();
            builder.append(level, cmp.level);
            builder.append(column, cmp.column);
            return builder.isEquals();
        }

        @Override
        public int hashCode() {
            HashCodeBuilder builder = new HashCodeBuilder();
            builder.append(level);
            builder.append(column);
            return builder.hashCode();
        }
    }

    /* getter */

    public HashMap<TruthTable, Circuit> getResultCircuits() {
        return resultCircuits;
    }
}
