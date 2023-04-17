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
import org.logicng.formulas.Formula;
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
    private boolean isGuided = false;
    private final GateLibrary gateLib;
    private final int feasibility;

    private int currentMinSize;

    private final int numGroupsInLib;
    final List<LogicType> gateTypes;

    private HashMap<TruthTable, Set<Circuit>> resultCircuits;

    /* enumeration buffers */
    private final List<String> gateInputNames;
    private final List<String> inputVars;
    private final List<String> intermediateVariables;
    private final List<Gate> inputGates;

    /* constructor for guided enumeration */

    public EnumeratorFast(GateLibrary gateLib, TruthTable targetTruthTable, SynthesisConfiguration synConfig) {

        this(gateLib, targetTruthTable.getSupportSize(), synConfig);
        currentMinSize = Integer.MAX_VALUE - synConfig.getWeightRelaxation();

        this.targetTruthTable = targetTruthTable;
        this.isGuided = true;
    }

    /* constructor for free enumeration */

    public EnumeratorFast(GateLibrary gateLib, int feasibility, SynthesisConfiguration synConfig) {

        this.gateLib = gateLib;
        this.feasibility = feasibility;
        this.synConfig = synConfig;

        this.gateTypes = new ArrayList<>();
        this.gateTypes.add(LogicType.EMPTY);
        this.gateTypes.addAll(gateLib.getRealizations().keySet());
        this.gateTypes.remove(LogicType.INPUT);
        this.gateTypes.remove(LogicType.OUTPUT_BUFFER);

        this.numGroupsInLib = gateLib.getGroups().size();

        /* init final fields */
        gateInputNames = new ArrayList<>();
        gateInputNames.add("x");
        gateInputNames.add("y");

        inputVars = new ArrayList<>();
        inputGates = new ArrayList<>();
        char varName = 'a';
        for (int i = 0; i < feasibility; i ++) {
            inputVars.add(""+varName);
            inputGates.add(new Gate(Character.toString(varName), LogicType.INPUT));
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

        List<TreeCircuit> rawCircuits = new ArrayList<>();
        resultCircuits = new HashMap<>();

        /* call recursive circuit build function */

        List<BuildWorker> buildWorkers = new ArrayList<>();

        /* initialize build workers with differing start gates */
        for (LogicType type : gateTypes) {

            if (type == LogicType.EMPTY)
                continue;

            TreeCircuit newCircuit = new TreeCircuit(type);
            buildWorkers.add(new BuildWorker(this, gateLib, newCircuit, Arrays.asList(LogicType.NOT, LogicType.NOR2, LogicType.EMPTY), 1, synConfig.getMaxDepth(), synConfig.getMaxWeight()));
        }

        int numThreads = synConfig.getSynLimitThreadsNum() == 0 ? availableProcessors : Math.min(synConfig.getSynLimitThreadsNum(), availableProcessors);
        logger.info("possible parallelism during building: " + buildWorkers.size());
        logger.info("number of threads resulting from hardware and synthesis settings: " + numThreads);

        ExecutorService executor = Executors.newFixedThreadPool(numThreads);
        List<Future<List<TreeCircuit>>> buildResults = Collections.emptyList();

        try {
            buildResults = executor.invokeAll(buildWorkers);
        } catch (InterruptedException e) {
            logger.error(e.getMessage());
        }

        for (Future<List<TreeCircuit>> resultList : buildResults) {

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
        for (TreeCircuit circuit : rawCircuits){
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
                //logger.error(e.getMessage());
            }
        }

        executor.shutdownNow();
    }

    boolean circuitAllowedByPreFilter(TreeCircuit circuit) {

        // check if circuit is implementable with library
        if (!isCoveredByLibrary(circuit.serializePreOrder().stream().map(n -> n.type).collect(Collectors.toList())))
            return false;

        // check is feasibility is met
        if (circuit.getNumOpenInputs(true) < feasibility)
            return false;

        // check circuit weight. if it exceeds maximum --> continue
        if (circuit.getWeight() > synConfig.getMaxWeight())
            return false;

        return true;
    }

    private void buildCircuitGraph(TreeCircuit candidate, String inputMapping) {

        Circuit circuit = new Circuit();
        int gateCounter = 0;

        List<TreeNode> gates = candidate.serializePostOrder();
        Map<TreeNode, Gate> transformation = new HashMap<>();
        int inputCount = 0;

        for (TreeNode gate : gates) {

            if (gate.type == LogicType.EMPTY)
                continue;

            gateCounter++;
            Gate newGate = new Gate(gate.type.name() + "_" + gateCounter, gate.type);
            circuit.addVertex(newGate);
            transformation.put(gate, newGate);

            if (gate.child0 != null && gate.child0.type != LogicType.EMPTY) {
                circuit.addEdge(transformation.get(gate.child0), newGate, new Wire(ExpressionParser.parse("x").variables().first()));
            }

            if (gate.child1 != null && gate.child1.type != LogicType.EMPTY) {
                circuit.addEdge(transformation.get(gate.child1), newGate, new Wire(ExpressionParser.parse("y").variables().first()));
            }

            if (gate.child0 == null || gate.child0.type == LogicType.EMPTY) {
                Gate inputGate = inputGates.get(Integer.parseInt(String.valueOf(inputMapping.charAt(inputCount))));

                if (!circuit.containsVertex(inputGate))
                    circuit.addVertex(inputGate);

                if (circuit.containsEdge(inputGate, transformation.get(gate)))
                    return;

                circuit.addEdge(inputGate, newGate, new Wire(ExpressionParser.parse("x").variables().first()));

                inputCount ++;
            }

            if (gate.type != LogicType.NOT && (gate.child1 == null || gate.child1.type == LogicType.EMPTY)) {
                Gate inputGate = inputGates.get(Integer.parseInt(String.valueOf(inputMapping.charAt(inputCount))));

                if (!circuit.containsVertex(inputGate))
                    circuit.addVertex(inputGate);

                if (circuit.containsEdge(inputGate, transformation.get(gate)))
                    return;

                circuit.addEdge(inputGate, newGate, new Wire(ExpressionParser.parse("y").variables().first()));

                inputCount ++;
            }
        }

        // if or gate (which serves as buffer) is not the last gate, add buffer
        if (candidate.getRootNode().type != LogicType.OUTPUT_OR2) {
            gateCounter ++;
            Gate outputBuffer = new Gate(LogicType.OUTPUT_BUFFER.name() + "_" + gateCounter, LogicType.OUTPUT_BUFFER);
            circuit.addVertex(outputBuffer);
            circuit.addEdge(transformation.get(candidate.getRootNode()), outputBuffer, new Wire(outputBuffer.getExpression().variables().first()));
        }

        // remove redundant gates
        if (!circuit.removeRedundantGates())
            return;

        int weight = circuit.getWeight();

        // if guided, compare to current min size
        if (isGuided) {
            if (weight > currentMinSize + synConfig.getWeightRelaxation())
                return;
        }

        Formula expression = circuit.getExpression();
        TruthTable truthTable = new TruthTable(expression);

        /* filter circuits with reduced support after redundancy removal */
        if (!isGuided) {
            if (expression.variables().size() < feasibility) {
                return;
            }
        }

        if (isGuided) {
            if (!truthTable.equals(targetTruthTable)) {
                logger.error("Synthesis error: " + truthTable);
            }
        }

        // filter out structurally equivalent circuits
        if (resultCircuits.containsKey(truthTable)) {
            for (Circuit existingCircuit : resultCircuits.get(truthTable)) {
                if (circuit.isEquivalent(existingCircuit))
                    return;
            }
        }

        if (isGuided) {
            if (weight < currentMinSize) {
                currentMinSize = weight;

                // update result circuits
                if (resultCircuits.containsKey(truthTable)) {
                    resultCircuits.get(truthTable).removeIf(c -> c.getWeight() > (currentMinSize + synConfig.getWeightRelaxation()));
                }
            }
        }

        resultCircuits.putIfAbsent(truthTable, new HashSet<>());
        resultCircuits.get(truthTable).add(circuit);
    }

    /* helper functions for handling of primitive circuits */

    boolean isCoveredByLibrary(List<LogicType> gates) {

        // check if enough gate groups are available
        int numGates = 0;

        for (LogicType element : gates) {
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

            for (LogicType gateType : gates) {
                if (gateType == type)
                    occurrences++;
            }

            if (occurrences > gateLib.getNumAvailableGates(type))
                return false;
        }

        return true;
    }

    boolean fulfillsMaxWeight(List<LogicType> row) {

        int weight = 0;

        for (LogicType type : row) {
            weight += type.getWeight();
        }

        return weight <= synConfig.getMaxWeight();
    }

    boolean isNotEmpty(List<LogicType> row) {

        for (LogicType type : row) {
            if (type != LogicType.EMPTY)
                return true;
        }

        return false;
    }

    String evaluateTreeCircuit(TreeCircuit circuit, String inputMapping) {

        List<TreeNode> postOrder = circuit.serializePostOrder();

        Map<String, TreeNode> inputVariableMapping = new HashMap<>();

        int i = 0;

        for (TreeNode node : postOrder) {

            if (node.type == LogicType.EMPTY)
                continue;

            if (node.type == LogicType.NOT) {
                if (node.child0 == null || node.child0.type == LogicType.EMPTY) {
                    node.expression = "~" + intermediateVariables.get(i++);
                } else {
                    node.expression = "~(" + node.child0.expression + ")";
                }
            } else {

                String left;
                String right;

                if (node.child0 == null || node.child0.type == LogicType.EMPTY) {
                    left = intermediateVariables.get(i++);
                    inputVariableMapping.put(left, node);
                } else {
                    left = node.child0.expression;
                }

                if (node.child1 == null || node.child1.type == LogicType.EMPTY) {
                    right = intermediateVariables.get(i++);
                    inputVariableMapping.put(right, node);
                } else {
                    right = node.child1.expression;
                }

                node.expression = (node.type == LogicType.NOR2 ? "~(" : "(") + left + "|" + right + ")";
            }
        }

        String expression = circuit.getRootNode().expression;

        Map<TreeNode, Set<String>> primaryInputMapping = new HashMap<>();

        for (int j = 0; j < i; j++) {

            String varToReplace = intermediateVariables.get(j);

            String substitution = inputVars.get(Character.getNumericValue(inputMapping.charAt(j)));

            TreeNode gate = inputVariableMapping.get(varToReplace);
            if (gate != null) {
                primaryInputMapping.putIfAbsent(gate, new HashSet<>());
                if (!primaryInputMapping.get(gate).add(substitution))
                    return null;
            }

            expression = expression.replaceAll(varToReplace, substitution);
        }

        return expression;
    }

    /* structural equivalence check */

    Map<List<Integer>, Integer> equivalenceClasses;
    int equivalenceClassCount;

    private List<TreeCircuit> filterRedundantCircuits(List<TreeCircuit> inputCircuits) {

        equivalenceClasses = new HashMap<>();
        equivalenceClassCount = 0;

        Map<Integer, TreeCircuit> results = new HashMap<>();

        for (TreeCircuit circuit : inputCircuits) {
            int equiv = getEquivalenceClass(circuit);
            results.putIfAbsent(equiv, circuit);
        }

        return new ArrayList<>(results.values());
    }

    private int getEquivalenceClass(TreeCircuit circuit) {

        List<TreeNode> postOrder = circuit.serializePostOrder();

        for (TreeNode node : postOrder) {

            List<Integer> list = new LinkedList<>();

            if (node.child0 == null && node.child1 == null)
                list.add(0);
            else if (node.child1 == null)
                list.add(node.child0.equivalenceClass);
            else {
                list.add(node.child0.equivalenceClass);
                list.add(node.child1.equivalenceClass);
            }

            list.sort(Integer::compareTo);

            list.add(0, node.type.ordinal());

            if (equivalenceClasses.containsKey(list))
                node.equivalenceClass = equivalenceClasses.get(list);
            else {
                equivalenceClasses.put(list, equivalenceClassCount);
                node.equivalenceClass = equivalenceClassCount;
                equivalenceClassCount++;
            }
        }

        return circuit.getRootNode().equivalenceClass;
    }

    /* getter */

    public HashMap<TruthTable, Set<Circuit>> getResultCircuits() {
        return resultCircuits;
    }
}
