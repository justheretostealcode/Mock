package de.tu_darmstadt.rs.synbio.mapping.util;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ObjectNode;
import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.common.circuit.*;
import de.tu_darmstadt.rs.synbio.synthesis.util.ExpressionParser;
import de.tu_darmstadt.rs.synbio.common.TruthTable;
import org.jgrapht.Graphs;
import org.jgrapht.traverse.TopologicalOrderIterator;
import org.logicng.datastructures.Assignment;
import org.logicng.formulas.FormulaFactory;
import org.logicng.formulas.Variable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class BranchAndBoundUtil {
    private static final Logger logger = LoggerFactory.getLogger(BranchAndBoundUtil.class);


    public static final String DEFAULT_INPUT_SPECIFICATION_ENTRY_FORMAT_STRING = "\"%s\" : {\"0\" : %s, \"1\" : %s}, ";
    public static final Map<String, Map<Boolean, Double>> CELLO_INPUT_SPECIFICATION = Map.ofEntries(
            new AbstractMap.SimpleEntry<>("a", Map.of(false, 0.0034, true, 2.8)),
            new AbstractMap.SimpleEntry<>("b", Map.of(false, 0.0013, true, 4.4)),
            new AbstractMap.SimpleEntry<>("c", Map.of(false, 0.0082, true, 2.5)),
            new AbstractMap.SimpleEntry<>("d", Map.of(false, 0.025, true, 0.31))
    );


    public static Circuit getCircuit(String path) throws IOException {
        File file = new File(path);
        final ObjectNode node;
        ObjectMapper mapper = new ObjectMapper();
        CircuitDeserializer deserializer = new CircuitDeserializer(Circuit.class);

        node = mapper.readValue(file, ObjectNode.class);

        return deserializer.deserializeString(node.get("graph").toString());
    }

    /**
     * This method creates the subproblems belonging to the provided circuit.   <br>
     * Each subproblem thereby contains one gate less than the previous subproblem.
     *
     * @param structure The circuit to from which the subproblems shall be created
     * @return A ordered list of subproblems (the first index is the original structure and the last is a circuit containing only one logic gate)
     */
    public static List<Circuit> getSubProblems(Circuit structure) {

        // Establish the reverse topological order of the logic gates
        ArrayList<Gate> gatesInReversedOrder = new ArrayList<>();
        Iterator<Gate> iterator = new TopologicalOrderIterator<>(structure);
        while (iterator.hasNext()) {
            Gate g = iterator.next();

            if (g.getType() == Gate.Type.LOGIC)
                gatesInReversedOrder.add(g);
        }

        Gate[] gateOrder = new Gate[gatesInReversedOrder.size()];
        gatesInReversedOrder.toArray(gateOrder);

        List<Gate> originalInputBuffers = structure.getInputBuffers();
        List<Map<Gate, List<Gate>>> substitutionsList;
        Map<Gate, String> substitutionTruthTables;
        // Create sub problems
        ArrayList<Circuit> subProblems = new ArrayList<>(gateOrder.length - 1);
        Circuit subProblem = new Circuit("Subproblem");
        Graphs.addGraph(subProblem, structure);
        FormulaFactory factory = new FormulaFactory();
        String structureIdentifier = structure.getIdentifier();
        String subProblemIdentifier = structureIdentifier + "_subproblem_";

        subProblems.add(subProblem.copy(subProblemIdentifier + 0));

        Map<Gate, Gate> insertedInputBuffers = new HashMap<>();
        for (int iX = 0; iX < gateOrder.length - 1; iX++) {
            Gate g = gateOrder[iX];

            // Remove Vertex
            Set<Wire> wires = subProblem.outgoingEdgesOf(g);

            // For each gate, which is removed, we need to add a new input to the circuit and connect this input to the following gates
            // Next to this, we mark, that the notes truthtable values result from a gate which was previous in the circuit.
            String inputIdentifier = "OUT_" + g.getIdentifier();    // Named "OUT", since it represents the output of the referenced gate
            Gate newInputBuffer = new InputGate(ExpressionParser.parse(inputIdentifier), inputIdentifier);
            subProblem.addVertex(newInputBuffer);
            insertedInputBuffers.put(newInputBuffer, g);
            for (Wire w : wires) {
                Gate target = subProblem.getEdgeTarget(w);
                subProblem.addEdge(newInputBuffer, target, new Wire(factory.variable(w.getVariable().name())));
            }

            subProblem.removeVertex(g);
            // Remove all Gates which do not contribute to the result
            cleanCircuitFromNonContributingGates(subProblem);
            subProblem.removeRedundantGates();

            // Get the Whitelist for the resulting structure
            String whiteList = determineWhitelist(structure, subProblem, insertedInputBuffers);
            subProblem.setWhitelist(whiteList);

            // Obtain substitutions list for the subproblem.
            substitutionsList = determineSubstitutionList(subProblem, originalInputBuffers);
            subProblem.setSubstitutionsList(substitutionsList);

            // Substitution Truthtables are not used anymore
            substitutionTruthTables = determineSubstitutionTruthTables(subProblem, substitutionsList, originalInputBuffers);
            subProblem.setSubstitutionTruthTables(substitutionTruthTables);

            subProblems.add(subProblem.copy(subProblemIdentifier + (iX + 1)));
        }

        return subProblems;
    }

    /**
     * Within this method, the whitelist of assignments is created. <br>
     * This is relevant, in order to only consider boolean input assignments, which are also present in the original circuit. <br>
     * A sub problem consisting of four input buffers has 16 possible input assignments, however, not all can be reached given the original structure with only 3 input buffers and thus 8 resulting logic input assignments.
     *
     * @param structure            The structure from which the sub problem is created
     * @param subProblem           The sub problem to determine the whitelist for
     * @param insertedInputBuffers A map of input buffers added to the circuit to the original logic gates replaced by them.
     * @return A whitelist string indicating which of the boolean input assignments of the sub problem are achieved by the original structure
     */
    private static String determineWhitelist(Circuit structure, Circuit subProblem, Map<Gate, Gate> insertedInputBuffers) {

        List<Assignment> inputAssignmentsCircuit = new ArrayList<>(TruthTable.getAllAssignments(structure.getExpression().variables()));

        // Generate the truth tables for every input
        HashMap<Gate, TruthTable> truthTables = new HashMap<>();
        TruthTable truthTable;

        for (Gate inputBuffer : subProblem.getInputBuffers()) {
            truthTable = new TruthTable(structure.getExpression(insertedInputBuffers.getOrDefault(inputBuffer, inputBuffer)), inputAssignmentsCircuit);
            truthTables.put(inputBuffer, truthTable);
        }

        // Obtain the available assignments from the truth tables
        List<Variable> inputVariablesSubProblem = new ArrayList<>(subProblem.getExpression().variables());
        int lastUsedBit = inputAssignmentsCircuit.size();
        Collections.reverse(inputVariablesSubProblem);  // Reversion is necessary, since the class TruthTable again reverts the order in the output String

        String[] prefix = {", ~", ", "};    // The two prefixes. The first indicates negative and the second one positive literals
        HashSet<String> availableCombinations = new HashSet<>(/*lastUsedBit*/);     // We use a set, since an assignment should only be included once.
        List<Gate> inputBuffers = subProblem.getInputBuffers();
        // Iterate over each entry of the truth table and then combine the values of the different truth tables into an assignment
        for (int iX = 0; iX < lastUsedBit; iX++) {
            StringBuilder[] stringBuilders = new StringBuilder[2];  // A StringBuilder for the positive and one for the negative Literals
            stringBuilders[0] = new StringBuilder();
            stringBuilders[1] = new StringBuilder();
            StringBuilder stringBuilder;

            // Iterate over the sub problem's input variables to combine the truth tables
            for (Variable v : inputVariablesSubProblem) {
                Gate vGate = inputBuffers.stream().filter(gate -> gate.getExpression().variables().first() == v).findFirst().get();
                TruthTable vT = truthTables.get(vGate);
                int vVal = vT.getTruthTable();
                vVal = (vVal & (1 << iX)) >> iX;

                stringBuilder = stringBuilders[vVal];
                stringBuilder.append(prefix[vVal]);
                stringBuilder.append(v.name());
            }

            // Create a string which would be equal to the .toString() reference of the assignment with the same values.
            String positiveLiterals = (stringBuilders[1].length() > 2) ? stringBuilders[1].substring(2) : "";
            String negativeLiterals = (stringBuilders[0].length() > 2) ? stringBuilders[0].substring(2) : "";

            String assig = "Assignment{pos=[" + positiveLiterals + "], neg=[" + negativeLiterals + "]}";

            // Add the assignment to the set of available assignments
            availableCombinations.add(assig);
        }

        // Iterate over all assignments for the sub problem in order to decide, whether the assignment is applicable or not
        StringBuilder whiteListBuilder = new StringBuilder();
        List<Assignment> inputAssignmentsSubProblem = new ArrayList<>(TruthTable.getAllAssignments(subProblem.getExpression().variables()));
        for (Assignment assignment : inputAssignmentsSubProblem) {

            if (availableCombinations.contains(assignment.toString()))
                whiteListBuilder.append("1");
            else
                whiteListBuilder.append("0");
        }

        whiteListBuilder.reverse(); // Reverse to be consistent with the reversed truth table reported by TruthTable class
        return whiteListBuilder.toString();
    }

    /**
     * This method removes gates from the structure, which do not take place in the evaluation of the circuit.
     *
     * @param structure The structure to clean.
     */
    private static void cleanCircuitFromNonContributingGates(Circuit structure) {
        Gate output = structure.getOutputBuffer();
        ArrayList<Gate> visitedNodes = new ArrayList<>(structure.vertexSet().size());
        traverseCircuit(structure, output, visitedNodes);

        structure.removeAllVertices(structure.vertexSet().stream().filter(gate -> !visitedNodes.contains(gate)).collect(Collectors.toList()));
    }

    /**
     * This method traverses a structure recursively from the given currentGate to the inputs, while visitedNodes contains all nodes visited.
     *
     * @param structure    The circuit to traverse
     * @param currentGate  The starting point of the traversation
     * @param visitedNodes The list of nodes visited
     */
    private static void traverseCircuit(Circuit structure, Gate currentGate, List<Gate> visitedNodes) {

        visitedNodes.add(currentGate);

        Gate source;

        for (Wire w : structure.incomingEdgesOf(currentGate)) {
            source = structure.getEdgeSource(w);
            traverseCircuit(structure, source, visitedNodes);
        }
    }

    /**
     * Determines the NOR gates and the assignments at which the values need to be replaced. <br>
     * Since the resulting string can be interpreted as a list of sets, the result can be used for every sub problem.
     *
     * @param structure The original circuit to assing and not any subproblem
     * @return A string which can be directly added
     */
    private static List<Map<Gate, List<Gate>>> determineSubstitutionList(Circuit structure, List<Gate> originalInputBuffers) {
        List<Assignment> inputAssignmentsCircuit = new ArrayList<>(TruthTable.getAllAssignments(structure.getExpression().variables()));

        List<LogicGate> logicGates = structure.getLogicGates();

        // Determine the ancestors of the relevant logic gates
        HashMap<LogicGate, Set<Gate>> ancestorsMap = new HashMap<>();
        for (LogicGate g : logicGates) {
            if (g.getLogicType() == LogicType.NOR2 || g.getLogicType() == LogicType.OR2) {
                Set<Gate> ancestorsSet = structure.edgesOf(g).stream().filter(edge -> structure.getEdgeTarget(edge) == g).map(structure::getEdgeSource).collect(Collectors.toSet());
                ancestorsMap.put(g, ancestorsSet);
            }
        }

        // Identify positions where a substitution is necessary
        List<Map<Gate, List<Gate>>> substitutionList = new ArrayList<>();
        Map<Gate, List<Gate>> pairsToSubstitute;
        for (Assignment assignment : inputAssignmentsCircuit) {
            pairsToSubstitute = new HashMap<>();
            for (LogicGate g : ancestorsMap.keySet()) {
                LogicType type = g.getLogicType();

                /*
                If the current NOR gates value is true, all inputs need to be zero -> correct bounds
                For the OR, it is only zero if both inputs are zero -> correct bounds
                This Filtering step simplifies the following determination of substitution, since the output of the gate is known and thus at least one input needs to be HIGH.
                 */
                if (type == LogicType.NOR2 && structure.getExpression(g).evaluate(assignment)
                        || type == LogicType.OR2 && !structure.getExpression(g).evaluate(assignment))
                    continue;

                Set<Gate> ancestors = ancestorsMap.get(g);
                for (Gate anc : ancestors) {
                    // TODO Watch out, whether this fits to the data to simulate. (If input distributions are considered, this could become a problem)
                    if (originalInputBuffers.contains(anc))
                        continue;       // Since the values of the real input Buffers are fixed, it is not necessary to substitute them.

                    boolean bInput = structure.getExpression(anc).evaluate(assignment);
                    if (!bInput) {   // If Input is false, the provided bound is not adequate and thus needs to be substituted.

                        if (!pairsToSubstitute.containsKey(anc))    // If anc not already included add the list
                            pairsToSubstitute.put(anc, new ArrayList<>());

                        pairsToSubstitute.get(anc).add(g);
                    }
                }
            }
            substitutionList.add(pairsToSubstitute);
        }
        return substitutionList;
    }

    /**
     * Converts the substitutionList in an equivalent representation containing only identifiers.
     *
     * @param substitutionList
     * @return
     */
    public static List<Map<String, List<String>>> substitutionListToString(List<Map<Gate, List<Gate>>> substitutionList) {
        return substitutionList.stream().map(map -> {
            Map<String, List<String>> newMap = new HashMap<>();
            map.forEach((key, value) -> newMap.put(key.getIdentifier(), value.stream().map(Gate::getIdentifier).collect(Collectors.toList())));
            return newMap;
        }).collect(Collectors.toList());
    }

    /**
     * Returns the truth tables of the gates for which substitution is to perform. <br>
     * This feature is currently not used by the simulator.
     *
     * @param structure            The sub problem for which the substitution truthtables shall be determined
     * @param substitutionList     The list of substitutions for each assignment
     * @param originalInputBuffers The input buffers of the original structure
     * @return Truth tables for the gates, whose output values get substituted.
     */
    private static Map<Gate, String> determineSubstitutionTruthTables(Circuit structure, List<Map<Gate, List<Gate>>> substitutionList, List<Gate> originalInputBuffers) {
        List<Assignment> assignments = TruthTable.getAllAssignments(structure.getExpression().variables());
        Set<Gate> gates = substitutionList.stream().map(Map::keySet).reduce(new HashSet<Gate>(), (subtotal, element) -> {
            subtotal.addAll(element);
            return subtotal;
        });
        Map<Gate, String> truthTables = new HashMap<>();
        gates
                .forEach(g -> {
                    String truthTable = new TruthTable(structure.getExpression(g), assignments).toString();
                    // Input Buffers added for the sub problem can not be bounded and shall be ignored.
                    if (structure.getInputBuffers().contains(g) && !originalInputBuffers.contains(g))
                        truthTable = "1".repeat(truthTable.length());
                    truthTables.put(g, truthTable);
                });

        return truthTables;
    }

    /**
     * Converts the substitutionTruthTables in an equivalent representation containing only identifiers instead of gates.
     *
     * @param substitutionTruthTables
     * @return
     */
    public static Map<String, String> substitutionTruthTablesToString(Map<Gate, String> substitutionTruthTables) {
        Map<String, String> truthTables = new HashMap<>();
        substitutionTruthTables.forEach((key, value) -> truthTables.put(key.getIdentifier(), value));
        return truthTables;
    }


    /**
     * Creates a custom input specification as required by the simulator.
     *
     * @param missingInputBufferIDs The elements to add beside the standard inputs ("A", "B", "C", "D")
     * @param minVal                The value used for "0"
     * @param maxVal                The value used for "1"
     * @return A string representing the input specification.
     */
    public static String createCustomInputSpecification(Set<String> missingInputBufferIDs, double minVal, double maxVal) {
        StringBuilder cis = new StringBuilder();
        cis.append("{");

        CELLO_INPUT_SPECIFICATION.entrySet().stream().filter(stringMapEntry -> !missingInputBufferIDs.contains(stringMapEntry.getKey())).forEach(entry -> {
            Map<Boolean, Double> value = entry.getValue();
            cis.append(String.format(DEFAULT_INPUT_SPECIFICATION_ENTRY_FORMAT_STRING, entry.getKey(), value.get(Boolean.FALSE), value.get(Boolean.TRUE)));
        });

        missingInputBufferIDs.forEach(entry -> {
            cis.append(String.format(DEFAULT_INPUT_SPECIFICATION_ENTRY_FORMAT_STRING, entry, minVal, maxVal));
        });

        cis.append("}");
        return cis.toString();
    }

    //TODO: implement creation of particle input specification
    /*public static String createCustomInputSpecificationParticles(Set<String> missingInputBufferIDs, List<Double> minVals, List<Double> maxVals) {
        StringBuilder cis = new StringBuilder();
        cis.append("{");

        CELLO_INPUT_SPECIFICATION.entrySet().stream().filter(stringMapEntry -> !missingInputBufferIDs.contains(stringMapEntry.getKey())).forEach(entry -> {
            Map<Boolean, Double> value = entry.getValue();
            cis.append(String.format(DEFAULT_INPUT_SPECIFICATION_ENTRY_FORMAT_STRING, entry.getKey(), value.get(Boolean.FALSE), value.get(Boolean.TRUE)));
        });

        missingInputBufferIDs.stream().forEach(entry -> {
            cis.append(String.format(DEFAULT_INPUT_SPECIFICATION_ENTRY_FORMAT_STRING, entry, minVal, maxVal));
        });

        cis.append("}");
        return cis.toString();
    }*/
}
