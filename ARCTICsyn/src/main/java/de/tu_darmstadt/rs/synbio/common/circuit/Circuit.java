package de.tu_darmstadt.rs.synbio.common.circuit;

import com.fasterxml.jackson.core.Version;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.module.SimpleModule;
import de.tu_darmstadt.rs.synbio.common.library.GateRealization;
import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.mapping.util.BranchAndBoundUtil;
import de.tu_darmstadt.rs.synbio.common.TruthTable;
import org.jgrapht.GraphPath;
import org.jgrapht.Graphs;
import org.jgrapht.alg.isomorphism.VF2GraphIsomorphismInspector;
import org.jgrapht.alg.shortestpath.AllDirectedPaths;
import org.jgrapht.graph.DirectedAcyclicGraph;
import org.jgrapht.io.*;
import org.jgrapht.traverse.TopologicalOrderIterator;
import org.logicng.formulas.Formula;
import org.logicng.formulas.Variable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.FileWriter;
import java.io.Writer;
import java.util.*;
import java.util.stream.Collectors;

public class Circuit extends DirectedAcyclicGraph<Gate, Wire> implements Comparable<Circuit> {

    private static final Logger logger = LoggerFactory.getLogger(Circuit.class);

    private String identifier;

    // Additional parameters only required by Branch and Bound
    private String whitelist;
    private List<Map<Gate, List<Gate>>> substitutionsList;
    private Map<Gate, String> substitutionTruthTables;

    public Circuit() {
        super(Wire.class);
    }

    public Circuit(String identifier) {
        super(Wire.class);
        this.identifier = identifier;
        this.whitelist = null;
        this.substitutionsList = null;
        this.substitutionTruthTables = null;
    }

    public Circuit(GateRealization realization, String identifier) {
        super(Wire.class);
        this.identifier = identifier;
        this.whitelist = null;
        this.substitutionsList = null;
        this.substitutionTruthTables = null;

        Gate gate = new Gate(realization.getIdentifier(), realization.getLogicType());

        Gate outputBuffer = new Gate("F", LogicType.OUTPUT);
        this.addVertex(gate);
        this.addVertex(outputBuffer);
        this.addEdge(gate, outputBuffer, new Wire(outputBuffer.getExpression().variables().first()));

        Gate inputBuffer;

        for (Variable inputVar : gate.getExpression().variables()) {
            inputBuffer = new Gate(inputVar.name(), LogicType.INPUT);
            this.addVertex(inputBuffer);
            this.addEdge(inputBuffer, gate, new Wire(inputVar));
        }
    }

    public Gate getOutputBuffer() {
        return vertexSet().stream().filter(gate -> gate.getLogicType() == LogicType.OUTPUT).findFirst().get();
    }

    public List<Gate> getInputBuffers() {
        return vertexSet().stream().filter(gate -> (gate.getLogicType() == LogicType.INPUT)).collect(Collectors.toList());
    }

    public Gate getInputBuffer(Variable variable) {
        return vertexSet().stream().filter(gate -> (gate.getLogicType() == LogicType.INPUT && gate.getExpression().containsVariable(variable))).findFirst().get();
    }

    public List<Gate> getLogicGates() {
        return vertexSet().stream().filter(Gate::isLogicGate).collect(Collectors.toList());
    }

    public Formula getExpression() {
        return getExpression(getOutputBuffer());
    }

    public Formula getExpression(Gate startNode) {

        if (startNode.getLogicType() == LogicType.INPUT)
            return startNode.getExpression();

        Formula expression = startNode.getExpression();

        for (Variable var : expression.variables()) {

            Gate connectedGate = getEdgeSource(edgesOf(startNode).stream().filter(w -> w.getVariable().equals(var)).findFirst().get());
            expression = expression.substitute(var, getExpression(connectedGate));
        }

        return expression;
    }

    public TruthTable getTruthTable() {
        return new TruthTable(getExpression());
    }

    public String getIdentifier() {
        return identifier;
    }

    public Integer getWeight() {
        return vertexSet().stream().mapToInt(Gate::getWeight).sum();
    }

    public int getNumberLogicGates() {
        return (int) vertexSet().stream().filter(Gate::isLogicGate).count();
    }

    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }

    public void setWhitelist(String whitelist) {
        this.whitelist = whitelist;
    }

    public void setSubstitutionsList(List<Map<Gate, List<Gate>>> substitutionsList) {
        this.substitutionsList = substitutionsList;
    }

    public void setSubstitutionTruthTables(Map<Gate, String> substitutionTruthTables) {
        this.substitutionTruthTables = substitutionTruthTables;
    }

    public void replaceGate(Gate node, Gate replacement) {
        this.addVertex(replacement);
        for (Wire wire : this.outgoingEdgesOf(node))
            this.addEdge(replacement, this.getEdgeTarget(wire), new Wire(wire.getVariable()));
        for (Wire wire : this.incomingEdgesOf(node))
            this.addEdge(this.getEdgeSource(wire), replacement, new Wire(wire.getVariable()));
        this.removeVertex(node);
    }

    private boolean isValid() {

        for (Gate logicGate : vertexSet().stream().filter(Gate::isLogicGate).collect(Collectors.toList())) {

            // test if number of input wires equals support size
            if (logicGate.getExpression().variables().size() > incomingEdgesOf(logicGate).size())
                return false;

        }

        for (Gate gate : vertexSet().stream().filter(gate -> gate.getLogicType() != LogicType.OUTPUT).collect(Collectors.toList())) {

            // test if gate output is connected
            if (outgoingEdgesOf(gate).isEmpty())
                return false;
        }

        return true;
    }

    public boolean isEquivalent(Circuit cmp) {

        if (vertexSet().size() != cmp.vertexSet().size())
            return false;

        if (edgeSet().size() != cmp.edgeSet().size())
            return false;

        if (!usesSameSetOfGates(cmp))
            return false;

        return new VF2GraphIsomorphismInspector<>(this, cmp, new GateComparator(), new WireComparator()).isomorphismExists();
    }

    public boolean removeRedundantGates() {

        TruthTable beforeTT = new TruthTable(getExpression());

        // build redundancy map
        HashMap<String, List<Gate>> redMap = new HashMap<>();

        Iterator<Gate> iterator = new TopologicalOrderIterator<>(this);

        while (iterator.hasNext()) {

            Gate gate = iterator.next();

            if (gate.isLogicGate()) {

                Formula expr = getExpression(gate);
                TruthTable gateTT = new TruthTable(getExpression(gate));
                String gateClassIdentifier = gateTT.toString() + expr.variables().toString();

                redMap.putIfAbsent(gateClassIdentifier, new ArrayList<>());
                redMap.get(gateClassIdentifier).add(gate);
            }
        }

        // remove redundancies
        for (List<Gate> gateClass : redMap.values()) {

            if (gateClass.size() > 1) {

                Gate keeper = gateClass.get(0);

                for (int i = 1; i < gateClass.size(); i++) {

                    // collect information on redundant gate
                    Gate redGate = gateClass.get(i);
                    Set<Wire> inWires = new HashSet<>(incomingEdgesOf(redGate));
                    Optional<Wire> outWire = outgoingEdgesOf(redGate).stream().findFirst();

                    if (outWire.isEmpty())
                        continue;

                    Gate targetGate = getEdgeTarget(outWire.get());
                    Variable outVar = outWire.get().getVariable();

                    // remove redundant gate and edges
                    removeEdge(outWire.get());
                    inWires.forEach(this::removeEdge);
                    removeVertex(redGate);

                    // add new edge
                    if (!containsEdge(keeper, targetGate))
                        addEdge(keeper, targetGate, new Wire(outVar));
                }
            }
        }

        // remove dangling gates

        List<Gate> danglingGates = new ArrayList<>();
        do {

            for (Gate danglingGate : danglingGates) {

                Set<Wire> inWires = new HashSet<>(incomingEdgesOf(danglingGate));
                inWires.forEach(this::removeEdge);
                removeVertex(danglingGate);
            }

            danglingGates.clear();

            for (Gate logicGate : vertexSet().stream().filter(gate -> gate.isLogicGate()).collect(Collectors.toList())) {

                Optional<Wire> outWire = outgoingEdgesOf(logicGate).stream().findFirst();

                if (outWire.isEmpty())
                    danglingGates.add(logicGate);
            }

        } while (!danglingGates.isEmpty());

        if (!isValid())
            return false;

        TruthTable afterTT = new TruthTable(getExpression());

        if (!beforeTT.equalsLogically(afterTT))
            logger.error("Circuit with redundancy removed does not equal input circuit!");

        return true;
    }

    public boolean usesSameSetOfGates(Circuit cmp) {

        HashMap<LogicType, Integer> logicTypeCount = new HashMap<>();
        HashMap<LogicType, Integer> cmpLogicTypeCount = new HashMap<>();

        vertexSet().stream().filter(Gate::isLogicGate).forEach(g -> {
            LogicType logicType = g.getLogicType();
            logicTypeCount.putIfAbsent(logicType, 0);
            logicTypeCount.put(logicType, logicTypeCount.get(logicType) + 1);
        });

        cmp.vertexSet().stream().filter(Gate::isLogicGate).forEach(g -> {
            LogicType logicType = g.getLogicType();
            cmpLogicTypeCount.putIfAbsent(logicType, 0);
            cmpLogicTypeCount.put(logicType, cmpLogicTypeCount.get(logicType) + 1);
        });

        return logicTypeCount.equals(cmpLogicTypeCount);
    }

    public HashMap<Gate, List<Variable>> getUnconnectedGateInputs() {

        HashMap<Gate, List<Variable>> gateInputs = new HashMap<>();

        for (Gate gate : vertexSet().stream().filter(Gate::isLogicGate).collect(Collectors.toList())) {

            ArrayList<Variable> inputs = new ArrayList<>(gate.getExpression().variables());
            inputs.removeAll(incomingEdgesOf(gate).stream().map(Wire::getVariable).collect(Collectors.toList()));

            if (inputs.size() > 0)
                gateInputs.put(gate, inputs);
        }

        return gateInputs;
    }

    public int getDepth() {

        AllDirectedPaths<Gate, Wire> allPaths = new AllDirectedPaths<>(this);

        Gate output = getOutputBuffer();
        List<Gate> inputs = getInputBuffers();

        int depth = 0;

        for (Gate input : inputs) {

            List<GraphPath<Gate, Wire>> paths = allPaths.getAllPaths(input, output, true, null);
            OptionalInt pathLength = paths.stream().mapToInt(GraphPath::getLength).max();

            if (pathLength.isPresent()) {
                int length = pathLength.getAsInt();

                if (length > depth)
                    depth = length;
            }
        }

        return depth - 1;
    }

    /**
     * Estimates the circuit complexity based on the size and the used gates.
     * Gates with n inputs are weighted with n.
     *
     * @return
     */
    public double estimateCircuitComplexity() {
        double circuitComplexity = 0;

        Iterator<Gate> iterator = new TopologicalOrderIterator<Gate, Wire>(this);

        double complexity;

        int iX = 1;
        while (iterator.hasNext()) {
            Gate g = iterator.next();
            if (g.isLogicGate()) {
                complexity = this.getAncestors(g).stream().filter(g2 -> this.getEdge(g2, g) != null).count();
                circuitComplexity += complexity * iX;

                iX++;
            }
        }
        return circuitComplexity;
    }

    /**
     * @return A copy of the circuit on which this method was invoked
     */
    public Circuit copy() {
        return copy(this.getIdentifier());
    }

    /**
     * @param identifier The identifier chosen for the copied circuit
     * @return A copy of the circuit on which this method was invoked
     */
    public Circuit copy(String identifier) {
        Circuit newCircuit = new Circuit(identifier);
        Graphs.addGraph(newCircuit, this);
        newCircuit.setWhitelist(this.whitelist);
        if (substitutionsList != null)
            newCircuit.setSubstitutionsList(substitutionsList.stream().map(map -> {
                Map<Gate, List<Gate>> newMap = new HashMap<Gate, List<Gate>>();
                map.forEach((key, value) -> newMap.put(key, new ArrayList<>(value)));
                return newMap;
            }).collect(Collectors.toList()));  // Copy list and maps. Gates are equal in the maps.

        if (substitutionTruthTables != null)
            newCircuit.setSubstitutionTruthTables(new HashMap<>(substitutionTruthTables));
        return newCircuit;
    }

    private Map<Gate, TruthTable> getGateTruthTables() {

        Map<Gate, TruthTable> truthTables = new HashMap<>();

        List<org.logicng.datastructures.Assignment> assignments = TruthTable.getAllAssignments(this.getExpression().variables());

        for (Gate gate : this.vertexSet()) {
            truthTables.put(gate, new TruthTable(this.getExpression(gate), assignments));
        }

        return truthTables;
    }


    /* comparators */

    private static class GateComparator implements Comparator<Gate> {
        @Override
        public int compare(Gate first, Gate second) {

            if (first.getLogicType() != second.getLogicType())
                return 1;

            switch (first.getLogicType()) {
                case INPUT:
                case OUTPUT:
                    return 0;
                default:
                    return first.getLogicType().compareTo(second.getLogicType());
            }
        }
    }

    private static class WireComparator implements Comparator<Wire> {
        @Override
        public int compare(Wire first, Wire second) {
            return 0;
        }
    }

    @Override
    public int compareTo(Circuit cmp) {
        return cmp.getWeight() - getWeight();
    }

    /* IO */

    public void saveGml(File outputFile, Assignment assignment) {

        ComponentNameProvider<Gate> vertexIdProvider = gate -> "\"" + gate.getIdentifier() + "\"";

        ComponentNameProvider<Gate> vertexLabelProvider = gate -> {
            if (gate.isLogicGate()) {
                String altIdentifier = assignment == null ? "UNASSIGNED" : assignment.get(gate).getIdentifier();

                if (!altIdentifier.equals(""))
                    return gate.getLogicType() + "," + altIdentifier;
                else
                    return gate.getLogicType().toString();
            }
            return gate.getIdentifier().equals("O") ? "X" : gate.getIdentifier().toUpperCase();
        };

        ComponentNameProvider<Wire> edgeLabelProvider = wire -> "";

        ComponentNameProvider<Wire> edgeIdProvider = new IntegerComponentNameProvider<>();//wire -> String.valueOf(wire.hashCode());

        GmlExporter<Gate, Wire> exporter = new GmlExporter<>(vertexIdProvider, vertexLabelProvider, edgeIdProvider, edgeLabelProvider);
        exporter.setParameter(GmlExporter.Parameter.EXPORT_VERTEX_LABELS, true);

        try {
            Writer writer = new FileWriter(outputFile);
            exporter.exportGraph(this, writer);
        } catch(Exception e) {
            logger.error(e.getMessage());
        }

    }

    public void print(File outputFile) {

        ComponentNameProvider<Gate> vertexIdProvider = Gate::getIdentifier;

        ComponentNameProvider<Gate> vertexLabelProvider = gate -> {

            if (gate.isLogicGate()) {
                return gate.getIdentifier();
            }
            return gate.getLogicType() + " " + gate.getIdentifier();
        };

        ComponentNameProvider<Wire> edgeLabelProvider = wire -> "";

        GraphExporter<Gate, Wire> exporter = new DOTExporter<>(vertexIdProvider, vertexLabelProvider, edgeLabelProvider);

        try {
            Writer writer = new FileWriter(outputFile);
            exporter.exportGraph(this, writer);
            writer.close();
        } catch(Exception e) {
            logger.error(e.getMessage());
        }

    }

    public void print(File outputFile, Assignment assignment) {

        ComponentNameProvider<Gate> vertexIdProvider = Gate::getIdentifier;

        ComponentNameProvider<Gate> vertexLabelProvider = gate -> gate.getIdentifier() + "\n" + assignment.get(gate).getIdentifier();

        ComponentNameProvider<Wire> edgeLabelProvider = wire -> "";

        ComponentAttributeProvider<Gate> vertexAttrProvider = gate -> {
            Map<String, Attribute> map = new LinkedHashMap<>();
            if (gate.isLogicGate()) {
                map.put("shape", DefaultAttribute.createAttribute("box"));
                map.put("fixedsize", DefaultAttribute.createAttribute(true));
                map.put("width", DefaultAttribute.createAttribute(1.0));
                map.put("height", DefaultAttribute.createAttribute(0.6));
            } else {
                map.put("shape", DefaultAttribute.createAttribute("circle"));
            }
            return map;
        };

        GraphExporter<Gate, Wire> exporter = new DOTExporter<>(vertexIdProvider, vertexLabelProvider, edgeLabelProvider, vertexAttrProvider, null);

        try {
            Writer writer = new FileWriter(outputFile);
            exporter.exportGraph(this, writer);
            writer.close();
        } catch(Exception e) {
            logger.error(e.getMessage());
        }
    }

    public void save(File file) {

        ObjectMapper mapper = new ObjectMapper();
        CircuitSerializer circuitSerializer = new CircuitSerializer(Circuit.class);
        SimpleModule module = new SimpleModule("CircuitSerializer", new Version(1, 0, 0, null, null, null));
        module.addSerializer(Circuit.class, circuitSerializer);
        mapper.registerModule(module);

        try {
            Map<String, Object> jsonMap = new HashMap<>();
            jsonMap.put("graph", this);
            jsonMap.put("truthtable", getTruthTable().toString());
            if (whitelist != null) {
                jsonMap.put("whitelist", whitelist);
            }
            if (substitutionsList != null) {
                jsonMap.put("substitutions_list", substitutionsList.stream().map(map -> {
                    Map<String, List<String>> newMap = new HashMap<>();
                    map.forEach((key, value) -> newMap.put(key.getIdentifier(), value.stream().map(Gate::getIdentifier).collect(Collectors.toList())));
                    return newMap;
                }).collect(Collectors.toList()));
            }
            if (substitutionTruthTables != null) {
                jsonMap.put("substitution_truthtables", substitutionTruthTables);
            }

            jsonMap.put("gate_truthtables", getGateTruthTables());

            mapper.writerWithDefaultPrettyPrinter().writeValue(file, jsonMap);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
