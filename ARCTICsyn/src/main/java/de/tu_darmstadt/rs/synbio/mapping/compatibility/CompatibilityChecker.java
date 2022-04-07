package de.tu_darmstadt.rs.synbio.mapping.compatibility;

import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.circuit.Gate;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.common.library.GateRealization;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.logicng.datastructures.Tristate;
import org.logicng.formulas.Formula;
import org.logicng.formulas.FormulaFactory;
import org.logicng.formulas.Variable;
import org.logicng.io.parsers.PropositionalParser;
import org.logicng.solvers.MiniSat;
import org.logicng.solvers.SATSolver;
import org.logicng.solvers.SolverState;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.stream.Collectors;

public class CompatibilityChecker {

    private static final Logger logger = LoggerFactory.getLogger(CompatibilityChecker.class);

    /* library information */

    final GateLibrary library;
    final Map<LogicType, List<GateRealization>> devices;

    final Map<String, String> deviceToTf;
    final Map<String, Set<String>> tfToDevice;

    /* circuit information */

    final Set<GateTriple> triples;
    final Set<Gate> gates;

    /* compatibility matrix */

    final CompatibilityMatrix<String> matrix;

    /* SAT utilities */

    final FormulaFactory factory;
    final PropositionalParser parser;
    final Formula constantTrue;
    final SATSolver miniSat;

    /* SAT */

    final SolverState state;

    public CompatibilityChecker(GateLibrary library, Circuit structure) {

        this.library = library;

        /* extract devices and groups from library */

        this.devices = library.getRealizations();

        this.deviceToTf = new HashMap<>();
        this.tfToDevice = new HashMap<>();

        for (List<GateRealization> devicesList : devices.values()) {
            for (GateRealization device : devicesList) {

                deviceToTf.putIfAbsent(device.getIdentifier(), device.getGroup());

                tfToDevice.putIfAbsent(device.getGroup(), new HashSet<>());
                tfToDevice.get(device.getGroup()).add(device.getIdentifier());
            }
        }

        /* compute triples of gates from the circuit topology */

        this.triples = extractTriples(structure);
        this.gates = structure.vertexSet().stream()
                .filter(g -> g.getLogicType() != LogicType.OUTPUT_BUFFER && g.getLogicType() != LogicType.OUTPUT_OR2)
                .collect(Collectors.toSet());

        /* initialize formula factory and MiniSAT */

        this.factory = new FormulaFactory();
        this.parser =  new PropositionalParser(factory);
        this.constantTrue = factory.constant(true);
        this.miniSat = MiniSat.miniSat(factory);

        // dummy matrix
        this.matrix = new CompatibilityMatrix<>();

        Random rand = new Random(123456789);
        double threshold = 0.25;

        for (String source : deviceToTf.keySet()) {
            for (String destination : deviceToTf.keySet()) {

                double randomNumber = rand.nextDouble();
                matrix.addEntry(source, destination, null, randomNumber > threshold);

                for (String secondSource : deviceToTf.keySet()) {

                    if (source.equals(destination) || secondSource.equals(destination) || source.equals(secondSource)) {
                        matrix.addEntry(source, destination, secondSource, false);
                    } else {
                        randomNumber = rand.nextDouble();
                        matrix.addEntry(source, destination, secondSource, randomNumber > threshold);
                    }
                }
            }
        }

        /* build SAT */

        List<Formula> clauses = new ArrayList<>();

        //long start = System.currentTimeMillis();

        /* Clause 1: Every gate has to be assigned one device */

        for (Gate gate : gates) {

            StringBuilder builder = new StringBuilder();

            for (String deviceId : deviceToTf.keySet()) {

                builder.append(gate + "_" + deviceId + "& ~(");

                for (String device2Id : deviceToTf.keySet()) {
                    if (device2Id.equals(deviceId))
                        continue;

                    builder.append(gate + "_" + device2Id + "|");
                }

                builder.deleteCharAt(builder.length() - 1);
                builder.append(")|");
            }

            builder.deleteCharAt(builder.length() - 1);

            try {
                Formula formula = parser.parse(builder.toString());
                clauses.add(formula);
            } catch (Exception e) {
                logger.error(e.getMessage());
            }
        }

        /* Clause 2: Every device and its group members must be assigned maximally once */

        for (String deviceId : deviceToTf.keySet()) {

            StringBuilder builder = new StringBuilder();

            for (Gate gate : gates) {

                builder.append("(" + gate + "_" + deviceId + "& ~(");

                for (Gate gate2 : gates) {

                    if (gate2.equals(gate))
                        continue;

                    for (String groupMember : tfToDevice.get(deviceToTf.get(deviceId))) {
                        builder.append(gate2 + "_" + groupMember + "|");
                    }
                }

                builder.deleteCharAt(builder.length() - 1);
                builder.append("))|");
            }

            builder.append("(");
            for (Gate gate : gates) {

                builder.append("~" + gate + "_" + deviceId + "&");
            }

            builder.deleteCharAt(builder.length() - 1);
            builder.append(")");

            try {
                Formula formula = parser.parse(builder.toString());
                clauses.add(formula);
            } catch (Exception e) {
                logger.error(e.getMessage());
            }

        }

        /* Clause 3: Every gate has to be assigned a device with the right type */

        for (Gate gate : gates) {

            StringBuilder builder = new StringBuilder();

            for (GateRealization device : devices.get(gate.getLogicType())) {
                builder.append(gate + "_" + device.getIdentifier() + "|");
            }

            builder.deleteCharAt(builder.length() - 1);

            try {
                Formula formula = parser.parse(builder.toString());
                clauses.add(formula);
            } catch (Exception e) {
                logger.error(e.getMessage());
            }
        }

        /* Clause 4: All triples of devices have to be compatible */

        for (GateTriple triple : triples) {

            StringBuilder builder = new StringBuilder();

            for (String destinationDevice : deviceToTf.keySet()) {

                builder.append("(" + triple.destination + "_" + destinationDevice + "&(");

                for (String sourceDevice : deviceToTf.keySet()) {

                    if (sourceDevice.equals(destinationDevice))
                        continue;

                    builder.append(triple.source + "_" + sourceDevice + "&(");

                    if (triple.secondSource == null) {
                        builder.append("(" + (matrix.isCompatible(sourceDevice, destinationDevice, null) ? "$true" : "$false") + ")|");
                    } else {

                        for (String secondSourceDevice : deviceToTf.keySet()) {

                            if (secondSourceDevice.equals(destinationDevice) || secondSourceDevice.equals(sourceDevice))
                                continue;

                            builder.append("(" + triple.secondSource + "_" + secondSourceDevice + "&" +
                                    (matrix.isCompatible(sourceDevice, destinationDevice, secondSourceDevice) ? "$true" : "$false") + ")|");
                        }
                    }

                    builder.deleteCharAt(builder.length() - 1);
                    builder.append(")|");
                }

                builder.deleteCharAt(builder.length() - 1);
                builder.append("))|");
            }

            builder.deleteCharAt(builder.length() - 1);

            try {
                Formula formula = parser.parse(builder.toString());
                clauses.add(formula);
            } catch (Exception e) {
                logger.error(e.getMessage());
            }
        }

        /* add clauses */
        for (Formula clause : clauses) {
            miniSat.add(clause);
        }

        state = miniSat.saveState();
    }

    public boolean isCompatible(Assignment incompleteAssigment) {

        miniSat.loadState(state);
        /* determine vars to substitute with constant TRUE to account for incomplete assignment */

        for (Gate gate : incompleteAssigment.keySet()) {
            if (gate.getLogicType() != LogicType.OUTPUT_BUFFER && gate.getLogicType() != LogicType.OUTPUT_OR2)
                miniSat.add(factory.variable(gate.getIdentifier() + "_" + incompleteAssigment.get(gate).getIdentifier()));
                //constants.add(factory.variable(gate.getIdentifier() + "_" + incompleteAssigment.get(gate).getIdentifier()));
        }

        /* solve SAT */

        Tristate result = miniSat.sat();


        /*logger.info(incompleteAssigment.getIdentifierMap()+"");
        org.logicng.datastructures.Assignment ass = miniSat.model();
        if (ass != null) {
            logger.info("SAT model:");
            for (Variable pos : ass.positiveLiterals()) {
                if (!pos.name().startsWith("@"))
                    logger.info(pos.name());
            }
        }*/

        return result == Tristate.TRUE;
    }

    public boolean verify(Assignment assignment) {

        for (GateTriple triple : triples) {

            GateRealization destination = assignment.get(triple.destination);
            GateRealization source = assignment.get(triple.source);
            GateRealization secondSource = triple.secondSource != null ? assignment.get(triple.secondSource) : null;

            if (destination == null || source == null || (triple.secondSource != null && assignment.get(triple.secondSource) == null))
                continue;

            if (!matrix.isCompatible(source.getIdentifier(), destination.getIdentifier(), secondSource != null ? secondSource.getIdentifier() : null))
                return false;
        }

        return true;
    }

    /* private helper functions */

    private Formula substituteVars(Formula formula, List<Variable> vars) {

        for (Variable var : vars) {
            formula = formula.substitute(var, constantTrue);
        }

        return formula;
    }

    private Set<GateTriple> extractTriples(Circuit circuit) {

        Set<GateTriple> triples = new HashSet<>();

        for (Gate destination : circuit.vertexSet()) {

            if (!destination.isLogicGate())
                continue;

            List<Gate> sources = new ArrayList<>();
            circuit.incomingEdgesOf(destination).forEach(w -> sources.add(circuit.getEdgeSource(w)));

            if (sources.size() > 1)
                triples.add(new GateTriple(sources.get(0), destination, sources.get(1)));
            else
                triples.add(new GateTriple(sources.get(0), destination));

        }
        return triples;
    }

    private class GateTriple {

        public final Gate source;
        public final Gate destination;
        public final Gate secondSource;

        public GateTriple(Gate source, Gate destination, Gate secondSource) {
            this.source = source;
            this.destination = destination;
            this.secondSource = secondSource;
        }

        public GateTriple(Gate source, Gate destination) {
            this(source, destination, null);
        }

        @Override
        public boolean equals(Object o) {

            if (o == this) {
                return true;
            }

            if (!(o instanceof GateTriple)) {
                return false;
            }

            GateTriple t = (GateTriple) o;

            return new EqualsBuilder()
                    .append(this.source, t.source)
                    .append(this.destination, t.destination)
                    .append(this.secondSource, t.secondSource)
                    .isEquals();
        }

        @Override
        public int hashCode() {
           return new HashCodeBuilder()
                   .append(source)
                   .append(destination)
                   .append(secondSource)
                   .toHashCode();
        }
    }
}
