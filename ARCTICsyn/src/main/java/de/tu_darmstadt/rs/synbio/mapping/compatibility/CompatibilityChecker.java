package de.tu_darmstadt.rs.synbio.mapping.compatibility;

import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.circuit.Gate;
import de.tu_darmstadt.rs.synbio.common.circuit.LogicGate;
import de.tu_darmstadt.rs.synbio.common.circuit.Wire;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.common.library.GateRealization;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import org.logicng.datastructures.Tristate;
import org.logicng.formulas.Formula;
import org.logicng.formulas.FormulaFactory;
import org.logicng.formulas.Variable;
import org.logicng.io.parsers.PropositionalParser;
import org.logicng.solvers.MiniSat;
import org.logicng.solvers.SATSolver;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;

public class CompatibilityChecker {

    private static final Logger logger = LoggerFactory.getLogger(CompatibilityChecker.class);

    /* library information */

    final GateLibrary library;
    final List<String> repressors;
    final Map<String, Set<String>> groups;

    /* compatibility matrix */

    final CompatibilityMatrix<String> matrix;

    /* SAT utilities */

    final FormulaFactory factory;
    final PropositionalParser parser;
    final Formula constantTrue;
    final SATSolver miniSat;

    public CompatibilityChecker(GateLibrary library) {

        this.library = library;

        /* extract repressors and groups from library */

        this.repressors = extractRepressors(library);

        this.groups = new HashMap<>();
        for (String repressor : repressors) {
            groups.putIfAbsent(getGroup(repressor), new HashSet<>());
            groups.get(getGroup(repressor)).add(repressor);
        }

        /* initialize formula factory and MiniSAT */

        this.factory = new FormulaFactory();
        this.parser =  new PropositionalParser(factory);
        this.constantTrue = factory.constant(true);
        this.miniSat = MiniSat.miniSat(factory);

        // dummy matrix
        this.matrix = new CompatibilityMatrix<>();

        Random rand = new Random();

        for (String source : repressors) {
            for (String destination : repressors) {

                if (source.equals(destination)) {
                    matrix.addEntry(source, destination, false);
                    continue;
                }

                double randomNumber = rand.nextDouble();
                matrix.addEntry(source, destination, randomNumber > 0.5);
            }
        }
    }

    public boolean isCompatible(Circuit structure, Assignment incompleteAssigment) {

        factory.clear();
        miniSat.reset();

        /* extract gates from the circuit */

        List<String> gates = extractGates(structure);

        /* compute pairs of gates from the circuit topology */

        Map<String, Set<String>> pairs = extractPairs(structure);

        /* determine vars to substitute with constant TRUE to account for incomplete assignment */

        List<Variable> constants = new ArrayList<>();

        if (incompleteAssigment != null) {
            for (LogicGate gate : incompleteAssigment.keySet()) {
                constants.add(factory.variable(gate.getIdentifier() + "_" + incompleteAssigment.get(gate).getIdentifier()));
            }
        }

        /* Clause 1: Every gate has to be assigned one repressor */

        for (String gate : gates) {

            StringBuilder builder = new StringBuilder();

            for (String repressor : repressors) {

                builder.append(gate + "_" + repressor + "& ~(");

                for (String repressor2 : repressors) {
                    if (repressor2.equals(repressor))
                        continue;

                    builder.append(gate + "_" + repressor2 + "|");
                }

                builder.deleteCharAt(builder.length() - 1);
                builder.append(")|");
            }

            builder.deleteCharAt(builder.length() - 1);

            try {
                Formula formula = parser.parse(builder.toString());
                miniSat.add(substituteVars(formula, constants));
            } catch (Exception e) {
                logger.error(e.getMessage());
            }

        }

        /* Clause 2: Every repressor and its group members must be assigned maximally once */

        for (String repressor : repressors) {

            StringBuilder builder = new StringBuilder();

            for (String gate : gates) {

                builder.append(gate + "_" + repressor + "& ~(");

                for (String gate2 : gates) {
                    if (gate2.equals(gate))
                        continue;

                    for (String groupMember : groups.get(getGroup(repressor))) {
                        builder.append(gate2 + "_" + groupMember + "|");
                    }
                }

                builder.deleteCharAt(builder.length() - 1);
                builder.append(")|");
            }

            builder.append("(");
            for (String gate : gates) {
                for (String groupMember : groups.get(getGroup(repressor))) {
                    builder.append("~" + gate + "_" + groupMember + "&");
                }
            }

            builder.deleteCharAt(builder.length() - 1);
            builder.append(")");

            try {
                Formula formula = parser.parse(builder.toString());
                miniSat.add(substituteVars(formula, constants));
            } catch (Exception e) {
                logger.error(e.getMessage());
            }

        }

        /* Clause 3: All pairs of source and destination repressors have to be compatible */

        for (String source : pairs.keySet()) {

            StringBuilder builder = new StringBuilder();

            for (String destination : pairs.get(source)) {

                builder.append("(");

                for (String repressor1 : repressors) {

                    builder.append("(" + source + "_" + repressor1 + "&(");

                    for (String repressor2 : repressors) {
                            builder.append("(" + destination + "_" + repressor2 + "&" +
                                    (matrix.isCompatible(repressor1, repressor2) ? "$true" : "$false") + ")|");
                    }

                    builder.deleteCharAt(builder.length() - 1);
                    builder.append("))|");
                }

                builder.deleteCharAt(builder.length() - 1);
                builder.append(")&");
            }

            builder.deleteCharAt(builder.length() - 1);

            try {
                Formula formula = parser.parse(builder.toString());
                miniSat.add(substituteVars(formula, constants));
            } catch (Exception e) {
                logger.error(e.getMessage());
            }
        }

        Tristate result = miniSat.sat();
        org.logicng.datastructures.Assignment ass = miniSat.model();

        if (ass != null) {
            logger.info("SAT model:");
            for (Variable pos : ass.positiveLiterals()) {
                if (!pos.name().startsWith("@"))
                    logger.info(pos.name());
            }
        }

        return result == Tristate.TRUE;
    }

    public boolean isCompatible(Circuit structure) {
        return isCompatible(structure, null);
    }

    /* private helper functions */

    private List<String> extractRepressors(GateLibrary library) {

        List<String> extractedRepressors = new ArrayList<>();

        for (List<GateRealization> realizations : library.getRealizations().values()) {
            for(GateRealization realization : realizations) {
                if (realization.isCharacterized()) {

                    String repressor = realization.getIdentifier();

                    if (!extractedRepressors.contains(repressor))
                        extractedRepressors.add(repressor);
                }
            }
        }
        return extractedRepressors;
    }

    private String getGroup(String altIdentifier) {
        String[] repressorString = altIdentifier.split("_");
        return repressorString[repressorString.length - 1];
    }

    private Formula substituteVars(Formula formula, List<Variable> vars) {

        for (Variable var : vars) {
            formula = formula.substitute(var, constantTrue);
        }

        return formula;
    }

    private List<String> extractGates(Circuit structure) {

        List<String> extractedGates = new ArrayList<>();

        for (Gate gate : structure.vertexSet()) {
            if (gate instanceof LogicGate && ((LogicGate) gate).getLogicType() != LogicType.OR2) {
                extractedGates.add(gate.getIdentifier());
            }
        }
        return extractedGates;
    }

    private Map<String, Set<String>> extractPairs(Circuit circuit) {

        Map<String, Set<String>> pairs = new HashMap<>();

        for (Gate gate : circuit.vertexSet()) {

            if (!(gate instanceof LogicGate))
                continue;

            LogicGate logicGate = (LogicGate) gate;

            if (logicGate.getLogicType() == LogicType.OR2)
                continue;

            for (Wire wire : circuit.outgoingEdgesOf(gate)) {

                Gate target = circuit.getEdgeTarget(wire);

                if (target instanceof LogicGate && ((LogicGate) target).getLogicType() != LogicType.OR2) {
                    pairs.putIfAbsent(gate.getIdentifier(), new HashSet<>());
                    pairs.get(gate.getIdentifier()).add(target.getIdentifier());
                }
            }
        }
        return pairs;
    }
}
