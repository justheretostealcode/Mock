package de.tu_darmstadt.rs.synbio.mapping.assigner;

import de.tu_darmstadt.rs.synbio.common.*;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.circuit.Gate;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.common.library.GateRealization;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.stream.Collectors;

public class RandomAssigner implements Assigner {

    private static final Logger logger = LoggerFactory.getLogger(RandomAssigner.class);

    final private Random randomGenerator;

    // gates
    private final Map<LogicType, List<GateRealization>> availableGates;
    private final List<Gate> circuitGates;

    private final Gate outputGate;
    private final GateRealization outputRealization;

    public RandomAssigner(GateLibrary gateLib, Circuit circuit) {

        // initialize gate library
        this.availableGates = gateLib.getRealizations();

        // initialize circuit gate list
        this.circuitGates = circuit.vertexSet().stream().filter(g -> g.getLogicType() != LogicType.OUTPUT).collect(Collectors.toList());

        outputGate = circuit.getOutputBuffer();
        outputRealization = gateLib.getOutputDevice();

        this.randomGenerator = new Random();
    }

    public Assignment getNextAssignment() {

        Assignment assignment = new Assignment(outputGate, outputRealization);

        do {
            for (Gate gate : circuitGates) {

                if (availableGates.get(gate.getLogicType()).size() <= 1)
                    continue;

                int rand = randomGenerator.nextInt(availableGates.get(gate.getLogicType()).size());
                assignment.put(gate, availableGates.get(gate.getLogicType()).get(rand));
            }
        } while (!assignment.isValid());

        return assignment;
    }
}
