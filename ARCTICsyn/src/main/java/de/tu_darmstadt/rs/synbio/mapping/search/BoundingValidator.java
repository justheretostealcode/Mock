package de.tu_darmstadt.rs.synbio.mapping.search;

import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.circuit.Gate;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.mapping.MappingConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.assigner.ExhaustiveAssigner;
import de.tu_darmstadt.rs.synbio.simulation.SimulationConfiguration;
import de.tu_darmstadt.rs.synbio.simulation.SimulationResult;
import de.tu_darmstadt.rs.synbio.simulation.SimulatorInterface;
import org.jgrapht.traverse.TopologicalOrderIterator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class BoundingValidator extends AssignmentSearchAlgorithm {

    private static final Logger logger = LoggerFactory.getLogger(BoundingValidator.class);

    private final ExhaustiveAssigner assigner;
    private final SimulatorInterface simulator;

    private final List<Gate> gates;

    public BoundingValidator(Circuit structure, GateLibrary lib, MappingConfiguration mapConfig, SimulationConfiguration simConfig) {
        super(structure, lib, mapConfig, simConfig);
        assigner = new ExhaustiveAssigner(lib, structure);
        simulator = new SimulatorInterface(simConfig, lib);

        gates = new ArrayList<>();
        Iterator<Gate> iterator = new TopologicalOrderIterator<>(structure);

        while (iterator.hasNext()) {
            Gate g = iterator.next();
            if (g.isLogicGate() || g.getLogicType() == LogicType.INPUT)
                gates.add(g);
        }
    }


    @Override
    public SimulationResult assign() {

        Assignment assignment;
        simulator.initSimulation(structure);
        long numSims = 0;

        do {
            assignment = assigner.getNextAssignment();
        } while (assignment != null && !assignment.fulfilsConstraints(structure));

        while (assignment != null) {

            Assignment completeAssignment = new Assignment(assignment);
            double exactScore = simulator.simulate(completeAssignment, SimulatorInterface.PropagationMode.NORMAL);

            for (int i = 0; i < gates.size() - 1; i++) {

                Gate gate = gates.get(i);

                assignment.remove(gate);
                double estimation = simulator.simulate(assignment, SimulatorInterface.PropagationMode.ITA_OPTIMAL);

                if (estimation < exactScore) {
                    logger.warn("Estimation is not optimistic! Aborting check for " + structure.getIdentifier() + ".");
                    logger.warn("Exact result: " + exactScore + ", " + completeAssignment.toString());
                    logger.warn("Estimation: " + estimation + ", " + assignment.toString());

                    simulator.shutdown();
                    return null;
                }

                numSims++;
            }

            do {
                assignment = assigner.getNextAssignment();
            } while (assignment != null && !assignment.fulfilsConstraints(structure));
        }

        simulator.shutdown();

        SimulationResult dummyResult = new SimulationResult(structure, null, 0.0);
        dummyResult.setNeededSimulations(numSims);

        logger.info("Bounding validated successfully for " + structure.getIdentifier() + ".");

        return dummyResult;
    }
}