package de.tu_darmstadt.rs.synbio.mapping.search;

import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.mapping.MappingConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.assigner.ExhaustiveAssigner;
import de.tu_darmstadt.rs.synbio.simulation.SimulationConfiguration;
import de.tu_darmstadt.rs.synbio.simulation.SimulationResult;

public class AssignmentCounter extends AssignmentSearchAlgorithm {

    private final ExhaustiveAssigner assigner;

    public AssignmentCounter(Circuit structure, GateLibrary lib, MappingConfiguration mapConfig, SimulationConfiguration simConfig) {
        super(structure, lib, mapConfig, simConfig);
        assigner = new ExhaustiveAssigner(lib, structure);
    }

    @Override
    public SimulationResult assign() {

        Assignment assignment;
        long assignments = 0;

        do {
            assignment = assigner.getNextAssignment();

            if (assignment != null && assignment.fulfilsConstraints(structure))
                assignments ++;

        } while(assignment != null);

        SimulationResult result = new SimulationResult(structure, null, 0.0);
        result.setNeededSimulations(assignments);

        return result;
    }
}