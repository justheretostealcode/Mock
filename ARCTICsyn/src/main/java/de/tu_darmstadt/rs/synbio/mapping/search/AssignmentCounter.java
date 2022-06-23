package de.tu_darmstadt.rs.synbio.mapping.search;

import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.mapping.MappingConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.assigner.ExhaustiveAssigner;
import de.tu_darmstadt.rs.synbio.mapping.compatibility.CompatibilityChecker;
import de.tu_darmstadt.rs.synbio.simulation.SimulationConfiguration;
import de.tu_darmstadt.rs.synbio.simulation.SimulationResult;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class AssignmentCounter extends AssignmentSearchAlgorithm {

    private static final Logger logger = LoggerFactory.getLogger(AssignmentCounter.class);

    private final ExhaustiveAssigner assigner;
    private final CompatibilityChecker checker;

    public AssignmentCounter(Circuit structure, GateLibrary lib, MappingConfiguration mapConfig, SimulationConfiguration simConfig) {
        super(structure, lib, mapConfig, simConfig);
        assigner = new ExhaustiveAssigner(lib, structure);
        checker = new CompatibilityChecker(lib, structure);
    }

    @Override
    public SimulationResult assign() {

        Assignment assignment;
        long assignments = 0;
        long valid = 0;

        do {
            assignment = assigner.getNextAssignment();

            if (assignment != null) {
                assignments ++;
                if (checker.checkSimple(assignment)) {
                    valid++;
                }
            }

        } while(assignment != null);

        logger.info(structure.getIdentifier() + "," + assignments + "," + valid);

        SimulationResult result = new SimulationResult(structure, null, 0.0);
        result.setNeededSimulations(assignments);

        return result;
    }
}
