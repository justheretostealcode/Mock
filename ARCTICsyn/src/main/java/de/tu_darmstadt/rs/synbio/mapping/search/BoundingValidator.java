package de.tu_darmstadt.rs.synbio.mapping.search;

import com.fasterxml.jackson.databind.ObjectMapper;
import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.circuit.Gate;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.mapping.MappingConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.assigner.ExhaustiveAssigner;
import de.tu_darmstadt.rs.synbio.mapping.compatibility.CompatibilityChecker;
import de.tu_darmstadt.rs.synbio.simulation.AssignmentCompiler;
import de.tu_darmstadt.rs.synbio.simulation.SimulationConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.MappingResult;
import de.tu_darmstadt.rs.synbio.simulation.SimulatorInterfaceThermo;
import org.jgrapht.traverse.TopologicalOrderIterator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class BoundingValidator extends AssignmentSearchAlgorithm {

    private static final Logger logger = LoggerFactory.getLogger(BoundingValidator.class);

    private final ExhaustiveAssigner assigner;
    private final SimulatorInterfaceThermo simulator;
    private final AssignmentCompiler compiler;
    private final CompatibilityChecker checker;
    private final ObjectMapper mapper = new ObjectMapper();


    private final List<Gate> gates;

    public BoundingValidator(Circuit structure, GateLibrary lib, MappingConfiguration mapConfig, SimulationConfiguration simConfig) {
        super(structure, lib, mapConfig, simConfig);
        assigner = new ExhaustiveAssigner(lib, structure);
        simulator = new SimulatorInterfaceThermo(simConfig, lib);

        gates = new ArrayList<>();
        Iterator<Gate> iterator = new TopologicalOrderIterator<>(structure);

        while (iterator.hasNext()) {
            Gate g = iterator.next();
            if (g.isLogicGate() || g.getLogicType() == LogicType.INPUT)
                gates.add(g);
        }

        compiler = new AssignmentCompiler(structure, lib);
        checker = new CompatibilityChecker(lib, structure);
    }


    @Override
    public MappingResult assign() {

        Assignment assignment;
        simulator.initSimulation(structure);
        long numSims = 0;

        double biggestError = 0.0;
        int num= 0;

        do {
            assignment = assigner.getNextAssignment();
        } while (assignment != null && !checker.checkSimple(assignment));

        while (assignment != null) {

            Assignment completeAssignment = new Assignment(assignment);
            Double exactScore = simulator.simulate(completeAssignment, SimulatorInterfaceThermo.PropagationMode.NORMAL);

            for (int i = 0; i < gates.size() - 1; i++) {

                Gate gate = gates.get(i);

                assignment.remove(gate);
                Double estimation = simulator.simulate(assignment, SimulatorInterfaceThermo.PropagationMode.OPTIMAL);

                if (exactScore == null || estimation == null)
                    continue;

                if (estimation < exactScore) {
                    /*logger.warn("Estimation is not optimistic! Aborting check for " + structure.getIdentifier() + ".");
                    logger.warn("Exact result: " + exactScore + ", " + completeAssignment);
                    logger.warn("Estimation: " + estimation + ", " + assignment);

                    simulator.shutdown();
                    return null;*/

                    double error = exactScore - estimation;

                    if (error > biggestError) {
                        biggestError = error;
                        logger.warn(biggestError + " | " + biggestError / exactScore);

                        /*Map<String, Map<?, ?>> completeAssignmentMap = compiler.compile(completeAssignment);
                        Map<String, Map<?, ?>> assignmentMap = compiler.compile(assignment);

                        try {
                            mapper.writerWithDefaultPrettyPrinter().writeValue(new File("assignments_5", "assignment_" + num + "_full_" + structure.getIdentifier() + "_" + error + ".json"), completeAssignmentMap);
                            mapper.writerWithDefaultPrettyPrinter().writeValue(new File("assignments_5", "assignment_" + num + "_partial_" + structure.getIdentifier() + "_" + error + ".json"), assignmentMap);
                            num++;
                        } catch (Exception e) {}*/
                    }
                }

                numSims++;
            }

            do {
                assignment = assigner.getNextAssignment();
            } while (assignment != null && !checker.checkSimple(assignment));
        }

        simulator.shutdown();

        MappingResult dummyResult = new MappingResult(structure, null, 0.0);
        dummyResult.setNeededSimulations(numSims);

        logger.info("Bounding validated successfully for " + structure.getIdentifier() + ".");

        return dummyResult;
    }
}