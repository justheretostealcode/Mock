package de.tu_darmstadt.rs.synbio.mapping.search;

import de.tu_darmstadt.rs.synbio.mapping.MappingConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.compatibility.CompatibilityChecker;
import de.tu_darmstadt.rs.synbio.simulation.SimulationConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.mapping.assigner.ExhaustiveAssigner;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.simulation.SimulationResult;
import de.tu_darmstadt.rs.synbio.simulation.SimulatorInterface;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.concurrent.Callable;

public class ExhaustiveSearchWorker implements Callable<SimulationResult> {

    private static final Logger logger = LoggerFactory.getLogger(ExhaustiveSearchWorker.class);

    private final SimulatorInterface simulator;
    private final ExhaustiveAssigner assigner;
    private final Circuit structure;
    private final MappingConfiguration mapConfig;

    private final CompatibilityChecker checker;

    public ExhaustiveSearchWorker(ExhaustiveAssigner assigner, Circuit structure, MappingConfiguration mapConfig,
                                  SimulationConfiguration simConfig, GateLibrary gateLibrary) {
        this.mapConfig = mapConfig;
        this.simulator = new SimulatorInterface(simConfig, gateLibrary);
        this.assigner = assigner;
        this.structure = structure;

        this.checker = new CompatibilityChecker(gateLibrary, structure);
    }

    @Override
    public SimulationResult call() {

        Assignment assignment;
        simulator.initSimulation(structure);

        do {
            assignment = assigner.getNextAssignment();
        } while(assignment != null && !checker.verify(assignment));

        SimulationResult bestRes = null;
        long numSims = 0;

        while (assignment != null && !Thread.interrupted()) {
            SimulationResult result = new SimulationResult(structure, assignment, simulator.simulate(assignment, SimulatorInterface.PropagationMode.EXACT));

            numSims ++;

            //if (simulator.getLastGrowth() >= 0.75) {
                if (bestRes == null || (mapConfig.getOptimizationType().compare(bestRes.getScore(), result.getScore()))) {
                    bestRes = result;
                }
            //}

            do {
                assignment = assigner.getNextAssignment();
            } while(assignment != null && !checker.verify(assignment));
        }

        simulator.shutdown();

        if (bestRes != null)
            bestRes.setNeededSimulations(numSims);

        return bestRes;
    }
}
