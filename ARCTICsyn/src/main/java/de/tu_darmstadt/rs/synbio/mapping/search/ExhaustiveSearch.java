package de.tu_darmstadt.rs.synbio.mapping.search;

import de.tu_darmstadt.rs.synbio.mapping.MappingConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.assigner.ExhaustiveAssigner;
import de.tu_darmstadt.rs.synbio.simulation.SimulationConfiguration;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.simulation.SimulationResult;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class ExhaustiveSearch extends AssignmentSearchAlgorithm {

    private static final Logger logger = LoggerFactory.getLogger(ExhaustiveSearch.class);

    private final ExhaustiveAssigner assigner;

    public ExhaustiveSearch(Circuit structure, GateLibrary lib, MappingConfiguration mapConfig, SimulationConfiguration simConfig) {
        super(structure, lib, mapConfig, simConfig);
        assigner = new ExhaustiveAssigner(lib, structure);
    }

    public SimulationResult assign() {

        List<ExhaustiveSearchWorker> workers = new ArrayList<>();
        int maxThreads = Runtime.getRuntime().availableProcessors() - 1;
        int availableProcessors = simConfig.getSimLimitThreadsNum() != 0 ? Math.min(simConfig.getSimLimitThreadsNum(), maxThreads) : maxThreads;

        //logger.info("Simulating \"" + structure.getIdentifier() + "\" (up to " + assigner.getNumTotalPermutations() + " assignments) with " + availableProcessors + " threads");

        for (int i = 0; i < availableProcessors; i ++) {
            workers.add(new ExhaustiveSearchWorker(assigner, structure, mapConfig, simConfig, gateLib));
        }

        ExecutorService executor = Executors.newFixedThreadPool(availableProcessors);

        List<Future<SimulationResult>> simResults = Collections.emptyList();

        try {
            simResults = executor.invokeAll(workers);
        } catch (InterruptedException e) {
            logger.error(e.getMessage());
        }

        SimulationResult bestRes = null;

        long numSims = 0;

        for (Future<SimulationResult> result : simResults) {

            try {
                SimulationResult res = result.get();

                if (res != null) {

                    numSims += res.getNeededSimulations();

                    if (bestRes == null || (mapConfig.getOptimizationType().compare(bestRes.getScore(), res.getScore()))) {
                        bestRes = res;
                    }
                }

            } catch (Exception e) {
                logger.error(e.getMessage());
            }
        }

        //logger.info("Finished simulating " + structure.getIdentifier() + ", score: " + (bestRes != null ? bestRes.getScore() : 0));

        executor.shutdownNow();

        if (bestRes != null)
            bestRes.setNeededSimulations(numSims);

        return bestRes;
    }
}
