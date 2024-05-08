package de.tu_darmstadt.rs.synbio.mapping.search;

import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.mapping.MappingConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.MappingResult;
import de.tu_darmstadt.rs.synbio.simulation.SimulationConfiguration;

import java.util.concurrent.Callable;

public class ParetoWorker implements Callable<MappingResult> {

    private final SimulatedAnnealingSearch annealer;

    public ParetoWorker(Circuit structure, GateLibrary lib, MappingConfiguration mapConfig, SimulationConfiguration simConfig, double lowerScore, double upperEnergy) {
        this.annealer = new SimulatedAnnealingSearch(structure, lib, mapConfig, simConfig, upperEnergy, lowerScore, SimulatedAnnealingSearch.Objective.SCORE);
    }

    @Override
    public MappingResult call() throws Exception {
        return annealer.assign();
    }

}
