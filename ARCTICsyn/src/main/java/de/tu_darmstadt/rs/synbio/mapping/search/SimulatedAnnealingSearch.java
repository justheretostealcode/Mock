package de.tu_darmstadt.rs.synbio.mapping.search;

import de.tu_darmstadt.rs.synbio.common.*;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.circuit.Gate;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.common.library.GateRealization;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.mapping.MappingConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.assigner.ExhaustiveAssigner;
import de.tu_darmstadt.rs.synbio.mapping.assigner.RandomAssigner;
import de.tu_darmstadt.rs.synbio.simulation.SimulationConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.MappingResult;
import de.tu_darmstadt.rs.synbio.simulation.SimulatorInterfaceEnergy;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;
import java.util.stream.Collectors;

public class SimulatedAnnealingSearch extends AssignmentSearchAlgorithm {

    private static final Logger logger = LoggerFactory.getLogger(SimulatedAnnealingSearch.class);

    private Map<LogicType, List<GateRealization>> realizations;

    public SimulatedAnnealingSearch(Circuit structure, GateLibrary lib, MappingConfiguration mapConfig, SimulationConfiguration simConfig) {
        super(structure, lib, mapConfig, simConfig);
    }

    private static final boolean printTrajectory = false;
    private static final boolean useRadius = true;

    private double maxDistance = 0.0;
    private double minDistance = Double.MAX_VALUE;

    private enum Objective {
        SCORE_ONLY,
        SCORE,
        ENERGY,
        MULTI
    }

    private double lowerScore = 200.0;
    private double upperEnergy = 65000;//Double.POSITIVE_INFINITY;
    private final Objective objective = Objective.ENERGY;

    private double energyScaler;

    public void pareto() {

        SimulatorInterfaceEnergy simulator = new SimulatorInterfaceEnergy(simConfig, gateLib);
        simulator.initSimulation(structure);
        getInitialTemperature(simulator);
        simulator.shutdown();

        int iterations = 100;

        for (int i = 0; i < iterations; i++) {
            lowerScore = minScoreSample + (((maxScoreSample - minScoreSample)/iterations) * i);
            upperEnergy = maxEnergySample - (((maxEnergySample - minEnergySample)/iterations) * i);
            assign();
        }
    }

    public MappingResult assign() {

        File output;
        PrintWriter out;
        if (printTrajectory) {
            output = new File("traj_" + System.currentTimeMillis() + ".txt");
            try {
                out = new PrintWriter(output);
            } catch (Exception e) {
                e.printStackTrace();
                return null;
            }
        }

        // initialize gate lib

        realizations = gateLib.getRealizations();

        for (GateRealization r1 : realizations.get(LogicType.NOR2)) {
            for (GateRealization r2 : realizations.get(LogicType.NOR2)) {
                if (!r1.equals(r2) && r1.isCharacterized() && r2.isCharacterized()) {
                    double distance = r1.getCharacterization().getEuclidean(r2.getCharacterization(), gateLib.getProxNormalization(), gateLib.getProxWeights());
                    maxDistance = Math.max(maxDistance, distance);
                    minDistance = Math.min(minDistance, distance);
                }
            }
        }

        // initialize simulator

        SimulatorInterfaceEnergy simulator = new SimulatorInterfaceEnergy(simConfig, gateLib);
        simulator.initSimulation(structure);

        // get initial assignment

        ExhaustiveAssigner exhaustiveAssigner = new ExhaustiveAssigner(gateLib, structure);
        RandomAssigner randomAssigner = new RandomAssigner(gateLib, structure);

        Assignment current;
        long problemSize;
        double currentScore;
        double currentEnergy;

        do {
            current = randomAssigner.getNextAssignment();
            problemSize = exhaustiveAssigner.getNumTotalPermutations();
            currentScore = simulator.simulate(current);
            currentEnergy = simulator.getEnergy();
        } while (!current.isValid() ||
                (objective == Objective.SCORE && currentEnergy > upperEnergy) ||
                (objective == Objective.ENERGY && currentScore < lowerScore));

        // initialize search

        int acceptCount = 0;
        int simCount = 0;
        double acceptanceRate = 1.0;
        double radius = 1.0;

        Assignment neighbor;

        // annealing parameters

        double temperature = getInitialTemperature(simulator);
        double coolingFactor;

        int movesPerTemp = 2 * (int) Math.pow(problemSize, 0.5);

        while(acceptanceRate >= 0.25) {

            neighbor = getNeighbor(current, radius);

            if (neighbor == null)
                break;

            double neighborScore = simulator.simulate(neighbor);
            double neighborEnergy = simulator.getEnergy();
            simCount ++;

            if (accept(currentScore, neighborScore, currentEnergy, neighborEnergy, temperature)) {
                current = neighbor;
                currentScore = neighborScore;
                currentEnergy = neighborEnergy;
                acceptCount ++;
                //logger.info("accepted score: " + neighborScore + ", energy: " + neighborEnergy + ", rate: " + acceptanceRate + ", temp: " + temperature);
            }

            if (simCount % movesPerTemp == 0) {

                acceptanceRate = (double) acceptCount / movesPerTemp;
                acceptCount = 0;

                if (acceptanceRate > 0.9) {
                    coolingFactor = 0.8;
                } else if (acceptanceRate > 0.2) {
                    coolingFactor = 0.99;
                } else {
                    coolingFactor = 0.9;
                }

                temperature *= coolingFactor;

                if (acceptanceRate > 0.4)
                    radius *= 1.05;
                else
                    radius *= 0.95;

                if (radius > 1.0)
                    radius = 1.0;
                else if (radius < 0.0)
                    radius = 0.0;
            }

            if (printTrajectory) {
                out.println(simCount + "," + neighborScore + "," + temperature + "," + radius + "," + acceptanceRate);
                out.flush();
            }
        }

        if (printTrajectory) {
            out.close();
            ProcessBuilder pb = new ProcessBuilder(simConfig.getPythonBinary(), "./python_plotting/sa_trajectory/sa_trajectory.py", output.getAbsolutePath());
            try {
                pb.start();
            } catch (Exception e) {
                return null;
            }
        }

        MappingResult result = new MappingResult(structure, current, simulator.simulate(current), simulator.getEnergy());

        //logger.info(upperEnergy + ", " + result.getScore() + ", " + simulator.getEnergy());
        //logger.info(result.getScore() + ", " + simulator.getEnergy());
         result.setNeededSimulations(simCount);
        simulator.shutdown();
        return result;
    }

    private double minScoreSample;
    private double maxScoreSample;
    private double minEnergySample;
    private double maxEnergySample;

    private double getInitialTemperature(SimulatorInterfaceEnergy simulator) {

        int numSamples = 10;

        double[] scoreSample = new double[numSamples];
        double[] energySample = new double[numSamples];

        RandomAssigner assigner = new RandomAssigner(gateLib, structure);

        for (int i = 0; i < numSamples; i++) {

            Assignment ass = assigner.getNextAssignment();

            scoreSample[i] = simulator.simulate(ass);
            energySample[i] = simulator.getEnergy();
        }

        minScoreSample = Arrays.stream(scoreSample).min().getAsDouble();
        maxScoreSample = Arrays.stream(scoreSample).max().getAsDouble();
        minEnergySample = Arrays.stream(energySample).min().getAsDouble();
        maxEnergySample = Arrays.stream(energySample).max().getAsDouble();

        double sDevScore = new StandardDeviation().evaluate(scoreSample);
        double sDevEnergy = new StandardDeviation().evaluate(energySample);

        double sDev;

        switch (objective) {
            case ENERGY:
                sDev = sDevEnergy;
                break;
            case MULTI:
                sDev = sDevEnergy;//Math.max(sDevEnergy, sDevScore);
                break;
            default:
                sDev = sDevScore;
        }

        energyScaler = Arrays.stream(scoreSample).average().orElse(Double.NaN) / Arrays.stream(energySample).average().orElse(Double.NaN);

        return sDev * 20;
    }

    private Assignment getNeighbor(Assignment current, double radius) {

        Random rand = new Random();

        Assignment neighbor;

        List<Gate> gates = new ArrayList<>(current.keySet());

        int tryCount = 0;

        do {
            tryCount ++;

            if (tryCount > 1000000) {
                logger.warn("no neighbor found, aborting.");
                return null;
            }

            neighbor = new Assignment(current);
            Gate selectedCircuitGate = gates.get(rand.nextInt(gates.size()));
            GateRealization currentRealization = neighbor.get(selectedCircuitGate);

            List<GateRealization> realizationsOfType = this.realizations.get(selectedCircuitGate.getLogicType());

            if (realizationsOfType.size() == 1) {
                continue;
            }

            long numUsedOfType = neighbor.values().stream()
                    .filter(r -> r.getLogicType() == selectedCircuitGate.getLogicType())
                    .count();

            double swapTargetSelect = Math.random();
            double swapTargetSelectThreshold = (double) numUsedOfType / realizationsOfType.size();

            /* all available realizations used in circuit --> swap in circuit */
            if (numUsedOfType == realizationsOfType.size() || swapTargetSelect < swapTargetSelectThreshold) {

                List<Gate> swapCandidates = gates.stream()
                        .filter(g -> g.getLogicType() == selectedCircuitGate.getLogicType())
                        .collect(Collectors.toList());

                Gate selectedForSwap = swapCandidates.get(rand.nextInt(swapCandidates.size()));

                if (selectedForSwap == selectedCircuitGate)
                    continue;

                GateRealization realizationToSwap = neighbor.get(selectedForSwap);

                neighbor.put(selectedCircuitGate, realizationToSwap);
                neighbor.put(selectedForSwap, currentRealization);

            } else { /* more realizations available in library --> swap external */

                List<GateRealization> relizationCandidates = new ArrayList<>(realizationsOfType);
                relizationCandidates.remove(currentRealization); // remove current realization from list of candidates
                List<GateRealization> selectedRealizations = new ArrayList<>();

                if (useRadius) {

                    Map<GateRealization, Double> distances = new HashMap<>();

                    for (GateRealization r : relizationCandidates) {
                        distances.put(r, r.getCharacterization().getEuclidean(currentRealization.getCharacterization(),
                                gateLib.getProxNormalization(),
                                gateLib.getProxWeights()));
                    }

                    List<Map.Entry<GateRealization, Double>> distanceList = new LinkedList<>(distances.entrySet());
                    distanceList.sort(Map.Entry.comparingByValue());

                    for (Map.Entry<GateRealization, Double> e : distanceList) {
                        if (selectedRealizations.size() < 5 || e.getValue() <= radius * maxDistance)
                            selectedRealizations.add(e.getKey());
                    }

                } else {
                    selectedRealizations = relizationCandidates;
                }

                if (selectedRealizations.size() < 1)
                    continue;

                neighbor.put(selectedCircuitGate, selectedRealizations.get(rand.nextInt(selectedRealizations.size())));
            }

        } while (!neighbor.isValid() || (neighbor.equals(current)));

        return neighbor;
    }

    private boolean accept(double currentScore, double neighborScore, double currentEnergy, double neighborEnergy, double temperature) {

        switch (objective) {
            case SCORE_ONLY:

                if (neighborScore > currentScore)
                    return true;
                else
                    return Math.random() < Math.exp(-1*((currentScore - neighborScore) / temperature));

            case SCORE:

                if (neighborEnergy > upperEnergy)
                    return false;
                else if (neighborScore > currentScore)
                    return true;
                else
                    return Math.random() < Math.exp(-1*((currentScore - neighborScore) / temperature));

            case ENERGY:

                if (neighborScore < lowerScore)
                    return false;
                else if (neighborEnergy < currentEnergy)
                    return true;
                else
                    return Math.random() < Math.exp(-1*((neighborEnergy - currentEnergy) / temperature));

            case MULTI:

                double omegaScore = 1.0;
                double omegaEnergy = 1.0 * energyScaler;

                double ps = Math.exp((omegaScore * (neighborScore - currentScore)) / temperature);
                double pe = Math.exp((omegaEnergy * (currentEnergy - neighborEnergy)) / temperature);

                double pMax = Math.min(ps, pe);
                double p = Math.min(1.0, pMax);

                return Math.random() < p;
        }

        return true;
    }
}
