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
import de.tu_darmstadt.rs.synbio.mapping.compatibility.CompatibilityChecker;
import de.tu_darmstadt.rs.synbio.simulation.SimulationConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.MappingResult;
import de.tu_darmstadt.rs.synbio.simulation.SimulatorInterface;
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

    private static final boolean printTrajectory = true;
    private static final boolean useRadius = true;

    private double maxDistance = 0.0;
    private double minDistance = Double.MAX_VALUE;

    private CompatibilityChecker checker;

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

        SimulatorInterface simulator = new SimulatorInterface(simConfig, gateLib);
        simulator.initSimulation(structure);

        // get initial assignment

        ExhaustiveAssigner exhaustiveAssigner = new ExhaustiveAssigner(gateLib, structure);
        RandomAssigner randomAssigner = new RandomAssigner(gateLib, structure);

        // compatibility checker
        checker = new CompatibilityChecker(gateLib, structure);

        Assignment current;
        long problemSize;
        double currentScore;
        double currentGrowth;

        do {
            current = randomAssigner.getNextAssignment();
            problemSize = exhaustiveAssigner.getNumTotalPermutations();
            currentScore = simulator.simulate(current, SimulatorInterface.PropagationMode.NORMAL);
            currentGrowth = 1.0;//simulator.getLastGrowth();
            currentScore = currentScore * (currentGrowth < 0.75 ? Math.pow(currentGrowth * 1.33, 1) : 1.0);
        } while (!checker.checkSimple(current) || !current.fulfilsConstraints(structure) || currentGrowth < 0.75);

        // initialize search

        int acceptCount = 0;
        int simCount = 0;
        double acceptanceRate = 1.0;
        double radius = 1.0;

        Assignment neighbor;

        Assignment best = current;
        double bestScore = currentScore;

        List<Double> scoreHistory = new LinkedList<>();
        final int historySize = 2000;

        // annealing parameters

        double temperature = getInitialTemperature(simulator);
        double coolingFactor;

        int movesPerTemp = 10 * (int) Math.pow(problemSize, 0.25);

        double histBestSDev = Math.abs(bestScore);

        while (histBestSDev > bestScore * 0.001 || acceptanceRate > 0.2) {

            neighbor = getNeighbor(current, radius);

            if (neighbor == null)
                break;

            double neighborScore = simulator.simulate(neighbor, SimulatorInterface.PropagationMode.NORMAL);
            double neighborGrowth = 1.0;//simulator.getLastGrowth();
            neighborScore =  neighborScore * (neighborGrowth < 0.75 ? Math.pow(neighborGrowth * 1.33, 1) : 1.0);
            simCount ++;

            if (accept(currentScore, neighborScore, temperature)) {
                current = neighbor;
                currentScore = neighborScore;
                currentGrowth = neighborGrowth;
                acceptCount ++;
            }

            scoreHistory.add(0, currentScore);
            if (scoreHistory.size() > historySize) {
                scoreHistory.remove(scoreHistory.size() - 1);
            }

            if (mapConfig.getOptimizationType().compare(bestScore, currentScore) && currentGrowth >= 0.75) {
                best = current;
                bestScore = currentScore;
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

                //radius *= (1 + ((acceptanceRate - 0.4)));

                if (radius > 1.0)
                    radius = 1.0;
                else if (radius < 0.0)
                    radius = 0.0;

                histBestSDev = getBestSD(scoreHistory);
            }

            if (printTrajectory) {
                out.println(simCount + "," + bestScore + "," + temperature + "," + radius);
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

        MappingResult result = new MappingResult(structure, best, simulator.simulate(best, SimulatorInterface.PropagationMode.NORMAL));
        result.setNeededSimulations(simCount);
        simulator.shutdown();
        return result;
    }

    private double getBestSD(List<Double> history) {

        StandardDeviation sDev = new StandardDeviation();

        List<Double> histBest = new ArrayList<>(history);
        histBest.sort(Double::compareTo);
        histBest = histBest.subList((int) (histBest.size() * 0.9), histBest.size() - 1);

        double[] histBestArray = new double[histBest.size()];
        for (int i = 0; i < histBest.size(); i++)
            histBestArray[i] = histBest.get(i);

        return sDev.evaluate(histBestArray);
    }

    private double getInitialTemperature(SimulatorInterface simulator) {

        int numSamples = 1000;

        double[] scores = new double[numSamples];
        RandomAssigner assigner = new RandomAssigner(gateLib, structure);

        for (int i = 0; i < numSamples; i++) {
            scores[i] = simulator.simulate(assigner.getNextAssignment(), SimulatorInterface.PropagationMode.NORMAL);
        }

        double sDev = new StandardDeviation().evaluate(scores);

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

            /* all available realizations used in circuit --> swap in circuit */
            if (numUsedOfType == realizationsOfType.size()) {

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

        } while (!neighbor.isValid() || !checker.checkSimple(neighbor) || !neighbor.fulfilsConstraints(structure) || (neighbor.equals(current)));

        return neighbor;
    }

    private boolean accept(double current, double neighbor, double temperature) {

        if (mapConfig.getOptimizationType().compare(current, neighbor))
            return true;

        if (mapConfig.getOptimizationType().equals(MappingConfiguration.OptimizationType.MAXIMIZE)) {
            return Math.random() < Math.exp(-1*((current - neighbor) / temperature));
        } else {
            return Math.random() < Math.exp(-1*((neighbor - current) / temperature));
        }
    }
}
