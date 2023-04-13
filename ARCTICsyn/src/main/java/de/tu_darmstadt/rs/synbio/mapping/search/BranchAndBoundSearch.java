package de.tu_darmstadt.rs.synbio.mapping.search;

import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.circuit.Gate;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.common.library.GateRealization;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.mapping.MappingConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.compatibility.CompatibilityChecker;
import de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.QueueItem;
import de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.SearchStatsLogger;
import de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.SearchTreeVisualizer;
import de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.searchstrategies.*;
import de.tu_darmstadt.rs.synbio.simulation.SimulationConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.MappingResult;
import de.tu_darmstadt.rs.synbio.simulation.SimulatorInterface;
import org.jgrapht.traverse.TopologicalOrderIterator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.stream.Collectors;

public class BranchAndBoundSearch extends AssignmentSearchAlgorithm {

    private static final Logger logger = LoggerFactory.getLogger(BranchAndBoundSearch.class);

    private final Map<LogicType, List<GateRealization>> realizations;
    private final SimulatorInterface simulator;
    private final CompatibilityChecker checker;
    private final Gate outputGate;
    private final GateRealization outputRealization;
    private final Gate[] logicGates;
    private final Gate[] reversedLogicGates;
    private Assignment initialAssignment;
    private double initialBestScore;
    private int iNeededSimulations;
    private final SearchTreeVisualizer searchTreeVisualizer;
    private final boolean bEagerBranchAndBound;
    private final boolean bVisualize;

    public BranchAndBoundSearch(Circuit structure, GateLibrary lib, MappingConfiguration mapConfig, SimulationConfiguration simConfig) {
        super(structure, lib, mapConfig, simConfig);

        if (mapConfig.getOptimizationType() == MappingConfiguration.OptimizationType.MINIMIZE)
            throw new Error("The Branch and Bound search is only capable of performing a maximization and not to minimize.");

        /*
        Realizations are sorted to improve performance of lazy branch and bound.
        Since the search strategy is lacking of any scoring information to guide which child to consider first, the sorting aids in exploiting beneficial areas first.
        For eager branch and bound the relevance is less, whereby the relevance of the sorting of the childs is much higher. Nevertheless, this sorting should lead to a sort of presorting.
         */
        realizations = gateLib.getRealizations();
        sortRealizations(realizations);

        checker = new CompatibilityChecker(gateLib, structure);

        simulator = new SimulatorInterface(simConfig, gateLib);
        simulator.initSimulation(structure);

        logicGates = getLogicGatesInTopologicalOrder(structure);
        reversedLogicGates = getReversedLogicGates(logicGates);

        outputGate = structure.getOutputGate();
        outputRealization = gateLib.getOutputDevice(outputGate.getLogicType());

        initialAssignment = null;
        initialBestScore = Double.NEGATIVE_INFINITY;

        bEagerBranchAndBound = mapConfig.getBabType() == MappingConfiguration.BAB_Type.EAGER;

        boolean bVisualize = mapConfig.getBabVisualization();

        searchTreeVisualizer = new SearchTreeVisualizer(structure.getIdentifier(), mapConfig, reversedLogicGates, bVisualize);
        this.bVisualize = searchTreeVisualizer.getbVisualize();
    }

    private Gate[] getLogicGatesInTopologicalOrder(Circuit structure) {
        ArrayList<Gate> gates = new ArrayList<>();
        Iterator<Gate> iterator = new TopologicalOrderIterator<>(structure);

        while (iterator.hasNext()) {
            Gate g = iterator.next();
            if (g.isLogicGate() || g.getLogicType() == LogicType.INPUT)
                gates.add(g);
        }
        Gate[] logicGates = new Gate[gates.size()];
        gates.toArray(logicGates);
        return logicGates;
    }

    private Gate[] getReversedLogicGates(Gate[] logicGates) {
        Gate[] reversedLogicGates = new Gate[logicGates.length];
        for (int iX = 0; iX < logicGates.length; iX++) {
            reversedLogicGates[reversedLogicGates.length - 1 - iX] = logicGates[iX];
        }
        return reversedLogicGates;
    }

    /**
     * Sorts the elements in the provided gate library. <br>
     * Ascending for DepthFirstSearch   <br>
     * Descending for BreadthFirstSearch <br>
     * By this, the gates having the largest ymax to ymin ratio are used first
     *
     * @param set
     */
    private void sortRealizations(Map<LogicType, List<GateRealization>> set) {
        MappingConfiguration.BAB_Sort_Order libraryOrder = mapConfig.getBabLibraryOrder();
        if (libraryOrder == MappingConfiguration.BAB_Sort_Order.UNSORTED)
            return;

        Comparator<GateRealization> comparator = (o1, o2) -> {
            if (o1.isCharacterized() && o2.isCharacterized())
                return (int) Math.signum(o2.getCharacterization().getMaxCelloScore() - o1.getCharacterization().getMaxCelloScore());    // Sort Descending
            else
                return 0;
        };

        List<GateRealization> realizations;
        for (LogicType type : LogicType.values()) {
            realizations = set.get(type);
            if (realizations == null)
                continue;

            if (libraryOrder == MappingConfiguration.BAB_Sort_Order.SORTED || libraryOrder == MappingConfiguration.BAB_Sort_Order.REVERSED) {
                realizations.sort(comparator);
                if (mapConfig.getBabSearchStrategy() == MappingConfiguration.BAB_SearchStrategy.DEPTH_FIRST_SEARCH)
                    // Reversion is necessary for depth first search, since the gate library is sorted descending based on the realizations possible scores (better first)
                    Collections.reverse(realizations);

                if (libraryOrder == MappingConfiguration.BAB_Sort_Order.REVERSED)
                    Collections.reverse(realizations);
            } else if (libraryOrder == MappingConfiguration.BAB_Sort_Order.SHUFFLED) {
                Collections.shuffle(realizations);
            }
        }
    }


    /**
     * A method for providing an initial assignment to the optimization algorithm.
     *
     * @param assignment The initial assignment
     */
    public void provideInitialAssignment(Assignment assignment) {
        if (assignment == null || !assignment.isValid() || assignment.size() != reversedLogicGates.length)
            return;

        for (Gate g : reversedLogicGates) {
            if (assignment.get(g) == null)
                return;
        }

        // Only set as initial assignment if the assignment is a valid assignment for the provided circuit
        this.initialAssignment = assignment;
    }

    /**
     * Enables to provide an initial best score to the algorithm enablign to prune branches worse than this score already at the start.
     *
     * @param score The score to use as initialBestScore
     */
    public void provideInitialBestScore(double score) {
        initialAssignment = null;
        initialBestScore = score;
    }

    @Override
    public MappingResult assign() {
        SearchStatsLogger searchStatsLogger = new SearchStatsLogger(structure, mapConfig, simConfig, reversedLogicGates.length);

        Comparator<QueueItem> comparator;
        if (mapConfig.getBabSearchStrategy() == MappingConfiguration.BAB_SearchStrategy.DEPTH_FIRST_SEARCH)
            // Sort ascending to improve Depth First Search
            comparator = (o1, o2) -> (int) Math.signum(o1.val - o2.val);
        else
            // Sort descending to improve Breath First Search
            comparator = (o1, o2) -> (int) Math.signum(o2.val - o1.val);

        // non-standard sorted cases:
        switch (mapConfig.getBabChildrenOrder()) {
            case SHUFFLED:
                // In the shuffled case, the branched assignments are shuffled and then not sorted later on
                // Therefore no break here since the following comparator does not apply any order after shuffling
            case UNSORTED:
                comparator = (o1, o2) -> 0;
                break;
            case REVERSED:
                comparator = comparator.reversed();
                break;
        }

        iNeededSimulations = 0;

        Assignment bestAssignment = initialAssignment;
        double bestScore = initialBestScore;

        if (bestAssignment != null)
            bestScore = bound(bestAssignment);

        SearchStrategy strategy = getSearchStrategy();

        strategy.addInitialItemToQueue(new Assignment(outputGate, outputRealization), Double.POSITIVE_INFINITY);

        long numberOfItemsAddedAndSkipped = 0;

        double errorThreshold = mapConfig.getBabFast() ? 0.0 : 0.001;

        QueueItem currentItem;
        Assignment currentAssignment;

        long iteration = 0;
        while ((currentItem = strategy.getNext()) != null) {    // Get the next element from queue and repeat until the queue is empty

            searchTreeVisualizer.add(currentItem, bestScore, false);

            //if (iteration % 100 == 0)  // Only update every hundred iterations
            //    System.out.print("\rIteration: " + iteration + " (" + iNeededSimulations + ")");

            if (currentItem.val <= bestScore * (1.0 - errorThreshold)) {
                // Removes item if the best score has changed after this queue item has been added to the queue
                numberOfItemsAddedAndSkipped++;
                continue;
            }

            currentAssignment = currentItem.assignment;
            List<Assignment> children;

            double val;

            if (bEagerBranchAndBound) {
                /*
                Eager Branch and Bound
                */
                children = branch(currentAssignment);

                if (children.size() == 0)
                    continue;

                int childSize = children.get(0).size();

                var ref = new Object() {
                    double highestScore = Double.NEGATIVE_INFINITY;
                };

                if (childSize < logicGates.length - 1) {    // The child assignments are intermediate nodes (no leaves)

                    double finalBestScore = bestScore  * (1.0 - errorThreshold);
                    List<QueueItem> queueItems = children.stream()
                            .map(assignment -> {                                    // Map to QueueItem
                                double dVal = bound(assignment);
                                ref.highestScore = Math.max(ref.highestScore, dVal);
                                return QueueItem.getQueueItem(assignment, dVal);})
                            .filter(item -> bVisualize || item.val > finalBestScore)              // Filter the assignments which do not have sufficient score (Filtering is skipped if visualisation is turned on, in order to result obtain complete Search Trees)
                            .sorted(comparator)// Ensures, that the search strategy visits the best node first
                            .collect(Collectors.toList());
                    strategy.addToQueue(queueItems);
                } else {   // The child assignments are leaves

                    for (Assignment assignment : children) {
                        val = bound(assignment);
                        ref.highestScore = Math.max(ref.highestScore, val);
                        searchTreeVisualizer.addLeafNode(assignment, val, bestScore, true);

                        if (val > bestScore) {   // Since the node is a terminal node (leave), one needs to check if it is better than the current best solution
                            bestScore = val;
                            bestAssignment = assignment;
                            //logger.info("new best score: " + bestScore);
                            searchStatsLogger.notifyNewBestAssignment(iteration, iNeededSimulations);
                        }
                    }
                }

                if (currentAssignment.size() > 0)
                    searchStatsLogger.addBoundingFunctionScore(currentAssignment.size(), currentItem.val, ref.highestScore);
                searchStatsLogger.addEntryToOverallSimulations(currentItem.val, children.size());

            } else {
                /*
                Lazy Branch and Bound
                */
                val = (currentAssignment.size() != 0) ? bound(currentAssignment) : Double.POSITIVE_INFINITY;
                if (val <= bestScore * (1.0 - errorThreshold))       // This check is required, if the best score has changed after this queue item has been added to the queue
                    continue;

                if (currentAssignment.size() < logicGates.length) {
                    // The assignment is an intermediate node and needs to be branched
                    children = branch(currentAssignment);

                    /*
                    The gate library is sorted ascending for depth first search and otherwise descending, based on the possible maximum scores of each gate.

                    Since depth first search adds the last to the front, the best gate will be used first, when it is added last.
                    For Breadth First Search the order is appropriate. However, Breadth First Search itself is not appropriate either.
                     */
                    double finalVal = val;
                    List<QueueItem> queueItems = children.stream()
                            .map(assignment -> QueueItem.getQueueItem(assignment, finalVal))
                            .collect(Collectors.toList());
                    strategy.addToQueue(queueItems);

                } else {
                    // The assignment is a leave node and needs to be compared to the current best solution
                    // Since the node is a terminal node (leave), one needs to check if it is better than the current best solution
                    if (val > bestScore) { // This check is complementary to the one above und should always be true in case the else is entered
                        bestScore = val;
                        bestAssignment = currentAssignment;
                        searchStatsLogger.notifyNewBestAssignment(iteration, iNeededSimulations);
                    }
                }
            }
            iteration++;
        }

        if (bestAssignment == null)
            return null;

        MappingResult result = new MappingResult(structure, bestAssignment, bestScore);
        result.setNeededSimulations(iNeededSimulations);
        result.setMinimumBranchAndBoundSimulations(minimalNumberOfSimulations(bestAssignment));

        searchStatsLogger.setNumberOfItemsAddedAndSkipped(numberOfItemsAddedAndSkipped);
        searchStatsLogger.setOptimalNumberOfSimulations(result.getMinimumBranchAndBoundSimulations());
        searchStatsLogger.setResult(result);
        searchStatsLogger.setMaximumNumberOfQueueItems(strategy.getMaximumNumberOfQueueEntries());
        searchStatsLogger.setAverageNumberOfQueueItems(strategy.getAverageNumberOfQueueEntries());

        if (mapConfig.getBabStatistics())
            searchStatsLogger.save("statistics/statistics_" + structure.getIdentifier() + "_" + mapConfig.toString().hashCode() + "_" + System.currentTimeMillis() + ".json");

        shutDown();

        return result;
    }


    /**
     * This function branches on a given assignment. <br>
     * This means, that all childs of the given assignment in reference to the search tree are created, while the group constraint is respected.
     *
     * @param assignment
     * @return
     */
    private List<Assignment> branch(Assignment assignment) {

        Gate logicGate = reversedLogicGates[assignment.size()];

        boolean checkNaive = false;

        long leftGatesOfType = Arrays.stream(reversedLogicGates)
                .filter(g -> !assignment.keySet().contains(g))
                .filter(g -> g.getLogicType() == logicGate.getLogicType()).count();

        Set<String> usedGroups = assignment.values().stream().map(GateRealization::getGroup).collect(Collectors.toSet());

        List<GateRealization> availableRealizations = realizations.get(logicGate.getLogicType()).stream()
                .filter(r -> !usedGroups.contains(r.getGroup())).collect(Collectors.toList());

        List<Assignment> assignments = new ArrayList<>();

        for (GateRealization realization : availableRealizations) {

            Assignment a = new Assignment(assignment);

            a.put(logicGate, realization);

            if (checkNaive ? !checker.checkSimple(a) : !checker.checkSat(a)) {
                //logger.info("suppressed branch of assignment with size " + assignment.size());
                continue;
            }

            /* if second last mapping determines mapping of last one --> jump straight to leaf nodes/complete assignments */
            if ((leftGatesOfType == 2) && (availableRealizations.size() == 2) && (assignment.size() == reversedLogicGates.length - 2)) {
                assignments.addAll(branch(a));
            } else {
                if (a.fulfilsConstraints(structure)) {
                    assignments.add(a);
                }
            }
        }

        if (mapConfig.getBabChildrenOrder() == MappingConfiguration.BAB_Sort_Order.SHUFFLED)
            Collections.shuffle(assignments);

        return assignments;
    }

    /**
     * The bounding function for the branch and bound algorithm. <br>
     * Thereby, the core functionality and design of the bounding function is given by the precomputed subproblems and intialised simulator interfaces.
     *
     * @param assignment The assignment to compute the bounding value for.
     * @return The bounding value.
     */
    private double bound(Assignment assignment) {
        iNeededSimulations++;

        SimulatorInterface.PropagationMode mode;

        if (assignment.size() == logicGates.length)
             mode = SimulatorInterface.PropagationMode.NORMAL;
        else
            mode = mapConfig.getBabFast() ? SimulatorInterface.PropagationMode.ITA : SimulatorInterface.PropagationMode.OPTIMAL;

        return simulator.simulate(assignment, mode);
    }

    /**
     * Returns the specified search strategy
     *
     * @return
     */
    private SearchStrategy getSearchStrategy() {
        // Get corresponding information from mapConfig
        SearchStrategy strategy;
        switch (mapConfig.getBabSearchStrategy()) {
            case DEPTH_FIRST_SEARCH:
                strategy = new DepthFirstSearch();
                break;
            case BREADTH_FIRST_SEARCH:
                strategy = new BreadthFirstSearch();
                break;
            case CYCLIC_BEST_FIRST_SEARCH:
                strategy = new CyclicBestFirstSearch(logicGates.length);
                break;
            default:
                strategy = new BestFirstSearch();
                break;
        }
        return strategy;
    }

    /**
     * Determines the number of simulations if the bounding function would be equivalent to max(childs.score) <br>
     * This value would be achieved with a perfect bounding function.
     *
     * @param assignment The assignment to determine the value for.
     * @return The minimum number of simulations.
     */
    private long minimalNumberOfSimulations(Assignment assignment) {
        long min = 0;
        long numberOfAvailableRealizations = 0;
        HashSet<String> usedGroups = new HashSet<>();

        for (Gate gate : reversedLogicGates) {
            numberOfAvailableRealizations = realizations.get(gate.getLogicType()).stream().filter(realization -> !usedGroups.contains(realization.getGroup())).count();
            min += numberOfAvailableRealizations;

            usedGroups.add(assignment.get(gate).getGroup());
        }

        return min;
    }

    /**
     * Stops the simulators and finalizes the searchtree visualisation.
     */
    private void shutDown() {
        searchTreeVisualizer.finish();
        shutDownInterfaces();
    }

    /**
     * Shuts down the used interfaces
     */
    private void shutDownInterfaces() {
        simulator.shutdown();
    }
}