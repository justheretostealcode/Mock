package de.tu_darmstadt.rs.synbio.mapping.search;

import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.circuit.Gate;
import de.tu_darmstadt.rs.synbio.common.circuit.LogicGate;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.common.library.GateRealization;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.mapping.MappingConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.QueueItem;
import de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.SearchStatsLogger;
import de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.SearchTreeVisualizer;
import de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.searchstrategies.*;
import de.tu_darmstadt.rs.synbio.mapping.util.BranchAndBoundUtil;
import de.tu_darmstadt.rs.synbio.simulation.SimulationConfiguration;
import de.tu_darmstadt.rs.synbio.simulation.SimulationResult;
import de.tu_darmstadt.rs.synbio.simulation.SimulatorInterface;
import org.jgrapht.traverse.TopologicalOrderIterator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.stream.Collectors;

public class BranchAndBoundSearch extends AssignmentSearchAlgorithm {

    private static final Logger logger = LoggerFactory.getLogger(BranchAndBoundSearch.class);

    private HashMap<LogicType, List<GateRealization>> realizations;

    private Circuit[] subproblems;
    private SimulatorInterface[] interfaces;

    private LogicGate[] logicGates;
    private LogicGate[] reversedLogicGates;

    private Assignment initialAssignment;
    private double initialBestScore;

    private int iNeededSimulations;

    private SearchTreeVisualizer searchTreeVisualizer;

    private boolean bEagerBranchAndBound;
    private boolean bVisualize;

    private boolean bFastMode;
    private String substitutionMode;

    private MappingConfiguration.BAB_INPUT_SPECIFICATION_TYPE babInputSpecificationType;


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


        // Bei n logic Gattern werden n-1 Teilprobleme erzeugt
        subproblems = getSubproblems(this.structure);
        interfaces = getSubproblemInterfaces(this.structure, this.subproblems);


        logicGates = getLogicGatesInTopologicalOrder(structure);
        reversedLogicGates = getReversedLogicGates(logicGates);

        initialAssignment = null;
        initialBestScore = Double.NEGATIVE_INFINITY;


        if (mapConfig.getBabType() != MappingConfiguration.BAB_Type.EAGER && mapConfig.getBabType() != MappingConfiguration.BAB_Type.LAZY)
            throw new Error(mapConfig.getBabType().name() + " is unknown to BranchAndBound");

        bEagerBranchAndBound = (mapConfig.getBabType() == MappingConfiguration.BAB_Type.EAGER) ? true : false;


        boolean bVisualize = mapConfig.getBabVisualization();

        searchTreeVisualizer = new SearchTreeVisualizer(structure.getIdentifier(), mapConfig, reversedLogicGates, bVisualize);
        this.bVisualize = searchTreeVisualizer.getbVisualize();

        bFastMode = mapConfig.getBabFast();

        /*
        If bFastMode == true -> No substitution is performed, leading to a bounding function which is not strictly greater equal
        If bFastMode == false -> Substitution is performed leading to a mathematical correct description but also decreasing the performance
         */
        substitutionMode = (bFastMode) ? "0" : "1";

        babInputSpecificationType = mapConfig.getBabInputSpecificationType();

        int iX = 0;
    }

    /**
     * Derives the subproblems of the provided structure
     *
     * @param structure The structure to consider
     * @return An array of circuits
     */
    private Circuit[] getSubproblems(Circuit structure) {
        List<Circuit> subproblems = BranchAndBoundUtil.getSubproblems(structure);
        Collections.reverse(subproblems);   // Reverse the order to achieve, that the index corresponds to the distance from the source in the search tree

        Circuit[] subproblemsA = new Circuit[subproblems.size()];
        subproblems.toArray(subproblemsA);
        return subproblemsA;
    }

    /**
     * Creates the subproblems and the corresponding simulator interfaces and returns the interfaces. <br>
     * The first element in the returned array belongs to the smallest subproblem. <br>
     * The last element in the returned array is the interface for the provided attribute @structure.
     *
     * @param structure The structure to crate the subproblems interfaces for
     * @return Simulator Interfaces to all subproblem of the structure and the structure itself
     */
    private SimulatorInterface[] getSubproblemInterfaces(Circuit structure, Circuit[] subproblems) {
        ArrayList<SimulatorInterface> interfaces = new ArrayList<>(structure.edgeSet().size() - 2);

        SimulatorInterface simulator = null;
        for (Circuit subproblem : subproblems) {
            simulator = new SimulatorInterface(simConfig, gateLib.getSourceFile());
            simulator.initSimulation(subproblem);
            interfaces.add(simulator);
        }

        SimulatorInterface[] array = new SimulatorInterface[interfaces.size()];
        interfaces.toArray(array);
        return array;
    }

    private LogicGate[] getLogicGatesInTopologicalOrder(Circuit structure) {
        ArrayList<Gate> gates = new ArrayList<>();
        Iterator<Gate> iterator = new TopologicalOrderIterator(structure);

        for (Iterator it = iterator; it.hasNext(); ) {
            Gate g = (Gate) it.next();
            if (g.getType() == Gate.Type.LOGIC)
                gates.add(g);
        }
        LogicGate[] logicGates = new LogicGate[gates.size()];
        gates.toArray(logicGates);
        return logicGates;
    }

    private LogicGate[] getReversedLogicGates(LogicGate[] logicGates) {
        LogicGate[] reversedLogicGates = new LogicGate[logicGates.length];
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
    private void sortRealizations(HashMap<LogicType, List<GateRealization>> set) {
        MappingConfiguration.BAB_Sort_Order libraryOrder = mapConfig.getBabLibraryOrder();
        if (libraryOrder == MappingConfiguration.BAB_Sort_Order.UNSORTED)
            return;

        Comparator<GateRealization> comparator = new Comparator<GateRealization>() {
            @Override
            public int compare(GateRealization o1, GateRealization o2) {
                if (o1.isCharacterized() && o2.isCharacterized())
                    return (int) Math.signum(o2.getCharacterization().getMaxCelloScore() - o1.getCharacterization().getMaxCelloScore());    // Sort Descending
                else
                    return 0;
            }
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

        for (LogicGate g : reversedLogicGates) {
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
    public SimulationResult assign() {
        SearchStatsLogger searchStatsLogger = new SearchStatsLogger(structure, mapConfig, simConfig, reversedLogicGates.length);

        Comparator<QueueItem> comparator;
        if (mapConfig.getBabSearchStrategy() == MappingConfiguration.BAB_SearchStrategy.DEPTH_FIRST_SEARCH)
            // Sort ascending to improve Depth First Search
            comparator = new Comparator<QueueItem>() {
                @Override
                public int compare(QueueItem o1, QueueItem o2) {
                    return (int) Math.signum(o1.val - o2.val);
                }
            };
        else
            // Sort descending to improve Breath First Search
            comparator = new Comparator<QueueItem>() {
                @Override
                public int compare(QueueItem o1, QueueItem o2) {
                    return (int) Math.signum(o2.val - o1.val);
                }
            };


        switch (mapConfig.getBabChildrenOrder()) {
            case SORTED:
                comparator = comparator;
                break;
            case SHUFFLED:
                // In the shuffled case, the branched assignments are shuffled and then not sorted later on
                // Therefore no break here since the following comparator does not apply any order after shuffling
            case UNSORTED:
                comparator = new Comparator<QueueItem>() {
                    @Override
                    public int compare(QueueItem o1, QueueItem o2) {
                        return 0;
                    }
                };
                break;
            case REVERSED:
                comparator = comparator.reversed();
                break;
            default:
                comparator = comparator;
                break;
        }


        iNeededSimulations = 0;

        Assignment bestAssignment = initialAssignment;
        double bestScore = initialBestScore;

        if (bestAssignment != null)
            bestScore = bound(bestAssignment);

        SearchStrategy strategy = getSearchStrategy();

        strategy.addInitialItemToQueue(Double.POSITIVE_INFINITY);

        long numberOfItemsAddedAndSkipped = 0;

        QueueItem currentItem;
        Assignment currentAssignment;
        System.out.print("Iteration: " + 0);
        long iteration = 0;
        while ((currentItem = strategy.getNext()) != null) {    // Get the next element from queue and repeat until the queue is empty

            searchTreeVisualizer.add(currentItem, bestScore);

            if (iteration % 100 == 0)  // Only update every hundred iterations
                System.out.print("\rIteration: " + iteration + " (" + iNeededSimulations + ")");

            if (currentItem.val <= bestScore) {
                // Removes item if the best score has changed after this queue item has been added to the queue
                numberOfItemsAddedAndSkipped++;
                continue;
            }

            currentAssignment = currentItem.assignment;
            List<Assignment> childs;

            double val;


            if (bEagerBranchAndBound) {
                /*
                Eager Branch and Bound
                */
                childs = branch(currentAssignment);

                var ref = new Object() {
                    double highestScore = Double.NEGATIVE_INFINITY;
                };
                double boundingFunctionScore = Double.NEGATIVE_INFINITY;


                if (currentAssignment.size() < logicGates.length - 1) {    // The child assignments are intermediate nodes (no leaves)


                    double finalBestScore = bestScore;
                    List<QueueItem> queueItems = childs.stream()
                            .map(assignment -> {                                    // Map to QueueItem

                                double dVal = bound(assignment);

                                ref.highestScore = Math.max(ref.highestScore, dVal);
                                return QueueItem.getQueueItem(assignment, dVal);
                            })
                            .filter(item -> bVisualize || item.val > finalBestScore)              // Filter the assignments which do not have sufficient score (Filtering is skipped if visualisation is turned on, in order to result obtain complete Search Trees)
                            .sorted(comparator)// Ensures, that the search strategy visits the best node first
                            .collect(Collectors.toList());
                    strategy.addToQueue(queueItems);
                } else {   // The child assignments are leaves

                    for (Assignment assignment : childs) {
                        val = bound(assignment);
                        ref.highestScore = Math.max(ref.highestScore, val);
                        searchTreeVisualizer.addLeafNode(assignment, val, bestScore);
                        if (val > bestScore) {   // Since the node is a terminal node (leave), one needs to check if it is better than the current best solution
                            bestScore = val;
                            bestAssignment = assignment;
                            searchStatsLogger.notifyNewBestAssignment(iteration, iNeededSimulations);
                        }
                    }
                }

                boundingFunctionScore = currentItem.val / ref.highestScore;
                if (currentAssignment.size() > 0)
                    searchStatsLogger.addBoundingFunctionScore(currentAssignment.size(), currentItem.val, ref.highestScore);
                searchStatsLogger.addEntryToOverallSimulations(currentItem.val, childs.size());

                // Activate Code in case the error is required
//                if (bFastMode && boundingFunctionScore < 1)
//                    logger.error(String.format("ERROR %s: Boundingfunction score is below 1. Score: %f, Parent Item Score: %f, Child Item Score: %f.", structure.getIdentifier(), boundingFunctionScore, currentItem.val, ref.highestScore));

            } else {
                /*
                Lazy Branch and Bound
                */
                val = (currentAssignment.size() != 0) ? bound(currentAssignment) : Double.POSITIVE_INFINITY;
                if (val <= bestScore)       // This check is required, if the best score has changed after this queue item has been added to the queue
                    continue;

                if (currentAssignment.size() < logicGates.length) {
                    // The assignment is an intermediate node and needs to be branched
                    childs = branch(currentAssignment);

                    /*
                    The gate library is sorted ascending for depth first search and otherwise descending, based on the possible maximum scores of each gate.

                    Since depth first search adds the last to the front, the best gate will be used first, when it is added last.
                    For Breadth First Search the order is appropriate. However, Breadth First Search itself is not appropriate either.
                     */
                    double finalVal = val;
                    List<QueueItem> queueItems = childs.stream()
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

        SimulationResult result = new SimulationResult(structure, bestAssignment, bestScore);
        result.setNeededSimulations(iNeededSimulations);
        result.setMinimumBranchAndBoundSimulations(minimalNumberOfSimulations(bestAssignment));

        searchStatsLogger.setNumberOfItemsAddedAndSkipped(numberOfItemsAddedAndSkipped);
        searchStatsLogger.setOptimalNumberOfSimulations(result.getMinimumBranchAndBoundSimulations());
        searchStatsLogger.setResult(result);
        searchStatsLogger.setMaximumNumberOfQueueItems(strategy.getMaximumNumberOfQueueEntries());
        searchStatsLogger.setAverageNumberOfQueueItems(strategy.getAverageNumberOfQueueEntries());
        int randVal = (int) Math.random() * 1000;
        if (mapConfig.getBabStatistics())
            searchStatsLogger.save("statistics/statistics_" + structure.getIdentifier() + "_" + mapConfig.toString().hashCode() + "_" + randVal + "_" + System.currentTimeMillis() + ".json");



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
        int index = assignment.size();
        LogicGate logicGate = reversedLogicGates[index];
        List<Assignment> assignments = new ArrayList<>();

        Set<String> usedGroups = assignment.values().stream().map(gateRealization -> gateRealization.getGroup()).collect(Collectors.toSet());

        List<GateRealization> availableRealizations = realizations.get(logicGate.getLogicType());

        availableRealizations.stream()
                .filter(gateRealization -> !usedGroups.contains(gateRealization.getGroup()))
                .forEach(gateRealization -> {
                    Assignment a = new Assignment(assignment);
                    a.put(logicGate, gateRealization);
                    assignments.add(a);
                });

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
        int assignmentSize = assignment.size();
        String additionalArgs;


        if (assignmentSize != interfaces.length)    // The assignment to simulate does not belong to a leaf node
            additionalArgs = getAdditionalSimArgs(assignment);  // Use optimistic values if not leaf node
        else
            additionalArgs = " custom_input_specification=0 substitute=0 "; // Use default values if leaf node

        double score = interfaces[assignmentSize - 1].simulate(assignment, additionalArgs);
        return score;
    }

    /**
     * This method generates a custom input specification to enable the bounding function to be as tightly as currently possible.
     *
     * @param assignment The assignment for which the custom input specification needs to be determined
     * @return A string consisting of the custom input specification
     */
    private String getAdditionalSimArgs(Assignment assignment) {
        String additionalArgs = null;

        int subproblemPointer = assignment.size() - 1;
        int problemPointer = subproblems.length - 1;
        if (subproblemPointer == problemPointer)
            return " cis=0 ";

        /*
        Determine which input gates of the original circuit are present in the actual circuit.
         */
        Circuit subproblem = subproblems[subproblemPointer]; // assignment.size() sollte immer größer 0 sein
        Circuit problem = subproblems[problemPointer];

        // The original input buffer ids
        Set<String> problemInputBufferIDs = problem.getInputBuffers()
                .stream()
                .map(Gate::getIdentifier)
                .collect(Collectors.toSet());

        List<Gate> inputGates = subproblem.getInputBuffers();
        // The input buffer ids not present in the original problem
        Set<String> artificialInputBufferIDs = inputGates.stream()
                .map(Gate::getIdentifier)
                .filter(ident -> !problemInputBufferIDs.contains(ident))
                .collect(Collectors.toSet());


        double maxVal = Double.NEGATIVE_INFINITY;
        double minVal = Double.POSITIVE_INFINITY;
        // Case differentiation for Input Specification Type
        if (babInputSpecificationType == MappingConfiguration.BAB_INPUT_SPECIFICATION_TYPE.INPRECISE) {
            maxVal = Math.pow(10, 12);  // Double.POSITIVE_INFINITY can not be transferred to the Simulator -> 10^12
            minVal = 0;
            artificialInputBufferIDs.addAll(BranchAndBoundUtil.CELLO_INPUT_SPECIFICATION.keySet());
        } else {
            if (babInputSpecificationType == MappingConfiguration.BAB_INPUT_SPECIFICATION_TYPE.EQUAL) {
                Collection<Map<Boolean, Double>> inputSpec = BranchAndBoundUtil.CELLO_INPUT_SPECIFICATION.values();

                maxVal = inputSpec.stream().mapToDouble(booleanDoubleMap -> booleanDoubleMap.get(true)).max().getAsDouble();
                minVal = inputSpec.stream().mapToDouble(booleanDoubleMap -> booleanDoubleMap.get(false)).min().getAsDouble();
                artificialInputBufferIDs.addAll(BranchAndBoundUtil.CELLO_INPUT_SPECIFICATION.keySet());
            }
            /*
            Only determine widest output values if there is an artificial input present.
            */
            if (artificialInputBufferIDs.size() > 0) {
            /*
            Determining which groups of repressors are used within the assignment.
            */
                Set<String> usedGroups = assignment.values().stream().map(GateRealization::getGroup).collect(Collectors.toSet());


            /*
            Set up the prerequisites for the maximization of input interval.
             */
                double yMax;
                double yMin;

                List<Double> yMaxParticles = new ArrayList<>();
                List<Double> yMinParticles = new ArrayList<>();


                /*
                Determine all realizations which belong to groups not already used.
                 */
                List<GateRealization.GateCharacterization> availableRealizations = realizations.values().stream().flatMap(Collection::stream)
                        .filter(realization -> !usedGroups.contains(realization.getGroup()) && realization.isCharacterized())
                        .map(GateRealization::getCharacterization)
                        .collect(Collectors.toList());

                /*
                Determine largest ymax and smallest ymin
                 */
                yMax = availableRealizations.stream().mapToDouble(GateRealization.GateCharacterization::getYmax).max().orElse(Math.pow(10, 15));
                yMin = availableRealizations.stream().mapToDouble(GateRealization.GateCharacterization::getYmin).min().orElse(0);
                maxVal = Math.max(maxVal, yMax);
                minVal = Math.min(minVal, yMin);

                // TODO: remove quick fix
                if (true) {
                    for (int i = 0; i < 5000; i++) {
                        int finalI = i;
                        yMax = availableRealizations.stream().mapToDouble(r -> r.getParticles().getYmax(finalI)).max().orElse(Math.pow(10, 15));
                        yMin = availableRealizations.stream().mapToDouble(r -> r.getParticles().getYmin(finalI)).min().orElse(0);

                        yMaxParticles.add(i, yMax);
                        yMinParticles.add(i, yMin);
                    }
                }
            }
        }

        /*
        substitutionMode.equals("f") -> No substitution is performed -> calculating with inexact bounds
        substitutionMode.equals("t") -> Substitution is performed -> Leads to guaranteed optimal result but also to less tightness
         */
        additionalArgs = String.format(" substitute=%s use_custom_input_specification=1 cis=%s", substitutionMode, BranchAndBoundUtil.createCustomInputSpecification(artificialInputBufferIDs, minVal, maxVal));
        return additionalArgs;
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
            case BEST_FIRST_SEARCH:
                strategy = new BestFirstSearch();
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

        for (LogicGate gate : reversedLogicGates) {
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
     * Shutsdown the used interfaces
     */
    private void shutDownInterfaces() {
        for (int iX = 0; iX < interfaces.length; iX++) {
            if (interfaces[iX] != null)
                interfaces[iX].shutdown();
        }
    }
}