package de.tu_darmstadt.rs.synbio.mapping.search;

import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.circuit.Gate;
import de.tu_darmstadt.rs.synbio.common.circuit.Wire;
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
import org.jgrapht.Graphs;
import org.jgrapht.traverse.TopologicalOrderIterator;
import org.logicng.formulas.FormulaFactory;
import org.logicng.util.Pair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.stream.Collectors;

public class BranchAndBoundSearch extends AssignmentSearchAlgorithm {

    private static final Logger logger = LoggerFactory.getLogger(BranchAndBoundSearch.class);

    private final Map<LogicType, List<GateRealization>> realizations;

    private final Circuit[] subProblems;
    private final Map<String, LogicType> replacedLogicTypes;
    //private final SimulatorInterface[] interfaces;
    private final SimulatorInterface simulator;

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

    private final boolean bFastMode;

    private final MappingConfiguration.BAB_INPUT_SPECIFICATION_TYPE babInputSpecificationType;


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
        Pair<Circuit[], Map<String, LogicType>> subProblemResult = getSubProblems(this.structure);
        subProblems = subProblemResult.first();
        replacedLogicTypes = subProblemResult.second();

        /*structure.print(new File("structure.dot"));
        structure.save(new File("structure.json"));
        int i = 0;
        for(Circuit sp : subProblems) {
            sp.print(new File("sb_" + i + ".dot"));
            sp.save(new File("sb_" + i + ".json"));
            i++;
        }*/

        //interfaces = getSubProblemInterfaces(this.structure, this.subProblems);
        simulator = new SimulatorInterface(simConfig, gateLib);
        simulator.initSimulation(structure);

        logicGates = getLogicGatesInTopologicalOrder(structure);
        reversedLogicGates = getReversedLogicGates(logicGates);

        outputGate = structure.getOutputGate();
        outputRealization = gateLib.getOutputDevice(outputGate.getLogicType());

        initialAssignment = null;
        initialBestScore = Double.NEGATIVE_INFINITY;


        if (mapConfig.getBabType() != MappingConfiguration.BAB_Type.EAGER && mapConfig.getBabType() != MappingConfiguration.BAB_Type.LAZY)
            throw new Error(mapConfig.getBabType().name() + " is unknown to BranchAndBound");

        bEagerBranchAndBound = mapConfig.getBabType() == MappingConfiguration.BAB_Type.EAGER;


        boolean bVisualize = mapConfig.getBabVisualization();

        searchTreeVisualizer = new SearchTreeVisualizer(structure.getIdentifier(), mapConfig, reversedLogicGates, bVisualize);
        this.bVisualize = searchTreeVisualizer.getbVisualize();

        /*
        If bFastMode == true -> No substitution is performed, leading to a bounding function which is not strictly greater equal
        If bFastMode == false -> Substitution is performed leading to a mathematical correct description but also decreasing the performance
         */
        bFastMode = mapConfig.getBabFast();

        babInputSpecificationType = mapConfig.getBabInputSpecificationType();
    }

    /**
     * This method creates the subproblems belonging to the provided circuit.   <br>
     * Each subproblem thereby contains one gate less than the previous subproblem.
     *
     * @param structure The circuit to from which the subproblems shall be created
     * @return A pair of: An array of subproblem circuits and a map of IDs of gate replacements and their replaced logic type
     */
    private Pair<Circuit[], Map<String, LogicType>> getSubProblems(Circuit structure) {

        Map<String, LogicType> replacedLogicTypes = new HashMap<>();

        // Establish the reverse topological order of the logic gates
        ArrayList<Gate> gatesInReversedOrder = new ArrayList<>();
        Iterator<Gate> iterator = new TopologicalOrderIterator<>(structure);
        while (iterator.hasNext()) {
            Gate g = iterator.next();

            if (g.isLogicGate() || g.getLogicType() == LogicType.INPUT)
                gatesInReversedOrder.add(g);
        }

        Gate[] gateOrder = new Gate[gatesInReversedOrder.size()];
        gatesInReversedOrder.toArray(gateOrder);

        List<Gate> originalInputBuffers = structure.getInputBuffers();
        List<Map<Gate, List<Gate>>> substitutionsList;
        Map<Gate, String> substitutionTruthTables;
        // Create sub problems
        ArrayList<Circuit> subProblems = new ArrayList<>(gateOrder.length - 1);
        Circuit subProblem = new Circuit("Subproblem");
        Graphs.addGraph(subProblem, structure);
        FormulaFactory factory = new FormulaFactory();
        String structureIdentifier = structure.getIdentifier();
        String subProblemIdentifier = structureIdentifier + "_subproblem_";

        subProblems.add(subProblem.copy(subProblemIdentifier + 0));

        Map<Gate, Gate> insertedInputBuffers = new HashMap<>();
        for (int iX = 0; iX < gateOrder.length - 1; iX++) {
            Gate g = gateOrder[iX];

            // Remove Vertex
            Set<Wire> wires = subProblem.outgoingEdgesOf(g);

            // For each gate, which is removed, we need to add a new input to the circuit and connect this input to the following gates
            // Next to this, we mark, that the node's truthtable values result from a gate which was previously in the circuit.
            String inputIdentifier = "OUT_" + g.getIdentifier();    // Named "OUT", since it represents the output of the referenced gate
            Gate newInputBuffer = new Gate(inputIdentifier, LogicType.INPUT);
            replacedLogicTypes.putIfAbsent(inputIdentifier, g.getLogicType());
            subProblem.addVertex(newInputBuffer);
            insertedInputBuffers.put(newInputBuffer, g);
            for (Wire w : wires) {
                Gate target = subProblem.getEdgeTarget(w);
                subProblem.addEdge(newInputBuffer, target, new Wire(factory.variable(w.getVariable().name())));
            }

            subProblem.removeVertex(g);
            // Remove all Gates which do not contribute to the result
            BranchAndBoundUtil.cleanCircuitFromNonContributingGates(subProblem);
            subProblem.removeRedundantGates();

            // Get the Whitelist for the resulting structure
            String whiteList = BranchAndBoundUtil.determineWhitelist(structure, subProblem, insertedInputBuffers);
            subProblem.setWhitelist(whiteList);

            // Obtain substitutions list for the subproblem. //TODO remove
            substitutionsList = BranchAndBoundUtil.determineSubstitutionList(subProblem, originalInputBuffers);
            subProblem.setSubstitutionsList(substitutionsList);

            // Substitution Truthtables are not used anymore //TODO remove
            substitutionTruthTables = BranchAndBoundUtil.determineSubstitutionTruthTables(subProblem, substitutionsList, originalInputBuffers);
            subProblem.setSubstitutionTruthTables(substitutionTruthTables);

            subProblems.add(subProblem.copy(subProblemIdentifier + (iX + 1)));
        }

        Collections.reverse(subProblems);   // Reverse the order to achieve, that the index corresponds to the distance from the source in the search tree

        Circuit[] subProblemsA = new Circuit[subProblems.size()];
        subProblems.toArray(subProblemsA);

        return new Pair<>(subProblemsA, replacedLogicTypes);
    }

    /**
     * Creates the sub problems and the corresponding simulator interfaces and returns the interfaces. <br>
     * The first element in the returned array belongs to the smallest sub problem. <br>
     * The last element in the returned array is the interface for the provided attribute @structure.
     *
     * @param structure The structure to crate the sub problems interfaces for
     * @return Simulator Interfaces to all sub problem of the structure and the structure itself
     */
    private SimulatorInterface[] getSubProblemInterfaces(Circuit structure, Circuit[] subProblems) {
        ArrayList<SimulatorInterface> interfaces = new ArrayList<>(structure.edgeSet().size() - 2);

        SimulatorInterface simulator;
        for (Circuit subProblem : subProblems) {
            simulator = new SimulatorInterface(simConfig, gateLib);
            simulator.initSimulation(subProblem);
            interfaces.add(simulator);
        }

        SimulatorInterface[] array = new SimulatorInterface[interfaces.size()];
        interfaces.toArray(array);
        return array;
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
    public SimulationResult assign() {
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

        double errorThreshold = 0.0;

        QueueItem currentItem;
        Assignment currentAssignment;
        //System.out.print("Iteration: " + 0);
        long iteration = 0;
        while ((currentItem = strategy.getNext()) != null) {    // Get the next element from queue and repeat until the queue is empty

            searchTreeVisualizer.add(currentItem, bestScore, false);

            if (iteration % 100 == 0)  // Only update every hundred iterations
                System.out.print("\rIteration: " + iteration + " (" + iNeededSimulations + ")");

            if (currentItem.val <= bestScore - errorThreshold) {
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
                int childSize = childs.get(0).size();

                var ref = new Object() {
                    double highestScore = Double.NEGATIVE_INFINITY;
                };

                if (childSize < logicGates.length - 1) {    // The child assignments are intermediate nodes (no leaves)

                    double finalBestScore = bestScore;
                    List<QueueItem> queueItems = childs.stream()
                            .map(assignment -> {                                    // Map to QueueItem

                                double dVal = bound(assignment);

                                ref.highestScore = Math.max(ref.highestScore, dVal);
                                return QueueItem.getQueueItem(assignment, dVal);
                            })
                            .filter(item -> bVisualize || item.val > finalBestScore - errorThreshold)              // Filter the assignments which do not have sufficient score (Filtering is skipped if visualisation is turned on, in order to result obtain complete Search Trees)
                            .sorted(comparator)// Ensures, that the search strategy visits the best node first
                            .collect(Collectors.toList());
                    strategy.addToQueue(queueItems);
                } else {   // The child assignments are leaves

                    for (Assignment assignment : childs) {
                        val = bound(assignment);
                        ref.highestScore = Math.max(ref.highestScore, val);
                        searchTreeVisualizer.addLeafNode(assignment, val, bestScore, true);

                        double growth = 1.0;//interfaces[assignment.size() - 1].getLastGrowth();
                        if (val > bestScore && growth >= 0.75) {   // Since the node is a terminal node (leave), one needs to check if it is better than the current best solution
                            bestScore = val;
                            bestAssignment = assignment;
                            searchStatsLogger.notifyNewBestAssignment(iteration, iNeededSimulations);
                        }
                    }
                }

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
             mode = SimulatorInterface.PropagationMode.EXACT;
        else
            mode = mapConfig.getBabFast() ? SimulatorInterface.PropagationMode.HEURISTIC : SimulatorInterface.PropagationMode.BOUNDING;

        return simulator.simulate(assignment, mode);
    }

    /**
     * This method generates a custom input specification to enable the bounding function to be as tightly as currently possible.
     *
     * @param assignment The assignment for which the custom input specification needs to be determined
     * @return A string consisting of the custom input specification
     */
    private String getAdditionalSimArgs(Assignment assignment) {
        String additionalArgs;

        int subProblemPointer = assignment.size() - 1;
        int problemPointer = subProblems.length - 1;
        if (subProblemPointer == problemPointer)
            return " --cis=0 ";

        Map<String, Pair<Double, Double>> inputIntervals = new HashMap<>();

        /*
        Determine which input gates of the original circuit are present in the actual circuit.
         */
        Circuit subProblem = subProblems[subProblemPointer]; // assignment.size() sollte immer größer 0 sein
        Circuit problem = subProblems[problemPointer];

        // The original input buffer ids
        Set<String> problemInputBufferIDs = problem.getInputBuffers()
                .stream()
                .map(Gate::getIdentifier)
                .collect(Collectors.toSet());

        List<Gate> inputGates = subProblem.getInputBuffers();
        // The input buffer ids not present in the original problem
        Set<String> artificialInputBufferIDs = inputGates.stream()
                .map(Gate::getIdentifier)
                .filter(ident -> !problemInputBufferIDs.contains(ident))
                .collect(Collectors.toSet());


        double maxVal = Double.NEGATIVE_INFINITY;
        double minVal = Double.POSITIVE_INFINITY;
        // Case differentiation for Input Specification Type //TODO: Unused --> tidy up
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

                for (String inputId : artificialInputBufferIDs) {

                    LogicType replacedLogicType = replacedLogicTypes.get(inputId);
                /*
                Set up the prerequisites for the maximization of input interval.
                 */
                    double yMax;
                    double yMin;

                /*
                Determine all realizations which belong to groups not already used.
                 */
                    List<GateRealization.GateCharacterization> availableRealizations = realizations.get(replacedLogicType).stream()
                            .filter(realization -> !usedGroups.contains(realization.getGroup()) && realization.isCharacterized())
                            .map(GateRealization::getCharacterization)
                            .collect(Collectors.toList());

                /*
                Determine largest ymax and smallest ymin
                 */
                    yMax = availableRealizations.stream().mapToDouble(GateRealization.GateCharacterization::getYmax).max().orElse(Math.pow(10, 15));
                    yMin = availableRealizations.stream().mapToDouble(GateRealization.GateCharacterization::getYmin).min().orElse(0);
                    //maxVal = Math.max(maxVal, yMax);
                    //minVal = Math.min(minVal, yMin);

                    inputIntervals.put(inputId, new Pair<>(yMin, yMax));

                    // TODO: uncomment when implementing particle case
                /*List<Double> yMaxParticles = new ArrayList<>();
                List<Double> yMinParticles = new ArrayList<>();
                for (int i = 0; i < 5000; i++) {
                    int finalI = i;
                    yMax = availableRealizations.stream().mapToDouble(r -> r.getParticles().getYmax(finalI)).max().orElse(Math.pow(10, 15));
                    yMin = availableRealizations.stream().mapToDouble(r -> r.getParticles().getYmin(finalI)).min().orElse(0);

                    yMaxParticles.add(i, yMax);
                    yMinParticles.add(i, yMin);
                }*/
                }
            }
        }

        //TODO: implement additional args for particle case
        additionalArgs = String.format(" --substitute=%s --use_custom_input_specification=1 --cis=%s", bFastMode ? "0" : "1", BranchAndBoundUtil.createCustomInputSpecification(inputIntervals));
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
        /*for (SimulatorInterface anInterface : interfaces) {
            if (anInterface != null)
                anInterface.shutdown();
        }*/
    }
}