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
        // TODO Remove fixed circuit later


        // TODO Beachten, dass die Bounding Function für Minimierung angepasst werden müsste

        // TODO Möglichkeit geben eine Initiale Lösung zu übernehmen

        // TODO Möglichkeit implementieren die lower bound iterativ zu decrementieren (am besten als Mögliche Optimierungsart)

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
            //System.out.println(g.getIdentifier());
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
                    Collections.reverse(realizations);    // Reversion is necessary for depth first search, since the gate library is sorted descending based on the realizations possible scores (better first)

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
        // TODO: Check whether initial Assignment exists. If not create one or start blank.


        // TODO: Generate Subproblems (inclusive Whitelists etc.) (+)
        //       (+) Erzeugen der Simulatoren für jedes Subproblem
        //       (+)   Ordnung über Gatter aufstellen
        //       (+)   Whitelists Berechnen
        //       (+) Strukturen erzeugen
        //       (?) (Heuristik für gute initiale Lösung)

        // TODO Optimierung Implementieren
        //      (+)Bounding Function implementieren
        //      (+) Selection Strategy implementieren (Allgemein, sodass die vier genannten Strategien möglich sind)
        //           Überprüfen wie das modifizierte Depth First Search eingesetzt werden kann, welches den besten Kindsknoten präveriert
        //           Kann durch Selection Strategy anhand der Ebene des Subproblems erkannt werden und dem Score (oder Beihilfe durch die Branching Strategy)
        //      (+) Branching Function implementieren
        //           Nur valide Assignments erzeugen -> Auf Gruppen Zugehörigkeit prüfen
        //           Der Wert der Bounding Function wird nicht in der Branching Function evaluiert
        //           Hier werden wirklich nur alle Kindsknoten erzeugt
        //           List aller verfügbaren Gatter nehmen und diese auf die Group Constraint filtern
        //      Selection Strategy Switch implementieren

        // TODO Logging implementieren
        //      (+) Anzahl der benötigten Simulationen
        //      Score in Abhängigkeit der Iteration
        //      Wie viele Knoten wurden besucht
        //      Wie viele Knoten wurden auf jeder Ebene besucht
        //      Auf welcher Ebene wurden die meisten Knoten gepruned
        //      Durchschnittlicher Score auf einer Ebene (Auch Median und Varianz betrachten)

        SearchStatsLogger searchStatsLogger = new SearchStatsLogger(structure, mapConfig, simConfig, reversedLogicGates.length);

        Comparator<QueueItem> comparator;
        if (mapConfig.getBabSearchStrategy() == MappingConfiguration.BAB_SearchStrategy.DEPTH_FIRST_SEARCH)
            comparator = new Comparator<QueueItem>() {                   // Sort ascending to improve Depth First Search
                @Override
                public int compare(QueueItem o1, QueueItem o2) {
                    return (int) Math.signum(o1.val - o2.val);
                }
            };
        else
            comparator = new Comparator<QueueItem>() {                   // Sort descending to improve Breath First Search
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
                // Therefore no break here since the following comparator achieves
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


        // TODO Optionally determine a value based on the library
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

            if (currentItem.val <= bestScore) { // This check is required, if the best score has changed after this queue item has been added to the queue
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

//                    QueueItem finalCurrentItem = currentItem;

                    double finalBestScore = bestScore;
                    List<QueueItem> queueItems = childs.stream()
                            .map(assignment -> {                                    // Map to QueueItem

//                                if (assignment.toString().equals("{\"NOR2_4 : NOR2_4\", \"NOR2_7 : NOR2_1\", \"NOT_0 : NOT_17\", \"NOR2_3 : NOR2_13\", \"NOR2_2 : NOR2_10\"}")) {
//                                    logger.info("Assignment of interest");
//                                }
                                double dVal = bound(assignment);

//                                if (dVal > finalCurrentItem.val) {
//                                    logger.info("Error regarding optimality of bounding function.");
//                                }

                                ref.highestScore = Math.max(ref.highestScore, dVal);
                                return QueueItem.getQueueItem(assignment, dVal);
                            })
                            .filter(item -> bVisualize || item.val > finalBestScore)              // Filter the assignments which do not have sufficient score (Filtering is skipped if visualisation is turned on, in order to result obtain complete Search Trees)
                            .sorted(comparator)// Ensures, that the search strategy visits the best node first
                            .collect(Collectors.toList());
                    //.forEach(item -> strategy.addToQueue(item));            // Add to queue
                    strategy.addToQueue(queueItems);
                } else {   // The child assignments are leaves

                    for (Assignment assignment : childs) {
//                        if (currentItem.val < 200) {
//                            logger.info("Assignment of interest.");     // TODO Remove!
//                        }

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

//                ToDo Remove
                if (boundingFunctionScore < 1)
                    logger.error(String.format("ERROR %s: Boundingfunction score is below 1. Score: %f, Parent Item Score: %f, Child Item Score: %f.", structure.getIdentifier(), boundingFunctionScore, currentItem.val, ref.highestScore));


            } else {
                /*
                Lazy Branch and Bound
                Verifizieren ob numberofOverallSimulations identisch ist.
                Optional die Stats wie BoundingFunctionScores hier hinzufügen.
                */
                val = (currentAssignment.size() != 0) ? bound(currentAssignment) : Double.POSITIVE_INFINITY;
                if (val <= bestScore)       // This check is required, if the best score has changed after this queue item has been added to the queue
                    continue;

                if (currentAssignment.size() < logicGates.length) {
                    // The assignment is an intermediate node and needs to be branched
                    childs = branch(currentAssignment);


//                    if (mapConfig.getBabSearchStrategy() == MappingConfiguration.BAB_SearchStrategy.DEPTH_FIRST_SEARCH)
//                        Collections.reverse(childs);    // Reversion is necessary for depth first search, since the gate library is sorted descending based on the realizations possible scores (better first)
                    /*
                    The gate library is sorted ascending for depth first search and otherwise descending, based on the possible maximum scores of each gate.

                    Since depth first search adds the last to the front, the best gate will be used first, when it is added last.
                    For Breadth First Search the order is appropriate. However, Breadth First Search itself is not appropriate either.
                     */
//                    // TODO Adding elements to queue twice is a huge mistake!
//                    for (Assignment assignment : childs) {
//                        strategy.addToQueue(QueueItem.getQueueItem(assignment, val));   // At child with parents score
//                    }
                    double finalVal = val;
                    List<QueueItem> queueItems = childs.stream()
                            .map(assignment -> QueueItem.getQueueItem(assignment, finalVal))
                            .collect(Collectors.toList());
                    //.forEach(item -> strategy.addToQueue(item));
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

//        // TODO bestScore needs to be evaluated for non custom input specification
//        //      Also therefore it would help, that the simulator is capable of changing the input specification for any simulation
//        double score = interfaces[bestAssignment.size() - 1].simulate(bestAssignment, " cis=0 substitute=f ");
//        if (score != bestScore)
//            throw new Error(String.format("Mismatch in scores. Check correctness of bounding. Current best score: %f vs. obtained score: %f", bestScore, score));
//        bestScore = score;  // TODO Evaluieren ob noch notwendig, da die Bound Funktion ja

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


//        String searchStatsEval = searchStatsLogger.evaluate(bestScore);
//        logger.info(searchStatsEval);

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
        //double score = interfaces[assignment.size() - 1].simulate(assignment, " custom_input_specification=1 custom_input_low=0.00001 custom_input_high=100 ");
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
                .map(gate -> gate.getIdentifier())
                .collect(Collectors.toSet());

        List<Gate> inputGates = subproblem.getInputBuffers();
        // The input buffer ids not present in the original problem
        Set<String> artificialInputBufferIDs = inputGates.stream()
                .map(gate -> gate.getIdentifier())
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
            Determining which groups of repressesors are used within the assignment.
            */
                Set<String> usedGroups = assignment.values().stream().map(realization -> realization.getGroup()).collect(Collectors.toSet());


            /*
            Set up the prerequisites for the maximization of input intervall.
             */
                double yMax;
                double yMin;
                for (LogicType type : LogicType.values()) {
                    if (realizations.get(type) == null)
                        continue;

                /*
                Determine all realizations which belong to groups not already used.
                 */
                    List<GateRealization.GateCharacterization> availableRealizations = realizations.get(type).stream()
                            .filter(realization -> !usedGroups.contains(realization.getGroup()) && realization.isCharacterized())
                            .map(realization -> realization.getCharacterization())
                            .collect(Collectors.toList());

                    if (availableRealizations.isEmpty())
                        continue;

                /*
                Determine largest ymax and smallest ymin
                 */
                    yMax = availableRealizations.stream().mapToDouble(gateCharacterization -> gateCharacterization.getYmax()).max().orElse(Math.pow(10, 15));
                    yMin = availableRealizations.stream().mapToDouble(gateCharacterization -> gateCharacterization.getYmin()).min().orElse(0);
                    maxVal = Math.max(maxVal, yMax);
                    minVal = Math.min(minVal, yMin);
                }
            }
        }

        /*
        substitutionMode.equals("f") -> No substitution is performed -> calculating with inexact bounds
        substitutionMode.equals("t") -> Substitution is performed -> Leads to guaranteed optimal result but also to less tightness
         */
        additionalArgs = String.format(" substitute=%s use_custom_input_specification=1 cis=%s", substitutionMode, BranchAndBoundUtil.createCustomInputSpecification(artificialInputBufferIDs, minVal, maxVal));
        //additionalArgs = additionalArgs.replace(",", ".");
        return additionalArgs;
    }

    /**
     * Returns the specified search strategy
     *
     * @return
     */
    private SearchStrategy getSearchStrategy() {
        // Get corresponding informaiton from mapConfig
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
        //LogicGate gate;
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