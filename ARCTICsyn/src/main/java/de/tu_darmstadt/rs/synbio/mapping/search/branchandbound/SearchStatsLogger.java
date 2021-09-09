package de.tu_darmstadt.rs.synbio.mapping.search.branchandbound;

import com.fasterxml.jackson.annotation.JsonAutoDetect;
import com.fasterxml.jackson.annotation.PropertyAccessor;
import com.fasterxml.jackson.databind.ObjectMapper;
import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.mapping.MappingConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.search.BranchAndBoundSearch;
import de.tu_darmstadt.rs.synbio.simulation.SimulationConfiguration;
import de.tu_darmstadt.rs.synbio.simulation.SimulationResult;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.IntSummaryStatistics;

public class SearchStatsLogger {

    private static final Logger logger = LoggerFactory.getLogger(BranchAndBoundSearch.class);

    private ArrayList<Double>[] boundingFunctionScores;
    private ArrayList<Double[]>[] boundingFunctionValues;
    private ArrayList<SimulationsPerScore> numberOfOverallSimulations;
    private long numberOfSimulationAtBestAssignmentDiscovery;
    private long iterationOfBestAssignmentDiscovery;
    private long numberOfSimulationsRequired;
    private long numberOfItemsAddedAndSkipped;
    private long optimalNumberOfSimulations;

    private double averageNumberOfQueueItems;
    private int maximumNumberOfQueueItems;

    private Assignment bestAssignment;
    private double bestScore;

    private long duration = -1;

    private Circuit structure;

    private MappingConfiguration mappingConfiguration;
    private SimulationConfiguration simulationConfiguration;

    public SearchStatsLogger(Circuit structure, MappingConfiguration mappingConfiguration, SimulationConfiguration simulationConfiguration) {
        constructor(structure, mappingConfiguration, simulationConfiguration, 1);
    }

    public SearchStatsLogger(Circuit structure, MappingConfiguration mappingConfiguration, SimulationConfiguration simulationConfiguration, int iNumberOfLevels) {
        constructor(structure, mappingConfiguration, simulationConfiguration, iNumberOfLevels);
    }

    private void constructor(Circuit structure, MappingConfiguration mappingConfiguration, SimulationConfiguration simulationConfiguration, int iNumberOfLevels) {
        this.boundingFunctionScores = new ArrayList[iNumberOfLevels - 1];
        this.boundingFunctionValues = new ArrayList[iNumberOfLevels - 1];
        for (int iX = 0; iX < boundingFunctionScores.length; iX++) {
            this.boundingFunctionScores[iX] = new ArrayList<>();
            this.boundingFunctionValues[iX] = new ArrayList<>();
        }

        this.numberOfOverallSimulations = new ArrayList<>();
        //iterationOfBestAssignmentDiscovery = -1;
        //bestScore = -1;
        this.structure = structure;
        this.mappingConfiguration = mappingConfiguration;
        this.simulationConfiguration = simulationConfiguration;
    }

    public void addBoundingFunctionScore(int assignmentSize, double parentVal, double childVal) {
        boundingFunctionScores[assignmentSize - 1].add(parentVal / childVal);
        boundingFunctionValues[assignmentSize - 1].add(new Double[]{parentVal, childVal});
    }

    public void addEntryToOverallSimulations(double score, int childCount) {
        numberOfOverallSimulations.add(new SimulationsPerScore(score, childCount));
    }

    public void notifyNewBestAssignment(long iteration, long neededSimulations) {
        iterationOfBestAssignmentDiscovery = iteration;
        this.numberOfSimulationAtBestAssignmentDiscovery = neededSimulations;
    }

    public void setNumberOfItemsAddedAndSkipped(long numberOfItemsAddedAndSkipped) {
        this.numberOfItemsAddedAndSkipped = numberOfItemsAddedAndSkipped;
    }

    public void setOptimalNumberOfSimulations(long optimalNumberOfSimulations) {
        this.optimalNumberOfSimulations = optimalNumberOfSimulations;
    }

    public void setResult(SimulationResult result) {
        this.bestScore = result.getScore();
        this.bestAssignment = result.getAssignment();
        this.numberOfSimulationsRequired = result.getNeededSimulations();
    }

    public void setDuration(long duration) {
        this.duration = duration;
    }

    public void setAverageNumberOfQueueItems(double averageNumberOfQueueItems) {
        this.averageNumberOfQueueItems = averageNumberOfQueueItems;
    }

    public void setMaximumNumberOfQueueItems(int maximumNumberOfQueueItems) {
        this.maximumNumberOfQueueItems = maximumNumberOfQueueItems;
    }

    private long getMinimalNumberOfSimulations() {
        IntSummaryStatistics simulationStats = numberOfOverallSimulations.stream()
                .filter(item -> item.score > bestScore)
                .mapToInt(item -> item.childCount)
                .summaryStatistics();
        long minimalNumberOfSimulations = simulationStats.getSum();
        return minimalNumberOfSimulations;
    }

    public String evaluate(double overallBestScore) {


        long minimalNumberOfSimulations = getMinimalNumberOfSimulations();

        DescriptiveStatistics boundingFunctionScoreStatistics = new DescriptiveStatistics();
        for (ArrayList<Double> list : boundingFunctionScores) {
            list.stream().forEach(boundingFunctionScoreStatistics::addValue);
        }

        DescriptiveStatistics[] boundingFunctionScoreLevelStatistics = new DescriptiveStatistics[boundingFunctionScores.length];
        for (int iX = 0; iX < boundingFunctionScores.length; iX++) {
            boundingFunctionScoreLevelStatistics[iX] = new DescriptiveStatistics();
            int finalIX = iX;
            boundingFunctionScores[iX].stream().forEach(item -> boundingFunctionScoreLevelStatistics[finalIX].addValue(item));
        }

        SearchStats stats = new SearchStats(minimalNumberOfSimulations, boundingFunctionScoreStatistics, boundingFunctionScoreLevelStatistics, iterationOfBestAssignmentDiscovery);

        return stats.toString();
    }

    /*
    TODO Tracking der Anzahl von NOT/NOR und OR
    Überlegen ob ein Score möglich ist der Beschreibt ob die NOR/OR eher vorne oder hinten liegen.
    Tracken was die maximale Anzahl an Elementen in der Queue ist.
    Tracken was die durchschnittliche Anzahl an Elementen in der Queue ist.

     */

    public void save(String path) {

        ObjectMapper mapper = new ObjectMapper();
        mapper.setVisibility(PropertyAccessor.FIELD, JsonAutoDetect.Visibility.ANY);

        HashMap<String, Object> jsonMap = new HashMap<>();
        File file = new File(path);
        file.getParentFile().mkdirs();
        try {
            jsonMap.put("TIMESTAMP", new Timestamp(System.currentTimeMillis()).toString());
            jsonMap.put("BEST_SCORE", bestScore);
            jsonMap.put("BEST_ASSIGNMENT", bestAssignment.getIdentifierMap());
            jsonMap.put("NUMBER_OF_SIMULATIONS_REQUIRED", numberOfSimulationsRequired);

            jsonMap.put("MAPPING_CONFIGURATION", mappingConfiguration);
            jsonMap.put("SIMULATOR_CONFIGURATION", simulationConfiguration);
            jsonMap.put("STRUCTURE_STATS", new HashMap<>() {{
                put("NUMBER_OF_INPUT_BUFFERS", structure.getInputBuffers().size());
                put("NUMBER_OF_LOGIC_GATES", structure.getLogicGates().size());
                put("NUMBER_OF_OUTPUT_BUFFERS", 1);
                put("OR_GATE_COUNT", structure.getLogicGates().stream().filter(logicGate -> logicGate.getLogicType() == LogicType.OR2).count());
                put("NOR_GATE_COUNT", structure.getLogicGates().stream().filter(logicGate -> logicGate.getLogicType() == LogicType.NOR2).count());
                put("NOT_GATE_COUNT", structure.getLogicGates().stream().filter(logicGate -> logicGate.getLogicType() == LogicType.NOT).count());
                put("ESTIMATED_CIRCUIT_COMPLEXITY", structure.estimateCircuitComplexity());
            }});
            jsonMap.put("EXECUTION_TIME", duration);

            if (mappingConfiguration.getBabStatistics()) {
                jsonMap.put("OPTIMAL_NUMBER_OF_SIMULATIONS", optimalNumberOfSimulations); // An theoretic construct which depends on the quality of the bounding function. Could be achieved if the bounding function is equal to max(childScores)
                jsonMap.put("MINIMAL_NUMBER_OF_SIMULATIONS", getMinimalNumberOfSimulations());
                jsonMap.put("NUMBER_OF_ITEMS_ADDED_AND_SKIPPED", numberOfItemsAddedAndSkipped);
                jsonMap.put("BOUNDING_FUNCTION_SCORES", boundingFunctionScores);
                jsonMap.put("BOUNDING_FUNCTION_VALUES", boundingFunctionValues);

                jsonMap.put("NUMBER_OF_SIMULATIONS_PER_SCORE", numberOfOverallSimulations);
                jsonMap.put("NUMBER_OF_SIMULATION_AT_BEST_ASSIGNMENT_DISCOVERY", numberOfSimulationAtBestAssignmentDiscovery);
                jsonMap.put("NUMBER_OF_ITERATION_OF_BEST_ASSIGNMENT_DISCOVERY", iterationOfBestAssignmentDiscovery);
                jsonMap.put("AVERAGE_NUMBER_OF_QUEUE_ITEMS", averageNumberOfQueueItems);
                jsonMap.put("MAXIMUM_NUMBER_OF_QUEUE_ITEMS", maximumNumberOfQueueItems);
            }

            mapper.writerWithDefaultPrettyPrinter().writeValue(file, jsonMap);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public class SimulationsPerScore {
        double score;
        int childCount;

        public SimulationsPerScore(double score, int childCount) {
            this.score = score;
            this.childCount = childCount;
        }
    }


    public static class SearchStats {
        private long minimumNumberOfSimulations;
        private DescriptiveStatistics boundingFunctionScoresStatistics;
        private DescriptiveStatistics[] boundingFunctionScoresLevelStatistics;
        private long iterationOfBestAssignmentDiscovery;

        private SearchStats(long minimumNumberOfSimulations, DescriptiveStatistics boundingFunctionScoresStatistics, DescriptiveStatistics[] boundingFunctionScoresLevelStatistics, long iterationOfBestAssignmentDiscovery) {
            this.minimumNumberOfSimulations = minimumNumberOfSimulations;
            this.boundingFunctionScoresStatistics = boundingFunctionScoresStatistics;
            this.boundingFunctionScoresLevelStatistics = boundingFunctionScoresLevelStatistics;
            this.iterationOfBestAssignmentDiscovery = iterationOfBestAssignmentDiscovery;
        }

        public long getMinimumNumberOfSimulations() {
            return minimumNumberOfSimulations;
        }

        public DescriptiveStatistics getBoundingFunctionScoresStatistics() {
            return boundingFunctionScoresStatistics;
        }

        public long getIterationOfBestAssignmentDiscovery() {
            return iterationOfBestAssignmentDiscovery;
        }

        public String toString() {
            StringBuilder builder = new StringBuilder();
            builder.append("Search Statistics: \n");
            builder.append(String.format("Minimum Number of Simulations: %d\n", minimumNumberOfSimulations));
            builder.append(String.format("Iteration of best Assignment Discovery: %d\n", iterationOfBestAssignmentDiscovery));
            builder.append(descriptiveStatisticsToString(boundingFunctionScoresStatistics));

            for (int iX = 0; iX < boundingFunctionScoresLevelStatistics.length; iX++) {
                builder.append(String.format("Bounding Function Statistics for Level %d to %d\n", iX + 1, iX + 2));
                builder.append(descriptiveStatisticsToString(boundingFunctionScoresLevelStatistics[iX]));
            }

            return builder.toString();
        }


        private String descriptiveStatisticsToString(DescriptiveStatistics statistics) {
            StringBuilder builder = new StringBuilder();

            builder.append(String.format("Median Bounding Function Score: %f\n", statistics.getPercentile(50)));
            builder.append(String.format("Average Bounding Function Score: %f\n", statistics.getMean()));
            builder.append(String.format("Variance Bounding Function Scores: %f\n", statistics.getVariance()));
            builder.append(String.format("[Min, Max]: [%f, %f]\n", statistics.getMin(), statistics.getMax()));
            builder.append(String.format("Quantiles: Q001 = %f\n", statistics.getPercentile(1)));
            builder.append(String.format("           Q005 = %f\n", statistics.getPercentile(5)));
            builder.append(String.format("           Q025 = %f\n", statistics.getPercentile(25)));
            builder.append(String.format("           Q045 = %f\n", statistics.getPercentile(45)));
            builder.append(String.format("           Q055 = %f\n", statistics.getPercentile(55)));
            builder.append(String.format("           Q075 = %f\n", statistics.getPercentile(75)));
            builder.append(String.format("           Q095 = %f\n", statistics.getPercentile(95)));
            builder.append(String.format("           Q099 = %f\n", statistics.getPercentile(99)));

            return builder.toString();
        }
    }
}