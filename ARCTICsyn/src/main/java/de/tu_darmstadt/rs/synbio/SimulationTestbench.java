package de.tu_darmstadt.rs.synbio;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ObjectNode;
import de.tu_darmstadt.rs.synbio.mapping.MappingConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.search.AssignmentSearchAlgorithm;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.circuit.CircuitDeserializer;
import de.tu_darmstadt.rs.synbio.simulation.SimulationConfiguration;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.SearchStatsLogger;
import de.tu_darmstadt.rs.synbio.simulation.SimulationResult;
import org.apache.commons.cli.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;

public class SimulationTestbench {

    private static final Logger logger = LoggerFactory.getLogger(SimulationTestbench.class);

    private static String mappingConfigFile = "map.config";
    private static String simulationConfigFile = "sim.config";
    private static Double[] proxWeights = {0.9, 0.1, 0.0};

    private static int numRepetitions = 1;

    public static void main(String[] args) throws Exception {

        Options options = new Options();

        /* parameters */
        options.addRequiredOption("i", "inputPath", true, "path to the input directory or file");

        /* optional parameters */
        options.addOption("w", "proxWeights", true, "weights for the gate proximity measure");
        options.addOption("mc", "mappingConfig", true, "path to the mapping configuration");
        options.addOption("sc", "simulationConfig", true, "path to the simulation configuration");
        options.addOption("n", "numRepetitions", true, "iteration number for the test bench");

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            logger.error(e.getMessage());
            formatter.printHelp("ARCTIC SimulationTestbench", options);
            System.exit(1);
            return;
        }

        if (!cmd.hasOption("inputPath")) {
            logger.error("Input directory or file not given!");
            formatter.printHelp("ARCTIC SimulationTestbench", options);
            System.exit(1);
            return;
        }

        if (cmd.hasOption("proxWeights")) {
            String[] pw = cmd.getOptionValue("proxWeights").split(",");
            proxWeights = Arrays.stream(pw).map(Double::valueOf).toArray(Double[]::new);
        }

        if (cmd.hasOption("numRepetitions")) {
            numRepetitions = Integer.parseInt(cmd.getOptionValue("numRepetitions"));
        }

        // Override mapping config file if provided
        if (cmd.hasOption("mappingConfig")) {
            mappingConfigFile = cmd.getOptionValue("mappingConfig");
        }

        // Override simulation config file if provided
        if (cmd.hasOption("simulationConfig")) {
            simulationConfigFile = cmd.getOptionValue("simulationConfig");
        }

        MappingConfiguration mapConfig = new MappingConfiguration(mappingConfigFile);
        SimulationConfiguration simConfig = new SimulationConfiguration(simulationConfigFile);

        //GateLibrary gateLib = new GateLibrary(mapConfig.getLibrary(), proxWeights);
        GateLibrary gateLib = new GateLibrary(mapConfig.getLibrary(), true);

        //CompatibilityChecker checker = new CompatibilityChecker(gateLib);

        File inputPath = new File(cmd.getOptionValue("inputPath"));

        File[] directoryListing;

        if (inputPath.isDirectory()) {

            directoryListing = inputPath.listFiles();

            if (directoryListing == null) {
                logger.info("Empty input directory.");
                return;
            }

            Arrays.sort(directoryListing);
        } else {
            directoryListing = new File[1];
            directoryListing[0] = inputPath;
        }

        File output = new File(inputPath.isDirectory() ? inputPath : inputPath.getParentFile(), "results_" + System.currentTimeMillis() + ".txt");
        PrintWriter out;
        try {
            out = new PrintWriter(output);
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }

        long startTime = System.currentTimeMillis();

        for (File child : directoryListing) {

            // test if json
            if (!child.getName().endsWith(".json"))
                continue;

            Circuit structure = null;

            final ObjectNode node;
            ObjectMapper mapper = new ObjectMapper();
            CircuitDeserializer circuitDeserializer = new CircuitDeserializer(Circuit.class);

            try {
                node = mapper.readValue(child, ObjectNode.class);

                if (node.has("graph")) {
                    structure = circuitDeserializer.deserializeString(node.get("graph").toString());
                    structure.setIdentifier(child.getName().replace(".json", ""));
                }

            } catch (Exception e) {
                e.printStackTrace();
            }

            if (structure != null) {

                /*Assignment ass = new Assignment();

                LogicGate k = (LogicGate) structure.vertexSet().stream().filter(g -> g instanceof LogicGate).findAny().get();
                GateRealization l = gateLib.getRealizations().get(k.getLogicType()).get(0);
                ass.put(k,l);
                checker.isCompatible(structure, ass);*/
                
                try {
                    out.print(child.getName());
                } catch (Exception e) {
                    e.printStackTrace();
                }

                int neededSims = 0;

                for (int i = 0; i < numRepetitions; i ++) {
                    AssignmentSearchAlgorithm search = mapConfig.getSearchAlgorithm(structure, gateLib, simConfig);
                    long start = System.nanoTime(); // Added for measuring the execution time.

                    SimulationResult result = search.assign();

                    long stop = System.nanoTime();
                    long duration = stop - start;

                    if (result != null) {

                        neededSims += result.getNeededSimulations();

                        if (result.getStructure() != null && result.getAssignment() != null) {
                            logger.info(child.getName() + "," + result.getScore() + "," + result.getStructure().getWeight() + "," + result.getNeededSimulations() + "," + duration + "," + result.getAssignment().getIdentifierMap());

                    /*AssignmentCounter counter = new AssignmentCounter(structure, gateLib, mapConfig, simConfig);
                    long maxAssignments = counter.assign().getNeededSimulations();
                    logger.info("Simulations: " + (double) result.getNeededSimulations()/maxAssignments*100 + "% (of " + maxAssignments + ")");*/

                            //mapper.writerWithDefaultPrettyPrinter().writeValue(new File("11101100_exhaustive_assignment.json"), result.getAssignment());

                            //result.getStructure().print(new File(inputPath.isDirectory() ? inputPath : inputPath.getParentFile(),
                            //        "result_" + child.getName() + ".dot"), result.getAssignment());

                            //result.getStructure().saveGml(new File(inputPath.isDirectory() ? inputPath : inputPath.getParentFile(),
                            //        "result_" + child.getName() + "_" + i + ".gml"), result.getAssignment());

                    /*try {
                        mapper.writerWithDefaultPrettyPrinter().writeValue(new File(inputPath.isDirectory() ? inputPath : inputPath.getParentFile(), "result_" + result.getStructure().getTruthTable() + "_" + i + "_assignment.json"), result.getAssignment());
                    } catch (Exception e) {
                        e.printStackTrace();
                    }*/

                            try {
                                out.print("," + result.getStructure().getNumberLogicGates() + "," + result.getScore() + "," + result.getAssignment().getIdentifierMap());
                                out.flush();
                            } catch (Exception e) {
                                e.printStackTrace();
                            }

                            if (mapConfig.getStatistics()) {
                                SearchStatsLogger logger = new SearchStatsLogger(structure, mapConfig, simConfig);
                                logger.setResult(result);
                                logger.setDuration(duration);
                                logger.save("testbench_statistics/statistics_" + structure.getIdentifier() + "_" + mapConfig.toString().hashCode() + "_" + i + "_" + System.currentTimeMillis() + ".json");
                            }
                        }
                    }
                }

                try {
                    out.print("," + neededSims + "\n");
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }

        out.println("total time: " + (System.currentTimeMillis() - startTime) + " ms for " + numRepetitions + " repetitions");
        out.close();

        System.exit(0);
    }
}
