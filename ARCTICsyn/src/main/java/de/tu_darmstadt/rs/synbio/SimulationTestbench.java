package de.tu_darmstadt.rs.synbio;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ObjectNode;
import de.tu_darmstadt.rs.synbio.common.circuit.LogicGate;
import de.tu_darmstadt.rs.synbio.common.library.GateRealization;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.mapping.MappingConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.compatibility.CompatibilityChecker;
import de.tu_darmstadt.rs.synbio.mapping.search.AssignmentSearchAlgorithm;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.circuit.CircuitDeserializer;
import de.tu_darmstadt.rs.synbio.simulation.SimulationConfiguration;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.mapping.MappingConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.search.AssignmentSearchAlgorithm;
import de.tu_darmstadt.rs.synbio.mapping.search.BranchAndBoundSearch;
import de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.SearchStatsLogger;
import de.tu_darmstadt.rs.synbio.simulation.SimulationConfiguration;
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

    private static int numRepetitions = 1;

    public static void main(String[] args) throws Exception {

        Options options = new Options();

        Option inputDirString = new Option("i", "inputPath", true, "path to the input directory or file");
        options.addOption(inputDirString);
        Option gateLibraryFile = new Option("l", "library", true, "path of the gate library file");
        options.addOption(gateLibraryFile);
        Option proxWeightsOpt = new Option("w", "proxWeights", true, "weights for the gate proximity measure");
        options.addOption(proxWeightsOpt);


        Option mappingConfigString = new Option("mc", "mappingConfig", true, "path to the mapping configuration");
        options.addOption(mappingConfigString);
        Option simulationConfigString = new Option("sc", "simulationConfig", true, "path to the simulation configuration");
        options.addOption(simulationConfigString);

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            logger.error(e.getMessage());
            formatter.printHelp("AssignmentBenchmark", options);
            System.exit(1);
            return;
        }

        if (!cmd.hasOption("inputPath") || !cmd.hasOption("library")) {
            logger.error("Input directory or gate library file not given!");
            formatter.printHelp("enuMap", options);
            System.exit(1);
            return;
        }

        Double[] proxWeights = {0.9, 0.1, 0.0};

        if (cmd.hasOption("proxWeights")) {
            String[] pw = cmd.getOptionValue("proxWeights").split(",");
            proxWeights = Arrays.stream(pw).map(Double::valueOf).toArray(Double[]::new);
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

        if (mapConfig.getNumRepetitions() > 0)
            // TODO Dirty coding since mapping config ideally should not change the simulation testbench
            numRepetitions = mapConfig.getNumRepetitions();

        GateLibrary gateLib = new GateLibrary(new File(cmd.getOptionValue("library")), proxWeights);

        CompatibilityChecker checker = new CompatibilityChecker(gateLib);

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

                    neededSims += result.getNeededSimulations();

                    logger.info(child.getName() + "," + result.getScore() + "," + result.getStructure().getWeight() + "," + result.getNeededSimulations() + "," + result.getAssignment().getIdentifierMap());

                    //result.getStructure().print(new File(inputPath.isDirectory() ? inputPath : inputPath.getParentFile(),
                    //        "result_" + child.getName() + ".dot"), result.getAssignment());

                    //result.getStructure().saveGml(new File(inputPath.isDirectory() ? inputPath : inputPath.getParentFile(),
                    //        "result_" + child.getName() + ".gml"), result.getAssignment());

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
