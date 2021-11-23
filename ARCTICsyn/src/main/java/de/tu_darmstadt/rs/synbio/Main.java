package de.tu_darmstadt.rs.synbio;

import de.tu_darmstadt.rs.synbio.mapping.MappingConfiguration;
import de.tu_darmstadt.rs.synbio.simulation.SimulationConfiguration;
import de.tu_darmstadt.rs.synbio.synthesis.SynthesisConfiguration;
import de.tu_darmstadt.rs.synbio.synthesis.util.ExpressionParser;
import de.tu_darmstadt.rs.synbio.common.TruthTable;
import org.apache.commons.cli.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Main {

    private static final Logger logger = LoggerFactory.getLogger(Main.class);

    private static final Options options = new Options();
    private static final CommandLineParser parser = new DefaultParser();
    private static final HelpFormatter formatter = new HelpFormatter();

    private static String synthesisConfigFile = "syn.config";
    private static String mappingConfigFile = "map.config";
    private static String simulationConfigFile = "sim.config";

    public static void main(String[] args) throws Exception {

        // parse command line arguments

        /* input */
        options.addOption("f", "function", true, "input function");
        options.addOption("t", "truthtable", true, "input truth table");

        /* config files*/
        options.addOption("mc", "mappingConfig", true, "path to the mapping configuration");
        options.addOption("sc", "simulationConfig", true, "path to the simulation configuration");
        options.addOption("synconf", "synthesisConfig", true, "path to the synthesis configuration");

        CommandLine cmd;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            exit(e.getMessage());
            return;
        }

        // sanity check arguments

        if (!(cmd.hasOption("function")  || cmd.hasOption("truthtable"))) {
            exit("Input function or truth table not given!");
        }

        if (cmd.hasOption("function")  && cmd.hasOption("truthtable")) {
            exit("Both input function and truth table given!");
        }

        // Override config files if provided
        if (cmd.hasOption("mappingConfig")) {
            mappingConfigFile = cmd.getOptionValue("mappingConfig");
        }
        if (cmd.hasOption("simulationConfig")) {
            simulationConfigFile = cmd.getOptionValue("simulationConfig");
        }
        if (cmd.hasOption("synthesisConfig")) {
            synthesisConfigFile = cmd.getOptionValue("synthesisConfig");
        }

        SynthesisConfiguration synConfig = new SynthesisConfiguration(synthesisConfigFile);
        SimulationConfiguration simConfig = new SimulationConfiguration(simulationConfigFile);
        MappingConfiguration mapConfig = new MappingConfiguration(mappingConfigFile);

        /* input handling */

        TruthTable inputTruthTable;

        if (cmd.hasOption("function")) {
            inputTruthTable = new TruthTable(ExpressionParser.parse(cmd.getOptionValue("function")));
        } else {
            int ttLength = cmd.getOptionValue("truthtable").length();
            if ((ttLength & (ttLength - 1)) != 0) {
                exit("Length of truth table has to be a power of two.");
            }
            inputTruthTable = new TruthTable(cmd.getOptionValue("truthtable"));
        }

        // call main program
        ARCTICsyn syn = new ARCTICsyn(inputTruthTable, synConfig, mapConfig, simConfig);
        syn.synthesize();
    }

    private static void exit(String message) {
        logger.error("Error: " + message);
        formatter.printHelp("ARCTICsyn", options);
        System.exit(1);
    }
}
