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

        Option function = new Option("f", "function", true, "input function");
        options.addOption(function);

        Option truthTable = new Option("t", "truthtable", true, "input truth table");
        options.addOption(truthTable);

        /* library */

        Option gateLibraryFile = new Option("l", "library", true, "path of the gate library file");
        options.addOption(gateLibraryFile);


        Option mappingConfigString = new Option("mc", "mappingConfig", true, "path to the mapping configuration");
        options.addOption(mappingConfigString);
        Option simulationConfigString = new Option("sc", "simulationConfig", true, "path to the simulation configuration");
        options.addOption(simulationConfigString);
        Option synthesisConfigString = new Option("synconf", "synthesisConfig", true, "path to the synthesis configuration");
        options.addOption(synthesisConfigString);

        CommandLine cmd;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            exit(e.getMessage());
            return;
        }

        // sanity check arguments

        if (!(cmd.hasOption("function")  || cmd.hasOption("truthtable")) || !cmd.hasOption("library")) {
            exit("Input function or gate library file not given!");
        }

        if (cmd.hasOption("function")  && cmd.hasOption("truthtable")) {
            exit("Input function and truth table given!");
        }

        // Override mapping config file if provided
        if (cmd.hasOption("mappingConfig")) {
            mappingConfigFile = cmd.getOptionValue("mappingConfig");
        }

        // Override simulation config file if provided
        if (cmd.hasOption("simulationConfig")) {
            simulationConfigFile = cmd.getOptionValue("simulationConfig");
        }

        // Override synthesis config file if provided
        // Override simulation config file if provided
        if (cmd.hasOption("synthesisConfig")) {
            synthesisConfigFile = cmd.getOptionValue("synthesisConfig");
        }

        SynthesisConfiguration synConfig = new SynthesisConfiguration(synthesisConfigFile);
        //synConfig.print();
        SimulationConfiguration simConfig = new SimulationConfiguration(simulationConfigFile);
        //simConfig.print();
        MappingConfiguration mapConfig = new MappingConfiguration(mappingConfigFile);
        //mapConfig.print();

        /* input handling */

        TruthTable inputTruthTable;

        if (cmd.hasOption("function")) {
            inputTruthTable = new TruthTable(ExpressionParser.parse(cmd.getOptionValue("function")));
        } else {
            int ttLength = cmd.getOptionValue("truthtable").length();
            if ((ttLength & (ttLength - 1)) != 0) {
                exit("Length of truth table has to be power of two.");
            }
            inputTruthTable = new TruthTable(cmd.getOptionValue("truthtable"));
        }

        // call main program

        ARCTICsyn syn = new ARCTICsyn(inputTruthTable, cmd.getOptionValue("library"), synConfig, mapConfig, simConfig);
        syn.synthesize();
    }

    private static void exit(String message) {
        logger.error("Error: " + message);
        formatter.printHelp("ARCTICsyn", options);
        System.exit(1);
    }
}
