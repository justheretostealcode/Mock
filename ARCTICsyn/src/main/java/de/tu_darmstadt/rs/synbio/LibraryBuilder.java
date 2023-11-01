package de.tu_darmstadt.rs.synbio;

import com.fasterxml.jackson.databind.ObjectMapper;
import de.tu_darmstadt.rs.synbio.common.TruthTable;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.mapping.MappingConfiguration;
import de.tu_darmstadt.rs.synbio.synthesis.SynthesisConfiguration;
import de.tu_darmstadt.rs.synbio.synthesis.enumeration.EnumeratorFast;
import org.apache.commons.cli.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.util.*;

public class LibraryBuilder {

    private static final Logger logger = LoggerFactory.getLogger(LibraryBuilder.class);

    public static <SyntheisConfiguration> void main(String[] args) throws Exception {

        Options options = new Options();

        options.addRequiredOption("o", "outputFile", true, "path of the output file");

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            logger.error(e.getMessage());
            formatter.printHelp("ARCTIC LibraryBuilder", options);
            System.exit(1);
            return;
        }

        if (!cmd.hasOption("outputFile")) {
            logger.error("Output file not given!");
            formatter.printHelp("ARCTIC LibraryBuilder", options);
            System.exit(1);
            return;
        }

        String mappingConfigFile = "map.config";
        MappingConfiguration mapConfig = new MappingConfiguration(mappingConfigFile);
        String synthesisConfigFile = "syn.config";
        SynthesisConfiguration synConfig = new SynthesisConfiguration(synthesisConfigFile);

        GateLibrary gateLib = new GateLibrary(mapConfig.getLibrary(), mapConfig.getCompatibilityLibrary(), GateLibrary.Type.ENERGY);

        EnumeratorFast enumerator = new EnumeratorFast(gateLib, 3, synConfig);
        enumerator.enumerate();
        Map<TruthTable, Set<Circuit>> circuits = enumerator.getResultCircuits();

        Map<TruthTable, Set<String>> output = new HashMap<>();
        circuits.keySet().forEach(t -> output.put(t, new HashSet<>()));

        for (TruthTable tt : circuits.keySet()) {
            output.put(tt, new HashSet<>());
            for (Circuit circuit : circuits.get(tt)) {
                output.get(tt).add(circuit.print());
            }
        }

        ObjectMapper mapper = new ObjectMapper();
        mapper.writerWithDefaultPrettyPrinter().writeValue(new File(cmd.getOptionValue("outputFile")), output);
    }
}
