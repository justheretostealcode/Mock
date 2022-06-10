package de.tu_darmstadt.rs.synbio;

import com.fasterxml.jackson.databind.ObjectMapper;
import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.mapping.MappingConfiguration;
import de.tu_darmstadt.rs.synbio.mapping.search.AssignmentSearchAlgorithm;
import de.tu_darmstadt.rs.synbio.simulation.SimulationConfiguration;
import de.tu_darmstadt.rs.synbio.simulation.SimulationResult;
import de.tu_darmstadt.rs.synbio.synthesis.SynthesisConfiguration;
import de.tu_darmstadt.rs.synbio.common.TruthTable;
import de.tu_darmstadt.rs.synbio.synthesis.enumeration.EnumeratorFast;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class ARCTICsyn {

    private static final Logger logger = LoggerFactory.getLogger(ARCTICsyn.class);

    private final TruthTable inputTruthTable;
    private final File outputDir;

    private final SynthesisConfiguration synConfig;
    private final MappingConfiguration mapConfig;
    private final SimulationConfiguration simConfig;

    public ARCTICsyn(TruthTable inputTruthTable, SynthesisConfiguration synConfig,
                     MappingConfiguration mapConfig, SimulationConfiguration simConfig) throws IOException {

        this.inputTruthTable = inputTruthTable;
        this.synConfig = synConfig;
        this.mapConfig = mapConfig;
        this.simConfig = simConfig;

        /* output dir handling */
        outputDir = (synConfig.getOutputDir() == null) ? new File("run_" + System.currentTimeMillis() + "_" + inputTruthTable.toString()) :
                new File(synConfig.getOutputDir(), "run_" + System.currentTimeMillis() + "_" + inputTruthTable.toString());

        if (!this.outputDir.mkdirs())
            throw new IOException("Error creating output directory " + outputDir.getAbsolutePath());
    }

    public void synthesize() {

        /* initialize gate library */

        //final GateLibrary gateLib = new GateLibrary(mapConfig.getLibrary(), new Double[]{0.9,0.1,0.0});
        final GateLibrary gateLib = new GateLibrary(mapConfig.getLibrary(), mapConfig.getCompatibilityLibrary(),true);
        logger.info("Loaded gate library " + gateLib.getSourceFile() + ".");

        if (inputTruthTable.getSupportSize() > gateLib.getNumAvailableGates(LogicType.INPUT)) {
            logger.error("Number of available input gates (" + gateLib.getNumAvailableGates(LogicType.INPUT) +
                    ") < number of input function support (" +inputTruthTable.getSupportSize() + ")." );
            return;
        }

        /* circuit enumeration */

        logger.info("Enumeration of circuit variants...");
        EnumeratorFast enumerator = new EnumeratorFast(gateLib, inputTruthTable, synConfig);

        enumerator.enumerate();
        List<Circuit> circuits = new ArrayList<>(enumerator.getResultCircuits().get(inputTruthTable));

        if (circuits.isEmpty()) {
            logger.info("No circuit structures found.");
            return;
        }

        logger.info("Found " + circuits.size() + " circuits.");

        circuits.sort(Circuit::compareTo);

        if (circuits.size() > synConfig.getLimitStructuresNum()) {
            circuits = circuits.subList(circuits.size() - synConfig.getLimitStructuresNum(), circuits.size());
        }

        for (int i = 0; i < circuits.size(); i ++) {
            Circuit circ = circuits.get(i);
            circ.setIdentifier("structure_" + i);
            circ.print(new File(outputDir, inputTruthTable.toString() + "_" + circ.getIdentifier() + ".dot"));
            circ.save(new File(outputDir, inputTruthTable.toString() + "_" + circ.getIdentifier() + ".json"));
        }

        /* technology mapping */

        if (simConfig.isSimEnabled()) {

            logger.info("Technology mapping...");

            SimulationResult[] bestResults = new SimulationResult[synConfig.getWeightRelaxation() + 1];
            int minSize = circuits.get(circuits.size() - 1).getWeight();

            for (int i = circuits.size() - 1; i >= 0; i--) {

                AssignmentSearchAlgorithm sim = mapConfig.getSearchAlgorithm(circuits.get(i), gateLib, simConfig);

                try {
                    SimulationResult result = sim.assign();
                    int relSize = circuits.get(i).getWeight() - minSize;

                    if (bestResults[relSize] == null || result.getScore() > bestResults[relSize].getScore())
                        bestResults[relSize] = result;

                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

            double currentBestScore = 0.0;

            for(int i = 0; i < bestResults.length; i++) {
                SimulationResult result = bestResults[i];

                if (result != null && result.getScore() > currentBestScore) {
                    result.getStructure().save(new File(outputDir, "result_" + result.getStructure().getTruthTable() + "_" + i + "_relax.json"));
                    result.getStructure().print(new File(outputDir, "result_" + result.getStructure().getTruthTable() + "_" + i + "_relax.dot"));

                    try {
                        ObjectMapper mapper = new ObjectMapper();
                        mapper.writerWithDefaultPrettyPrinter().writeValue(new File(outputDir, "result_" + result.getStructure().getTruthTable() + "_" + i + "_relax_assignment.json"), result.getAssignment());
                    } catch (Exception e) {
                        e.printStackTrace();
                    }

                    result.getStructure().saveGml(new File(outputDir.getParent(), "result_" + result.getStructure().getTruthTable() + "_" + i + "_relax.gml"), result.getAssignment());

                    logger.info(result.getStructure().getTruthTable() + "," + i + "," + result.getStructure().getNumberLogicGates() + "," + result.getScore());

                    currentBestScore = result.getScore();

                    logger.info("Finished. Result:");
                    logger.info(result.getStructure().getIdentifier() + "," + result.getScore() + "," + result.getStructure().getWeight() + "," + result.getAssignment().getIdentifierMap().toString());
                }
            }
        }
        System.exit(0);
    }
}
