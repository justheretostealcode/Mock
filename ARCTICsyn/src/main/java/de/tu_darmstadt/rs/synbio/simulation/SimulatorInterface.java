package de.tu_darmstadt.rs.synbio.simulation;

import com.fasterxml.jackson.databind.ObjectMapper;
import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.common.circuit.Gate;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.mapping.util.BranchAndBoundUtil;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class SimulatorInterface {

    private static final Logger logger = LoggerFactory.getLogger(SimulatorInterface.class);

    private final String pythonBinary;
    private final File simulatorPath;
    private final String simScript;
    private final String simInitArgs;
    private final String simArgs;
    private final GateLibrary library;

    private Process simProcess;
    private BufferedReader reader;
    private BufferedWriter writer;
    private BufferedReader errorReader;
    private final ObjectMapper mapper = new ObjectMapper();

    private Set<Gate> circuitGates;

    public SimulatorInterface(SimulationConfiguration config, GateLibrary gateLibrary) {
        pythonBinary = config.getPythonBinary();
        simulatorPath = config.getSimPath();
        simScript = config.getSimScript();
        simInitArgs = config.getSimInitArgs();
        simArgs = config.getSimArgs();
        library = gateLibrary;
    }

    public SimulatorInterface(String pythonBinary, File simPath, String simScript, String simInitArgs, String simArgs, GateLibrary gateLibrary) {
        this.pythonBinary = pythonBinary;
        simulatorPath = simPath;
        this.simScript = simScript;
        this.simInitArgs = simInitArgs;
        this.simArgs = simArgs;
        this.library = gateLibrary;
    }

    public boolean initSimulation(Circuit circuit) { //TODO: handle return value

        if (simProcess!= null && simProcess.isAlive())
            simProcess.destroy();

        //TODO: remove filter when output gates become part of the assignments
        circuitGates = circuit.vertexSet().stream().filter(g -> g.getLogicType() != LogicType.OUTPUT).collect(Collectors.toSet());

        /*try {
            String structureFileName = "structure_" + circuit.getIdentifier() + "_tid" + Thread.currentThread().getId() + "_" + System.nanoTime() + ".json";
            File structureFile = new File(simulatorPath, structureFileName);
            circuit.save(structureFile);

            ProcessBuilder pb = new ProcessBuilder(pythonBinary, simScript,
                    "--structure=" + structureFileName + " --lib_path=" + library.getSourceFile().getAbsolutePath() + " " + simInitArgs);
            pb.directory(simulatorPath);

            simProcess = pb.start();

            reader = new BufferedReader(new InputStreamReader(simProcess.getInputStream()));
            writer = new BufferedWriter(new OutputStreamWriter(simProcess.getOutputStream()));
            errorReader = new BufferedReader(new InputStreamReader(simProcess.getErrorStream()));

            String error = getError();

            String output;
            do {
                output = reader.readLine();

                if (!simProcess.isAlive()) {
                    logger.error("Simulator exited on initialization:\n" + error);
                    structureFile.delete();
                    return false;
                }

            } while (output == null || !output.startsWith("ready:"));

            structureFile.delete();

            return true;

        } catch (Exception e) {
            e.printStackTrace();
            return false;
        }*/ return true;
    }

    public Double simulate(Assignment assignment, String additionalArgs) { //TODO: handle null return value

        Map<String, String> assignmentIdentifiers = assignment.getIdentifierMap();

        /* handle incomplete assignments of B&B */
        if (assignment.keySet().size() < circuitGates.size()) {
            Map<String, String> dummyInfos = BranchAndBoundUtil.compileDummyInfos(library, circuitGates, assignment);
            if (dummyInfos != null)
                assignmentIdentifiers.putAll(dummyInfos);
        }

        assignmentIdentifiers.put("O", "output_1"); // work around until output gate is part of library

        double score = 0.0;

        try {

            if (!simProcess.isAlive()) {
                logger.error("Simulator exited before simulation start:\n" + getError());
                return null;
            }

            String assignmentStr = mapper.writeValueAsString(assignmentIdentifiers);

            writer.write("update_settings " + simArgs + additionalArgs + " --assignment=" + assignmentStr);
            writer.newLine();
            writer.flush();
            writer.write("start");
            writer.newLine();
            writer.flush();

            String scoreStr = reader.readLine();
            //String growthString = reader.readLine();

            if (scoreStr == null || !scoreStr.startsWith("O ")) {
                logger.error("Simulator did not return legal score value.\n" + getError());
                return null;
            }

            scoreStr = scoreStr.substring(2);

            // relevant to correctly parse infinity as returned score
            score = scoreStr.equals("inf") ? Double.POSITIVE_INFINITY : Double.parseDouble(scoreStr);

            //growthString = growthString.substring(7);
            //growth = Double.parseDouble(growthString);

        } catch (Exception e) {
            e.printStackTrace();
        }

        return score;
    }

    Double growth;

    public Double getLastGrowth() {
        return growth;
    }

    public Double simulate(Assignment assignment) {
        return simulate(assignment, "");
    }

    public void shutdown() {
        simProcess.destroy();
    }

    private String getError() throws IOException {
        String line;
        StringBuilder error = new StringBuilder();
        while (errorReader.ready() && (line = errorReader.readLine()) != null) {
            error.append(line).append("\n");
        }
        return error.toString();
    }
}
