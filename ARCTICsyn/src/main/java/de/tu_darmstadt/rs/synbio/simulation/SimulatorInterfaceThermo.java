package de.tu_darmstadt.rs.synbio.simulation;

import com.fasterxml.jackson.databind.ObjectMapper;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;

public class SimulatorInterfaceThermo {

    private static final Logger logger = LoggerFactory.getLogger(SimulatorInterfaceThermo.class);

    private static final String scorePrefix = "score: ";

    private final String pythonBinary;
    private final File simulatorPath;
    private final String simScript;
    private final String simInitArgs;
    private final String simArgs;
    private final GateLibrary library;

    private ProcessBuilder pb;
    private Process simProcess;
    private File structureFile;
    private BufferedReader reader;
    private BufferedWriter writer;
    private BufferedReader errorReader;
    private final ObjectMapper mapper = new ObjectMapper();
    private AssignmentCompiler compiler;

    public enum PropagationMode {
        NORMAL,
        _NAIVE_OPTIMAL,
        _NAIVE_HEURISTIC,
        ITA,
        _ITA_HEURISTIC,
        _FULL_HEURISTIC,
        OPTIMAL
    }

    public SimulatorInterfaceThermo(SimulationConfiguration config, GateLibrary gateLibrary) {
        pythonBinary = config.getPythonBinary();
        simulatorPath = config.getSimPath();
        simScript = config.getSimScript();
        simInitArgs = config.getSimInitArgs();
        simArgs = config.getSimArgs();
        library = gateLibrary;
    }

    public boolean initSimulation(Circuit circuit) { //TODO: handle return value

        if (simProcess!= null && simProcess.isAlive())
            simProcess.destroy();

        this.compiler = new AssignmentCompiler(circuit, library);

        try {
            String structureFileName = "structure_" + circuit.getIdentifier() + "_tid" + Thread.currentThread().getId() + "_" + System.nanoTime() + ".json";
            structureFile = new File(simulatorPath, structureFileName);
            circuit.save(structureFile);

            List<String> arguments = new ArrayList<>();
            arguments.addAll(Arrays.asList(pythonBinary, simScript, "-s=" + structureFileName, "-l=" + library.getSourceFile().getAbsolutePath()));
            arguments.addAll(Arrays.asList(simInitArgs.split(" ")));

            pb = new ProcessBuilder(arguments.toArray(new String[0]));
            pb.directory(simulatorPath);

            simProcess = pb.start();

            reader = new BufferedReader(new InputStreamReader(simProcess.getInputStream()));
            writer = new BufferedWriter(new OutputStreamWriter(simProcess.getOutputStream()));
            errorReader = new BufferedReader(new InputStreamReader(simProcess.getErrorStream()));

            String output;
            do {
                output = reader.readLine();

                if (!simProcess.isAlive()) {
                    logger.error("Simulator exited on initialization:\n" + getError());
                    shutdown();
                    return false;
                }

            } while (output == null || !output.startsWith("ready"));

            return true;

        } catch (Exception e) {
            e.printStackTrace();
            return false;
        }
    }



    public Double simulate(Assignment assignment, PropagationMode mode) {

        SimulationResult result;

        result = trySimulation(assignment, mode);

        switch (result.code) {
            case ERROR:
                return null;
            case NO_CONVERGENCE:
                return 0.0;
            default:
                return result.score;
        }
    }

    private static class SimulationResult {

        public Double score;
        public ReturnCode code;

        public enum ReturnCode {
            OK,
            NO_CONVERGENCE,
            ERROR
        }

        public SimulationResult(Double score, ReturnCode code) {
            this.score = score;
            this.code = code;
        }
    }
    private SimulationResult trySimulation(Assignment assignment, PropagationMode mode) {

        String additionalArgs = "--propagation_mode=" + mode.ordinal();

        Map<String, Map<?, ?>> assignmentMap = compiler.compile(assignment);

        double score = 0.0;

        try {

            if (!simProcess.isAlive()) {
                logger.error("Simulator exited before simulation start:\n" + getError());
                shutdown();
                return new SimulationResult(null, SimulationResult.ReturnCode.ERROR);
            }

            String assignmentStr = mapper.writeValueAsString(assignmentMap);

            writer.write("update_settings " + simArgs + additionalArgs + " --assignment=" + assignmentStr);
            writer.newLine();
            writer.flush();
            writer.write("start");
            writer.newLine();
            writer.flush();

            String output;
            do {
                output = reader.readLine();

                if (!simProcess.isAlive()) {
                    logger.error("Simulator exited during simulation: " + getError());
                    shutdown();
                    return new SimulationResult(null, SimulationResult.ReturnCode.ERROR);
                }

                if (output != null && output.equals("error above tolerance, score dismissed.")) {
                    return new SimulationResult(null, SimulationResult.ReturnCode.NO_CONVERGENCE);
                }
            } while (output == null || !output.startsWith(scorePrefix));

            output = output.substring(scorePrefix.length());

            // relevant to correctly parse infinity as returned score
            try {
                score = output.equals("inf") ? Double.POSITIVE_INFINITY : Double.parseDouble(output);
            } catch (Exception e) {
                logger.error(e.getMessage());
            }

        } catch (Exception e) {
            e.printStackTrace();
            shutdown();
        }

        /*try {
            mapper.writerWithDefaultPrettyPrinter().writeValue(new File("assignments", "assignment_00000110_"+ num + "_" + score + ".json"), assignmentMap);
            num++;
        } catch (Exception e) {
        }*/

        return new SimulationResult(score, SimulationResult.ReturnCode.OK);
    }

    public void shutdown() {
        structureFile.delete();
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
