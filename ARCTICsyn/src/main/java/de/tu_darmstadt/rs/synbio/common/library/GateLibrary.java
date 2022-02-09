package de.tu_darmstadt.rs.synbio.common.library;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ObjectNode;
import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.synthesis.util.ExpressionParser;
import org.logicng.formulas.Formula;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.util.*;

public class GateLibrary {

    private static final Logger logger = LoggerFactory.getLogger(GateLibrary.class);

    private final File sourceFile;

    private final HashMap<LogicType, List<GateRealization>> gateRealizations = new HashMap<>();

    private final Double[] proxNormalization;
    private final Double[] proxWeights;

    public GateLibrary(File libraryFile, boolean isThermo) {

        this.sourceFile = libraryFile;
        this.proxNormalization = new Double[]{1.0, 1.0, 1.0};
        this.proxWeights = new Double[]{1.0, 1.0, 1.0};

        if (isThermo) {
            loadThermoLibrary(libraryFile);
        } else {
            loadConvLibrary(libraryFile);
        }
    }

    public GateLibrary(File libraryFile, Double[] proxWeights) {

        this.sourceFile = libraryFile;
        this.proxWeights = proxWeights;

        loadConvLibrary(libraryFile);

        proxNormalization = calcProxNormalization();
    }

    /* library file handling */

    private void loadThermoLibrary(File libraryFile) {

        ObjectMapper mapper = new ObjectMapper();

        JsonNode content;

        try {
            content = mapper.readTree(libraryFile);
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }

        JsonNode devices = null;
        JsonNode tfs = null;

        for (int i = 0; i < content.size(); i++) {

            if (content.get(i).get("class").textValue().equals("devices"))
                devices = content.get(i).get("members");

            if (content.get(i).get("class").textValue().equals("transcription_factors"))
                tfs = content.get(i).get("members");
        }

        if (devices == null || tfs == null) {
            logger.error("Invalid thermo library.");
            return;
        }

        Map<String, List<LogicType>> deviceFunctions = new HashMap<>();

        for (JsonNode device : devices) {

            String name = device.get("name").textValue();

            for (JsonNode function : device.get("functions")) {

                if (!deviceFunctions.containsKey(name))
                    deviceFunctions.put(name, new ArrayList<>());

                deviceFunctions.get(name).add(LogicType.valueOf(function.textValue()));
            }
        }

        for (JsonNode tf : tfs) {

            String tfName = tf.get("name").textValue();

            for (JsonNode promoter : tf.get("associated_devices")) {

                String promoterName = promoter.textValue();

                for (LogicType function : deviceFunctions.get(promoterName)) {

                    GateRealization newGate = new GateRealization(promoterName, function, tfName);
                    addGateRealization(newGate);
                }
            }
        }
    }

    private void loadConvLibrary(File libraryFile) {

        HashMap<String, Object>[] parsedRealizations;

        ObjectMapper mapper = new ObjectMapper();

        try {
            parsedRealizations = mapper.readValue(libraryFile, HashMap[].class);
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }

        for (HashMap<String, Object> realization : parsedRealizations) {

            String primitiveIdentifier = (String) Optional.ofNullable(realization.get("primitiveIdentifier"))
                    .orElseThrow(() -> new RuntimeException("Invalid gate library: Key \"primitiveIdentifier\" not found!"));

            String identifier = (String) Optional.ofNullable(realization.get("identifier"))
                    .orElseThrow(() -> new RuntimeException("Invalid gate library: Key \"identifier\" not found!"));

            String group = (String) Optional.ofNullable(realization.get("group"))
                    .orElseThrow(() -> new RuntimeException("Invalid gate library: Key \"group\" not found!"));

            LinkedHashMap biorep = (LinkedHashMap) Optional.ofNullable(realization.get("biorep"))
                    .orElseThrow(() -> new RuntimeException("Invalid gate library: Key \"biorep\" not found!"));

            LinkedHashMap responseFunction = (LinkedHashMap) Optional.ofNullable(biorep.get("response_function"))
                    .orElseThrow(() -> new RuntimeException("Invalid gate library: Key \"response_function\" not found!"));

            LinkedHashMap parameters = (LinkedHashMap) Optional.ofNullable(responseFunction.get("parameters"))
                    .orElseThrow(() -> new RuntimeException("Invalid gate library: Key \"parameters\" not found!"));

            Optional ymaxOpt = Optional.ofNullable(parameters.get("ymax"));
            Optional yminOpt = Optional.ofNullable(parameters.get("ymin"));
            Optional kOpt = Optional.ofNullable(parameters.get("K"));
            Optional nOpt = Optional.ofNullable(parameters.get("n"));

            LinkedHashMap particles = (LinkedHashMap) Optional.ofNullable(biorep.get("particles"))
                    .orElseThrow(() -> new RuntimeException("Invalid gate library: Key \"particles\" not found!"));

            Optional ymaxParticles = Optional.ofNullable(particles.get("ymax"));
            Optional yminParticles = Optional.ofNullable(particles.get("ymin"));

            GateRealization newRealization;

            if (ymaxOpt.isPresent() && yminOpt.isPresent() && kOpt.isPresent() && nOpt.isPresent()) {

                double ymax = (double) ymaxOpt.get();
                double ymin = (double) yminOpt.get();
                double k = (double) kOpt.get();
                double n = (double) nOpt.get();

                GateRealization.Particles gateParticles = null;

                if (ymaxParticles.isPresent() && yminParticles.isPresent()) {
                    List<Double> ymaxList = (ArrayList<Double>) ymaxParticles.get();
                    List<Double> yminList = (ArrayList<Double>) yminParticles.get();
                    gateParticles = new GateRealization.Particles(ymaxList, yminList);
                }

                newRealization = new GateRealization(identifier, LogicType.valueOf(primitiveIdentifier), group,
                        new GateRealization.GateCharacterization(ymax, ymin, k ,n, gateParticles));

            } else {
                newRealization = new GateRealization(identifier, LogicType.valueOf(primitiveIdentifier), group);
            }

            addGateRealization(newRealization);
        }
    }

    private void addGateRealization(GateRealization newGate) {
        if (!gateRealizations.containsKey(newGate.getLogicType())) {
            gateRealizations.put(newGate.getLogicType(), new ArrayList<>());
        }
        gateRealizations.get(newGate.getLogicType()).add(newGate);
    }

    private Double[] calcProxNormalization() {

        Double[] normalizationValues = new Double[]{1.0, 1.0, 1.0};

        double ym_max = Double.NEGATIVE_INFINITY;
        double xm_max = Double.NEGATIVE_INFINITY;
        double grad_max = Double.NEGATIVE_INFINITY;
        double ym_min = Double.POSITIVE_INFINITY;
        double xm_min = Double.POSITIVE_INFINITY;
        double grad_min = Double.POSITIVE_INFINITY;

        for (List<GateRealization> allRealizations : gateRealizations.values()) {

            if (allRealizations.size() < 2)
                continue;

            for (GateRealization r1 : allRealizations) {

                ym_max = Math.max(ym_max, r1.getCharacterization().getYm());
                xm_max = Math.max(xm_max, r1.getCharacterization().getXm());
                grad_max = Math.max(grad_max, r1.getCharacterization().getGrad());

                ym_min = Math.min(ym_min, r1.getCharacterization().getYm());
                xm_min = Math.min(xm_min, r1.getCharacterization().getXm());
                grad_min = Math.min(grad_min, r1.getCharacterization().getGrad());
            }
        }

        normalizationValues[0] = xm_max - xm_min;
        normalizationValues[1] = ym_max - ym_min;
        normalizationValues[2] = grad_max - grad_min;

        return normalizationValues;
    }

    /* getter and utility functions */

    public File getSourceFile() {
        return sourceFile;
    }

    public HashMap<LogicType, List<GateRealization>> getRealizations() {
        return gateRealizations;
    }

    public void print() {

        logger.info("Circuit library:");

        logger.info("Gate number constraints:");

        for (LogicType type : gateRealizations.keySet()) {
            logger.info(type.name() + ": " + gateRealizations.get(type).size());
        }
    }

    public Integer getNumAvailableGates(LogicType type) {

        if (gateRealizations.get(type) == null)
            return 0;

        return gateRealizations.get(type).size();
    }

    public List<LogicType> getGateTypes() {
        return new ArrayList<>(gateRealizations.keySet());
    }

    public Double[] getProxWeights() {
        return proxWeights;
    }

    public Double[] getProxNormalization() { return proxNormalization; }

    public Set<String> getGroups() {

        Set<String> groupSet = new HashSet<>();

        for (List<GateRealization> realizationsOfType : gateRealizations.values()) {
            for (GateRealization realization : realizationsOfType) {
                groupSet.add(realization.getGroup());
            }
        }

        return groupSet;
    }

}