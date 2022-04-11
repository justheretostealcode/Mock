package de.tu_darmstadt.rs.synbio.common.library;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.mapping.compatibility.CompatibilityMatrix;
import org.logicng.util.Pair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.util.*;

public class GateLibrary {

    private static final Logger logger = LoggerFactory.getLogger(GateLibrary.class);

    private final File sourceFile;

    private final File compatibilityFile;

    private final Map<LogicType, List<GateRealization>> gateRealizations = new HashMap<>();
    private final Map<String, Map<String, Double>> devicePromoterFactors = new HashMap<>();

    private Double[] proxNormalization = new Double[]{1.0, 1.0, 1.0};
    private Double[] proxWeights = new Double[]{1.0, 1.0, 1.0};

    public GateLibrary(File libraryFile, File compatibilityFile, boolean isThermo) {

        this.sourceFile = libraryFile;
        this.compatibilityFile = compatibilityFile;

        if (isThermo) {
            if (compatibilityFile != null)
                loadCompatibility(compatibilityFile);
            loadThermoLibrary(libraryFile);
        } else {
            loadConvLibrary(libraryFile);
        }
    }

    public GateLibrary(File libraryFile, Double[] proxWeights) {

        this(libraryFile, null,false);

        this.proxWeights = proxWeights;
        proxNormalization = calcProxNormalization();
    }

    /* library file handling */

    private Map<String, List<Double>> inputs3dB;
    private Map<String, List<Double>> outputs3dB;

    private CompatibilityMatrix<String> compatibilityMatrix;

    private void loadCompatibility(File compatibilityFile) {

        ObjectMapper mapper = new ObjectMapper();

        JsonNode content;

        try {
            content = mapper.readTree(compatibilityFile);
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }

        /* load 3dB in-/outputs */

        inputs3dB = new HashMap<>();
        outputs3dB = new HashMap<>();

        Iterator<Map.Entry<String, JsonNode>> deviceIterator = content.get("critical_points").fields();

        deviceIterator.forEachRemaining(e -> {
                    inputs3dB.put(e.getKey(), Arrays.asList(e.getValue().get("x_c").get(0).doubleValue(), e.getValue().get("x_c").get(1).doubleValue()));
                    outputs3dB.put(e.getKey(), Arrays.asList(e.getValue().get("y_c").get(0).doubleValue(), e.getValue().get("y_c").get(1).doubleValue()));
                });

        /* load compatibility info */

        compatibilityMatrix = new CompatibilityMatrix<>();

        List<String> deviceOrder = new ArrayList<>();
        Iterator<JsonNode> orderIterator = content.get("compatibility").get("order").elements();
        orderIterator.forEachRemaining(d -> deviceOrder.add(d.asText()));

        JsonNode compatPairs = content.get("compatibility").get("table_2");
        JsonNode compatTriples = content.get("compatibility").get("table_3");

        for (int s = 0; s < deviceOrder.size(); s ++) {
            for (int t = 0; t < deviceOrder.size(); t ++) {

                boolean isCompatible = compatPairs.get(s).get(t).doubleValue() == 1.0;
                compatibilityMatrix.addEntry(deviceOrder.get(s), deviceOrder.get(t), null, isCompatible);

                for (int ss = 0; ss < deviceOrder.size(); ss ++) {

                    isCompatible = compatTriples.get(s).get(ss).get(t).doubleValue() == 1.0;
                    compatibilityMatrix.addEntry(deviceOrder.get(s), deviceOrder.get(t), deviceOrder.get(ss), isCompatible);
                }
            }
        }
    }

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
        JsonNode promoters = null;

        /* get and check root nodes */

        for (int i = 0; i < content.size(); i++) {

            if (content.get(i).get("class").textValue().equals("devices"))
                devices = content.get(i).get("members");

            if (content.get(i).get("class").textValue().equals("transcription_factors"))
                tfs = content.get(i).get("members");

            if (content.get(i).get("class").textValue().equals("promoters"))
                promoters = content.get(i).get("members");
        }

        if (devices == null || tfs == null || promoters == null) {
            logger.error("Invalid thermo library.");
            return;
        }

        /* get logical functions of devices */

        Map<String, List<LogicType>> deviceFunctions = new HashMap<>();

        for (JsonNode device : devices) {

            String name = device.get("name").textValue();

            for (JsonNode function : device.get("functions")) {

                if (!deviceFunctions.containsKey(name))
                    deviceFunctions.put(name, new ArrayList<>());

                deviceFunctions.get(name).add(LogicType.valueOf(function.textValue()));
            }
        }

        /* get mapping of devices to promoters and maximum promoter levels */

        Map<String, String> deviceToPromoter = new HashMap<>();
        Map<String, Pair<Double, Double>> promoterLevels = new HashMap<>();

        for (JsonNode promoter : promoters) {

            String name = promoter.get("name").textValue();

            if (!promoterLevels.containsKey(name))
                promoterLevels.put(name, new Pair<>(promoter.get("levels").get("off").asDouble(), promoter.get("levels").get("on").asDouble()));

            if (promoter.get("associated_devices").size() != 1) {
                logger.error("Thermo library invalid: Promoter " + name + " has invalid number (!= 1) of associated devices.");
                return;
            }

            String device = promoter.get("associated_devices").get(0).textValue();

            if (!deviceToPromoter.containsKey(device))
                deviceToPromoter.put(device, name);

            /* fill device <-> promoter factors */
            Iterator<Map.Entry<String, JsonNode>> it = promoter.get("factors").get("tf_only").fields();
            devicePromoterFactors.put(device, new HashMap<>());
            it.forEachRemaining(e -> devicePromoterFactors.get(device).put(e.getKey(), e.getValue().get(e.getValue().size() - 1).doubleValue()));
        }

        /* add gate for each tf by getting its device and promoter */

        for (JsonNode tf : tfs) {

            String tfName = tf.get("name").textValue();

            for (JsonNode device : tf.get("associated_devices")) {

                String promoterName = deviceToPromoter.get(device.textValue());
                String deviceName = device.textValue();

                for (LogicType function : deviceFunctions.get(device.textValue())) {

                    GateRealization newGate;

                    if (function == LogicType.OUTPUT_BUFFER || function == LogicType.OUTPUT_OR2) {
                        newGate = new GateRealization(deviceName, function, tfName);
                    } else {
                        newGate = new GateRealization(deviceName, function, tfName,
                                new GateRealization.GateCharacterization(
                                        promoterLevels.get(promoterName).second(),
                                        promoterLevels.get(promoterName).first(),
                                        function == LogicType.INPUT ? 0.0 : outputs3dB.get(deviceName).get(0),
                                        function == LogicType.INPUT ? 0.0 : outputs3dB.get(deviceName).get(1),
                                        function == LogicType.INPUT ? 0.0 : inputs3dB.get(deviceName).get(0),
                                        function == LogicType.INPUT ? 0.0 : inputs3dB.get(deviceName).get(1),
                                        0,
                                        0,
                                        null));
                    }
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
                        new GateRealization.GateCharacterization(ymax, ymin, 0.0, 0.0, 0.0, 0.0, k ,n, gateParticles));

            } else {
                newRealization = new GateRealization(identifier, LogicType.valueOf(primitiveIdentifier), group);
            }

            addGateRealization(newRealization);
        }

        /* add input gates */

        addGateRealization(new GateRealization("pTac", LogicType.INPUT, "pTac",
                new GateRealization.GateCharacterization(2.8, 0.0034 , 0.0, 0.0, 0.0, 0.0, 0, 0, null)));

        addGateRealization(new GateRealization("pTet", LogicType.INPUT, "pTet",
                new GateRealization.GateCharacterization(4.4, 0.0013 , 0.0, 0.0, 0.0, 0.0, 0, 0, null)));

        addGateRealization(new GateRealization("pBAD", LogicType.INPUT, "pBAD",
                new GateRealization.GateCharacterization(2.5, 0.0082 , 0.0, 0.0, 0.0, 0.0, 0, 0, null)));
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

    public Map<LogicType, List<GateRealization>> getRealizations() {
        return gateRealizations;
    }

    public GateRealization getOutputDevice(LogicType type) {

        if (gateRealizations.get(type).size() != 1) {
            logger.error("Unsupported number of output devices in library != 1.");
            return null;
        }

        return gateRealizations.get(type).get(0);
    }

    public Map<String, Double> getTfFactorsForDevice(String promoter) {
        return devicePromoterFactors.get(promoter);
    }

    public CompatibilityMatrix<String> getCompatibilityMatrix() {
        return compatibilityMatrix;
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

    public List<LogicType> getLogicGateTypes() {

        List<LogicType> types = new ArrayList<>(gateRealizations.keySet());
        types.remove(LogicType.INPUT);
        types.remove(LogicType.OUTPUT_BUFFER);
        types.remove(LogicType.OUTPUT_OR2);

        return types;
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