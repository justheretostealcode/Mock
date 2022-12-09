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

    private final Map<LogicType, List<GateRealization>> gateRealizations = new HashMap<>();
    private final Map<String, Map<String, Double>> devicePromoterFactors = new HashMap<>();

    private Double[] proxNormalization = new Double[]{1.0, 1.0, 1.0};
    private Double[] proxWeights = new Double[]{1.0, 1.0, 1.0};

    public GateLibrary(File libraryFile, File compatibilityFile, boolean isThermo) {

        this.sourceFile = libraryFile;

        if (isThermo) {
            if (compatibilityFile != null)
                loadCompatibility(compatibilityFile);
            loadThermoLibrary(libraryFile);
        } else {
            loadConvLibrary(libraryFile);
        }
    }

    public GateLibrary(File libraryFile, File compatibilityFile, boolean isThermo, Double[] proxWeights) {

        this(libraryFile, compatibilityFile,isThermo);

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

        /* get mapping of devices to promoters, maximum promoter levels and hill parameters */

        Map<String, String> deviceToPromoter = new HashMap<>();
        Map<String, Pair<Double, Double>> promoterLevels = new HashMap<>();
        Map<String, Double> k = new HashMap<>(), n = new HashMap<>();
        
        for (JsonNode promoter : promoters) {

            String name = promoter.get("name").textValue();

            if (!promoterLevels.containsKey(name))
                promoterLevels.put(name, new Pair<>(promoter.get("levels").get("off").asDouble(), promoter.get("levels").get("on").asDouble()));

            if (promoter.get("hill_parameters") != null) {

                if (!k.containsKey(name))
                    k.put(name, promoter.get("hill_parameters").get("K").asDouble());

                if (!n.containsKey(name))
                    n.put(name, promoter.get("hill_parameters").get("n").asDouble());
            }

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

        Boolean compat = compatibilityMatrixLoaded();

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
                                        compat ? function == LogicType.INPUT ? promoterLevels.get(promoterName).first() : outputs3dB.get(deviceName).get(0) : 0.0,
                                        compat ? function == LogicType.INPUT ? promoterLevels.get(promoterName).second() : outputs3dB.get(deviceName).get(1) : 0.0,
                                        compat ? function == LogicType.INPUT ? 0.0 : inputs3dB.get(deviceName).get(0) : 0.0,
                                        compat ? function == LogicType.INPUT ? 0.0 : inputs3dB.get(deviceName).get(1) : 0.0,
                                        k.getOrDefault(deviceName, 0.0),
                                        n.getOrDefault(deviceName, 0.0),
                                        null));
                    }
                    addGateRealization(newGate);
                }
            }
        }
    }

    private void loadConvLibrary(File libraryFile) {

        ObjectMapper mapper = new ObjectMapper();

        JsonNode content;

        try {
            content = mapper.readTree(libraryFile);
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }

        for (JsonNode r : content) {

            String primitiveIdentifier = r.get("primitiveIdentifier").textValue();

            if (primitiveIdentifier.equals("OR2"))
                primitiveIdentifier = "OUTPUT_OR2";

            String identifier = r.get("identifier").textValue();
            String group = r.get("group").textValue();

            JsonNode parameters = r.get("biorep").get("response_function").get("parameters");

            if (r.get("biorep").get("response_function").get("parameters").has("ymax")) {

                double yMax = parameters.get("ymax").asDouble();
                double yMin = parameters.get("ymin").asDouble();
                double k = parameters.get("K").asDouble();
                double n = parameters.get("n").asDouble();

                addGateRealization(new GateRealization(identifier, LogicType.valueOf(primitiveIdentifier), group,
                        new GateRealization.GateCharacterization(yMax, yMin, 0.0, 0.0, 0.0, 0.0, k ,n, null)));
            } else {
                addGateRealization(new GateRealization(identifier, LogicType.valueOf(primitiveIdentifier), group));
            }
        }

        /* add input gates */

        addGateRealization(new GateRealization("pTac", LogicType.INPUT, "pTac",
                new GateRealization.GateCharacterization(2.8, 0.0034 , 0.0, 0.0, 0.0, 0.0, 0, 0, null)));

        addGateRealization(new GateRealization("pTet", LogicType.INPUT, "pTet",
                new GateRealization.GateCharacterization(4.4, 0.0013 , 0.0, 0.0, 0.0, 0.0, 0, 0, null)));

        addGateRealization(new GateRealization("pBAD", LogicType.INPUT, "pBAD",
                new GateRealization.GateCharacterization(2.5, 0.0082 , 0.0, 0.0, 0.0, 0.0, 0, 0, null)));

        /* add output buffer */

        addGateRealization(new GateRealization("output_1", LogicType.OUTPUT_BUFFER, "YFP"));
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

        for (GateRealization realization : gateRealizations.get(LogicType.NOR2)) {
            ym_max = Math.max(ym_max, realization.getCharacterization().getYm());
            xm_max = Math.max(xm_max, realization.getCharacterization().getXm());
            grad_max = Math.max(grad_max, realization.getCharacterization().getGrad());

            ym_min = Math.min(ym_min, realization.getCharacterization().getYm());
            xm_min = Math.min(xm_min, realization.getCharacterization().getXm());
            grad_min = Math.min(grad_min, realization.getCharacterization().getGrad());
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

    public Boolean compatibilityMatrixLoaded() {
        return compatibilityMatrix != null;
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