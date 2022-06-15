package de.tu_darmstadt.rs.synbio.simulation;

import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.circuit.Gate;
import de.tu_darmstadt.rs.synbio.common.library.GateLibrary;
import de.tu_darmstadt.rs.synbio.common.library.GateRealization;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.stream.Collectors;

public class AssignmentCompiler {

    private static final Logger logger = LoggerFactory.getLogger(AssignmentCompiler.class);

    private final Circuit circuit;
    private final GateLibrary library;

    private final Set<Gate> circuitGates;

    public AssignmentCompiler(Circuit circuit, GateLibrary library) {

        this.circuit = circuit;
        this.library = library;

        this.circuitGates = circuit.vertexSet().stream()
                .filter(g -> g.getLogicType() != LogicType.OUTPUT_BUFFER)
                .filter(g -> g.getLogicType() != LogicType.OUTPUT_OR2)
                .collect(Collectors.toSet());
    }

    public Map<String, Map<?, ?>> compile(Assignment assignment) {

        /* possible refinement: treat possible dummies exclusively when compiling c and k lists */

        /* initialize */

        Map<String, Map<?, ?>> output = new HashMap<>();

        Set<Gate> dummyGates = circuitGates.stream().filter(g -> !assignment.keySet().contains(g)).collect(Collectors.toSet());

        Set<String> usedGroups = new HashSet<>();
        assignment.values().forEach(r -> usedGroups.add(r.getGroup()));

        /* write assigned gates */

        for (Gate gate : assignment.keySet()) {

            Map<String, Object> gateMap = new HashMap<>();

            if (assignment.get(gate).isCharacterized()) {
                List<Double> iList = Arrays.asList(assignment.get(gate).getCharacterization().getILower(), assignment.get(gate).getCharacterization().getIUpper());
                gateMap.put("i", iList);

                List<Double> jList = Arrays.asList(assignment.get(gate).getCharacterization().getJLower(), assignment.get(gate).getCharacterization().getJUpper());
                gateMap.put("j", jList);
            }

            /* compose c and k lists */

            Set<Gate> sourceGates = circuit.incomingEdgesOf(gate).stream().map(circuit::getEdgeSource).collect(Collectors.toSet());

            List<Double> c0Vals = new ArrayList<>();
            List<Double> c1Vals = new ArrayList<>();
            List<Double> k0Vals = new ArrayList<>();
            List<Double> k1Vals = new ArrayList<>();

            for (Gate sourceGate : sourceGates) {

                if (assignment.keySet().contains(sourceGate)) { /* source gate is assigned */

                    GateRealization sourceRealization = assignment.get(sourceGate);

                    if (!sourceRealization.isCharacterized())
                        continue;

                    c0Vals.add(sourceRealization.getCharacterization().getYmin());
                    c1Vals.add(sourceRealization.getCharacterization().getYmax());
                    k0Vals.add(sourceRealization.getCharacterization().getILower());
                    k1Vals.add(sourceRealization.getCharacterization().getIUpper());

                } else { /* source gate is dummy */

                    Set<GateRealization> possibleRealizations = library.getRealizations().get(sourceGate.getLogicType()).stream()
                            .filter(r -> !usedGroups.contains(r.getGroup()))
                            .filter(GateRealization::isCharacterized)
                            .collect(Collectors.toSet());

                    Optional<Double> c0 = possibleRealizations.stream()
                            .map(GateRealization::getCharacterization)
                            .map(GateRealization.GateCharacterization::getYmin)
                            .min(Double::compareTo);

                    Optional<Double> c1 = possibleRealizations.stream()
                            .map(GateRealization::getCharacterization)
                            .map(GateRealization.GateCharacterization::getYmax)
                            .max(Double::compareTo);

                    Optional<Double> k0 = possibleRealizations.stream()
                            .map(GateRealization::getCharacterization)
                            .map(GateRealization.GateCharacterization::getILower)
                            .max(Double::compareTo);

                    Optional<Double> k1 = possibleRealizations.stream()
                            .map(GateRealization::getCharacterization)
                            .map(GateRealization.GateCharacterization::getIUpper)
                            .min(Double::compareTo);

                    if (c0.isPresent()) {
                        c0Vals.add(c0.get());
                        c1Vals.add(c1.get());
                    }

                    if (k0.isPresent()) {
                        k0Vals.add(k0.get());
                        k1Vals.add(k1.get());
                    }
                }
            }

            gateMap.put("c", Arrays.asList(c0Vals.stream().reduce(0.0, Double::sum), c1Vals.stream().reduce(0.0, Double::sum)));

            if (gate.getLogicType() == LogicType.OUTPUT_OR2)
                gateMap.put("k", Arrays.asList(k0Vals.stream().reduce(0.0, Double::sum), c0Vals.stream().min(Double::compareTo).get() + k1Vals.stream().min(Double::compareTo).get()));

            gateMap.put("d", assignment.get(gate).getIdentifier());

            output.put(gate.getIdentifier(), gateMap);
        }

        /* write dummy gates */

        /* sets of used TFs for dummy gates for each gate */
        Map<Gate, Set<String>> setMinFactors = new HashMap<>();
        Map<Gate, Set<String>> setMaxFactors = new HashMap<>();

        for (Gate dummy : dummyGates) {

            Map<String, Map<?, ?>> dummyMap = new HashMap<>();
            String dummyName = dummy.getIdentifier();
            LogicType dummyType = dummy.getLogicType();

            Set<GateRealization> possibleRealizations = library.getRealizations().get(dummyType).stream()
                    .filter(r -> !usedGroups.contains(r.getGroup()))
                    .collect(Collectors.toSet());

            /* list a contains widest promoter activity interval for given dummy logic type */

            Optional<Double> minActivity = possibleRealizations.stream()
                    .map(r -> r.getCharacterization().getYmin())
                    .min(Double::compareTo);

            Optional<Double> maxActivity = possibleRealizations.stream()
                    .map(r -> r.getCharacterization().getYmax())
                    .max(Double::compareTo);

            if (minActivity.isEmpty()) {
                logger.error("Unable to obtain activities for dummy gate.");
                return null;
            }

            List<Double> a = Arrays.asList(minActivity.get(), maxActivity.get());

            /* list c */

            List<LogicType> sourceTypes = circuit.incomingEdgesOf(dummy).stream()
                    .map(circuit::getEdgeSource)
                    .map(Gate::getLogicType)
                    .collect(Collectors.toList());

            List<Double> c0Vals = new ArrayList<>();
            List<Double> c1Vals = new ArrayList<>();

            for (LogicType sourceDummy : sourceTypes) {

                Set<GateRealization> possibleSourceRealizations = library.getRealizations().get(sourceDummy).stream()
                        .filter(r -> !usedGroups.contains(r.getGroup()))
                        .filter(GateRealization::isCharacterized)
                        .collect(Collectors.toSet());

                Optional<Double> c0 = possibleSourceRealizations.stream()
                        .map(GateRealization::getCharacterization)
                        .map(GateRealization.GateCharacterization::getYmin)
                        .min(Double::compareTo);

                Optional<Double> c1 = possibleSourceRealizations.stream()
                        .map(GateRealization::getCharacterization)
                        .map(GateRealization.GateCharacterization::getYmax)
                        .max(Double::compareTo);

                if (c0.isPresent()) {
                    c0Vals.add(c0.get());
                    c1Vals.add(c1.get());
                }
            }

            List<Double> c = Arrays.asList(c0Vals.stream().reduce(0.0, Double::sum), c1Vals.stream().reduce(0.0, Double::sum));

            /* list i contains narrowest output promoter activity interval for given dummy logic type */

            Optional<Double> i0 = possibleRealizations.stream()
                    .map(r -> r.getCharacterization().getILower())
                    .max(Double::compareTo);

            Optional<Double> i1 = possibleRealizations.stream()
                    .map(r -> r.getCharacterization().getIUpper())
                    .min(Double::compareTo);

            List<Double> i = Arrays.asList(i0.get(), i1.get());

            /* list j contains narrowest input promoter activity interval for given dummy logic type */

            Optional<Double> j0 = possibleRealizations.stream()
                    .map(r -> r.getCharacterization().getJLower())
                    .max(Double::compareTo);

            Optional<Double> j1 = possibleRealizations.stream()
                    .map(r -> r.getCharacterization().getJUpper())
                    .min(Double::compareTo);

            List<Double> j = Arrays.asList(j0.get(), j1.get());

            /* relations to assigned gates */

            for (Gate gate : assignment.keySet()) {

                String gateName = gate.getIdentifier();

                setMinFactors.putIfAbsent(gate, new HashSet<>());
                setMaxFactors.putIfAbsent(gate, new HashSet<>());

                /*
                    list b contains TFs with least/biggest binding factor to given gate
                    - factors corresponding to given gate promoter are streamed and their TFs are filtered to match dummy gate type
                    - TFs used in assignment are excluded
                    - TFs used in other dummy gates for current gate are excluded
                 */

                List<String> b = new ArrayList<>();

                /* if gate is YFP -> add YFP as TFs dummy-wise */
                if (gate.getLogicType() == LogicType.OUTPUT_OR2 || gate.getLogicType() == LogicType.OUTPUT_BUFFER) {

                    b.add(assignment.get(gate).getGroup());
                    b.add(assignment.get(gate).getGroup());

                } else {

                    Set<Map.Entry<String, Double>> factors = library.getTfFactorsForDevice(assignment.get(gate).getIdentifier()).entrySet().stream()
                            .filter(e -> library.getRealizations().get(dummyType).stream().map(GateRealization::getGroup).collect(Collectors.toList()).contains(e.getKey()))
                            .filter(e -> !usedGroups.contains(e.getKey()))
                            .collect(Collectors.toSet());

                    Optional<Map.Entry<String, Double>> minFactor = factors.stream()
                            .filter(e -> !setMinFactors.get(gate).contains(e.getKey()))
                            .min(Map.Entry.comparingByValue());

                    Optional<Map.Entry<String, Double>> maxFactor = factors.stream()
                            .filter(e -> !setMaxFactors.get(gate).contains(e.getKey()))
                            .max(Map.Entry.comparingByValue());

                    if (minFactor.isEmpty() || maxFactor.isEmpty()) {
                        logger.error("Unable to obtain TF factors for dummy gate.");
                        return null;
                    }

                    b.add(minFactor.get().getKey());
                    b.add(maxFactor.get().getKey());

                    setMinFactors.get(gate).add(minFactor.get().getKey());
                    setMaxFactors.get(gate).add(maxFactor.get().getKey());
                }

                Map<String, List<?>> gateMap = new HashMap<>();
                gateMap.put("a", a);
                gateMap.put("b", b);
                gateMap.put("c", c);
                gateMap.put("i", i);
                gateMap.put("j", j);

                dummyMap.put(gateName, gateMap);
            }

            try {
                output.put(dummyName, dummyMap);
            } catch (Exception e) {
                logger.error(e.getMessage());
                return null;
            }
        }

        return output;
    }
}
