package de.tu_darmstadt.rs.synbio.mapping;

import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;
import de.tu_darmstadt.rs.synbio.common.circuit.Gate;
import de.tu_darmstadt.rs.synbio.common.circuit.Wire;
import de.tu_darmstadt.rs.synbio.common.library.GateRealization;
import de.tu_darmstadt.rs.synbio.common.circuit.LogicGate;

import java.util.*;

public class Assignment {

    private final Map<LogicGate, GateRealization> map;

    public Assignment() {
        this.map = new HashMap<>();
    }

    public Assignment(Assignment assignment) {
        this.map = new HashMap<>(assignment.map);
    }

    public void put(LogicGate circuitGate, GateRealization realization) {
        map.put(circuitGate, realization);
    }

    public Set<LogicGate> keySet() {
        return map.keySet();
    }

    public GateRealization get(LogicGate circuitGate) {
        return map.get(circuitGate);
    }

    public Collection<GateRealization> values() {
        return map.values();
    }

    public int size() {
        return map.size();
    }

    @Override
    public boolean equals(Object o) {

        if (o == this)
            return true;

        if (!(o instanceof Assignment))
            return false;

        return map.equals(((Assignment) o).map);
    }

    @Override
    public int hashCode() {
        return map.hashCode();
    }

    public Map<String, String> getIdentifierMap() {

        Map<String, String> stringMap = new HashMap<>();

        for (LogicGate circuitGate : map.keySet()) {
            stringMap.put(circuitGate.getIdentifier(), map.get(circuitGate).getIdentifier());
        }

        return stringMap;
    }

    public boolean isValid() {

        // check gate instance redundancy
        List<GateRealization> realizationList = new ArrayList<>(map.values());
        Set<GateRealization> realizationSet = new HashSet<>(map.values());

        if (realizationList.size() != realizationSet.size()) {
            return false;
        }

        // check group constraints
        List<String> usedGroups = new ArrayList<>();

        for (GateRealization realization : map.values()) {
            if (usedGroups.contains(realization.getGroup())) {
                return false;
            } else {
                usedGroups.add(realization.getGroup());
            }
        }

        return true;
    }

    public boolean fulfilsConstraints(Circuit structure) {

        for (LogicGate dest : keySet()) {

            if (!dest.getLogicType().equals(LogicType.NOR2))
                continue;

            Set<Wire> wires = structure.incomingEdgesOf(dest);

            List<LogicGate> gates = new ArrayList<>();
            for(Wire wire : wires) {
                Gate gateSrc = structure.getEdgeSource(wire);

                if (gateSrc instanceof LogicGate)
                    gates.add((LogicGate)gateSrc);
            }

            if (gates.size() == 2) {

                GateRealization gr0 = get(gates.get(0));
                GateRealization gr1 = get(gates.get(1));

                if (gr0 == null || gr1 == null)
                    continue;

                if (gr0.getGroup().equals("PhlF") || gr0.getGroup().equals("SrpR") || gr0.getGroup().equals("BM3R1") || gr0.getGroup().equals("QacR"))
                    if (gr1.getGroup().equals("PhlF") || gr1.getGroup().equals("SrpR") || gr1.getGroup().equals("BM3R1") || gr1.getGroup().equals("QacR"))
                        return false;
            }

        }
        return true;
    }

    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();
        builder.append("{");
        map.entrySet().stream().forEach(logicGateGateRealizationEntry -> {
            builder.append("\"");
            builder.append(logicGateGateRealizationEntry.getKey().getIdentifier());
            builder.append(" : ");
            builder.append(logicGateGateRealizationEntry.getValue().getIdentifier());
            builder.append("\", ");
        });
        builder.append("}");
        return builder.toString().replace(", }", "}");
    }
}
