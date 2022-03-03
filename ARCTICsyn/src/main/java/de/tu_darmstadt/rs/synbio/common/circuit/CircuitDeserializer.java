package de.tu_darmstadt.rs.synbio.common.circuit;

import com.fasterxml.jackson.core.JsonParser;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.DeserializationContext;
import com.fasterxml.jackson.databind.deser.std.StdDeserializer;
import de.tu_darmstadt.rs.synbio.common.LogicType;
import org.jgrapht.io.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.Optional;
import java.util.Set;

public class CircuitDeserializer extends StdDeserializer<Circuit> {

    private static final Logger logger = LoggerFactory.getLogger(CircuitDeserializer.class);

    private int supergateCounter = 0;

    public CircuitDeserializer(Class<Circuit> t) {
        super(t);
    }

    @Override
    public Circuit deserialize(JsonParser jsonParser, DeserializationContext deserializationContext) throws IOException, JsonProcessingException {

        Circuit circuit = new Circuit("supergate_" + supergateCounter);
        supergateCounter ++;

        Gate.GateProvider gateProv = new Gate.GateProvider();
        Wire.WireProvider wireProv = new Wire.WireProvider();

        JSONImporter<Gate, Wire> importer = new JSONImporter<Gate, Wire>(gateProv, wireProv);

        String jsonContent = jsonParser.readValueAsTree().toString();
        StringReader reader = new StringReader(jsonContent);

        try {
            importer.importGraph(circuit, reader);
        } catch (ImportException e) {
            e.printStackTrace();
        }

        return circuit;
    }

    public Circuit deserializeString(String content) throws IOException, JsonProcessingException {

        Circuit circuit = new Circuit("circuit_" + supergateCounter);
        supergateCounter ++;

        Gate.GateProvider gateProv = new Gate.GateProvider();
        Wire.WireProvider wireProv = new Wire.WireProvider();

        JSONImporter<Gate, Wire> importer = new JSONImporter<Gate, Wire>(gateProv, wireProv);

        try {
            importer.importGraph(circuit, new StringReader(content));
        } catch (ImportException e) {
            e.printStackTrace();
        }

        // for backwards compatibility: remove output buffer, of output OR gate is present
        if (circuit.vertexSet().stream().filter(g -> g.getLogicType() == LogicType.OUTPUT_OR2 || g.getLogicType() == LogicType.OUTPUT_BUFFER).count() >= 2) {

            Optional<Gate> exessBuffer = circuit.vertexSet().stream().filter(g -> g.getLogicType() == LogicType.OUTPUT_BUFFER).findFirst();

            if (exessBuffer.isEmpty())
                return circuit;

            Set<Wire> edges = circuit.incomingEdgesOf(exessBuffer.get());
            circuit.removeVertex(exessBuffer.get());
            circuit.removeAllEdges(edges);
        }

        return circuit;
    }
}
