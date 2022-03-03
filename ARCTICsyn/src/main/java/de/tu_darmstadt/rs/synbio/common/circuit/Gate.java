package de.tu_darmstadt.rs.synbio.common.circuit;

import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.common.TruthTable;
import de.tu_darmstadt.rs.synbio.synthesis.util.ExpressionParser;
import org.jgrapht.io.Attribute;
import org.jgrapht.io.VertexProvider;
import org.logicng.formulas.Formula;

import java.util.Map;

public class Gate {

    private final String identifier;
    private final LogicType type;

    public Gate(String identifier, LogicType type) {
        this.identifier = identifier;
        this.type = type;
    }

    private static final Map<String, Formula> inputExpressions = Map.of(
            "a", ExpressionParser.parse("a"),
            "b", ExpressionParser.parse("b"),
            "c", ExpressionParser.parse("c"));

    public Formula getExpression() {
        if (type == LogicType.INPUT) {
            if (inputExpressions.containsKey(identifier))
                return inputExpressions.get(identifier);
            else
                return ExpressionParser.parse(identifier);
        } else {
            return type.getExpression();
        }
    }

    public TruthTable getTruthTable() {
        return new TruthTable(getExpression());
    }

    public String getIdentifier() {
        return identifier;
    }

    public int getWeight() {
        return type.getWeight();
    };

    public LogicType getLogicType() {
        return type;
    }

    public boolean isLogicGate() {
        return type != LogicType.INPUT && type != LogicType.OUTPUT_BUFFER && type != LogicType.OUTPUT_OR2;
    }

    @Override
    public String toString() {
        return identifier;
    }

    static class GateProvider implements VertexProvider<Gate> {

        @Override
        public Gate buildVertex(String s, Map<String, Attribute> map) {

            /* this is for backwards compatibility with old structure files */
            if (map.containsKey("expression")) {

                String primitiveIdentifier = map.get("primitiveIdentifier").getValue();

                switch (map.get("type").getValue()) {
                    case "INPUT":
                        return new Gate(s, LogicType.INPUT);
                    case "OUTPUT":
                        return new Gate(s, LogicType.OUTPUT_BUFFER);
                    default:
                        if (primitiveIdentifier.equals("OR2"))
                            return new Gate(s, LogicType.OUTPUT_OR2);
                        else
                            return new Gate(s, LogicType.valueOf(primitiveIdentifier));
                }
            } else {
                return new Gate(s, LogicType.valueOf(map.get("type").getValue()));
            }
        }
    }
}
