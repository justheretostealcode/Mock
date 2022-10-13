package de.tu_darmstadt.rs.synbio.synthesis.enumeration;

public class EnumerationResult {

    private final PrimitiveCircuit circuit;
    private final String inputMapping;

    public EnumerationResult(PrimitiveCircuit circuit, String inputMapping) {
        this.circuit = circuit;
        this.inputMapping = inputMapping;
    }

    public String getInputMapping() {
        return inputMapping;
    }

    public PrimitiveCircuit getCircuit() {
        return circuit;
    }
}
