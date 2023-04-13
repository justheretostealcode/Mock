package de.tu_darmstadt.rs.synbio.synthesis.enumeration;

public class EnumerationResult {

    private final TreeCircuit circuit;
    private final String inputMapping;

    public EnumerationResult(TreeCircuit circuit, String inputMapping) {
        this.circuit = circuit;
        this.inputMapping = inputMapping;
    }

    public String getInputMapping() {
        return inputMapping;
    }

    public TreeCircuit getCircuit() {
        return circuit;
    }
}
