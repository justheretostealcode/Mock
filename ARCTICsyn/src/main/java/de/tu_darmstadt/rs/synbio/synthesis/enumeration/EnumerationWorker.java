package de.tu_darmstadt.rs.synbio.synthesis.enumeration;

import de.tu_darmstadt.rs.synbio.common.TruthTable;
import org.apache.commons.lang3.StringUtils;
import org.logicng.formulas.Formula;
import org.logicng.formulas.FormulaFactory;
import org.logicng.io.parsers.ParserException;
import org.logicng.io.parsers.PropositionalParser;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

public class EnumerationWorker implements Callable<List<EnumerationResult>> {

    private final EnumeratorFast enumerator;
    private final TreeCircuit circuit;
    private final int feasibility;
    private final TruthTable targetTruthTable;

    private final List<EnumerationResult> results;

    public EnumerationWorker(EnumeratorFast enumerator, TreeCircuit circuit, int feasibility, TruthTable targetTruthTable) {
        this.enumerator = enumerator;
        this.circuit = circuit;
        this.feasibility = feasibility;
        this.targetTruthTable = targetTruthTable;
        this.results = new ArrayList<>();
    }

    @Override
    public List<EnumerationResult> call() throws Exception {

        FormulaFactory factory =  new FormulaFactory();
        PropositionalParser parser = new PropositionalParser(factory);

        int numUnboundInputs = circuit.getNumOpenInputs(true);

        for (int i = 0; i < (int) Math.pow(feasibility, numUnboundInputs); i++) {

            String inputMapping = Integer.toString(i, feasibility);
            inputMapping = StringUtils.leftPad(inputMapping, numUnboundInputs, '0');

            // test if mapping contains all primary inputs
            boolean valid = true;
            for (int j = 0; j < feasibility; j++) {
                if (inputMapping.indexOf(Character.forDigit(j, feasibility)) == -1) {
                    valid = false;
                    break;
                }
            }
            if (!valid)
                continue;

            // evaluate primitive circuit with input mapping
            String circuitExpression = enumerator.evaluateTreeCircuit(circuit, inputMapping);

            if (circuitExpression == null)
                continue;

            Formula f = null;
            try {
                f = parser.parse(circuitExpression);
            } catch(ParserException e) {
                e.printStackTrace();
            }

            if (f == null)
                continue;

            TruthTable circuitTT = new TruthTable(f);

            // match
            if(targetTruthTable == null || targetTruthTable.equals(circuitTT)) {
                results.add(new EnumerationResult(circuit, inputMapping));
            }
        }
        return results;
    }
}
