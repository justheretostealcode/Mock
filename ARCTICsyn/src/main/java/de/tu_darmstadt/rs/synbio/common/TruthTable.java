package de.tu_darmstadt.rs.synbio.common;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonValue;
import org.apache.commons.lang3.builder.EqualsBuilder;
import org.logicng.datastructures.Assignment;
import org.logicng.formulas.Formula;
import org.logicng.formulas.Literal;
import org.logicng.formulas.Variable;

import java.util.*;

public class TruthTable {

    private final BitSet truthTable = new BitSet();
    private final int firstUnusedBit;

    public TruthTable(Formula f) {
        int i = 0;
        for (Assignment assignment : getAllAssignments(f.variables())) {
            truthTable.set(i, f.evaluate(assignment));
            i++;
        }
        firstUnusedBit = i;
    }

    public TruthTable(Formula f, List<Assignment> assignments)  {
        int i = 0;
        for (Assignment assignment : assignments) {
            truthTable.set(i, f.evaluate(assignment));
            i++;
        }
        firstUnusedBit = i;
    }

    public BitSet getBitSet() {
        return truthTable;
    }

    public int getSupportSize() {

        int support = 0;

        while((1 << support) < firstUnusedBit) {
            support ++;
        }

        return support;
    }

    public int getLength() {
        return firstUnusedBit;
    }

    public static LinkedList<Assignment> getAllAssignments(SortedSet<Variable> variables) {

        LinkedList<Assignment> assignments = new LinkedList<>();
        NavigableSet<Variable> variableSet = new TreeSet<>(variables);
        int pos;

        // iterate over all possible assignments
        for (int i = 0; i < (1 << variables.size()); i ++) {

            LinkedList<Literal> literals = new LinkedList<>();
            pos = 0;

            // iterate over all variables per assignment in descending order (like in abc)
            Iterator<Variable> variableIterator = variableSet.descendingIterator();
            while (variableIterator.hasNext()) {

                Variable var = variableIterator.next();

                // add non-negated or negated variable to assigment
                boolean positiveVar = (i & (1 << pos)) != 0;
                literals.add(positiveVar ? var : var.negate());
                pos ++;

            }
            assignments.add(new Assignment(literals));
        }
        return assignments;
    }

    public String toString() {

        StringBuilder builder = new StringBuilder();

        for (int i = 0; i < firstUnusedBit; i ++) {
            builder.insert(0, truthTable.get(i) ? "1" : "0");
        }

        return builder.toString();
    }

    public boolean equalsLogically(TruthTable cmp) {
        if (this.firstUnusedBit != cmp.firstUnusedBit)
            return false;

        return this.truthTable.equals(cmp.truthTable);
    }

    @JsonValue
    public String toJson() {
        return toString();
    }

    @JsonCreator
    public TruthTable(String jsonValue) {

        int i = jsonValue.length() - 1;
        for (Character character : jsonValue.toCharArray()) {
            truthTable.set(i, (character == '1'));
            i--;
        }
        firstUnusedBit = jsonValue.length();
    }
}
