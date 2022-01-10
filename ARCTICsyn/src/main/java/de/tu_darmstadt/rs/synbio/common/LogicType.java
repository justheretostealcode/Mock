package de.tu_darmstadt.rs.synbio.common;

import de.tu_darmstadt.rs.synbio.synthesis.util.ExpressionParser;
import org.logicng.formulas.Formula;

public enum LogicType {

    EMPTY,
    NOT,
    NOR2,
    OR2;

    private static final Formula notExpr = ExpressionParser.parse("~x");
    private static final Formula nor2Expr = ExpressionParser.parse("~(x|y)");
    private static final Formula or2Expr = ExpressionParser.parse("(x|y)");

    public Formula getExpression() {
        switch (this) {
            case NOT:
                return notExpr;
            case NOR2:
                return nor2Expr;
            default:
                return or2Expr;
        }
    }

    public int getNumInputs() {
        switch (this) {
            case EMPTY: return 0;
            case NOT: return 1;
            default: return 2;
        }
    }

    public int getWeight() {
        if (this.equals(EMPTY)/* || this.equals(OR2)*/)
            return 0;
        else
            return 1;
    }
}
