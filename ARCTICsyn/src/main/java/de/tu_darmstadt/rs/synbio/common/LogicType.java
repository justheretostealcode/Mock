package de.tu_darmstadt.rs.synbio.common;

import de.tu_darmstadt.rs.synbio.synthesis.util.ExpressionParser;
import org.logicng.formulas.Formula;

public enum LogicType {

    EMPTY,
    NOT,
    NOR2,
    NAND2,
    OR2;

    private static final Formula notExpr = ExpressionParser.parse("~x");
    private static final Formula nor2Expr = ExpressionParser.parse("~(x|y)");
    private static final Formula nand2Expr = ExpressionParser.parse("~(x&y)");
    private static final Formula or2Expr = ExpressionParser.parse("(x|y)");

    private static final String notStr = notExpr.toString();
    private static final String nor2Str = nor2Expr.toString();
    private static final String nand2Str = nand2Expr.toString();
    private static final String or2Str = or2Expr.toString();

    public Formula getExpression() {
        switch (this) {
            case NOT:
                return notExpr;
            case NOR2:
                return nor2Expr;
            case NAND2:
                return nand2Expr;
            default:
                return or2Expr;
        }
    }

    public String getExpressionString() {
        switch (this) {
            case NOT:
                return notStr;
            case NOR2:
                return nor2Str;
            case NAND2:
                return nand2Str;
            default:
                return or2Str;
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
