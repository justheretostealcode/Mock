package de.tu_darmstadt.rs.synbio.common;

import de.tu_darmstadt.rs.synbio.synthesis.util.ExpressionParser;
import org.logicng.formulas.Formula;

public enum LogicType {

    EMPTY,
    INPUT,
    OUTPUT_BUFFER,
    OUTPUT_OR2,
    NOT,
    NOR2;

    private static final Formula inputExpr = ExpressionParser.parse("x");
    private static final Formula outputBufExpr = ExpressionParser.parse("x");
    private static final Formula outputOr2Expr = ExpressionParser.parse("(x|y)");
    private static final Formula notExpr = ExpressionParser.parse("~x");
    private static final Formula nor2Expr = ExpressionParser.parse("~(x|y)");

    private static final String inputStr = inputExpr.toString();
    private static final String outputBufStr = outputBufExpr.toString();
    private static final String outputOr2Str = outputOr2Expr.toString();
    private static final String notStr = notExpr.toString();
    private static final String nor2Str = nor2Expr.toString();

    public Formula getExpression() {
        switch (this) {
            case INPUT:
                return inputExpr;
            case OUTPUT_BUFFER:
                return outputBufExpr;
            case OUTPUT_OR2:
                return outputOr2Expr;
            case NOT:
                return notExpr;
            case NOR2:
                return nor2Expr;
            default:
                return outputBufExpr;
        }
    }

    public String getExpressionString() {
        switch (this) {
            case INPUT:
                return inputStr;
            case OUTPUT_BUFFER:
                return outputBufStr;
            case OUTPUT_OR2:
                return outputOr2Str;
            case NOT:
                return notStr;
            case NOR2:
                return nor2Str;
            default:
                return outputBufStr;
        }
    }

    public int getNumInputs() {
        switch (this) {
            case EMPTY:
            case INPUT:
                return 0;
            case OUTPUT_BUFFER:
            case NOT: return 1;
            case OUTPUT_OR2:
            default: return 2;
        }
    }

    public int getWeight() {
        if (this.equals(EMPTY) || this.equals(INPUT) || this.equals(OUTPUT_BUFFER))
            return 0;
        else if (this.equals(OUTPUT_OR2))
            return 0;
        else
            return 1;
    }
}
