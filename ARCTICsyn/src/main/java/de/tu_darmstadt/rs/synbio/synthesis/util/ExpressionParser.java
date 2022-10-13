package de.tu_darmstadt.rs.synbio.synthesis.util;

import org.logicng.formulas.Formula;
import org.logicng.formulas.FormulaFactory;
import org.logicng.io.parsers.ParserException;
import org.logicng.io.parsers.PropositionalParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class ExpressionParser {

    private static final Logger logger = LoggerFactory.getLogger(ExpressionParser.class);
    private static final FormulaFactory factory = new FormulaFactory();
    private static final PropositionalParser parser = new PropositionalParser(factory);

    public static Formula parse(String expression) {

        Formula f = null;
        try {
            f = parser.parse(expression);
        } catch(ParserException e) {
            e.printStackTrace();
        }

        return f;
    }
}
