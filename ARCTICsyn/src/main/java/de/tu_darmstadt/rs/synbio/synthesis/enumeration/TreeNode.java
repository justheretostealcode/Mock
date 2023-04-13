package de.tu_darmstadt.rs.synbio.synthesis.enumeration;

import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;

public class TreeNode {

    public LogicType type;
    public TreeNode child0;
    public TreeNode child1;

    public int equivalenceClass;

    public String expression;

    public TreeNode(LogicType type) {
        this.type = type;
    }

    public TreeNode(LogicType type, TreeNode child0, TreeNode child1) {
        this.type = type;
        this.child0 = child0;
        this.child1 = child1;
    }

    TreeNode copy() {

        TreeNode child0 = null;
        TreeNode child1 = null;

        if (this.child0 != null) {
            child0 = this.child0.copy();
        }
        if (this.child1 != null) {
            child1 = this.child1.copy();
        }

        return new TreeNode(type, child0, child1);
    }
}
