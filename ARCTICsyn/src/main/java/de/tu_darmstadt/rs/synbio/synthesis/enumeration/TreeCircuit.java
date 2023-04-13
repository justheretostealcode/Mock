package de.tu_darmstadt.rs.synbio.synthesis.enumeration;

import de.tu_darmstadt.rs.synbio.common.LogicType;

import java.util.ArrayList;
import java.util.List;

public class TreeCircuit {

    private final TreeNode rootNode;
    public TreeCircuit(LogicType rootType) {
        this.rootNode = new TreeNode(rootType);
    }

    public TreeCircuit(TreeCircuit copyCircuit) {
        this.rootNode = copyCircuit.rootNode.copy();
    }

    public boolean addLevel(List<LogicType> gates) {

        List<TreeNode> unconnected = unconnectedNodes();

        int i = 0;

        for (TreeNode node : unconnected) {

            if (node.type == LogicType.EMPTY)
                continue;

            if (node.type == LogicType.NOT) {

                /* avoid subsequent inverters */
                if (gates.get(i) == LogicType.NOT)
                    return false;

                node.child0 = new TreeNode(gates.get(i));
                i++;
            } else {

                if (node.child0 == null) {
                    node.child0 = new TreeNode(gates.get(i));
                    i++;
                }

                if (node.child1 == null) {
                    node.child1 = new TreeNode(gates.get(i));
                    i++;
                }
            }
        }

        return true;
    }

    public int getNumOpenInputs(boolean countEmpty) {

        List<TreeNode> unconnected = unconnectedNodes();

        int numOpenInputs = 0;

        for (TreeNode node : unconnected) {

            if (node.type == LogicType.EMPTY) {
                if (countEmpty)
                    numOpenInputs++;
            } else if (node.type == LogicType.NOT) {
                numOpenInputs++;
            } else {
                if (node.child0 == null)
                    numOpenInputs ++;
                if (node.child1 == null)
                    numOpenInputs++;
            }
        }

        return numOpenInputs;
    }

    public int getWeight() {

        List<TreeNode> nodes = new ArrayList<>();
        traversePreOrder(rootNode, nodes);

        int weight = 0;

        for (TreeNode node : nodes) {
            weight += node.type.getWeight();
        }

        return weight;
    }

    public List<TreeNode> serializePreOrder() {
        List<TreeNode> nodes = new ArrayList<>();
        traversePreOrder(rootNode, nodes);
        return nodes;
    }

    public List<TreeNode> serializePostOrder() {
        List<TreeNode> nodes = new ArrayList<>();
        traversePostOrder(rootNode, nodes);
        return nodes;
    }

    public TreeNode getRootNode() {
        return rootNode;
    }

    private List<TreeNode> unconnectedNodes() {

        List<TreeNode> nodes = new ArrayList<>();
        List<TreeNode> unconnected = new ArrayList<>();

        traversePreOrder(rootNode, nodes);

        for(TreeNode node : nodes) {

            if (node.type == LogicType.EMPTY)
                unconnected.add(node);
            else if (node.type == LogicType.NOT && node.child0 == null)
                unconnected.add(node);
            else if (node.type != LogicType.NOT && (node.child0 == null || node.child1 == null))
                unconnected.add(node);
        }

        return unconnected;
    }

    private void traversePreOrder(TreeNode node, List<TreeNode> nodes) {

        if (node == null)
            return;

        nodes.add(node);

        traversePreOrder(node.child0, nodes);
        traversePreOrder(node.child1, nodes);
    }

    private void traversePostOrder(TreeNode node, List<TreeNode> nodes) {

        if (node == null)
            return;

        traversePostOrder(node.child0, nodes);
        traversePostOrder(node.child1, nodes);

        nodes.add(node);
    }
}
