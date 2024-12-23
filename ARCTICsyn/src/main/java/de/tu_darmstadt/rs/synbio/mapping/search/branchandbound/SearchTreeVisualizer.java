package de.tu_darmstadt.rs.synbio.mapping.search.branchandbound;

import de.tu_darmstadt.rs.synbio.common.circuit.Gate;
import de.tu_darmstadt.rs.synbio.common.library.GateRealization;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.mapping.MappingConfiguration;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * A class to visualize the branch and bound search tree
 */
public class SearchTreeVisualizer {
    private final String initialNodeID = "ROOT";
    private BufferedWriter writer;
    private int expansionIndex;
    private final Gate[] reversedLogicGates;
    private boolean bVisualize;


    /**
     * Creates a search tree visualizer. <br>
     * A file is created to store the created graph.
     * @param structureName         The name of the structure considered
     * @param mappingConfiguration  The mapping configuration
     * @param reversedLogicGates    The gates of the structure in reversed topological order
     * @param bVisualize            Whether to visualize or not
     */
    public SearchTreeVisualizer(String structureName, MappingConfiguration mappingConfiguration, Gate[] reversedLogicGates, boolean bVisualize) {

        this.reversedLogicGates = reversedLogicGates;
        this.bVisualize = bVisualize;

        if (bVisualize) {
            try {
                File file = new File("visualisation/");
                file.mkdirs();
                file = new File(file.getAbsolutePath(), structureName + "_" + System.currentTimeMillis() + ".dot");
                file.createNewFile();
                writer = new BufferedWriter(new FileWriter(file));

                writer.write("digraph \"" + structureName + "_" + System.currentTimeMillis() + "\" {\n");

                writer.write("ranksep=\"1.2 equally\";\n");
                writer.write(String.format("%s [label=\"root\"];\n", initialNodeID));
            } catch (IOException e) {
                this.bVisualize = false;
                e.printStackTrace();
            }
        }
        expansionIndex = 0;
    }

    public boolean getbVisualize()  {
        return this.bVisualize;
    }

    public void addLeafNode(Assignment assignment, double val, double currentBound, boolean skipParent) {
        if (!bVisualize)
            return;

        QueueItem item = QueueItem.getQueueItem(assignment, val);
        add(item, currentBound, skipParent);
    }

    public void addIntermediateNode(QueueItem item, double currentBound) {
        if (!bVisualize)
            return;

        add(item, currentBound, false);
    }

    public void add(QueueItem item, double currentBound, boolean skipParent) {
        if (!bVisualize)
            return;

        Assignment assignment = item.assignment;
        double value = item.val;
        int size = assignment.size();
        if (size == 0)
            return;

        Gate gate = reversedLogicGates[size - 1];
        GateRealization realization = assignment.get(gate);
        String sourceIdentifier = "";
        String targetIdentifier = initialNodeID;
        for (int iX = 0; iX < assignment.size(); iX++) {
            Gate g = reversedLogicGates[iX];
            if (!(skipParent && iX == (assignment.size() - 1)))
                sourceIdentifier = targetIdentifier;
            targetIdentifier += "__" + getNodeIdentifier(g, assignment.get(g));
        }

        String targetLabel = String.format("%s (%d)\n[%f | %f]", realization.getIdentifier(), expansionIndex, value, currentBound);
        try {
            writer.write(String.format("%s [label=\"%s\", shape=rect];\n", targetIdentifier, targetLabel));
            writer.write(String.format("%s -> %s;\n", sourceIdentifier, targetIdentifier));
        } catch (IOException e) {
            e.printStackTrace();
        }

        expansionIndex++;
    }

    private String getNodeIdentifier(Gate gate, GateRealization realization) {
        return gate.getIdentifier() + "_" + gate.getLogicType().toString() + "_" + realization.getIdentifier() + "_" + realization.getGroup();
    }

    /**
     * Completes the created graph and finalizes its export to the file.
     */
    public void finish() {
        if (!bVisualize)
            return;

        try {
            writer.write("}");
            writer.flush();
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
