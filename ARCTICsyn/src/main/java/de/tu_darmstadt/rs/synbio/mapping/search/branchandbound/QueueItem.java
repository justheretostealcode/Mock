package de.tu_darmstadt.rs.synbio.mapping.search.branchandbound;

import de.tu_darmstadt.rs.synbio.mapping.Assignment;

import java.util.ArrayList;
import java.util.List;

/**
 * This class represents an item for the queues realized by each particular search strategy
 */
public class QueueItem {
    /**
     * The represented assignment
     */
    public Assignment assignment;

    /**
     * The value assigned to this assignment. <br>
     * Whether it is its own or its parents score depends on whether eager or lazy branch and bound is used.
     */
    public double val;

    // A list for buffering items and thus trying to reduce number of object creations
    private static final ArrayList<QueueItem> items = new ArrayList<>();


    /**
     * A method for obtaining a queue item of the asked values. <br>
     * If there exists at least one non used item, this one is used while otherwise a new one is created. <br>
     * This enables the recycling of already created QueueItem Objects
     *
     * @param assignment
     * @param val
     * @return The queue item featuring the provided information
     */
    public static QueueItem getQueueItem(Assignment assignment, double val) {
        if (items.size() == 0) {
            return new QueueItem(assignment, val);
        }

        QueueItem item1 = items.remove(0);
        item1.assignment = assignment;
        item1.val = val;
        return item1;
    }

    /**
     * A method for obtaining a queue item of the asked values. <br>
     *
     * @param assignment
     * @param val
     * @return The queue item featuring the provided information
     */
    public QueueItem(Assignment assignment, double val) {
        this.assignment = assignment;
        this.val = val;
    }

    /**
     * This method needs to be called in order to add this item to the list of available items for reuse.
     */
    public void free() {
        items.add(this);
        this.assignment = null;
        this.val = Double.NEGATIVE_INFINITY;
    }
}
