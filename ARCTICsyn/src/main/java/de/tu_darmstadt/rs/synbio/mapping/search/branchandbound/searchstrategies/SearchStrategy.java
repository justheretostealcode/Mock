package de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.searchstrategies;

import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.QueueItem;

import java.util.List;
import java.util.Queue;

/**
 * The interface for any SearchStrategy
 */
public interface SearchStrategy {

    /**
     * Adds an initial empty assignment with the corresponding value to the queue
     * @param val
     */
    public void addInitialItemToQueue(Assignment assignment, double val) ;

    /**
     * Adds a new QueueItem to the queue
     * @param item
     */
    public void addToQueue(QueueItem item);

    /**
     * Add multiple new QueueItems to the queue
     * @param items
     */
    public void addToQueue(List<QueueItem> items);

    /**
     * Returns the next QueueItem based on its ordering strategy
     * @return
     */
    public QueueItem getNext();

    /**
     * No guarantee on the ordering is provided.
     * @return
     */
    public List<QueueItem> getItems();

    /**
     * Gives rise to the number of elements within the queue
     * @return
     */
    public int size();

    /**
     * Gives rise to the maximum number of elements that has been in the queue
     * @return
     */
    public int getMaximumNumberOfQueueEntries();

    /**
     * Gives rise to the average number of elements in the queue
     * @return
     */
    public double getAverageNumberOfQueueEntries();
}
