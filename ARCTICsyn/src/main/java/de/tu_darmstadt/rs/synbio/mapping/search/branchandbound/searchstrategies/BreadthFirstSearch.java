package de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.searchstrategies;

import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.QueueItem;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;

/**
 * Implements the search strategy BreadthFirstSearch (BFS)
 */
public class BreadthFirstSearch implements SearchStrategy{
    private final ArrayList<QueueItem> queue;

    public BreadthFirstSearch()   {
        this.queue = new ArrayList<>();
    }

    public BreadthFirstSearch(List<QueueItem> items)   {
        this.queue = new ArrayList<>(items);

        /*
         * The queueitems are sorted according to their assignment size.
         * By this, a very simple breadth first ordering is created.
         */
        this.queue.sort(new Comparator<QueueItem>() {
            @Override
            public int compare(QueueItem o1, QueueItem o2) {
                return o1.assignment.size() - o2.assignment.size();
            }
        });
    }

    @Override
    public void addInitialItemToQueue(Assignment assignment, double val) {
        QueueItem item = new QueueItem(assignment, val);
        queue.add(item);

        updateStatistics();
    }

    @Override
    public void addToQueue(QueueItem item) {
        queue.add(item);

        updateStatistics();
    }

    @Override
    public void addToQueue(List<QueueItem> items) {
        queue.addAll(items);
        updateStatistics();
    }

    @Override
    public QueueItem getNext() {
        if (queue.size() == 0)
            return null;

        QueueItem nextItem = queue.remove(0);

        updateStatistics();

        return nextItem;
    }

    @Override
    public List<QueueItem> getItems() {
        return queue;
    }

    @Override
    public int size() {
        return queue.size();
    }

    private int maximumNumberOfEntries = 0;
    private double averageNumberOfEntries = 0;
    private long divisor = 0;

    @Override
    public int getMaximumNumberOfQueueEntries() {
        return maximumNumberOfEntries;
    }
    @Override
    public double getAverageNumberOfQueueEntries() {
        return averageNumberOfEntries;
    }

    private void updateStatistics() {
        int size = queue.size();
        maximumNumberOfEntries = Math.max(size, maximumNumberOfEntries);

        divisor += 1;

        averageNumberOfEntries = averageNumberOfEntries + (size - averageNumberOfEntries) / divisor;
    }
}
