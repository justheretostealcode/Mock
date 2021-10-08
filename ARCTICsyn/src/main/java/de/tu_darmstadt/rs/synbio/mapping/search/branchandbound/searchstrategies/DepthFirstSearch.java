package de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.searchstrategies;

import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.QueueItem;

import java.lang.reflect.Array;
import java.util.*;

/**
 * This class implements the search strategy DepthFirstSearch (DFS)
 */
public class DepthFirstSearch implements SearchStrategy {

    private final ArrayList<QueueItem> queue;

    public DepthFirstSearch() {
        this.queue = new ArrayList<>();
    }

    public DepthFirstSearch(List<QueueItem> items) {
        this.queue = new ArrayList<>(items);

        /*
         * The queueitems are sorted descending according to their assignment size.
         * By this, a very simple depthfirst ordering is created from the view of the first to the last item.
         */
        this.queue.sort(new Comparator<QueueItem>() {
            @Override
            public int compare(QueueItem o1, QueueItem o2) {
                return o2.assignment.size() - o1.assignment.size();
            }
        });
    }

    @Override
    public void addInitialItemToQueue(double val) {
        QueueItem item = new QueueItem(new Assignment(), val);
        queue.add(0, item);

        updateStatistics();
    }

    @Override
    public void addToQueue(QueueItem item) {
        queue.add(0, item);

        updateStatistics();
    }

    @Override
    public void addToQueue(List<QueueItem> items) {
        for (QueueItem item : items)    {
            queue.add(0, item);
        }
        updateStatistics();
    }

    @Override
    public QueueItem getNext() {
        int queueSize = queue.size();
        if (queueSize == 0)
            return null;

        QueueItem nextItem =  queue.remove(0);

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
