package de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.searchstrategies;

import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.QueueItem;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

/**
 * Class implementing BestFirstSearch (BeFS) as search strategy.
 */
public class BestFirstSearch implements SearchStrategy {
    private final ArrayList<QueueItem> queue;

    /**
     * The queue items are sorted descending to their value. <br>
     * By this, the first item is the one with the highest score.
     */
    private static final Comparator<QueueItem> comparator = (o1, o2) -> (int) Math.signum(o2.val - o1.val);

    public BestFirstSearch() {
        this.queue = new ArrayList<>();
    }

    public BestFirstSearch(List<QueueItem> items) {
        this.queue = new ArrayList<>(items);


        this.queue.sort(comparator);
    }

    @Override
    public void addInitialItemToQueue(Assignment assignment, double val) {
        QueueItem item = new QueueItem(assignment, val);
        queue.add(item);

        updateStatistics();
    }

    @Override
    public void addToQueue(QueueItem item) {
        if (item == null)
            return;

        int index = SearchStrategiesUtil.binarySearch(queue, item.val, false);
        queue.add(index, item);

        updateStatistics();
    }

    @Override
    public void addToQueue(List<QueueItem> items) {
        if (items == null)
            return;

        HashMap<Double, Integer> indexCache = new HashMap<>(items.size());
        for (QueueItem item : items) {
            double val = item.val;
            Integer index = indexCache.get(val);
            if (index == null) {
                index = SearchStrategiesUtil.binarySearch(queue, item.val, false);
                indexCache.put(val, index);
            }
            queue.add(indexCache.get(val), item);

            for (Double key : indexCache.keySet()) {
                if (key <= val)
                    indexCache.put(key, indexCache.get(key) + 1); // Increment as an node has been added in front
            }
        }
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


    private void isQueueOrdered(List<QueueItem> list) {
        if (!SearchStrategiesUtil.isSorted(list, false))
            throw new Error("Queue is not sorted");
    }
}
