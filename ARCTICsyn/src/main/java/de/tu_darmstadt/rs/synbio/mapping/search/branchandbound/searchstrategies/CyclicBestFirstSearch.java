package de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.searchstrategies;

import de.tu_darmstadt.rs.synbio.mapping.Assignment;
import de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.QueueItem;

import java.util.*;

/**
 * Implementation of CyclicBestFirstSearch (CBeFS) as search strategy.
 */
public class CyclicBestFirstSearch implements SearchStrategy {
    private final ArrayList<QueueItem>[] queues;

    private int pointer;
    private final int numberOfGates;
    /**
     * The queueitems are sorted descending to their value. <br>
     * By this, the first item is the one with the highest score.
     */
    private static final Comparator<QueueItem> comparator = new Comparator<QueueItem>() {
        @Override
        public int compare(QueueItem o1, QueueItem o2) {
            return (int) Math.signum(o2.val - o1.val);
        }
    };


    public CyclicBestFirstSearch(int numberOfGates) {
        this.numberOfGates = numberOfGates;
        this.pointer = 0;
        this.queues = new ArrayList[numberOfGates];
        for (int iX = 0; iX < this.queues.length; iX++) {
            this.queues[iX] = new ArrayList<>();
        }
    }

    public CyclicBestFirstSearch(int numberOfGates, List<QueueItem> items) {
        this.queues = new ArrayList[numberOfGates];
        for (int iX = 0; iX < this.queues.length; iX++) {
            this.queues[iX] = new ArrayList<>();
        }

        items.parallelStream().forEach(item -> this.queues[item.assignment.size() - 1].add(item));

        sortAll();
        this.pointer = 0;
        this.numberOfGates = numberOfGates;
    }

    @Override
    public void addInitialItemToQueue(Assignment assignment, double val) {
        QueueItem item = new QueueItem(assignment, val);
        queues[0].add(item);

        updateStatistics();
    }

    @Override
    public void addToQueue(QueueItem item) {
        int index = item.assignment.size() - 1;
        List<QueueItem> queue = queues[index];

        int indexToAdd = SearchStrategiesUtil.binarySearch(queue, item.val, false);
        queue.add(indexToAdd, item);

        updateStatistics();
    }

    @Override
    public void addToQueue(List<QueueItem> items) {
        int queueIndex;
        List<QueueItem> queue;
        if (items != null && items.size() > 0) {
            queueIndex = items.get(0).assignment.size() - 1;
            queue = queues[queueIndex];


            HashMap<Double, Integer> indexes = new HashMap<>(items.size());
            for (QueueItem item : items) {
                double val = item.val;
                Integer index = null;
                boolean useCache = item.assignment.size() - 1 == queueIndex;
                if (useCache) {
                    index = indexes.get(val);
                }

                if (index == null) {
                    index = SearchStrategiesUtil.binarySearch(queue, item.val, false);

                    if (useCache) {
                        indexes.put(val, index);
                    }
                }
                queue.add(indexes.get(val), item);

                for (Double key : indexes.keySet()) {
                    if (key <= val)
                        indexes.put(key, indexes.get(key) + 1); // Increment as an node has been added in front
                }
            }
        }
        updateStatistics();
    }


    @Override
    public QueueItem getNext() {
        QueueItem item = null;
        int currentVal = pointer;

        do {
            if (queues[pointer].size() != 0) {
                item = queues[pointer].remove(0);
            }

            pointer = (pointer + 1) % numberOfGates;
        } while (item == null && currentVal != pointer);

        updateStatistics();

        return item;
    }

    @Override
    public List<QueueItem> getItems() {
        List<QueueItem> items = new ArrayList<>();

        for (ArrayList<QueueItem> queue : queues) {
            items.addAll(queue);
        }

        return items;
    }

    private void sortAll() {
        for (ArrayList<QueueItem> queue : this.queues) {
            queue.sort(comparator);
        }
    }

    @Override
    public int size() {
        int size = 0;
        for (ArrayList<QueueItem> queue : queues) {
            size += queue.size();
        }
        return size;
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
        int size = size();
        maximumNumberOfEntries = Math.max(size, maximumNumberOfEntries);

        divisor += 1;

        averageNumberOfEntries = averageNumberOfEntries + (size - averageNumberOfEntries) / divisor;
    }

    private void isQueueOrdered(List<QueueItem> list) {
        if (!SearchStrategiesUtil.isSorted(list, false))
            throw new Error("Queue is not sorted");
    }
}
