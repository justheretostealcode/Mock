package de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.searchstrategies;

import de.tu_darmstadt.rs.synbio.mapping.search.branchandbound.QueueItem;

import java.util.ArrayList;
import java.util.List;

/**
 * This class provides means for efficiently inserting an element into an ordered list.
 */
public class SearchStrategiesUtil {
    @Deprecated
    public static List<QueueItem> binaryInsert(List<QueueItem> list, QueueItem elem) {
        int l = 0;
        int h = list.size() - 1;
        int m = 0;

        double x = elem.val;
        double val = 0;

        if (x <= list.get(l).val)
            list.add(0, elem);
        else if (x >= list.get(h).val)
            list.add(h + 1, elem);
        else {
            while (l <= h) {
                m = l + (h - l) / 2;
                val = list.get(m).val;

                if (l == h) {
                    if (list.get(l).val < x)
                        list.add(l + 1, elem);
                    else
                        list.add(l, elem);
                    break;
                } else if (val == x) {
                    list.add(m + 1, elem);
                    break;
                } else if (val > x)
                    h = m - 1;
                else
                    l = m + 1;
            }
        }

        return list;
    }


    /**
     * Gives rise to the index at which an element with value x can be inserted into the ordered list,
     * whereby the order is preserved.
     *
     * @param list      An ordered list
     * @param x         The value to determine its position
     * @param ascending True if the list is ordered ascending otherwise false
     * @return The index at which value x could be inserted in list, while preserving the order
     */
    public static int binarySearch(List<QueueItem> list, double x, boolean ascending) {
        int l = 0;
        int h = list.size() - 1;
        int m = 0;

        if (h == -1)
            return 0;

        double val = 0;
        double firstVal = list.get(l).val;
        double lastVal = list.get(h).val;
        if (ascending && x <= firstVal
                || !ascending && x >= firstVal)
            return 0;
        else if (ascending && x >= lastVal
                || !ascending && x <= lastVal)
            return h + 1;
        else {
            while (l <= h) {
                m = l + (h - l) / 2;
                val = list.get(m).val;
                if (val == x)
                    return m;
                else if (ascending && val > x
                        || !ascending && val < x)
                    h = m - 1;      // m has already been considered
                else
                    l = m + 1;      // m has already been considered
            }
        }

        return l;
    }


    public static boolean isSorted(List<QueueItem> list, boolean ascending) {
        int factor = (ascending) ? 1 : -1;
        double diff = 0;
        for (int iX = 0; iX < list.size() - 1; iX++) {
            diff = list.get(iX + 1).val - list.get(iX).val;
            if (!(factor * diff >= 0))
                return false;
        }
        return true;
    }
}
