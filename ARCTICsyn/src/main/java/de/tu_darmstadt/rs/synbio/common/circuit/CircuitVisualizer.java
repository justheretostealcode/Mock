package de.tu_darmstadt.rs.synbio.common.circuit;

import de.tu_darmstadt.rs.synbio.common.LogicType;
import de.tu_darmstadt.rs.synbio.common.library.GateRealization;
import de.tu_darmstadt.rs.synbio.mapping.Assignment;

import java.awt.*;

import org.jfree.svg.SVGGraphics2D;
import org.logicng.util.Pair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileWriter;
import java.io.Writer;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

public class CircuitVisualizer {

    private static final Logger logger = LoggerFactory.getLogger(CircuitVisualizer.class);

    private final Circuit structure;
    private final Assignment assignment;

    public CircuitVisualizer(Circuit structure, Assignment assignment) {
        this.structure = structure;
        this.assignment = assignment;
    }

    static final int tileWidth = 150;
    static final int tileHeight = 100;

    public void visualize() {

        /* get gate levels and connections */

        List<Pair<Gate, Gate>> connections = new ArrayList<>();
        Map<Gate, Integer> gateLevels = new HashMap<>();

        Gate outputGate = structure.getOutputGate();
        gateLevels.put(outputGate, 0);
        Set<Wire> incomingEdges = new HashSet<>(structure.incomingEdgesOf(outputGate));

        Set<Gate> mappedGates = new HashSet<>();
        mappedGates.add(outputGate);

        int level = 1;

        while (!incomingEdges.isEmpty()) {

            Set<Gate> sourceGates = incomingEdges.stream().map(structure::getEdgeSource).collect(Collectors.toSet());

            Set<Gate> mappedGatesOnLevel = new HashSet<>();

            int finalLevel = level;

            for (Gate gate : sourceGates) {
                if (mappedGates.containsAll(structure.getDescendants(gate))) {
                    gateLevels.putIfAbsent(gate, finalLevel);
                    mappedGatesOnLevel.add(gate);
                }
            }

            mappedGates.addAll(mappedGatesOnLevel);

            incomingEdges.forEach(e -> connections.add(new Pair<>(structure.getEdgeSource(e), structure.getEdgeTarget(e))));

            incomingEdges.clear();
            incomingEdges.addAll(sourceGates.stream().map(structure::incomingEdgesOf).flatMap(Set::stream).collect(Collectors.toSet()));

            level++;
        }

        int depth = level;

        /* move all inputs to last level */

        gateLevels.entrySet().stream()
                .filter(e -> e.getKey().getLogicType().equals(LogicType.INPUT))
                .forEach(e -> e.setValue(depth - 1));

        /* gate rows */

        Map<Gate, Integer> gateYValues = getYValues(gateLevels, connections, depth);

        /* optimize y values of gates for straighter connections */

          int height = gateYValues.entrySet().stream().max(Map.Entry.comparingByValue()).get().getValue() / tileHeight;
//        Loader.loadNativeLibraries();
//        MPSolver solver = MPSolver.createSolver("GLOP");
//
//        Map<Gate, MPVariable> gateYValue = new HashMap<>();
//
//        /* add variables */
//        for (Gate gate : structure.vertexSet()) {
//            gateYValue.put(gate, solver.makeNumVar(0.0, height * tileHeight, gate.getIdentifier()));
//        }
//
//        /* add constraints */
//        for (level = 0; level < depth; level++) {
//
//            int finalLevel = level;
//            List<Gate> gatesOnLevel = gateLevels.entrySet().stream()
//                    .filter(e -> e.getValue() == finalLevel)
//                    .map(Map.Entry::getKey)
//                    .collect(Collectors.toList());
//
//            if (gatesOnLevel.size() == 1)
//                continue;
//
//            gatesOnLevel.sort(Comparator.comparing(gateRows::get));
//
//            for (int i = 0; i < gatesOnLevel.size() - 1; i++) {
//                MPConstraint constraint = solver.makeConstraint(Double.NEGATIVE_INFINITY, -1 * tileHeight, "level_" + finalLevel + "_row_" + i);
//                constraint.setCoefficient(gateYValue.get(gatesOnLevel.get(i)), 1);
//                constraint.setCoefficient(gateYValue.get(gatesOnLevel.get(i + 1)), -1);
//            }
//        }
//
//        /* objective */
//
//        Map<Gate, Double> coefficients = new HashMap<>();
//
//        for (Pair<Gate, Gate> pair : connections) {
//
//            coefficients.putIfAbsent(pair.first(), 0.0);
//            coefficients.putIfAbsent(pair.second(), 0.0);
//
//            coefficients.put(pair.first(), coefficients.get(pair.first()) + 1.0);
//            coefficients.put(pair.second(), coefficients.get(pair.second()) - 1.0);
//        }
//
//        MPObjective objective = solver.objective();
//
//        for (Map.Entry<Gate, Double> entry : coefficients.entrySet()) {
//            objective.setCoefficient(gateYValue.get(entry.getKey()), entry.getValue());
//        }
//
//        objective.setMinimization();
//
//        final MPSolver.ResultStatus resultStatus = solver.solve();
//
//        for (MPVariable var : gateYValue.values()) {
//            logger.info(var.name() + ": " + var.solutionValue());
//        }

        /* initialize image */

        SVGGraphics2D image = new SVGGraphics2D(tileWidth * depth, tileHeight * (height + 1));
        image.setColor(Color.WHITE);
        image.fillRect(0, 0, tileWidth * depth, tileHeight * (height + 1));

        image.setStroke(new BasicStroke(2, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        image.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 15));

        /* draw connections */

        Map<Gate, Map<Gate, List<Line2D>>> connectionLines = new HashMap<>();

        Map<Pair<Gate, Gate>, Integer> portAssignment = new HashMap<>();

        for (Pair<Gate, Gate> connection : connections) {

            Gate sourceGate1 = connection.first();
            Gate targetGate = connection.second();

            Optional<Pair<Gate, Gate>> secondConnection = connections.stream()
                    .filter(c -> c.second().equals(targetGate))
                    .findAny();

            if (secondConnection.isPresent()) {
                if (gateYValues.get(sourceGate1) > gateYValues.get(secondConnection.get().first())) {
                    portAssignment.put(connection, 1);
                    portAssignment.put(secondConnection.get(), 0);
                } else {
                    portAssignment.put(connection, 0);
                    portAssignment.put(secondConnection.get(), 1);
                }
            } else {
                portAssignment.put(connection, 0);
            }
        }

        Map<Pair<Gate, Gate>, Integer> beneficialMidpoints = new HashMap<>();

        for (Pair<Gate, Gate> connection : connections) {

            Gate sourceGate = connection.first();
            Gate targetGate = connection.second();

            double sourceY = gateYValues.get(sourceGate) + getOutputOffset(sourceGate.getLogicType()).getY();
            double sourceX = (depth - 1 - gateLevels.get(sourceGate)) * tileWidth + getOutputOffset(sourceGate.getLogicType()).getX();

            double targetY = gateYValues.get(targetGate) + getInputOffset(targetGate.getLogicType(), portAssignment.get(connection)).getY();
            double targetX = (depth - 1 - gateLevels.get(targetGate)) * tileWidth + getInputOffset(targetGate.getLogicType(), portAssignment.get(connection)).getX();

            connectionLines.putIfAbsent(sourceGate, new HashMap<>());
            connectionLines.get(sourceGate).putIfAbsent(targetGate, new ArrayList<>());


            /*Integer sourceYFreeUntil = gateLevels.get(targetGate);
            Integer targetYFreeUntil = gateLevels.get(sourceGate);

            if (gateLevels.get(sourceGate) - gateLevels.get(targetGate) > 1) {
                for (int l = gateLevels.get(sourceGate) - 1; l > gateLevels.get(targetGate); l --) {

                    int finalL = l;
                    List<Gate> gatesOnLevel = gateLevels.entrySet().stream()
                            .filter(e -> e.getValue() == finalL)
                            .map(Map.Entry::getKey)
                            .collect(Collectors.toList());

                    for (Gate gate : gatesOnLevel) {
                        int y = gateYValues.get(gate);

                        if (sourceYFreeUntil == gateLevels.get(targetGate) && (sourceY >= y && sourceY <= y + tileHeight))
                            sourceYFreeUntil = l + 1;

                        if (targetY >= y && targetY <= y + tileHeight)
                            targetYFreeUntil = l - 1;
                    }
                }
            }

            beneficialMidpoints.put(connection, depth - 1 - targetYFreeUntil * tileWidth);*/

            if (sourceY == targetY) {
                Line2D.Double directLine = new Line2D.Double(sourceX, sourceY, targetX, targetY);
                connectionLines.get(sourceGate).get(targetGate).add(directLine);
            } else {

                double midpointX = sourceX;//(sourceX + targetX) / 2;

                Line2D.Double startLine = new Line2D.Double(sourceX, sourceY, midpointX, sourceY);
                Line2D.Double midLine = new Line2D.Double(midpointX, sourceY, midpointX, targetY);
                Line2D.Double endLine = new Line2D.Double(midpointX, targetY, targetX, targetY);

                connectionLines.get(sourceGate).get(targetGate).add(startLine);
                connectionLines.get(sourceGate).get(targetGate).add(midLine);
                connectionLines.get(sourceGate).get(targetGate).add(endLine);
            }

            image.setPaint(new Color(0, 0, 0));
            connectionLines.get(sourceGate).get(targetGate).forEach(image::draw);
        }

        /* draw gates */

        for (Gate gate : structure.vertexSet()) {

            switch (gate.getLogicType()) {
                case INPUT:
                    InputGate.add(image, (depth - gateLevels.get(gate) - 1) * tileWidth, gateYValues.get(gate), gate.getIdentifier());
                    break;
                case NOT:
                    NotGate.add(image, (depth - gateLevels.get(gate) - 1) * tileWidth, gateYValues.get(gate), assignment.get(gate));
                    break;
                case NOR2:
                    OrNorGate.add(image, (depth - gateLevels.get(gate) - 1) * tileWidth, gateYValues.get(gate), assignment.get(gate), true);
                    break;
                case OUTPUT_OR2:
                    OrNorGate.add(image, (depth - gateLevels.get(gate) - 1) * tileWidth, gateYValues.get(gate), assignment.get(gate), false);
                    break;
                case OUTPUT_BUFFER:
                    break;
                default:
                    logger.error("Unsupported gate type for visualization: {}", gate.getLogicType());
            }
        }


        /* svg output */

        try {
            Writer writer = new FileWriter(new File("test.svg"));
            writer.write(image.getSVGDocument());
            writer.close();
        } catch (Exception e) {
            logger.error(e.getMessage());
        }


        int x = 0;
    }

    private static Point2D getOutputOffset(LogicType type) {

        switch (type) {
            case INPUT:
                return InputGate.outputOffset;
            case OUTPUT_OR2:
            case NOR2:
                return OrNorGate.outputOffset;
            case NOT:
                return NotGate.outputOffset;
            default:
                logger.error("Gate type has no output offset: {}", type);
                return null;
        }
    }

    private static Point2D getInputOffset(LogicType type, int port) {

        switch (type) {
            case OUTPUT_OR2:
            case NOR2:
                if (port == 0)
                    return OrNorGate.input0Offset;
                else
                    return OrNorGate.input1Offset;
            case NOT:
                return NotGate.inputOffset;
            case OUTPUT_BUFFER:
                return OutputBuffer.inputOffset;
            case INPUT:
            default:
                logger.error("Gate type has no input offset: {}", type);
                return null;
        }
    }

    private static class NotGate {

        static Point2D inputOffset = new Point2D.Double(0.25 * tileWidth, 0.5 * tileHeight);
        static Point2D outputOffset = new Point2D.Double(0.75 * tileWidth, 0.5 * tileHeight);

        static void add(Graphics2D image, double tileX, double tileY, GateRealization gate) {

            //TODO: get color

            image.setPaint(new Color(245, 189, 189));

            double a = 40;
            double height = a * 0.8;// * 0.721124785;
            double d = 6;

            double x = tileX + (tileWidth - (height + d)) / 2;
            double y = tileY + (tileHeight - a) / 2;

            Path2D.Double path = new Path2D.Double();
            path.moveTo(x, y);
            path.lineTo(x + height, y + a / 2);
            path.lineTo(x, y + a);
            path.closePath();

            Ellipse2D.Double circle = new Ellipse2D.Double(x + height, y + a / 2 - d / 2, d, d);

            Line2D connectors = new Line2D.Double(tileX + inputOffset.getX(), tileY + inputOffset.getY(),
                    tileX + outputOffset.getX(), tileY + outputOffset.getY());

            image.setPaint(new Color(0, 0, 0));
            image.draw(connectors);

            image.setPaint(new Color(245, 189, 189));
            image.fill(path);
            image.fill(circle);

            image.setPaint(new Color(0, 0, 0));
            image.draw(path);
            image.draw(circle);
        }
    }

    private static class OrNorGate {

        static Point2D input0Offset = new Point2D.Double(0.1 * tileWidth, 0.33 * tileHeight);
        static Point2D input1Offset = new Point2D.Double(0.1 * tileWidth, 0.66 * tileHeight);
        static Point2D outputOffset = new Point2D.Double(0.9 * tileWidth, 0.5 * tileHeight);

        static void add(Graphics2D image, double tileX, double tileY, GateRealization gate, boolean isNor) {

            //TODO: get color

            image.setPaint(new Color(245, 189, 189));

            double bodyLength = 80;
            double bodyHeight = 60;
            double d = 6;

            double x = tileX + (tileWidth - (bodyLength + d)) / 2;
            double y = tileY + (tileHeight - bodyHeight) / 2;

            Path2D.Double path = new Path2D.Double();
            path.moveTo(x, y);
            path.lineTo(x + 0.5 * bodyLength, y);
            path.curveTo(x + 0.75 * bodyLength, y, x + 0.9 * bodyLength, y + 0.25 * bodyHeight, x + bodyLength, y + bodyHeight / 2);
            path.curveTo(x + 0.9 * bodyLength, y + 0.75 * bodyHeight, x + 0.75 * bodyLength, y + bodyHeight, x + 0.5 * bodyLength, y + bodyHeight);
            path.lineTo(x + 0.5 * bodyLength, y + bodyHeight);
            path.lineTo(x, y + bodyHeight);
            path.curveTo(x + 0.15 * bodyLength, y + 0.6 * bodyHeight, x + 0.15 * bodyLength, y + 0.4 * bodyHeight, x, y);

            Ellipse2D.Double circle = new Ellipse2D.Double(x + bodyLength, y + bodyHeight / 2 - d / 2, d, d);

            Line2D input0 = new Line2D.Double(tileX + input0Offset.getX(), tileY + input0Offset.getY(),
                    tileX + input0Offset.getX() + 0.5 * tileWidth, tileY + input0Offset.getY());

            Line2D input1 = new Line2D.Double(tileX + input1Offset.getX(), tileY + input1Offset.getY(),
                    tileX + input1Offset.getX() + 0.5 * tileWidth, tileY + input1Offset.getY());

            Line2D output = new Line2D.Double(tileX + outputOffset.getX() - 0.5 * tileWidth, tileY + outputOffset.getY(),
                    tileX + outputOffset.getX(), tileY + outputOffset.getY());

            image.setPaint(new Color(0, 0, 0));
            image.draw(input0);
            image.draw(input1);
            image.draw(output);

            image.setPaint(new Color(245, 189, 189));
            image.fill(path);
            if (isNor)
                image.fill(circle);

            image.setPaint(new Color(0, 0, 0));
            image.draw(path);
            if (isNor)
                image.draw(circle);
        }

    }

    private static class InputGate {

        static Point2D outputOffset = new Point2D.Double(0.75 * tileWidth, 0.5 * tileHeight);

        static void add(Graphics2D image, double x, double y, String name) {

            int height = image.getFontMetrics().getHeight();
            int width = image.getFontMetrics().stringWidth(name);

            x = x + (double) (tileWidth - width) / 2;
            y = y + (double) (tileHeight - height) / 2 + height;

            image.drawString(name, (float) x, (float) y);
        }

    }

    private static class OutputBuffer {

        static Point2D inputOffset = new Point2D.Double(0.5 * tileWidth, 0.5 * tileHeight);

        static void add(Graphics2D image, double x, double y, String name) {
        }

    }

    private Map<Gate, Integer> getYValues(Map<Gate, Integer> gateLevels, List<Pair<Gate, Gate>> connections, int depth) {

        Map<Gate, Integer> yValues = new HashMap<>();

        /* naive approach: simply fill rows */

        int height = 0;

        for (int level = 0; level < depth; level++) {

            int finalLevel = level;
            Set<Gate> gates = gateLevels.entrySet().stream()
                    .filter(e -> e.getValue() == finalLevel)
                    .map(Map.Entry::getKey)
                    .collect(Collectors.toSet());

            int row = 0;

            for (Gate gate : gates) {
                yValues.put(gate, row * tileHeight);
                row++;
            }

            if (row > height)
                height = row;
        }

        logger.info("original crosspoints: " + numberOfCrosspoints(gateLevels, yValues, connections));

        /* greedy approach to minimize crosspoints */

        Map<Integer, List<Gate>> bestOrders = new HashMap<>();

        for (int iteration = 0; iteration < 5; iteration++) {
            for (int level = depth - 1; level >= 0; level--) {

                int originalCrosspoints = numberOfCrosspoints(gateLevels, yValues, connections);

                int finalLevel = level;
                List<Gate> gates = gateLevels.entrySet().stream()
                        .filter(e -> e.getValue() == finalLevel)
                        .map(Map.Entry::getKey)
                        .collect(Collectors.toList());

                for (int i = gates.size(); i < height; i++) {
                    gates.add(new Gate("dummy_" + i, LogicType.EMPTY));
                }

                bestOrders.putIfAbsent(level, new ArrayList<>(gates));

                org.apache.commons.collections4.iterators.PermutationIterator<Gate> iterator =
                        new org.apache.commons.collections4.iterators.PermutationIterator<>(gates);

                int minimumCrosspoints = originalCrosspoints;

                while (iterator.hasNext()) {

                    gates = iterator.next();

                    int row = 0;
                    for (Gate gate : gates) {

                        if (gate.getLogicType() != LogicType.EMPTY)
                            yValues.put(gate, row * tileHeight);

                        row++;
                    }

                    int crosspoints = numberOfCrosspoints(gateLevels, yValues, connections);

                    if (crosspoints < minimumCrosspoints) {
                        minimumCrosspoints = crosspoints;
                        bestOrders.put(level, new ArrayList<>(gates));
                    }

                    if (crosspoints == 0)
                        break;
                }

                int row = 0;
                for (Gate gate : bestOrders.get(level)) {

                    if (gate.getLogicType() != LogicType.EMPTY)
                        yValues.put(gate, row * tileHeight);

                    row++;
                }
            }
        }

        /* greedy approach to minimize y offset */

        for (int iteration = 0; iteration < 5; iteration ++) {
            for (int level = 0; level < depth - 1; level ++) {

                int finalLevel = level;
                List<Gate> gates = gateLevels.entrySet().stream()
                        .filter(e -> e.getValue() == finalLevel)
                        .map(Map.Entry::getKey)
                        .collect(Collectors.toList());

                for (Gate gate : gates) {

                    List<Integer> levelYValues = gates.stream()
                            .map(yValues::get)
                            .sorted(Comparator.naturalOrder())
                            .collect(Collectors.toList());

                    int gateIndex = levelYValues.indexOf(yValues.get(gate));

                    int lower;

                    if (gateIndex > 0)
                        lower = levelYValues.get(gateIndex - 1);
                    else
                        lower = -100;

                    int upper;

                    if (gateIndex < levelYValues.size() - 1)
                        upper = levelYValues.get(gateIndex + 1);
                    else
                        upper = (height - 1) * tileHeight;

                    int originalY = yValues.get(gate);
                    int yOffset = gateYOffset(gate, yValues, connections);

                    int newYValue = originalY + (yOffset / 2);
                    newYValue = Math.min(newYValue, upper - tileHeight);
                    newYValue = Math.max(newYValue, lower + tileHeight);

                    yValues.put(gate, newYValue);
                }
            }
        }

        logger.info("minimized crosspoints: " + numberOfCrosspoints(gateLevels, yValues, connections));

        return yValues;
    }

    private int gateYOffset(Gate gate, Map<Gate, Integer> gateYValues, List<Pair<Gate, Gate>> connections) {

        int yOffset = 0;

        for (Pair<Gate, Gate> connection : connections) {

            if (connection.second().equals(gate)) {

                Gate startGate = connection.first();
                Gate endGate = connection.second();

                yOffset += gateYValues.get(startGate) - gateYValues.get(endGate);
            }
        }

        return yOffset;
    }

    private int numberOfCrosspoints(Map<Gate, Integer> gateLevels, Map<Gate, Integer> yValues, List<Pair<Gate, Gate>> connections) {

        int crosspoints = 0;

        for (int i = 0; i < connections.size(); i++) {

            /* through-tile cross */

            Gate startGate1 = connections.get(i).first();
            Gate endGate1 = connections.get(i).second();

            //if (yValues.get(startGate1) + tileHeight >= yValues.get(endGate1) && yValues.get(startGate1) - tileHeight <= yValues.get(endGate1)) {

                int startLevel = gateLevels.get(startGate1);
                int endLevel = gateLevels.get(endGate1);

                crosspoints += (int) yValues.entrySet().stream()
                        .filter(e -> e.getValue() > yValues.get(startGate1) - tileHeight && e.getValue() < yValues.get(endGate1) + tileHeight)
                        .map(Map.Entry::getKey)
                        .map(gateLevels::get)
                        .filter(l -> (l < startLevel || l > endLevel))
                        .count();
            //}


            for (int j = i + 1; j < connections.size(); j++) {

                Gate startGate2 = connections.get(j).first();
                Gate endGate2 = connections.get(j).second();

                /* ignore unmapped gates */

                if (!yValues.containsKey(startGate1) || !yValues.containsKey(startGate2)
                        || !yValues.containsKey(endGate1) || !yValues.containsKey(endGate2))
                    continue;

                /* connections use (partly) the same levels (excl. endpoints) */

                if (gateLevels.get(startGate1) <= gateLevels.get(endGate2) || gateLevels.get(startGate2) <= gateLevels.get(endGate1))
                    continue;

                /* connections cross */

                if (yValues.get(startGate1) > yValues.get(startGate2) && yValues.get(endGate1) < yValues.get(endGate2))
                    crosspoints++;
                else if (yValues.get(startGate1) < yValues.get(startGate2) && yValues.get(endGate1) > yValues.get(endGate2))
                    crosspoints++;
            }
        }

        return crosspoints;
    }


}
