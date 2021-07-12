import matplotlib.pyplot as plt
import numpy as np
import os
import marshal
import types
import json
import sys
import base64

# Required for producing the combined plot
import svgutils as sg
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF

# import os

"""
@author: Erik Kubaczka
"""
"""
Constants
"""
VISUALISATION_TARGET_DIRECTORY_KEY = "CIRCUIT_VISUALISATION_LOCATION"
DEFAULT_TARGET_DIRECTORY = "visualisation/"

REPRESSOR_COLOR_ASSIGNMENT_KEY = "REPRESSOR_COLOR_ASSIGNMENT"

INCH_TO_CM = lambda inch: inch * 2.54
CM_TO_INCH = lambda cm: cm / 2.54

DIN_A4 = np.array((21.0, 29.7))
MATPLOTLIB_DEFAULT_FORMAT = np.array(plt.rcParams["figure.figsize"])
plt.rcParams["figure.dpi"] = 90
MATPLOTLIB_DEFAULT_DPI = plt.rcParams["figure.dpi"]

FONTSIZE = {"DEFAULT": 7,
            "BAR_LABEL": 4,
            }
LINE_INDICATOR_PROPERTIES = {"COLOR_AND_MARKER": "k--",
                             "ALPHA": 0.5,
                             "LINE_WIDTH": 1,
                             }


def SVG_VALUE_AND_UNIT(str, targetUnit=None):
    unit = sg.compose.Unit(str)
    if (targetUnit != None):
        unit = unit.to(targetUnit)

    return unit.value, unit.unit


def SVG_WIDTH_AND_HEIGHT(figure, targetUnit=None):
    width, unitW = SVG_VALUE_AND_UNIT(figure.width, targetUnit=targetUnit)
    height, unitH = SVG_VALUE_AND_UNIT(figure.height, targetUnit=targetUnit)

    return width, height


SVG_WIDTH_AND_HEIGHT_IN_CM = (lambda figure: SVG_WIDTH_AND_HEIGHT(figure, "cm"))
SVG_WIDTH_AND_HEIGHT_IN_PX = (lambda figure: SVG_WIDTH_AND_HEIGHT(figure, "px"))

AS_CM = (lambda val: "%fcm" % val)
AS_PX = (lambda val: "%fpx" % val)

"""
Helper Functions
"""


def escapePath(path, isDir=True):
    newPath = path.replace("\\", "/")
    newPath += "" if (newPath.endswith("/") or newPath == "" or (not isDir)) else "/"
    return newPath


def createDirectory(dirPath):
    if (not os.path.exists(dirPath)):
        os.makedirs(dirPath)


def getTargetDirectory(simContext):
    targetDirectory = DEFAULT_TARGET_DIRECTORY

    if (simContext != None and VISUALISATION_TARGET_DIRECTORY_KEY in simContext):
        targetDirectory = simContext[VISUALISATION_TARGET_DIRECTORY_KEY]
        targetDirectory = escapePath(targetDirectory)

    # Creates the directory if not already present
    createDirectory(targetDirectory)
    return targetDirectory


def loadJSON(filePath):
    with open(filePath, 'r') as jsonFile:
        data = json.load(jsonFile)
    return data


def writeJSON(filePath, data):
    with open(filePath, 'w') as jsonFile:
        # Sort keys needs to be false in order to preserve the node order in the circuit
        json.dump(data, jsonFile, indent=4, sort_keys=False)


def encodeAsBase64(data):
    return base64.standard_b64encode(data).decode("ascii")


def decodeAsBase64(data):
    return base64.standard_b64decode(data.encode("ascii"))


# onlyUpdate: If true, only keys included in targetDict are updated and no new are added
def updateDict(targetDict, sourceDict, onlyUpdate=True):
    if (targetDict != None and sourceDict != None):
        for key in sourceDict:
            if (not onlyUpdate or key in targetDict):  # Ensure only used values can be overwritten
                targetDict[key] = sourceDict[key]
        targetDict["UPDATED"] = True
    return targetDict


"""
Attention: Response Functions are only correct for the last provided parameter -> Particles not supported
"""


def visualiseCircuit(circuit, assignment, responseFunctions, circuitValues, simContext={}):
    visualizationType = simContext["visualise_circuit_visualisation_mode"]
    if (visualizationType == 1):
        visualiseCircuitWithValues(circuit=circuit, assignment=assignment, responseFunctions=responseFunctions,
                                   circuitValues=circuitValues)
    elif (visualizationType == 2):
        plotCircuitWithValueSummary(circuit=circuit, assignment=assignment, responseFunctions=responseFunctions,
                                    simContext=simContext,
                                    circuitVals=circuitValues)
    elif (visualizationType == 3):
        storeValuesForSubsequentVisualization(circuit=circuit, assignment=assignment,
                                              responseFunctions=responseFunctions, circuitVals=circuitValues,
                                              simContext=simContext)


# Calls plotCircuitWithValues for each possible input assignment and thus yields a visualisation of every combination possible
def visualiseCircuitWithValues(circuit, assignment, responseFunctions, circuitValues):
    for inputID in circuitValues:
        plotCircuitWithValues(inputID, circuit=circuit, assignment=assignment, responseFunctions=responseFunctions,
                              circuitVals=circuitValues[inputID])

    return


# Visualizes the transfer characteristic of the gates in combination with the concentrations present within the circuit
# The plot is saved to .svg and .png
# plotAsSubplots: If this is true, the plots of the single gates are grouped into on single figure,
# otherwise a separated folder is created for this input, which includes the single figures as images
def plotCircuitWithValues(inputID, circuit, assignment, responseFunctions, circuitVals, plotAsSubplots=False):
    # Define the definition domain of interest
    X = np.logspace(-4, 2, 100)
    # Determine the dimensions of the subplot to adequatly arange all plots
    plotDimensions = estimateDimensionsOfSubplots(len(assignment))

    if (plotAsSubplots):
        plt.figure()
        fig, ax = plt.subplots(int(plotDimensions[0]), int(plotDimensions[1]))
        fig.subplots_adjust(wspace=1, hspace=1)

    getCurrentAxis = (lambda i: ax[int(i / plotDimensions[1])][int(i % plotDimensions[1])]) if plotAsSubplots else (
        lambda i: plt.figure().gca())

    i = 0
    # Iterate over the single gates (input buffers, logic gates and output buffers) of a circuit
    for elem in circuit:
        # currentAxis = ax[int(i / plotDimensions[1])][int(i % plotDimensions[1])]
        currentAxis = getCurrentAxis(i=i)
        # plt.subplot(plotDimensions)
        if (elem not in assignment):  # non logic gates are skipped since there is no transfer characteristic to show
            continue

        # Plot the corresponding responsefunction
        responseFunction = responseFunctions[assignment[elem]]
        Y = plotResponseFunction(ax=currentAxis, X=X, responseFunction=responseFunction["equation"],
                                 parameters=responseFunction["parameters"])

        # Get the gates input values
        inputValues, inputVal, outputVal = getGateInputValues(elem, circuit, circuitVals)
        # Check if multiple input values exist and plot them if so
        if (len(inputValues) > 1):
            for inVal in inputValues:
                currentAxis.plot([inVal, inVal], [min(Y), outputVal], "gx")

        # Plot the mapping from cummulative input via transition point to output value
        currentAxis.plot([inputVal, inputVal, max(X)], [min(X), outputVal, outputVal], "k-", linewidth=1)

        currentAxis.set_title(elem)  # + " (" + assignment[elem] + ")")
        currentAxis.set_xscale("log")
        currentAxis.set_yscale("log")
        # currentAxis.set_ylim(bottom=min(Y) / 4)

        # Ensure quadratic domain and codomain
        currentAxis.set_xlim(left=min(X), right=max(X))
        currentAxis.set_ylim(bottom=min(X), top=max(X))

        i += 1
        if (not plotAsSubplots):
            targetDirectory = "%s/" % inputID
            createDirectory(targetDirectory)
            targetPath = targetDirectory + elem
            plt.savefig("%s.svg" % targetPath)
            plt.savefig("%s.png" % targetPath)
            # plt.close(fig)

    # print(os.path.abspath("visualisation/" + inputID + ".svg"))
    # Save the plot
    if (plotAsSubplots):
        fig.suptitle("Circuit for Input:" + inputID)
        plt.savefig("visualisation/" + inputID + ".svg")
        plt.savefig("visualisation/" + inputID + ".png")
        plt.close(fig)
    # plt.show()


def plotCircuitWithValueSummary(circuit, assignment, responseFunctions, circuitVals, plotAsSubplots=False,
                                simContext={},
                                updatedFigureProperties=None):
    figureProperties = {"fontsize": 11,
                        "figure_format": MATPLOTLIB_DEFAULT_FORMAT,
                        "X-LIM": [-4, 2],
                        "X-N": 100,
                        "INDICATOR_color_and_marker": "k-",
                        "INDICATOR_linewidth": 0.5,
                        REPRESSOR_COLOR_ASSIGNMENT_KEY: None}
    figureProperties["plot_as_subplots"] = plotAsSubplots

    updateDict(targetDict=figureProperties, sourceDict=updatedFigureProperties)
    figureProperties["figure_format"] = tuple(CM_TO_INCH(np.array(figureProperties["figure_format"])))

    targetDirectory = getTargetDirectory(simContext=simContext)

    repressorColorAssignment = None

    if (REPRESSOR_COLOR_ASSIGNMENT_KEY in figureProperties
            and figureProperties[REPRESSOR_COLOR_ASSIGNMENT_KEY] != None):
        repressorColorAssignment = figureProperties[REPRESSOR_COLOR_ASSIGNMENT_KEY]

    # Define the definition domain of interest
    X = np.logspace(figureProperties["X-LIM"][0], figureProperties["X-LIM"][1], figureProperties["X-N"])
    # Determine the dimensions of the subplot to adequatly arange all plots
    plotDimensions = estimateDimensionsOfSubplots(len(assignment))

    if (figureProperties["plot_as_subplots"]):
        plt.figure(figsize=figureProperties["figure_format"], frameon=False)
        fig, ax = plt.subplots(int(plotDimensions[0]), int(plotDimensions[1]))
        fig.subplots_adjust(wspace=1, hspace=1)

    getCurrentAxis = (lambda i: ax[int(i / plotDimensions[1])][int(i % plotDimensions[1])]) if plotAsSubplots else (
        lambda i: plt.figure(figsize=figureProperties["figure_format"], frameon=False).gca())

    i = 0
    # Iterate over the single gates (input buffers, logic gates and output buffers) of a circuit
    for elem in circuit:
        # currentAxis = ax[int(i / plotDimensions[1])][int(i % plotDimensions[1])]
        currentAxis = getCurrentAxis(i=i)
        # plt.subplot(plotDimensions)
        if (elem not in assignment):  # non logic gates are skipped since there is no transfer characteristic to show
            continue

        # Plot the corresponding responsefunction
        responseFunction = responseFunctions[assignment[elem]]

        properties = {"color_and_marker": None,
                      "color": None,
                      "marker":None}
        if (repressorColorAssignment != None
                and "group" in responseFunction
                and responseFunction["group"] in repressorColorAssignment):
            properties = repressorColorAssignment[responseFunction["group"]]

            # figureProperties["color"] = properties["color"]
            # figureProperties["marker"] = properties["marker"]
        elif (repressorColorAssignment != None):
            properties = repressorColorAssignment["DEFAULT"]


        updateDict(targetDict=figureProperties, sourceDict=properties, onlyUpdate=False)
        Y = plotResponseFunction(ax=currentAxis, X=X, responseFunction=responseFunction["equation"],
                                 parameters=responseFunction["parameters"], figureProperties=figureProperties)

        # Get the gates input values

        for inputID in circuitVals:
            inputValues, inputVal, outputVal = getGateInputValues(elem, circuit, circuitVals=circuitVals[inputID])

            # Plot the mapping from cummulative input via transition point to output value
            currentAxis.plot([inputVal, inputVal, max(X)], [min(X), outputVal, outputVal],
                             figureProperties["INDICATOR_color_and_marker"],
                             linewidth=figureProperties["INDICATOR_linewidth"])

        currentAxis.set_title(elem)  # + " (" + assignment[elem] + ")")
        currentAxis.set_xscale("log")
        currentAxis.set_yscale("log")
        # currentAxis.set_ylim(bottom=min(Y) / 4)

        # Ensure quadratic domain and codomain
        currentAxis.set_xlim(left=min(X), right=max(X))
        currentAxis.set_ylim(bottom=min(X), top=max(X))

        i += 1
        if (not plotAsSubplots):
            plt.tight_layout()
            targetPath = targetDirectory + elem
            plt.savefig("%s.svg" % targetPath)
            #plt.savefig("%s.png" % targetPath)
            # plt.close(fig)

    # print(os.path.abspath("visualisation/" + inputID + ".svg"))
    # Save the plot
    if (plotAsSubplots):
        plt.tight_layout()
        fig.suptitle("Circuit for Input:" + inputID)
        plt.savefig("%s/circuit.svg" % targetDirectory)
        plt.savefig("%s/circuit.png" % targetDirectory)
        plt.close(fig)
    # plt.show()

    return targetDirectory


# Determines the single input values of a gate.
# In case of a NOR, two input values exist, while they can not be reconstructed from the gates output values
def getGateInputValues(gate, circuit, circuitVals):
    inputValues = []
    for src in circuit:
        if (gate in circuit[src]):
            inputValues.append(circuitVals[src])

    gateValue = circuitVals[gate]

    # The input values as a list, their sum as well as the corresponding output of the gate
    return inputValues, sum(inputValues), gateValue


# Adds the response function (specified by responseFunction and parameters) sampled at positions X to the provided axes (ax)
def plotResponseFunction(ax, X, responseFunction, parameters, figureProperties=None):
    curveProperties = {
        "linewidth": 3,
        "color_and_marker": None,
        "color": None,
        "marker": None,
        "UPDATED": False
    }
    updateDict(targetDict=curveProperties, sourceDict=figureProperties)

    Y = np.zeros(len(X))
    for i in range(len(X)):
        Y[i] = responseFunction(X[i], parameters)

    if (curveProperties["color_and_marker"] != None):
        ax.plot(X, Y, curveProperties["color_and_marker"], linewidth=curveProperties["linewidth"])
    elif (curveProperties["color"] != None and curveProperties["marker"] != None):
        ax.plot(X, Y, c=curveProperties["color"], marker=curveProperties["marker"],
                linewidth=curveProperties["linewidth"])
    elif (curveProperties["color"] != None):
        ax.plot(X, Y, c=curveProperties["color"], linewidth=curveProperties["linewidth"])
    else:
        ax.plot(X, Y, linewidth=curveProperties["linewidth"])

    return Y


# Determines the number of subplots required for adequatly visualizing the circuit.
def estimateDimensionsOfSubplots(n):
    squareRoot = np.sqrt(n)
    dim = np.zeros(2)
    dim[0] = np.ceil(squareRoot)
    dim[1] = np.floor(n / dim[0])

    i = 1
    while (dim[0] * dim[1] < n):
        dim[i] = dim[i] + 1
        i = 1 - i
    return dim


# Stores the information relevant for visualizing a circuits values in a file
def storeValuesForSubsequentVisualization(circuit, assignment, responseFunctions, circuitVals, simContext={}):
    jsonMap = {}
    jsonMap["CIRCUIT"] = circuit
    jsonMap["ASSIGNMENT"] = assignment
    jsonMap["RESPONSE_FUNCTIONS"] = serializeDeserializeResponseFunctions(responseFunctions=responseFunctions,
                                                                          serialize=True)
    jsonMap["CIRCUIT_VALUES"] = serializeDeserializeCircuitValues(circuitVals=circuitVals, serialize=True)

    targetDirectory = getTargetDirectory(simContext=simContext)
    # Creates the directory if not already present
    createDirectory(targetDirectory)

    targetPath = targetDirectory + "circuit_visualisation_information.json"
    writeJSON(filePath=targetPath, data=jsonMap)


def serializeDeserializeResponseFunctions(responseFunctions, serialize=True):
    def serializeEquation(equation):
        equation = marshal.dumps(equation.__code__)
        return encodeAsBase64(equation)

    def deserializeEquation(equation):
        code = marshal.loads(decodeAsBase64(equation))
        return types.FunctionType(code, globals())

    operation = serializeEquation if serialize else deserializeEquation

    serializedResponseFunctions = {}
    for gate in responseFunctions:
        newEntry = {}
        for key in responseFunctions[gate]:
            newEntry[key] = responseFunctions[gate][key]

        # TODO implement serialization and deserialization of particles
        if ("equation" in newEntry):
            newEntry["equation"] = operation(newEntry["equation"])

        serializedResponseFunctions[gate] = newEntry
    return serializedResponseFunctions


def serializeDeserializeCircuitValues(circuitVals, serialize=True):
    operation = (lambda arg: list(arg)) if serialize else (lambda arg: np.array(arg))
    serializedCircuitValues = {}
    for inputID in circuitVals:
        serializedCircuitValues[inputID] = {}
        for key in circuitVals[inputID]:
            serializedCircuitValues[inputID][key] = operation(circuitVals[inputID][key])

    return serializedCircuitValues


# Plots two CDFs in the same graph
# Thereby CDF1 and CDF2 need to be defined for the values given in positions
# positions: The X-axes value
# CDF1: The first CDF
# CDF2: The second CDF
# [start, end]: The interval of the area between the two CDFs, which is taken into account for scoring purposes
def envelope_plotCDFs(positions, CDF1, CDF2, start, end):
    plt.figure()
    plt.plot(positions, CDF1, label="CDF1")
    plt.plot(positions, CDF2, label="CDF2")
    plt.plot(np.ones(2) * start, [0, 1], "--", label="Start")
    plt.plot(np.ones(2) * end, [0, 1], "--", label="End")
    plt.legend()
    plt.title("The CDFs")
    plt.show()


# Plots the difference between two CDFs
# Thereby CDF1 and CDF2 need to be defined for the values given in positions
# positions: The X-axes value
# CDF1: The first CDF
# CDF2: The second CDF
# [start, end]: The interval of the area between the two CDFs, which is taken into account for scoring purposes
def envelope_plotCDFsDiff(positions, CDF1, CDF2, start, end):
    cdf1 = np.array(CDF1)
    cdf2 = np.array(CDF2)

    cdf_diff = cdf1 - cdf2

    plt.figure()
    plt.plot(positions, cdf_diff, label="CDF Difference")
    plt.plot(positions, np.abs(cdf_diff), label="Absolute Difference")
    plt.plot(np.ones(2) * start, [0, 1], "--", label="Start")
    plt.plot(np.ones(2) * end, [0, 1], "--", label="End")
    plt.legend()
    plt.title("Difference of CDFs")
    plt.show()


"""
Methods for visualising a circuit independent of simulation
"""


def arrangeCombinedFigure(circuitFilePath, targetDirectory, config, circuitInformation, visualisationContext):
    def determineSubplotsLayout(svgFigures, margin):
        layout = {}
        subplotDimensions = estimateDimensionsOfSubplots(len(svgFigures.keys()))

        overallWidth = 0
        overallHeight = 0
        xOffset = 0
        yOffset = 0
        prevHeight = 0
        i = 0
        for key in svgFigures:
            iX = i % subplotDimensions[0]
            iY = int(i / subplotDimensions[0])
            if (iX == 0):
                xOffset = 0
                yOffset += prevHeight
            figWidth, figHeight = SVG_WIDTH_AND_HEIGHT_IN_CM(svgFigures[key])
            layout[key] = tuple([xOffset, yOffset])
            xOffset += figWidth

            overallWidth = max(overallWidth, xOffset)
            xOffset += margin

            if (iY + 1 != subplotDimensions[1]):
                prevHeight = max(prevHeight, figHeight + margin)
            else:
                prevHeight = max(prevHeight, figHeight)
            overallHeight = max(overallHeight, yOffset + prevHeight)
            i += 1

        return layout, overallWidth, overallHeight

    # Not the Roots of svgs
    # Horizontal Layout is default
    def determineLayout(circuitSVG, svgFigures, borderMargin=0.5, innerMargin=0.5, subplotsMargin=-0.4,
                        layoutDirection="HORIZONTAL"):
        # layout = {}
        # All measurements are in cm's
        circuitWidth = 0
        circuitHeight = 0
        if (circuitSVG != None):
            circuitWidth, circuitHeight = SVG_WIDTH_AND_HEIGHT_IN_CM(circuitSVG)

        layout, subplotsWidth, subplotsHeight = determineSubplotsLayout(svgFigures, margin=subplotsMargin)

        overallWidth = max(subplotsWidth, circuitWidth)
        overallHeight = max(subplotsHeight, circuitHeight)
        overallWidth += 2 * borderMargin
        overallHeight += 2 * borderMargin
        if (layoutDirection == "HORIZONTAL"):
            # maxHeight does not change
            # circuit must be in the middle of overallHeight
            # subLayout must be in the middle of overallHeight

            xOffsetCircuit = borderMargin
            yOffsetCircuit = int(np.round((overallHeight - circuitHeight) / 2))

            xOffsetSublayout = circuitWidth + innerMargin
            yOffsetSublayout = int(np.round((overallHeight - subplotsHeight) / 2))

            overallWidth = xOffsetSublayout + subplotsWidth

        elif (layoutDirection == "VERTICAL"):
            # maxWidth does not change
            xOffsetCircuit = int(np.round((overallWidth - circuitWidth) / 2))
            yOffsetCircuit = borderMargin

            xOffsetSublayout = int(np.round((overallWidth - subplotsWidth) / 2))
            yOffsetSublayout = circuitHeight + innerMargin

            overallHeight = yOffsetSublayout + subplotsHeight

        for key in layout:
            layout[key] = (layout[key][0] + xOffsetSublayout, layout[key][1] + yOffsetSublayout)

        layout["CIRCUIT"] = (xOffsetCircuit, yOffsetCircuit)

        return layout, overallWidth, overallHeight

    def applyLayout(circuitSVG, svgFigures, layout):
        def applyOffset(svgFigureRoot, offset):
            svgFigureRoot.moveto(x=SVG_VALUE_AND_UNIT(AS_CM(offset[0]), targetUnit="px")[0],
                                 y=SVG_VALUE_AND_UNIT(AS_CM(offset[1]), targetUnit="px")[0])

        # Caution, translation is performed based on pixels and not cm
        figureRoots = {}

        circuitOffset = layout["CIRCUIT"]
        if (circuitSVG != None):
            circuitSVGroot = circuitSVG.getroot()
            applyOffset(svgFigureRoot=circuitSVGroot, offset=circuitOffset)
            figureRoots["CIRCUIT"] = circuitSVGroot

        for key in svgFigures:
            currentOffset = layout[key]
            currentRoot = svgFigures[key].getroot()
            applyOffset(svgFigureRoot=currentRoot, offset=currentOffset)
            figureRoots[key] = currentRoot

        return figureRoots

    figureProperties = {"SVG-WIDTH": 20,
                        "SVG-HEIGHT": 20,
                        "SVG-LAYOUT-DIRECTION": "HORIZONTAL"}
    updateDict(targetDict=figureProperties, sourceDict=config, onlyUpdate=True)

    circuit = circuitInformation["CIRCUIT"]
    assignment = circuitInformation["ASSIGNMENT"]

    gateOrder = [gate for gate in circuit if (gate in assignment)]
    # filesList = ["%s%s.svg" % (targetDirectory, gate) for gate in gateOrder]

    circuitSVG = None
    if (circuitFilePath != None):
        circuitSVG = sg.transform.fromfile(circuitFilePath)
        print("Read circuit SVG")
    # fig = sg.transform.SVGFigure(figureProperties["SVG-WIDTH"], figureProperties["SVG-HEIGHT"])
    # fig = sg.transform.SVGFigure(20, 20)
    # fig = sg.transform.SVGFigure("%fpt" % figureProperties["SVG-WIDTH"], "%fpt" % figureProperties["SVG-HEIGHT"])

    svgFigures = {gate: sg.transform.fromfile("%s%s.svg" % (targetDirectory, gate)) for gate in gateOrder}
    print("Read SVG of plot figures")

    print("Determining Layout")
    layout, overallWidth, overallHeight = determineLayout(circuitSVG=circuitSVG, svgFigures=svgFigures,
                                                          layoutDirection=figureProperties["SVG-LAYOUT-DIRECTION"])

    print("Applying the layout")
    figureRoots = applyLayout(circuitSVG=circuitSVG, svgFigures=svgFigures, layout=layout)

    # svgFiguresRoots = {gate: svgFigures[gate].getroot() for gate in svgFigures.keys()}
    print("Creating the final figure")
    fig = sg.transform.SVGFigure(sg.compose.Unit("%fcm" % overallWidth).to("px"),
                                 sg.compose.Unit("%fcm" % overallHeight).to("px"))

    # if (circuitSVG != None):
    #     fig.append(circuitSVG.getroot())
    # fig.append(svgFiguresRoots.values())
    fig.append(figureRoots.values())


    combinedFigureName = "combined_figure.svg"
    combinedFigurePath = targetDirectory + combinedFigureName
    fig.save(combinedFigurePath)
    print("Saved the resulting figure to:\n %s" % combinedFigurePath)

    return combinedFigurePath, combinedFigureName


def visualiseCircuitAndTransferFunctions(config, circuitInformation):
    def cleanDirectory(dirPath, filesToPreserve):
        for file in os.listdir(dirPath):
            if (not file in filesToPreserve):
                filePath = dirPath + file
                if (os.path.exists(filePath)):
                    os.remove(filePath)

    visualisationContext = {}

    figurePropertiesToUpdate = {}
    updateDict(targetDict=figurePropertiesToUpdate, sourceDict=configFile, onlyUpdate=False)
    responseFunctions = serializeDeserializeResponseFunctions(
        responseFunctions=circuitInformation["RESPONSE_FUNCTIONS"], serialize=False)
    circuitValues = serializeDeserializeCircuitValues(circuitVals=circuitInformation["CIRCUIT_VALUES"], serialize=False)

    targetDirectory = getTargetDirectory(simContext=visualisationContext)
    filesInTargetDirectoryToPreserve = os.listdir(targetDirectory)

    # Plot the single transfer characteristics
    print("Plot circuit values and transfer characteristic")
    plotCircuitWithValueSummary(circuit=circuitInformation["CIRCUIT"], assignment=circuitInformation["ASSIGNMENT"],
                                responseFunctions=responseFunctions, circuitVals=circuitValues, plotAsSubplots=False,
                                updatedFigureProperties=figurePropertiesToUpdate, simContext=visualisationContext)

    circuitFilePath = None
    if ("CIRCUIT_FILE" in config):
        circuitFilePath = config["CIRCUIT_FILE"]
        escapePath(circuitFilePath, isDir=True)

    print("Arange combined figure")
    combinedFigurePath, combinedFigureName = arrangeCombinedFigure(circuitFilePath=circuitFilePath,
                                                                   targetDirectory=targetDirectory,
                                                                   config=config,
                                                                   circuitInformation=circuitInformation,
                                                                   visualisationContext=visualisationContext)

    print("Create PDF based on SVG")
    combinedFigure = svg2rlg(combinedFigurePath)
    renderPDF.drawToFile(combinedFigure, combinedFigurePath.replace(".svg", ".pdf"))

    # Ensure to preserve output file
    filesInTargetDirectoryToPreserve.append(combinedFigureName)
    filesInTargetDirectoryToPreserve.append(combinedFigureName.replace(".svg", ".pdf"))
    cleanDirectory(dirPath=targetDirectory, filesToPreserve=filesInTargetDirectoryToPreserve)


if __name__ == '__main__':
    argv = sys.argv[1:]
    if (len(argv) != 2):
        raise Exception(
            "The first argument needs to be the path to the .json config file and the second a path pointing to circuit_visualisation_information.json")

    configFilePath = argv[0]
    circuitInformationFilePath = argv[1]

    if (not os.path.exists(configFilePath)):
        raise Exception("Config File does not exist\n%s" % configFilePath)
    elif (not configFilePath.endswith(".json")):
        raise Exception("Config File is not a json file")

    if (not os.path.exists(circuitInformationFilePath)):
        raise Exception("Circuit Information File does not exist\n%s" % circuitInformationFilePath)
    elif (not circuitInformationFilePath.endswith(".json")):
        raise Exception("Circuit Information File is not a json file")

    configFile = loadJSON(filePath=configFilePath)
    circuitInformation = loadJSON(filePath=circuitInformationFilePath)

    visualiseCircuitAndTransferFunctions(config=configFile, circuitInformation=circuitInformation)
