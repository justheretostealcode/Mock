package de.tu_darmstadt.rs.synbio.mapping;

import de.tu_darmstadt.rs.synbio.common.circuit.Circuit;

public class MappingResult implements Comparable<MappingResult> {

    private final Circuit structure;
    private final Assignment assignment;
    private final double score;

    private double secondaryScore;

    private long neededSimulations = 0;
    private long minimumBranchAndBoundSimulations = -1;

    private double upperEnergyConstraint;
    private double lowerScoreConstraint;

    public MappingResult(Circuit structure, Assignment assignment, double score) {
        this.structure = structure;
        this.assignment = assignment;
        this.score = score;
    }

    public MappingResult(Circuit structure, Assignment assignment, double score, double secondaryScore) {
        this(structure, assignment, score);
        this.secondaryScore = secondaryScore;
    }

    public MappingResult(Circuit structure, Assignment assignment, double score, double secondaryScore,
                         double upperEnergyConstraint, double lowerScoreConstraint) {

        this(structure, assignment, score, secondaryScore);

        this.upperEnergyConstraint = upperEnergyConstraint;
        this.lowerScoreConstraint = lowerScoreConstraint;
    }

    public double getScore() {
        return score;
    }

    public Circuit getStructure() {
        return structure;
    }

    public Assignment getAssignment() {
        return assignment;
    }

    public double getSecondaryScore() {
        return secondaryScore;
    }

    public double getUpperEnergyConstraint() { return upperEnergyConstraint; }

    public double getLowerScoreConstraint() {
        return lowerScoreConstraint;
    }

    public void setNeededSimulations(long sims) {
        this.neededSimulations = sims;
    }
    public void setMinimumBranchAndBoundSimulations(long minimumBranchAndBoundSimulations)  {this.minimumBranchAndBoundSimulations = minimumBranchAndBoundSimulations;}
    public long getNeededSimulations() {
        return neededSimulations;
    }
    public long getMinimumBranchAndBoundSimulations() {return minimumBranchAndBoundSimulations;}

    @Override
    public int compareTo(MappingResult cmp) {
        if (this.score < cmp.score)
            return -1;
        else if (cmp.score < this.score)
            return 1;
        return 0;
    }
}
