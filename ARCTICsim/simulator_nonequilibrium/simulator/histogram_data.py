"""
Author: Erik Kubazka
"""
import numpy as np


class GateCytometryData:
    def __init__(self, data_dict, gate_name="", target_maximum_mean_var_ratio=None, scaling_factor=None,
                 cutoff_percentage=None):
        self._histograms = {}
        self._original_data_dict = data_dict
        self._input_val_scaling_mapping = {}
        self.gate_name = gate_name

        self.cutoff_percentage = cutoff_percentage
        self.histogram_is_cleaned = cutoff_percentage is not None

        for input_val in data_dict:
            cello_bins = data_dict[input_val]["cello_bins"]
            hist = data_dict[input_val]["hist"]

            self._histograms[input_val] = HistogramData(cello_bins=cello_bins,
                                                        histogram=hist,
                                                        scaling_factor=scaling_factor,
                                                        cutoff_percentage=cutoff_percentage)
            self._input_val_scaling_mapping[input_val] = input_val
        self.target_maximum_mean_var_ratio = target_maximum_mean_var_ratio
        self.scaling_factor = scaling_factor
        self.bins_are_rescaled = target_maximum_mean_var_ratio is not None or scaling_factor is not None
        if self.bins_are_rescaled:
            self.rescale()

    def derive_scaling_factor(self, scaling_factor):
        if scaling_factor is None:
            scaling_factor = self.scaling_factor

        mean_var_ratios = list(map(lambda input_val: self._histograms[input_val].mean_var_ratio(), self._histograms))
        max_mean_var_ratio = max(mean_var_ratios)

        if scaling_factor is None:
            gamma = max_mean_var_ratio
            if self.target_maximum_mean_var_ratio is not None and self.target_maximum_mean_var_ratio > 0:
                gamma /= self.target_maximum_mean_var_ratio

            scaling_factor = gamma

        target_mean_var_ratio = max_mean_var_ratio / scaling_factor

        return scaling_factor, target_mean_var_ratio

    def rescale(self, scaling_factor=None):
        scaling_factor, target_mean_var_ratio = self.derive_scaling_factor(scaling_factor=scaling_factor)
        old_input_val_scaling_mapping = self._input_val_scaling_mapping
        new_input_val_scaling_mapping = {}
        new_hist = {}
        for input_val in self._original_data_dict:
            new_input_val_scaling_mapping[input_val] = input_val * scaling_factor
            new_hist[new_input_val_scaling_mapping[input_val]] = self._histograms[
                old_input_val_scaling_mapping[input_val]].rescale_bins(scaling_factor)

        self.scaling_factor = scaling_factor
        self.target_maximum_mean_var_ratio = target_mean_var_ratio

        self._histograms = new_hist
        self._input_val_scaling_mapping = new_input_val_scaling_mapping
        return self

    def get_cytometry_data(self):
        return self._histograms

    def get_input_levels(self):
        return list(self._histograms.keys())

    def get_histogram(self, input_level):
        if input_level in self._histograms:
            return self._histograms[input_level]
        return None

    def get_max_mean_var_ratio(self):
        mean_var_ratios = [self._histograms[input_val].mean_var_ratio() for input_val in self._histograms]
        max_mean_var_ratio = max(mean_var_ratios)
        return max_mean_var_ratio

    @classmethod
    def cello_cytometry_data_to_dict(cls, cello_table):
        data_dict = {}
        for entry in cello_table:
            if "x" in entry:
                input_val = entry["x"]
                cello_bins = entry["bin"]
                hist = entry["output"]
            elif "input":
                input_val = entry["input"]
                cello_bins = entry["output_bins"]
                hist = entry["output_counts"]

            data_dict[input_val] = {"cello_bins": cello_bins, "hist": hist}

        return data_dict


class HistogramData:
    def __init__(self, bins=None, cello_bins=None, histogram=None, scaling_factor=None, cutoff_percentage=None):
        self.bin_centers = None
        self.bins = None

        if bins is None:
            if cello_bins is None:
                raise Exception("One of the bins must not be NONE")
            bins = self.cello_bins_to_bins(cello_bins=cello_bins)
        if histogram is None:
            raise Exception("The parameter histogram must not be NONE")

        self._original_bins = bins

        self.init_bins(bins)
        self.histogram = np.array(histogram)

        self.scaling_factor = scaling_factor
        self.bins_are_rescaled = scaling_factor is not None
        if self.bins_are_rescaled:
            self.rescale_bins(scaling_factor)

        self.probability_mass_preserved = 1
        self.cutoff_percentage = cutoff_percentage
        self.histogram_is_cleaned = cutoff_percentage is not None
        if self.histogram_is_cleaned:
            # Removes low mass regions to clean noise from measurements.
            # ToDo Think of alternative to preserve fixed amount of probability mass (e.g. 90 %)
            self.uncleaned_histogram = np.array(histogram)
            threshold = cutoff_percentage * np.max(self.histogram)
            self.histogram[self.histogram < threshold] = 0
            self.probability_mass_removed = np.sum(self.histogram)
            self.histogram = self.histogram / np.sum(self.histogram)
        pass
    def init_bins(self, bins):
        self.bins = np.array(bins)
        self.bin_centers = self._derive_bin_centers()

    def rescale_bins(self, scaling_factor):
        if scaling_factor is not None:
            rescaled_bins = self._original_bins * scaling_factor
            self.init_bins(rescaled_bins)
        return self

    def _derive_bin_centers(self):
        bin_to_bin_ratio = self.bins[1] / self.bins[0]
        bin_centers = self.bins[:-1] * bin_to_bin_ratio ** (0.5)  # Ist die Mitte der Bins im Logarithmischen
        return bin_centers

    def raw_moment(self, i):
        moment = np.power(self.bin_centers, i) @ self.histogram
        return moment

    def central_moment(self, i):
        mean = self.mean()
        moment = np.power((self.bin_centers - mean), i) @ self.histogram
        return moment

    def standardized_moment(self, i):
        standard_deviation = self.standard_deviation()

        central_moment = self.central_moment(i)
        standardized_moment = central_moment / np.power(standard_deviation, i)

        return standardized_moment

    def median(self):
        # Follow    https://www.statology.org/histogram-mean-median/
        # and       https://www2.microstrategy.com/producthelp/current/FunctionsRef/Content/FuncRef/Hist_Median.htm

        cdf = np.cumsum(self.histogram)
        median_mask = np.diff(cdf < 0.5)
        median_bins_index = np.where(median_mask)[0] + 1
        median_bin_lower_bound = self.bins[median_bins_index]
        median_bin_width = self.bins[median_bins_index + 1] - self.bins[median_bins_index]
        median = median_bin_lower_bound + median_bin_width * (0.5 - cdf[median_bins_index - 1]) / self.histogram[
            median_bins_index]
        median = median[0]
        return median

    def mean(self):
        mean = self.raw_moment(1)
        return mean

    def variance(self):
        variance = self.central_moment(2)
        return variance

    def standard_deviation(self):
        variance = self.variance()
        standard_deviation = np.sqrt(variance)
        return standard_deviation

    def mean_var_ratio(self):
        mean = self.mean()
        variance = self.variance()
        return mean / variance

    @classmethod
    def cello_bins_to_bins(cls, cello_bins):
        cello_bins = np.array(cello_bins)

        bin_to_bin_ratio = cello_bins[1] / cello_bins[0]
        bins = np.empty(cello_bins.shape[0] + 1)
        bins[:cello_bins.shape[0]] = cello_bins
        bins[-1] = cello_bins[-1] * bin_to_bin_ratio
        return bins
