import os
import pandas as pd
import numpy as np
from nilearn.maskers import NiftiLabelsMasker, NiftiSpheresMasker


def trlocked_events(events_path: str, onsets_column: str, trial_name: str,
                    bold_tr: float, bold_vols: int, separator: str = '\t'):
    """
    Loads behavior data, creates and merges into a TR (rounded) dataframe to match length of BOLD. Trial onsets are
    matched to nearby TR using rounding when acquisition is not locked to TR.

    Parameters:
        events_path (str): Path to the events data files for given subject/run.
        onsets_column (str): Name of the column containing onset times for the event/condition.
        trial_name (str): Name of the column containing condition/trial labels.
        bold_tr (int): TR acquisition of BOLD in seconds.
        bold_vols (int): Number of time points for BOLD acquisition
        separator (str): Separator used in the events data file, default = '\t'.
    Returns:
        pandas.DataFrame: Merged dataframe with time index and events data for each event + TR delays.
    """
    if not os.path.exists(events_path):
        raise FileNotFoundError(f"File '{events_path}' not found.")

    beh_df = pd.read_csv(events_path, sep=separator)

    missing_cols = [col for col in [onsets_column, trial_name] if col not in beh_df.columns]
    if missing_cols:
        raise KeyError(f"Missing columns: {', '.join(missing_cols)}")

    beh_df = beh_df[[onsets_column, trial_name]]
    beh_df["TimePoint"] = (beh_df[onsets_column] / bold_tr).round()

    time_index = pd.RangeIndex(start=0, stop=bold_vols, step=1)
    time_index_df = pd.DataFrame(index=time_index)
    # Merge behavior data with time index
    merged_df = pd.merge(time_index_df, beh_df, how='left', left_index=True, right_on='TimePoint')

    if len(merged_df) != bold_vols:
        raise ValueError(f"Merged data length ({len(merged_df)}) doesn't match volumes ({bold_vols}).")

    return merged_df


def extract_time_series_values(behave_df: pd.DataFrame, time_series_array: np.ndarray, delay: int):
    """
    Extracts time series data from the provided timeseries BOLD area for associated behavioral data
    that is acuiqred from trlocked_events w/ specified delay

    Parameters:
        behave_df (pandas.DataFrame): DataFrame containing behavioral data with a 'TimePoint' column
            indicating the starting point for each time series extraction.
        time_series_array (ndarray): Numpy Array containing time series data.
        delay (int): Number of data points to include in each extracted time series.

    Returns:
        np.ndarray: Array containing the extracted time series data for each time point in the behavioral DataFrame.
            Each row corresponds to a time point, and each column contains the extracted time series data.
    """
    extracted_series_list = []
    for row in behave_df['TimePoint']:
        start = int(row)
        end = start + delay
        extracted_series = time_series_array[start:end]
        extracted_series_list.append(extracted_series)
    return np.array(extracted_series_list)


def extract_time_series(bold_paths: list, roi_type: str, high_pass_sec: int, roi_mask: str = None,
                        roi_coords: tuple = None, radius_mm: int = None, bold_tr: float = None, detrend=True,
                        fwhm_smooth: float = None):
    """
    For each BOLD path, extract timeseries for either a specified mask or ROI coordinate. Mask and coordinate should be
    in same space/affine as BOLD data. Function leverages NiftiLabelsMasker (mask path) NiftiSpheresMasker (coordinates)
    to achieve this.

    Parameters:
        bold_paths (list): List of paths to subjects (list should match order of subs/runs/tasks for events file list).
        roi_type (str): Type of ROI ('mask' or 'coordinates').
        high_pass_sec (int): High-pass filter to use, in seconds. Used to convert to filter freq using 1/secs.
        roi_mask (str or None): Path to the ROI mask image (required if roi_type is 'mask').
        roi_coords (tuple or None): Coordinates (x,y,z) for the sphere center (required if roi_type is 'coordinates').
        radius_mm (int or None): Radius of the sphere in mm (required if roi_type is 'coordinates').
        bold_tr (float or None): TR value for acquisition of BOLD data.
        detrend: True/False, whether to use Nilearn's detrend function.
        fwhm_smooth (float or None): FWHM for spatial smoothing of data.

    Returns:
        list: List of time series for provided subjects/runs.
    """
    if roi_type not in ['mask', 'coordinates']:
        raise ValueError("Invalid ROI type. Choose 'mask' or 'coordinates'.")

    if roi_type == 'mask':
        masker = NiftiLabelsMasker(labels_img=roi_mask, standardize='psc', resampling_target='data',
                                   detrend=detrend, high_pass=1/high_pass_sec,
                                   t_r=bold_tr, smoothing_fwhm=fwhm_smooth)
        time_series = [masker.fit_transform(i) for i in bold_paths]
        return time_series

    else:
        masker_coord = NiftiSpheresMasker(seeds=[roi_coords], radius=radius_mm,
                                          standardize='psc', resampling_target='data',
                                          detrend=detrend, high_pass=1/high_pass_sec,
                                          t_r=bold_tr, smoothing_fwhm=fwhm_smooth)
        time_series_coord = [masker_coord.fit_transform(i) for i in bold_paths]
        return time_series_coord


def extract_postcue_trs_for_conditions(events_data: list, onset: str, trial_name: str,
                                       bold_tr: float, bold_vols: int, time_series: np.ndarray,
                                       conditions: list, tr_delay: int):
    """
    Extract TR coinciding with condition onset, plus TRs for specified delay for each file. Save this to a pandas
    dataframe (long) with associated Mean Signal value, for each subject, trial of condition and cue across the range
    of TRs (1 to TR delay)

    Parameters:
        events_data (list): List of paths to behavioral data files (list should match order for
        subs/runs/tasks as bold file list).
        onset (str): Name of the column containing onset values in the behavioral data.
        trial_name (str): Name of the column containing condition values in the behavioral data.
        bold_tr (int): TR for acquisiton of BOLD data.
        bold_vols (int): Number of volumes for BOLD.
        time_series (numpy.ndarray): numpy array of time series data.
        conditions (list): List of cue conditions to iterate over, min 1.
        tr_delay (int): Number of TRs to serve as delay (post onset)

    Returns:
        pd.DataFrame: DataFrame containing mean signal intensity values, subject labels,
            trial labels, TR values, and cue labels for all specified conditions.
    """
    dfs = []
    for cue in conditions:
        out_trs = []
        for index, beh_path in enumerate(events_data):
            subset_df = trlocked_events(events_path=beh_path, onsets_column=onset,
                                        trial_name=trial_name, bold_tr=bold_tr, bold_vols=bold_vols, separator='\t')
            trial_type = subset_df[subset_df[trial_name] == cue]

            out = extract_time_series_values(behave_df=trial_type, time_series_array=time_series[index], delay=tr_delay)
            out_trs.append(out)

        out_trs_array = np.array(out_trs)
        num_subs = out_trs_array.shape[0]
        num_trials = out_trs_array.shape[1]
        num_delay = out_trs_array.shape[2]
        reshaped_array = out_trs_array.reshape(-1, 1)

        df = pd.DataFrame(reshaped_array, columns=['Mean_Signal'])
        df['Subject'] = np.repeat(np.arange(1, num_subs + 1), num_trials * num_delay)
        df['Trial'] = np.tile(np.arange(1, num_trials + 1), num_subs * num_delay)

        tr_values = np.tile(np.arange(1, num_delay + 1), num_subs * num_trials)
        df['TR'] = np.tile(tr_values, 1)

        cue_values = np.repeat([cue], num_subs * num_trials * num_delay)
        df['Cue'] = np.tile(cue_values, 1)

        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)
