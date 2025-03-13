
The main script is "analyze_EODrecordings_testosterone.m" The rest of the ".m" files are called during the execution of this code.

This script takes raw EOD recordings and processes them. The main actions are: adjust baseline, normalize by peak-to-peak amplitude, center time (t = 0) at P1, determine the start and end points of the EOD, calculate EOD duration, averaged EODs, and matlab plots of EODs.

The R script "eodplots.R" plots the EOD duration data, tests for duration differences, and plots the days tested for duration differences.

