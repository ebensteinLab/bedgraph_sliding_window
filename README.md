sliding_window script:

This script has two parameters: 
1) Window size
2) Sample interval

The sliding window body is initially positioned at t_0 (represents the first position). 
When the first iteration starts, the sliding window heads forward one SI. 
Simultaneously, the original window decreases one SI at the bottom so that the total window size always remains the same length. 
In this way, the sliding window body moves to a new position t_1 (represents the second position). 
The script keeps running and after n-1 sample intervals, the sliding window will move to a new position t_n-1. 
In each interval, the average read counts are calculated and sets as the read count of the posotion in the middle of the interval.


#Input files:

chromozomes size txt file line 44: chromesizes_path

input bedgrah file line 47: input_bg_path

output destiny file path line 48: output_bg_path

#Parameters:

Window size in bp line 49: bin_size

Sample intrval i bp line 82: bin_start
