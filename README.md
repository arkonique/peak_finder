# Peak finding algorithm

This uses a couple of different techniques to find peaks in a time series signal with a very wide range of frequencies.

It uses two different techniques to find two different kinds of peaks (designed for seismic velocity signals of two kinds -

 - Very slow lasting for days aseismic peaks
 - Very fast lasting a few minutes seismic peaks

The two methods are -
- A statistical z-score based peak finder (https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data)
- A recursive peak removal algorithm

## Method
Final algorithm combines the strengths of the two algorithms (Z-score moving thresholds & recursive max-jumping) and a nearest neighbour clustering

For aseismic peaks

1. Use z-score thresholding to find a square wave signal. Each duration of the signal=1 is called a span
2. The troughs are labelled as -1, convert them to zero because trough information is not important
3. Now expand all of the spans to be at least as wide as the aseismic event window. They can be wider if detected by the ZST algoritnm.
4. Combine any overlapping spans
5. Steps 3 and 4 should remove any high frequency peaks that the algorithm finds
6. Find the max value in each span, that is one peak
7. Remove any peaks above the seismic event limit (i.e. only save aseismic events) because the aseismic event window is so long, there will be multiple earthquakes that happen within one window, and every earthquake except the largest one will be lost.
8. This method is good at finding aseismic peaks because of it's ability to remove any clusters by only taking the max value. This does not work for shorter events like earthquakes because a shorter window introduces chances of getting high frequency peaks, and since the window will be so small, all of those high frequency peaks will get counted.
 
 For seismic peaks
1. Use recursive max jumping to find all peaks using some specified seismic event  window
2. Remove aseismic peaks
3. This method is good at finding seismic peaks easily because above the seismic threshold, within a seismic window, the high frequency peaks are very few. So it can easily just keep all the peaks. Earthquakes have small enough windows that if separate peaks are detected, they are likely to be separate earthquakes. This doesn't work as well with aseismic peaks is because the window is very big, but the standard deviation of aseismic event duration is very large. This means that a big enough window that doesn't cause one 3 month long aseismic transient to appear as multiple will also cause two 5 day long transients to appear as one. ZST does not have that problem because apart from a minimum window size, it can adjust the width accordingly.
4. There might still be some clusters left which can then be removed by a NN filter using the formula -
    `((peak(n-1,time)/3.154e7-peak(n,time)/3.154e7)^2+(log10(peak(n-1,velocity))-log10(peak(n,velocity)))^2)^0.5`
5. This formula instead of using a regular x-y distance formula lets us get a value within 0-1 which is easier to deal with
6. The biggest weakness of this algorithm is that it needs a lot of tweaking to get the seismic and aseismic windows exactly right so that no peaks are missed and no extra peaks are detected. For my dataset, aseismic window works best at 10 days and seismic window works best at 300 seconds.


## Fail log
1. Using MATLAB findpeaks (problem - finds way too many peaks, a lot of which are useless)
2) Trying it by the logic that I find when it goes above my limit to where it comes below (problem - close enough peaks never come below the limit and they get counted as one)
3) Using gradients to find maxima (problem - higher frequency signal in there makes it useless)
4) Using a spline fit and then finding maxima (problem - too computationally expensive and a big loss of data at the peaks, making the max slip rate value wrong)
5) Using a lowpass filter to filter out the higher frequencies (problem - uneven sampling rate makes it impossible to successfully filter)
6) Even out sampling rate and then filter (problem - a lot, a. VERY computationally expensive, b. lowpass loses information exactly at peaks while highpass is not good enough to filter out just the peaks because of the Nyquist limit. Thus the tradeoff between even sampling and computational cost is very sharp and does not allow any filters or moving averages or envelopes)
7 ) Using MATLAB findpeaks on even sampling rate data with a useful minpeakdistance (problem - very expensive and still not a good enough filter {it found 728 peaks for a series with 3 events})
8 ) Use image processing to find lines in the plot saved as a picture (problems - if peaks are close enough together, this would count them as one peak, and also, if the peaks are slow enough, it will not count them as a peak at all)
9 ) Vertical top-down scanner (problem - no need to scan all y values, can just jump from max value to max value)
10) Recursive max jumping by removing peaks breaking signal into smaller signals before and after the peak. Best one yet, very fast. (problems - finds a lot of duplicates, but that can be fixed with a simple unique filter. Biggest problem is that slow slip events need a bigger window to be called one event {one aseismic transient can last 3-6 days} and using a large window slows down the code that many times {going from an event window of 150 to 5000}. Using a smaller window on the other hand leads to finding more than one peak during one event)
11) Wavelet transforms (problem- I had to learn hilbert spaces and some more linear algebra first, and second it seemed too time consuming to get to)
12) Z score based statistical peak detection (problem - similar to matlab findpeaks, this also finds too many peaks due to high frequency signals, although this method is significantly better)
13) Back to recursive max jumping. This algorithm barely takes any time, so it is easy to add another layer into it, i.e. run the function on the peaks that it finds. This way the first run sort of acts as an envelope, while the second run refines the peaks into only the reals ones that I need (problem - the problem of aseismic events having longer windows still remains, and one aseismic peak is detected as separate events if the window is not large enough)
14) Recursive max jumping with double pass can handle longer aseismic windows after I fixed the error of the infinite loop last night (problem - Double passing it introduces a new problem. If the aseismic event duration is long enough such that every event has just one peak from the first pass, the second pass creates an envelope and reduces the number of actual peaks)
15) Recursive max jumps with conditional double pass iff the number of peaks found is less than a predefined max number of peaks (problem - parameters need to be tweaked based on my dataset, mainly the aseismic event duration. Currently at 6 days. But the biggest issue is that if there are less than the maxn peaks in the first pass but more than one during the same peak, it counts them as multiple)
16) Combine recursive max jumps with z-score thresholding (clustered recursive max jumps), and remove peaks if there are multiple during one span (i.e. there is a cluster at one peak). Adjust window of aseismic event such that it does produce more than one peak during a transient, and then remove the extra ones using the thresholding.
