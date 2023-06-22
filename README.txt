File list:
1. ps_main: main code
2. ps_tran: identifying transitions
3. repackage: separate contours into 2 level sets and construct graph
4. drawBare: make the plots
5. plotByPieces: plot contours with breaks for NaNs
6. contoursep: change MATLAB contour format into NaN-separated
7. smush: bring PS list into format for tracking
8. drawGraph: display the graph at a certain frame
9. diff2, fitDiff, and fastDiff - not used for differentiation (at the moment)
10. drawOnly – not currently used
11. sideBySide – not currently used

Auxiliary functions
1. opread
2. intersections
3. track


Tunable variables:

ps_main:
1. opts.mode – which level sets to use
2. opts.exp – which kind of recording
3. opts.fou – which kind of smoothing of distance functions
4. st – start frame of analysis
5. L – length of analysis
6. de – number of frequencies thrown out in high-pass filter
7. deltat – minimum distance between peaks
8. MPP (mode 1) or MPP1, MPP2 (mode 2) – minimum peak prominence
9. sig – vector sigma for smoothing of high-passed data
10a. sig2 - vector sigma for smoothing of derivatives of smoothed data
10b. twind - size of time window for fit differentiation
11. rad (fou 1 only) – radius of Fourier mask
12. sigD (fou not 1) – sigma of 2D Gaussian
13. dis – minimum allowable distance of PS from edge/mask
14. record – whether movie of analysis is made
15. trQ – whether to rerun tracking
16. param.mem – how long a PS may disappear
17. sv - velocity scaling factor in tracking

ps_tran:
18. md – max distance between paired created/annihilated PSs
19. lm - minimum length of a PS track used
20. timeout – how far back we go to attempt to classify transition
21. l – ratio of ad to dist that is too "long"

repackage:
22. dmax – max jump along contour loop


