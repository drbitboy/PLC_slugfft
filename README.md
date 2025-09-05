# Characterizing FFTs of pulse data

* See https://www.plctalk.net/forums/threads/has-anyone-done-fast-fourier-transform-in-structured-text.147547/

* Sample images
  * 20s pulse every 100s
    * Lowest-frequency peak is 0.01Hz = 1/100s
    * ![](https://github.com/drbitboy/PLC_slugfft/blob/master/images/02.png?raw=true)
    * ![](https://github.com/drbitboy/PLC_slugfft/blob/master/images/02_zoom.png?raw=true)
  * 20s pulse every 33s
    * Lowest-frequency peak is 0.03Hz ~ 1/33s
    * ![](https://github.com/drbitboy/PLC_slugfft/blob/master/images/01.png?raw=true)
    * ![](https://github.com/drbitboy/PLC_slugfft/blob/master/images/01_zoom.png?raw=true)
  * 27s pulse every 25s (overlapping)
    * Lowest-frequency peak is 0.04Hz = 1/25s
    * ![](https://github.com/drbitboy/PLC_slugfft/blob/master/images/00.png?raw=true)
    * ![](https://github.com/drbitboy/PLC_slugfft/blob/master/images/00_zoom.png?raw=true)

