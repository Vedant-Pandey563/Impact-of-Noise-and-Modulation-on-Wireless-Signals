# Impact-of-Noise-and-Modulation-on-Wireless-Signals
This project demonstrates how different digital modulation schemes behave under Additive White Gaussian Noise (AWGN). Using Scilab, the simulation evaluates and visualizes the performance of BPSK, QPSK, and 16-QAM in terms of:

Bit Error Rate (BER)

Constellation distortion under noise

Waveform degradation

Analog carrier modulation effects

The simulation includes both static analysis and live animation, making it useful for learning, research, and wireless communication coursework.


------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Requirements

Scilab 6.1.1 or newer

Minimum 4GB RAM recommended (for animation and large symbol sets)

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

How to Run

1)Open Scilab
2)Load the .sci file into the editor
3)Press Run or execute:
    exec("Scilab_Simulation.sci", -1);
4)Watch the output figures update automatically.

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Output Generated

The script produces multiple live and static visualization windows, including:

  BER vs Eb/N0 curve
  
  Live QPSK and 16-QAM constellation animation
  
  Noisy signal visualization
  
  Eye diagram (BPSK)
  
  Analog carrier waveform (clean vs noisy)
  
  FFT spectrum plot
  
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Theory Summary

The performance under AWGN improves as Eb/N0 increases, with the following trade-offs:

| Modulation | Complexity | Bits/Symbol | Noise Sensitivity |
| ---------- | ---------- | ----------- | ----------------- |
| BPSK       | Low        | 1           | Best              |
| QPSK       | Medium     | 2           | Good              |
| 16-QAM     | High       | 4           | Worst             |

Higher-order modulation increases spectral efficiency but requires higher SNR to maintain the same BER.
