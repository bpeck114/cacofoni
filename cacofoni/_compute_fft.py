# FILE: _compute_fft.py

# Import packages
import numpy as np

def _compute_fft(signal_array, apply_hanning, hann_window):
    """Compute the fast fourier transform for a signal array."""
    n_steps = signal_array.shape[0]
    
    if apply_hanning:
        signal_array = signal_array * hann_window
        print("Applying a Hann window filter to centroids to reduce edge effects...\n")

    # IDL does: FFT(nsamp * signal), so multiply before
    scaled = signal_array * n_steps
    fft_result = np.fft.fft(scaled, axis=0).astype(np.complex64)
    
    return fft_result / n_steps


# # FILE: _compute_fft.py

# # Import packages
# import numpy as np

# def _compute_fft(signal_array, apply_hanning, hann_window):
#     """Compute the fast fourier transform for a signal array."""
#     n_steps = signal_array.shape[0]
#     if apply_hanning:
#         signal_array = signal_array * hann_window
#         print(f"Applying a Hann window filter to centroids to reduce edge effects...\n")
#     fft_result = np.fft.fft(signal_array, axis=0).astype(np.complex64)
    
#     return fft_result / n_steps

