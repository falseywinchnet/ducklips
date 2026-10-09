"""Original magma display mapping, with one output pixel per enhanced field cell."""
import numpy as np
import matplotlib


def enhanced_rgba(field, reference, db_range):
    if not (np.isfinite(reference) and reference > 0 and np.isfinite(db_range) and db_range > 0):
        raise ValueError('Display reference and range must be finite and positive')
    if field.ndim != 2 or not np.all(np.isfinite(field)):
        raise ValueError('Expected a finite two-dimensional field')
    db = 20 * np.log10(np.maximum(field / reference, 10 ** (-db_range / 20)))
    normalized = np.clip((db + db_range) / db_range, 0, 1)
    # Flip only for top-origin image coordinates; preserve every row and column.
    return matplotlib.colormaps['magma'](normalized[::-1], bytes=True)
