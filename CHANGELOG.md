# Changelog

## v.1.1

**Changes:**

- Sensitivity threshold added as an external parameter
- Default value for sensitivity threshold changed to 7
- Each HMM in contig counts only once (to mitigate issue with repeated domains)
- New blast parsing (outfmt6, no blast-based prediction, best hit reported)
- Python 3.6+ required
- Logging added
