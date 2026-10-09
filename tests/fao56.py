"""Worked examples from FAO Irrigation and Drainage Paper 56 (Allen et al., 1998)."""

# Example 18: Brussels, 6 July, 50°48'N, 100 m. Wind 10 km/h at 10 m and
# n = 9.25 h of sunshine; the published intermediate values are below.
EXAMPLE_18 = {
    "lat": 50.80,
    "elev": 100.0,
    "doy": 187,
    "T_max": 21.5,
    "T_min": 12.3,
    "T_mean": (21.5 + 12.3) / 2.0,
    "RH_max": 84.0,
    "RH_min": 63.0,
    "RH_mean": (84.0 + 63.0) / 2.0,
    "u2": 2.078,
    "Rs": 22.07,
    "P": 100.1,
}
EXAMPLE_18_RESULTS = {
    "es": 1.997,
    "ea": 1.409,
    "delta": 0.122,
    "gamma": 0.0666,
    "Ra": 41.09,
    "N": 16.1,
    "Rso": 30.90,
    "Rnl": 3.71,
    "Rn": 13.28,
    "ETo": 3.9,
}
