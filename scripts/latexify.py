import matplotlib
from math import sqrt

savefig_transparency = True
savefig_bbox = None
savefig_padinches = 0.0

def latexify(fig_width=None, fig_height=None, columns=1):
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}
    """
    # Ensure column input is valid
    assert columns in [1, 2], "Columns must be 1 or 2"

    # Set default figure width
    if fig_width is None:
        fig_width = 3.39 if columns == 2 else 6.9  # Width in inches

    # Set default figure height using golden ratio
    if fig_height is None:
        golden_mean = (sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
        fig_height = fig_width * golden_mean  # Height in inches

    # Limit maximum height to prevent oversized plots
    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        print(f"WARNING: fig_height too large: {fig_height}, reducing to {MAX_HEIGHT_INCHES} inches.")
        fig_height = MAX_HEIGHT_INCHES

    # Join preamble list into a single string
    preamble = r'\usepackage{gensymb},\usepackage{siunitx},\sisetup{detect-family=true},\usepackage{amsmath}'

    # Set rcParams for LaTeX compatibility
    params = {
        'text.usetex': True,  # Use LaTeX for rendering
        'pgf.preamble': preamble,  # LaTeX preamble as a single string
        'font.size': 10,  # Font size
        'legend.fontsize': 10,  # Legend font size
        'legend.handlelength': 2,
        'axes.labelsize': 10,  # Axis label size
        'axes.titlesize': 10,  # Title size
        'xtick.labelsize': 10,  # X-tick label size
        'ytick.labelsize': 10,  # Y-tick label size
        'figure.figsize': [fig_width, fig_height],  # Figure size
        'font.family': 'serif',  # Use serif fonts
    }

    matplotlib.rcParams.update(params)
