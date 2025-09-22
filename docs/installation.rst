Installation
~~~~~~~~~~~~

Pyrangeyes requires Python ≥ 3.12 and can be installed using pip.

Pyrangeyes supports two alternative graphical libraries ("engines"): plotly and matplotlib.
At least one must be installed. Use these commands to install Pyrangeyes together with
your engine of choice::

    pip install pyrangeyes[plotly]

    pip install pyrangeyes[matplotlib]

To install both engines, use instead::

    pip install pyrangeyes[all]

Note that the minimal installation by :code:`pip install pyrangeyes` is not able to produce
plots since the graphical dependencies are not installed.
