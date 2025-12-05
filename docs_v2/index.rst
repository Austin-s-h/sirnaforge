siRNAforge Documentation
========================

**Design effective siRNAs in minutes.**

siRNAforge is a comprehensive toolkit for designing siRNA candidates with thermodynamic scoring, off-target analysis, and chemical modification support.

.. note::

   This documentation is **self-generating**: CLI examples show live output from actual command execution, and API documentation is extracted directly from source code docstrings.

Quick Start
-----------

.. code-block:: bash

   # Install
   pip install sirnaforge

   # Design siRNAs for any gene
   sirnaforge workflow TP53 --output-dir results/

   # View your top candidates
   cat results/sirnaforge/TP53_pass.csv


.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   workflows

.. toctree::
   :maxdepth: 2
   :caption: Reference

   commands
   scoring
   models_and_scoring
   api_autodoc

.. toctree::
   :maxdepth: 2
   :caption: Advanced

   modifications
   api


Key Features
------------

* **Thermodynamic scoring** - Asymmetry analysis for guide strand selection
* **Multi-species off-target** - BWA-MEM2 alignment across genomes
* **Chemical modifications** - 2'-O-methyl, phosphorothioate patterns
* **Batch processing** - Analyze multiple genes efficiently
* **Docker/Nextflow** - Production-ready containerized workflows


Documentation Philosophy
------------------------

This documentation follows the principle of **single source of truth**:

* **CLI Reference** → Generated live from ``sirnaforge --help`` during build
* **API Reference** → Auto-generated from Python docstrings and type hints
* **User Guides** → Focus on concepts and workflows, not redundant details

The code documents itself through:

* Comprehensive type annotations (validated by mypy)
* Detailed docstrings following Google/NumPy style
* Rich CLI help text via Typer


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
