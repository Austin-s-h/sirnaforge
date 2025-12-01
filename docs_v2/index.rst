siRNAforge Documentation (v2)
=============================

**Design effective siRNAs in minutes.**

siRNAforge is a comprehensive toolkit for designing siRNA candidates with thermodynamic scoring, off-target analysis, and chemical modification support.

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

   commands_live
   commands
   scoring

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


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
