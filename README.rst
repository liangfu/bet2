# bet2
standalone Brain Extraction Tool (bet2), released by http://fsl.fmrib.ox.ac.uk/

The following features are added:

 * the program is reconfigured using CMake
 * ``vtkzlib`` is added to support generating `gzip`-ed nifti image files

to compile source code:

.. code-block:: bash

 mkdir build
 cd build
 cmake ..
 make -j6


.. image:: https://badges.gitter.im/liangfu/bet2.svg
   :alt: Join the chat at https://gitter.im/liangfu/bet2
   :target: https://gitter.im/liangfu/bet2?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge