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
