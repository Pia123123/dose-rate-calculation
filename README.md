# dose-rate-calculation

System requirements
A Version of python as well as the numpy-library has to be installed. The code was tested with python 3.14. Non-standard hardware is not required.

Installation guide
Make sure the input data files *.txt are in the same folder as the code doserate_caused_by_Am_or_Cm.py and doserate_caused_by_Pm.py.
Run the python code:
via console
Navigate to the scripts directory and run: python “doserate_caused_by_Am_or_Cm.py” for the Am- respectively Cm-solution or “doserate_caused_by_Pm.py” for the Pm-solution.
via editor
Open the script in an Integrated Development Environment (like Geany) and click the Run button.

Demo
At the head of the code, you can choose the isotope (only for the file “doserate_caused_by_Am_or_Cm.py”) and the cell geometry. 
Since every cells is unique in shape and size and none of them fits into a standard-geometry, the geometry can be approximated by two standard geometries. The results are in good agreement.
The possible geometries are a cylinder or a cuboid. If you choose cuboid, only the activity per volume will be calculated with this geometry, the doserate will be calculated using a cylinder.
Thereafter run the code.
The expected run time amounts 4 s for the code “doserate_caused_by_Am_or_Cm.py” and 20 s for the code “doserate_caused_by_Pm.py”.

Instructions for use
The results given in the paper were calculated by using the geometry “cuboid”.
