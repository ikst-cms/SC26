The code has three options:
(a) Calculate energy of any given configuration (in POSCAR format) or a directory
    of configurations.
(b) Perform random walk through configuration energy space. This would generate a 
    histogram of energies. This is Python-based and hence very slow. In future, this 
    will be discontinued.
(c) Perform finite-T Monte Carlo simulation.
    -Added CheckPoint. At the end of each temperature, a CheckPoint.mc file is created which can
     be used to restart from that point in case of any termination of the code.


Installation 
------------

(1) Unzip the package. 

(2) Copy _spglib.so and _xtalcomp.so into ./cogue. These two must be build separately 
    from cogue 
    NOTE::For making life easy, cogue-conform can be found in the package. Use the setup.py to build
    cogue. Transfer the above two files to ./cogue directory.

(3) In the directory C_Functions, run: ./compil (if required, change the permission using chmod)
    This will craete C_MCinterface.so in lib folder in the base-directory.

To start the calcualtion, create a separate run directory and issue command: 
  python <path>/CPyMonC/CPyMonC


NOTE on running: 
--------------

(1) Use following command at the prompt:
  python <path>/CPyMonC/CPyMonC

(2) In case of error:
    >>>>>>>  Error :: Unable to get hash-key. Exit.   Exit.  <<<<<<<<
   Do followin:
    - Open file "compil" in the directory C_Functions.
    - Set flag at the start OLD_HASH to True OLD_HASH=True
    - recompile using ./compil  
   This will force the code to use old hash method which is more robust.
