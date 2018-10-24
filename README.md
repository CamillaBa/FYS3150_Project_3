# FYS3150_Project_3

The file "source1.cpp" contains all the c++ code used in this project.


Different section of the main program generate different plot data.
For example, the section left uncommented now is:

//=========================================================================================================
//	 Study effect of Jupiters gravity on earth
//=========================================================================================================

When the program runs, it generates plot data for the jupiter-earth-sun system.
These must then manually be put into the folder "python_plotting/.." where
they can be used to make plots using "PythonApplication1.py".



The file source1.cpp contains the following functions:


- double ** sparse(int n)                             : creates a double array of doubles
- void clear_memory                                   : clears the memory from the above double array
- double* zeros(int n)                                : creates an array of doubles
- void initialize(body_list & bodies,                 : given a correctly formatted file, 
                  std::string filename)               : it initiates body_list with the values from that file
- void write_conserved_quantities_to_file(.....)      does at it says


- void unit_test_1(std::string method_name)
- void unit_test_converved_quantities(std::string method_name) 

The unit tests are discussed in the report.


Moreover, we have the class "body" and the class "planetary_system".
The contents of these classes are also discussed in the report.
