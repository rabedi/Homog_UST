Data From Dr Soghrati:
1-put them in a folder named (InputData)


2-open config folder


3-open "config_Read.txt" and set your targets for C++ code processor
in this folder you set which input files you want to precess (should type some properties of those files and their actual folder name which is put in "InputData")


4-open C++ code in visual studio and hit run or run it in g++ compiler
Be careful some parameters are hard-coded so you need to set them before hand. 
    num_sol = 3;//number of solution (loading)  it is always 3 since you need three independent solution for finding superposition coefficients
	double dTeta = 5; // it is the increment in angle sampling based on degree
	double dSpcAng = 30;//specific angle for presenting angle-sve  (not important)
	inclusionProp.density = 3210;//inclusion is SIC and its number is 1
	inclusionProp.tensileStrength = 359.0e6;
	inclusionProp.compresiveStrength = 2100.0e6;
	matrixProp.density = 6090;//matrix is Zr and its number is 0
	matrixProp.tensileStrength = 381.0e6;
	matrixProp.compresiveStrength = 2500.0e6;
	
	
5-The program prints output in the "output" folder (will be use for matlab ploting)



Matlab Processing:
1-need some data in "output" folder (done at fist task)
2-open "confing" folder 
3-open "config.txt" and set what data you ant to precess (plots)
-first you need to set how many data you want to process (first number & do not change the second number)
-then you need to set number of field for processing  (if you open one of the corresponding txt files you will see there are several fields are assigned to each SVE such as Elastic Modulus, density, fracture parameters and ....)
-then insert each of data set properties and name
4-open "MatlabFiles" folder and then main function.
run main function and comment out each part you want to see its results (same as Dr Acton Code)

5-The results will be generated in the folder "OutputPlot"